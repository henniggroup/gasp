package ga;

/*
 * Copyright 2011-2014 Will Tipton, Richard Hennig, Ben Revard, Stewart Wenner, Anna Yesypenko

This file is part of the Genetic Algorithm for Structure and Phase Prediction (GASP).

    GASP is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    GASP is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with GASP.  If not, see <http://www.gnu.org/licenses/>.
    
    
    */

import java.io.File;
import java.util.ArrayList;

import java.util.List;

import java.util.Arrays;   // Ben added this

import chemistry.Element;
import crystallography.Cell;
import crystallography.Site;
import pdvisual.PDAnalyzer;
import utility.Triplet;
import utility.Utility;
import utility.Vect;
import vasp.VaspIn;
import vasp.VaspOut;

public class SurfaceObjFcn extends ObjectiveFunction {
	
	List<String> objFcnArgs;
	StructureOrg org;
	double padding;
	
	public SurfaceObjFcn(List<String> subArray, Organism o) {
		padding = Double.parseDouble(subArray.get(0));
		objFcnArgs = Utility.subList(subArray, 1);
		org = (StructureOrg)o;
	}
	
	/**
	 * Rotates the cell into the principle directions and then adds padding to the magnitude 
	 * of the c lattice vector (which is now aligned vertically). Takes care of atomic sites 
	 * so that the relative positions of the atoms are the same in the padded cell as the original.
	 */
	void padOrg() {
		Cell oldCell = org.getCell();
		
		oldCell.rotatedIntoPrincDirs(); 
		
		List<Vect> basis = new ArrayList<Vect>();
		basis.add(oldCell.getLatticeVectors().get(0));
		basis.add(oldCell.getLatticeVectors().get(1));
		double c = oldCell.getCellLengths()[2];   
		basis.add(new Vect(0.0, 0.0, padding + c));
		
		List<Site> newSites = new ArrayList<Site>();
		for (Site s : oldCell.getSites())
			newSites.add(new Site(s.getElement(), s.getCoords().plus(new Vect(0.0, 0.0, (padding + c)/2))));	
		
		Cell newCell = new Cell(basis, newSites, oldCell.getLabel());
		org.setCell(newCell, false);
	}
	
	/**
	 * Rotates the cell into the principle directions and then removes vertical padding 
	 * from the cell without changing the relative positions of the atoms. The magnitude 
	 * of the vertical cell vector of the unpadded cell is equal to the vertical distance 
	 * between highest and lowest atoms in the cell plus the minimum interatomic distance 
	 * (mid). The atoms are arranged in the unpadded cell such that highest atom is mid/2 
	 * from the top of cell, and the lowest atom is mid/2 from the bottom of the cell.   
	 */
	void unpadOrg() {
		Cell oldCell = org.getCell();
		if (oldCell == null)
			return;
		
		oldCell.rotatedIntoPrincDirs();
		double mid; // minimum interatomic distance
		
		// get the largest minimum interatomic distance
		if (GAParameters.getParams().getMinInteratomicDistance() == -1) {
			mid =  maxMID();	
		} else {
			mid = GAParameters.getParams().getMinInteratomicDistance();
		}

		double [] bounds = getAtomBox(oldCell);
		double zlen = bounds[1] - bounds[0];
		
		// make new list of sites where we subtract (0,0,minz-mid/2) off all the old ones
		List<Site> newSites = new ArrayList<Site>();
		Vect minv = new Vect(0.0, 0.0, bounds[0] - mid/2);
		for (Site s : oldCell.getSites())
			newSites.add(new Site(s.getElement(), s.getCoords().subtract(minv)));
		
		// make a new box which is xlen x ylen x zlen plus the min interatomic distance
		List<Vect> newBasis = new ArrayList<Vect>();
		newBasis.add(oldCell.getLatticeVectors().get(0));
		newBasis.add(oldCell.getLatticeVectors().get(1));
		newBasis.add(new Vect(0.0, 0.0, zlen + mid));
		
		Cell newCell = new Cell(newBasis, newSites, oldCell.getLabel());
		org.setCell(newCell, false);
		
	}
	
	// returns a bounding box for the atoms (minz, maxz)
	double[] getAtomBox(Cell cell) {
		double minz = Double.MAX_VALUE;
		double maxz = Double.MIN_VALUE;
		for (Site s : cell.getSites()) {
			List<Double> cartComps = s.getCoords().getCartesianComponents();
			minz = Math.min(minz, cartComps.get(2));
			maxz = Math.max(maxz, cartComps.get(2));
		}
		return new double [] {minz, maxz};
	}

	/**
	 * Returns the largest per species minimum interatomic distance
	 * Precondition: the perSpeciesMID option has been used
	 */
	public double maxMID() {
		// Get an array of all the minimum interatomic distances
		List<Triplet<Element, Element, Double>> tripletList = GAParameters.getParams().getPerSpeciesMIDs();
		double[] MIDs = new double[tripletList.size()];
		for (int i = 0; i < tripletList.size(); i++) {
			MIDs[i] = tripletList.get(i).getThird();
		}
		// Searches through the array and returns the max value
		double mid = 0;
		for (int j = 0; j < MIDs.length; j++) {
			if (MIDs[j] > mid) {
				mid = MIDs[j];
			}
		}
		return mid;
	}
	
	
	public void run() {
		
		padOrg();
		
		ObjectiveFunction energyFcn = ObjFcnFactory.getObjectiveFunctionInstance(org, objFcnArgs);

		// relax the cell - have to wait for it to finish before using results
		Thread t = energyFcn.evaluate();
		try {
			if (t != null)
				t.join();
		} catch (InterruptedException x) {
			GAOut.out().stdout("InterruptedException in energy calc thread in SurfaceObjFcn: " + x.getMessage(), GAOut.WARNING, org.getID());
		} catch (Exception x) {
			GAOut.out().stdout("ERROR: exception in SurfaceObjFcn:run", GAOut.CRITICAL);
			x.printStackTrace();
		}
		
		// updating structure w/ relaxed version is this done in EPA already
		//org.setValue((new PDAnalyzer(pdbuilder.getPDData())).getEnergyPerAtomAboveHull(org));
		
		unpadOrg();
	}

	// dont update numCalculations when implementing this.. the underlying objfcn will do it
	public Thread evaluate() {		
		// short circuit here if we've done the calculation already
		if (org.knowsValue())
			return null;
		
		// start the calculation and return the Thread
		Thread t = new Thread(this);
		t.start();
		return t;
	}
	
	// an ObjectiveFunction should also overload toString();
	public String toString() {
		return "SurfaceObjFcn"; // TODO: more? 
	}
	
	
	// for testing
	public static void main(String args[]) {
		
		
		String arg[] = {"17", "hi"};
		
		StructureOrg c = new StructureOrg(VaspOut.getPOSCAR("/Users/benjaminrevard/Desktop/temp/9453.POSCAR"));

		
		// Convert array of strings to list of strings for the constructor
		List<String> larg = Arrays.asList(arg);
		
		SurfaceObjFcn cof = new SurfaceObjFcn(larg, c);
		
//		cof.padOrg();
		
//		c.getCell().writeCIF("/Users/benjaminrevard/GA/padding_testing/POSCAR.padded.cif");
		
//		cof.unpadOrg();
		
//		c.getCell().writeCIF("/Users/benjaminrevard/GA/padding_testing/POSCAR.unpadded.cif");
		
//		System.out.println(GAParameters.getParams().getMinInteratomicDistance());

	//	cof.unpadOrg();
		
//		c.getCell().writeCIF("/Users/benjaminrevard/GA/padding_testing/POSCAR.padded.cif");
		
		cof.unpadOrg();
		
		cof.padOrg();
		

		c.getCell().writeCIF("/Users/benjaminrevard/Desktop/temp/9453.cif");
				
		(new VaspIn(Cell.parseCif(new File("/Users/benjaminrevard/Desktop/temp/9453.cif")), null, null, null)).writePoscar("/Users/benjaminrevard/Desktop/temp/9453_padded.POSCAR", false);

		//c.getCell().writeCIF("/Users/benjaminrevard/GA/vasp/InP/ch5_initial_structs/8733.cif");
				
		//(new VaspIn(Cell.parseCif(new File("/Users/benjaminrevard/GA/vasp/InP/ch5_initial_structs/8733.cif")), null, null, null)).writePoscar("/Users/benjaminrevard/GA/vasp/InP/ch5_initial_structs/8733unpadded.POSCAR", false);
		
//		System.out.println(GAParameters.getParams().getMinInteratomicDistance());

	}

}
