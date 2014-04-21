<<<<<<< HEAD
=======
/*
 * Copyright 2011-2014 Will Tipton, Richard Hennig, Ben Revard, Stewart Wenner

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

>>>>>>> 0e3189c40547bbd59ea42c4f91890d7511fb7797
package ga;

import java.util.ArrayList;
import java.util.List;

import crystallography.Cell;
import crystallography.Site;
import pdvisual.PDAnalyzer;
import utility.Utility;
import utility.Vect;
import vasp.VaspOut;

public class ClusterObjFcn extends ObjectiveFunction {
	
	List<String> objFcnArgs;
	StructureOrg org;
	double padding;
	
	public ClusterObjFcn(List<String> subArray, Organism o) {
		padding = Double.parseDouble(subArray.get(0));
		objFcnArgs = Utility.subList(subArray, 1);
		org = (StructureOrg)o;
	}
	
	private void padOrg() {
		Cell oldCell = org.getCell();
		
		List<Vect> basis = new ArrayList<Vect>();
		basis.add(new Vect(padding, 0.0, 0.0));
		basis.add(new Vect(0.0, padding, 0.0));
		basis.add(new Vect(0.0, 0.0, padding));
		
		List<Site> newSites = new ArrayList<Site>();
		for (Site s : oldCell.getSites())
			newSites.add(new Site(s.getElement(), s.getCoords().plus(new Vect(padding/2, padding/2, padding/2))));	
		
		Cell newCell = new Cell(basis, newSites, oldCell.getLabel());
		org.setCell(newCell, false);
	}
	
	private void unpadOrg() {
		Cell oldCell = org.getCell();
		if (oldCell == null)
			return;
		double mid = GAParameters.getParams().getMinInteratomicDistance();

		// get a bounding box for the atoms
		double minx = Double.MAX_VALUE;
		double miny = Double.MAX_VALUE;
		double minz = Double.MAX_VALUE;
		double maxx = Double.MIN_VALUE;
		double maxy = Double.MIN_VALUE;
		double maxz = Double.MIN_VALUE;
		for (Site s : oldCell.getSites()) {
			List<Double> cartComps = s.getCoords().getCartesianComponents();
			minx = Math.min(minx, cartComps.get(0));
			miny = Math.min(miny, cartComps.get(1));
			minz = Math.min(minz, cartComps.get(2));
			maxx = Math.max(maxx, cartComps.get(0));
			maxy = Math.max(maxy, cartComps.get(1));
			maxz = Math.max(maxz, cartComps.get(2));
		}
		double xlen = maxx - minx;
		double ylen = maxy - miny;
		double zlen = maxz - minz;
		
		// make new list of sites where we subtract (minx,miny,minz) off all the old ones
		List<Site> newSites = new ArrayList<Site>();
		Vect minv = new Vect(minx - mid/2, miny - mid/2, minz - mid/2);
		for (Site s : oldCell.getSites())
			newSites.add(new Site(s.getElement(), s.getCoords().subtract(minv)));
		
		// make a new box which is xlen x ylen x zlen plus the min interatomic distance
		List<Vect> newBasis = new ArrayList<Vect>();
		newBasis.add(new Vect(xlen + mid, 0.0, 0.0));
		newBasis.add(new Vect(0.0, ylen + mid, 0.0));
		newBasis.add(new Vect(0.0, 0.0, zlen + mid));
		
		Cell newCell = new Cell(newBasis, newSites, oldCell.getLabel());
		org.setCell(newCell, false);
		
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
			GAOut.out().stdout("InterruptedException in energy calc thread in ClusterObjFcn: " + x.getMessage(), GAOut.NOTICE, org.getID());
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
		return "ClusterObjFcn"; // TODO: more? 
	}
	
	// for testing
	public static void main(String args[]) {
		
		/*
		String arg[] = {"20", "hi"};
		
		StructureOrg c = new StructureOrg(VaspOut.getPOSCAR("/home/wtipton/POSCAR"));
		
		ClusterObjFcn cof = new ClusterObjFcn(arg, c);
		
		cof.padOrg();
		
		c.getCell().writeCIF("/home/wtipton/POSCAR.padded.cif");
		
		cof.unpadOrg();
		
		c.getCell().writeCIF("/home/wtipton/POSCAR.unpadded.cif");
		*/

	}
}
