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

package ga;

import java.io.File;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;
import java.util.Arrays;   

import chemistry.Element;
import crystallography.Cell;
import crystallography.Site;
import pdvisual.PDAnalyzer;
import utility.Triplet;
import utility.Utility;
import utility.Vect;
import vasp.VaspIn;
import vasp.VaspOut;

public class IslandObjFcn extends ObjectiveFunction {
	
	List<String> objFcnArgs; 
	StructureOrg org;      // the island that gets sandwiched between larger sheet structures
	StructureOrg sandwich; // the org will be sandwiched between two layers with this structure
	
	public IslandObjFcn(List<String> subArray, Organism o) {
		sandwich = new StructureOrg(VaspOut.getPOSCAR(subArray.get(0))); // the path to the poscar file given in the input file (must end in .POSCAR)
		sandwich.getCell().getCellWithAllAtomsInCell().rotatedIntoPrincDirs(); // make sure the sandwich slab is lying in the x-y plane
		objFcnArgs = Utility.subList(subArray, 1);
		org = (StructureOrg)o;
		// TODO: maybe reading the sandwich structure from a file every time the class is instantiated isn't the best idea. Might be better to save it somewhere (GAParams, perhaps?)
		//   after reading it the first time and just get it from there in the future. Then we'd only have to rotate/move atoms into the cell once.
	}
	
	/**
	 * Places the island on top of a sheet of the sandwich structure (so it's sandwiched between them due to PCBs), such that the centroid of the island
	 * is located at the island's location variable, which is in fractional coordinates of the sandwich's a and b lattice vectors, and the vacuum spacing
	 * between the island and the sandwich structure is equal to the value of the island's interlayerDist variable.
	 */
	void sandwichOrg() {
		// remove the vertical vacuum from the island and sandwich sheet, leaving only 0.5*interlayerDist above and below
		sandwich.setCell(removeVerticalVacuum(sandwich.getCell(), org.getInterlayerDist()));
		org.setCell(removeVerticalVacuum(org.getCell(), org.getInterlayerDist()), false);
		
		// make sure all the atoms are in the cell
		sandwich.setCell(sandwich.getCell().getCellWithAllAtomsInCell());
		org.setCell(org.getCell().getCellWithAllAtomsInCell(), false);
		
		// make the lattice vectors of the sandwiched structure
		List<Vect> basis = new ArrayList<Vect>();
		basis.add(sandwich.getCell().getLatticeVectors().get(0));
		basis.add(sandwich.getCell().getLatticeVectors().get(1));
		double c = sandwich.getCell().getCellLengths()[2] + org.getCell().getCellLengths()[2];   
		basis.add(new Vect(0.0, 0.0, c));
		
		// add the sites of the sandwich sheet, and specify that they shouldn't be relaxed
		List<Site> sites = new ArrayList<Site>();
		for (Site s : sandwich.getCell().getSites())
			sites.add(new Site(s.getElement(), s.getCoords(), false));
		
		// compute the horizontal shift vector by taking the difference between the island's location and centroid
		Vect shift_vector = org.getCell().getCentroid().subtract(org.getLocation());
		double x_shift = shift_vector.getCartesianComponents().get(0);
		double y_shift = shift_vector.getCartesianComponents().get(1);
		
		// compute the vertical shift, which is the length of the c lattice vector of the sandwich sheet
		double z_shift = sandwich.getCell().getCellLengths()[2];
		
		// add the sites of the island, but shifted by the calculated amounts in each dimension
		for (Site s : org.getCell().getSites())
			sites.add(new Site(s.getElement(), s.getCoords().plus(new Vect(x_shift, y_shift, z_shift))));
		
		// assign the sandwiched cell to the island
		org.setCell(new Cell(basis, sites, org.getLabel()), false);
	}
	
	
	/**
	 * Returns the cell with all vertical vacuum removed, except for interlayerDist/2 of vacuum above and below
	 */
	private Cell removeVerticalVacuum(Cell cell, double interlayerDist) {
		// get the vertical bounds
		double zmin = cell.getAtomBoxZ()[0];
		double zmax = cell.getAtomBoxZ()[1];
		double zlen = zmax - zmin;
		
		// make new list of sites where we subtract (0, 0, minz - interlayerDist/2) off all the old ones
		List<Site> newSites = new ArrayList<Site>();
		Vect minv = new Vect(0.0, 0.0, zmin - interlayerDist/2);
		for (Site s : cell.getSites())
			newSites.add(new Site(s.getElement(), s.getCoords().subtract(minv)));
		
		// make a new box which is xlen x ylen x zlen
		List<Vect> newBasis = new ArrayList<Vect>();
		newBasis.add(cell.getLatticeVectors().get(0));
		newBasis.add(cell.getLatticeVectors().get(1));
		newBasis.add(new Vect(0.0, 0.0, zlen + interlayerDist));
				
		return new Cell(newBasis, newSites, cell.getLabel());
	}
	
	/**
	 * Removes the sandwich sheet from below the island, and also removes excess vacuum spacing from around the island.
	 * Places the island in an orthogonal cell, with a, b, and c parallel to Cartesian x, y and z, respectively.
	 */
	void unsandwichOrg() {
		// in case the sandwiched org rotated a little during relaxation
		org.getCell().rotatedIntoPrincDirs();
		
		// get the vertical bounds
		double total_zmin = org.getCell().getAtomBoxZ()[0];
		double total_zmax = org.getCell().getAtomBoxZ()[1];
		double total_zlen = total_zmax - total_zmin;
		
		// compute the vertical vacuum spacing between the island and the sandwich sheet
		double interlayerSpace = org.getCell().getLatticeVectors().get(2).getCartesianComponents().get(2) - total_zlen;
		
		// set the island's interlayerDist to the vertical vacuum between the island and the sandwich sheet
		org.setInterlayerDist(interlayerSpace);
		
		// get the vertical bounds of the sandwich structure
		double sandwich_zmin = sandwich.getCell().getAtomBoxZ()[0];
		double sandwich_zmax = sandwich.getCell().getAtomBoxZ()[1];
		
		// compute the vertical threshhold for deciding which sites belong to the island
		double sandwich_zlen = sandwich_zmax - sandwich_zmin;
		double verticalThreshhold = interlayerSpace + sandwich_zlen;
		
		// get the sites we think belong to the island based on their vertical coordinates
		List<Site> island_sites = new ArrayList<Site>();
		for (Site s : org.getCell().getSites()) {
			if (s.getCoords().getCartesianComponents().get(2) > verticalThreshhold) {
				island_sites.add(new Site(s.getElement(), s.getCoords()));
			}
		}
		
		// set the cell of the island so it only contains the island sites (and not the sites from the sandwich slab)
		org.setCell(new Cell(org.getCell().getLatticeVectors(), island_sites, org.getLabel()), false);
		
		// get the centroid of the island (in Cartesian coordinates)
		Vect centroid = org.getCell().getCentroid();
		
		// set the island's location equal the x and y components of its centroid
		Vect location = new Vect(centroid.getCartesianComponents().get(0), centroid.getCartesianComponents().get(1), 0.0);
		org.setLocation(location);
		
		// get the bounds of the island
		double xmin = org.getCell().getAtomBoxX()[0];
		double xmax = org.getCell().getAtomBoxX()[1];
		double ymin = org.getCell().getAtomBoxY()[0];
		double ymax = org.getCell().getAtomBoxY()[1];
		double zmin = org.getCell().getAtomBoxZ()[0];
		double zmax = org.getCell().getAtomBoxZ()[1];
		
		double xlen = xmax - xmin;
		double ylen = ymax - ymin;
		double zlen = zmax - zmin;
		
		// get the largest minimum interatomic distance
		double mid;
		if (GAParameters.getParams().getMinInteratomicDistance() == -1) {
			mid =  maxMID();	
		} else {
			mid = GAParameters.getParams().getMinInteratomicDistance();
		}
		
		// make the new lattice
		List<Vect> basis = new ArrayList<Vect>();
		basis.add(new Vect(xlen + mid, 0.0, 0.0));
		basis.add(new Vect(0.0, ylen + mid, 0.0));
		basis.add(new Vect(0.0, 0.0, zlen + mid));
		
		//  shift all the sites into the new lattice
		Vect shift_vector = new Vect(xmin - mid/2, ymin - mid/2, zmin - mid/2);
		List<Site> sites = new ArrayList<Site>();
		for (Site s : island_sites) {
			sites.add(new Site(s.getElement(), s.getCoords().subtract(shift_vector)));
		}
		
		// assign the unsandwiched cell to the island
		org.setCell(new Cell(basis, sites, org.getLabel()), false);
	}
	
	
	/**
	 * Returns the largest per species minimum interatomic distance
	 * Precondition: the perSpeciesMID option has been used
	 * 
	 * TODO: this method is copied directly from SurfaceObjFcn. Think about where to put this so all the objective functions can get it.
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
		
		sandwichOrg();
		
		ObjectiveFunction energyFcn = ObjFcnFactory.getObjectiveFunctionInstance(org, objFcnArgs);

		// relax the cell - have to wait for it to finish before using results
		Thread t = energyFcn.evaluate();
		try {
			if (t != null)
				t.join();
		} catch (InterruptedException x) {
			GAOut.out().stdout("InterruptedException in energy calc thread in IslandObjFcn: " + x.getMessage(), GAOut.WARNING, org.getID());
		} catch (Exception x) {
			GAOut.out().stdout("ERROR: exception in IslandObjFcn:run", GAOut.CRITICAL);
			x.printStackTrace();
		}
		
		// updating structure w/ relaxed version is this done in EPA already
		//org.setValue((new PDAnalyzer(pdbuilder.getPDData())).getEnergyPerAtomAboveHull(org));
		
		unsandwichOrg();
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
		return "IslandObjFcn"; // TODO: more? 
	}
	
	
	// for testing
	public static void main(String args[]) {
		
		String arg[] = {"/n/srv/brevard/structures/VSe2_8x8.POSCAR", "other_args"};
		
		StructureOrg islandOrg = new StructureOrg(VaspOut.getPOSCAR("/n/srv/brevard/structures/2_trial_run/use_in_next_search/12580_relaxed.POSCAR"));
		
	//	islandOrg.setInterlayerDist(4.0);
	//	double interlayerDist = islandOrg.getInterlayerDist();

		
		// Convert array of strings to list of strings for the constructor
		List<String> larg = Arrays.asList(arg);
		
		IslandObjFcn iof = new IslandObjFcn(larg, islandOrg);
		
		iof.unsandwichOrg();
		
		islandOrg.getCell().writeCIF("/n/srv/brevard/structures/2_trial_run/use_in_next_search/12580.cif");
		
		(new VaspIn(Cell.parseCif(new File("/n/srv/brevard/structures/2_trial_run/use_in_next_search/12580.cif")), null, null, null)).writePoscar("/n/srv/brevard/structures/2_trial_run/use_in_next_search/12580.POSCAR", false);
		
		
		
		
		
		// set the location of the island org. Maybe can set location earlier, and just set the basis here...
		//iof.org.setLocation(new Vect(0.5, 0.5, 0.0, iof.sandwich.getCell().getLatticeVectors()));
		
		//iof.sandwichOrg();
		
		//iof.org.getCell().writeCIF("/n/srv/brevard/structures/sandwiched_island.cif");
		
		//iof.unsandwichOrg();
		
		//iof.org.getCell().writeCIF("/n/srv/brevard/structures/unsandwiched_island.cif");
		
		
		// try removing all the vacuum padding from sandwich org
		// iof.sandwich.setCell(iof.removeVerticalVacuum(iof.sandwich.getCell(), interlayerDist));
		
		// iof.sandwich.getCell().writeCIF("/n/srv/brevard/structures/sandwich_removed_vac.cif");
		
		// try removing all the vacuum padding from the island org
		//islandOrg.setCell(iof.removeVerticalVacuum(islandOrg.getCell(), interlayerDist));
		
		//islandOrg.getCell().writeCIF("/n/srv/brevard/structures/island_removed_vac.cif");
		
	}
}
