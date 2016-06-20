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

package ga;

import java.util.*;

import utility.Vect;
import vasp.VaspOut;
import crystallography.*;
import chemistry.*;


// RandomSOCreator implements StructureOrgCreator.  It creates StructureOrgs with
// random lattice parameters and atomic positions within the hard constraints and 
// with a given number of stoichiometries.  Then, it optionally scales structures 
// to a particular initial volume.

public class RandomSOCreator implements StructureOrgCreator {
	
	private String newGenType;
	private double givenVolPerAtom;
	
	private boolean storedConstraints = false; // true if we're storing constraints, false if getting them from GAParameters
	private double maxll, minll, maxla, minla, maxh;
	private int maxna, minna;
	
	public RandomSOCreator(List<String> args, double maxll, double minll, double maxla, double minla, double maxh, int maxna, int minna) {
		init(args);
		storedConstraints = true;
		this.maxll = maxll;
		this.minll = minll;
		this.maxla = maxla;
		this.minla = minla;
		this.maxh = maxh;
		this.maxna = maxna;
		this.minna = minna;
		
	}
	
	public RandomSOCreator(List<String> args) {
		GAParameters params = GAParameters.getParams();
		
		init(args);
	}
	
	private void init(List<String> args) {
		// parse the options passed in
		if (args == null || args.size() < 1)
			GAParameters.usage("Not enough parameters given to RandomSOCreator", true);
		newGenType = args.get(0);
		
		if (newGenType.equalsIgnoreCase("randomVol")) {
			// nothing
		} else if (newGenType.equalsIgnoreCase("givenVol")) {
			if (args.size() != 2)
				GAParameters.usage("Incorrect number of parameters given to RandomSOCreator", true);
			givenVolPerAtom = Double.parseDouble(args.get(1));
		} else {
			GAParameters.usage("Unrecognized population type " + args.get(0), true);
		}
	}
		
	public String toString() {
		return "RandomSOCreator: " + newGenType + " with volume, " + givenVolPerAtom;
	}
	
	private void initCellConstraints() {
		GAParameters params = GAParameters.getParams();

		if (params.getMaxLatticeAngleDegrees() == -1 || params.getMinLatticeAngleDegrees() == -1
				|| params.getMaxLatticeLength() == -1 || params.getMinLatticeLength() == -1)
			GAParameters.usage("Error: Must set lattice parameter constraints to use makeRandomLattice().", true);
		
		maxll = params.getMaxLatticeLength();
		minll = params.getMinLatticeLength();
		maxla = params.getMaxLatticeAngleDegrees();
		minla = params.getMinLatticeAngleDegrees();		
		maxh = params.getMaxCellHeight();
		maxna = params.getMaxNumAtoms();
		minna = params.getMinNumAtoms();
	}

	// makes a basis which is random except that it satisfies the hard constraints
	public static List<Vect> makeRandomLattice() {
		GAParameters params = GAParameters.getParams();
		double maxll = params.getMaxLatticeLength();
		double minll = params.getMinLatticeLength();
		double maxla = params.getMaxLatticeAngleDegrees();
		double minla = params.getMinLatticeAngleDegrees();		
		double maxh = params.getMaxCellHeight();
		return makeRandomLattice(maxll, minll, maxla, minla, maxh);
	}
	
	public static List<Vect> makeRandomLattice(double maxll, double minll, double maxla, double minla, double maxh) {
		
		Random rand = GAParameters.getParams().getRandom();
		
		// sum of the angles needs to be < 360 degrees
		double adeg, bdeg, gdeg;
		do {
			adeg = (rand.nextDouble()*(maxla-minla)+minla);
			bdeg = (rand.nextDouble()*(maxla-minla)+minla);
			gdeg = (rand.nextDouble()*(maxla-minla)+minla);
		} while (!GAUtils.satisfiesTriangleInequality(adeg, bdeg, gdeg) || adeg + bdeg + gdeg >= 360);
		
		double la = rand.nextDouble()*(maxll-minll)+minll;
		double lb = rand.nextDouble()*(maxll-minll)+minll;
		
		//TODO: FIXME: maximum z lattice length
		// should be correctly constrained by maxCellHeight
		double maxZll = Math.min(maxll, 2*maxh);
		
		double lc = rand.nextDouble()*(maxZll-minll)+minll;
		
		return (new Cell(la,lb,lc,adeg,bdeg,gdeg,null,null)).getLatticeVectors();
	}
	
	private List<Vect> makeSubstrateLattice(Cell substrate, double maxll, double minll, double maxla, double minla, double maxh) {
		Random rand = GAParameters.getParams().getRandom();
		
		// sum of the angles needs to be < 360 degrees
		double adeg, bdeg;
		double gdeg = substrate.getCellAnglesDegrees()[2];
		do {
			adeg = (rand.nextDouble()*(maxla-minla)+minla);
			bdeg = (rand.nextDouble()*(maxla-minla)+minla);
		} while (!GAUtils.satisfiesTriangleInequality(adeg, bdeg, gdeg) || adeg + bdeg + gdeg >= 360);
		
		double la = substrate.getCellLengths()[0];
		double lb = substrate.getCellLengths()[1];
		
		//TODO: FIXME: maximum z lattice length
		// should be correctly constrained by maxCellHeight
		double maxZll = Math.min(maxll,2*maxh);
		
		double lc = rand.nextDouble()*(maxZll-minll)+minll;
		
		return (new Cell(la,lb,lc,adeg,bdeg,gdeg,null,null)).getLatticeVectors();
	}
	
	// fills up the population with random organisms (provided they satisfy the
	// hard constraints)
	protected StructureOrg makeRandomOrg() {
		GAParameters params = GAParameters.getParams();
		Random rand = params.getRandom();
		
		if (!storedConstraints)
			initCellConstraints();
		
		//HashMap<String,Integer> constituents = params.getConstituents();
		Composition comp = params.getCompSpace().getRandomIntegerCompInSpace(minna, maxna);

		// make random lattice parameters satisfying hard constraints
		// the Basis for our new organism
		List<Vect> latVects;
		if (params.usingSubstrate())
			latVects = makeSubstrateLattice(params.getSubstrate(), maxll, minll, maxla, minla, maxh);
		else
			latVects = makeRandomLattice(maxll, minll, maxla, minla, maxh);

		// make lists to hold the points and species which will make our new Structure
		ArrayList<Site> sitesList = new ArrayList<Site>();

		// loop through all the types of species
		final int maxFails = 100;
		int failCount = 0;
		for (Element e : comp.getElements())
			// add a stoichiometry's worth of that species
			for (int k = 0; k < comp.getOrigAmount(e); k++) {
				Vect potentialLocation = new Vect(rand.nextDouble(),rand.nextDouble(),rand.nextDouble(), latVects);
				if ((new Cell(latVects,sitesList)).getAtomsInSphereSorted(potentialLocation, params.getMinInteratomicDistance()).size() == 0)
					sitesList.add(new Site(e,potentialLocation));
				else if (failCount < maxFails) {
					failCount++;
					k--;
				}
			}
				
		Cell newStructure = new Cell(latVects, sitesList);
				
		// in the case of given volume, scale the cell to the desired volume
		if (newGenType.equalsIgnoreCase("givenVol"))
			newStructure = newStructure.scaleTo(givenVolPerAtom * newStructure.getBasisSize());

		StructureOrg org = new StructureOrg(newStructure);
		
		// If the island objective function is being used, need to assign random location and interlayer distance values within the constraints
		if (params.getObjFcnArgs().get(0) == "island") {
			// the constraints
			double maxloc = params.getMaxLocation();
			double minloc = params.getMinLocation();
			double maxinterlayer = params.getMaxInterlayerDist();
			double mininterlayer = params.getMinInterlayerDist();
			
			// get random values within the constraints
			double aloc = (rand.nextDouble()*(maxloc - minloc) + minloc);
			double bloc = (rand.nextDouble()*(maxloc - minloc) + minloc);
			double intelayer = (rand.nextDouble()*(maxinterlayer - mininterlayer) + mininterlayer);
			
			// get the sandwich slab from the poscar file and make sure it's lying in the x-y plane
			// TODO: this is also done every time IslandObjFcn is instantiated. Maybe reading it in from the poscar every time isn't the best approach...
			StructureOrg sandwich = new StructureOrg(VaspOut.getPOSCAR(params.getObjFcnArgs().get(1))); 
			sandwich.getCell().getCellWithAllAtomsInCell().rotatedIntoPrincDirs(); 
			
			// set the org's location using the random fractional coordinates and the basis of the sandwich sheet
			org.setLocation(new Vect(aloc, bloc, 0.0, sandwich.getCell().getLatticeVectors()));
			
			// set the org's interlayer distance to the randomly generated one
			org.setInterlayerDist(intelayer);
		}
		
		return org;
	}

	public StructureOrg makeOrganism(Generation g) {

		StructureOrg o = makeRandomOrg();
		
		return o;
	}
	
	// for testing
	public static void main(String args[]) {
		GAParameters params = GAParameters.getParams();
		
		Generation g = GAParameters.getParams().makeEmptyGeneration();
		
		String arg[] = {"--f", "/home/wtipton/test"};
		GAParameters.getParams().setArgs(arg);

	//	String socArgs[] = {"givenVol", "4", "40"};
		
	//	RandomSOCreator soc = new RandomSOCreator(socArgs);
		
	//	System.out.println(soc.makeOrganism(g));
		
		Composition comp = GAParameters.getParams().getCompSpace().getRandomIntegerCompInSpace(params.getMinNumAtoms(), params.getMaxNumAtoms());
		System.out.println(comp);
	}

}
