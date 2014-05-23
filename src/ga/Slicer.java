/* Genetic algorithm for crystal structure prediction.  Will Tipton.  Ceder Lab at MIT. Summer 2007 */
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

import java.util.List;
import java.util.Random;
import java.util.ArrayList;

import crystallography.Cell;
import crystallography.Site;
import crystallography.SupercellOptimizer;

import java.io.*;

import utility.Constants;
import utility.Vect;

// Slicer is a Variation operation.  It selects two organisms
// from the parent generation and creates an offspring organism.  The lattice
// parameters of the offspring are averages of those of its parents.  To choose
// the species (atom type and location) from the parents which will be passed
// to the child, we choose an axis along which to make a 'cut'.  For each site,
// we compute a function of the other two coordinates which is either a 
// constant function or cell-periodic. All atoms from one parent with fractional
// coordinate on the chosen axis greater than a certain distance away from the 
// function's value are copied to the child.  Atoms from the other parent with
// fractional coordinate on the chosen axis greater than that distance away are
// also copied to the child.

public final class Slicer implements Variation {
	
	// fraction of times we shift cells along the main and minor axes before crossover
	private double mainShiftFrac;
	private double minorShiftFrac;
	
	// used to determine the thickness of the slice we take from one parent.
	// (the thickness of the other slice is 1 minus that of the first, of course)
	private double thicknessMean;
	private double thicknessSigma;
	
	// if maxAmplitude is nonzero, we are doing a periodic type slicing
	private double maxAmplitude;
	private int maxFreq;
	
	// these constants are recalculated upon _each_ call to doVariation
	private double amplitude;
	private double split;
	private double etaFreq;
	private double zetaFreq;
	private double thickness;
	private int axis;
	private boolean growParents;
	private double doublingProb;
	
	public Slicer(List<String> args) {
		if (args == null || args.size() < 6)
			GAParameters.usage("Not enough parameters given to Slicer", true);
		
		thicknessMean = Double.parseDouble(args.get(0));
		thicknessSigma = Double.parseDouble(args.get(1));
		
		mainShiftFrac = Double.parseDouble(args.get(2));
		minorShiftFrac = Double.parseDouble(args.get(3));
		
		maxAmplitude = Double.parseDouble(args.get(4));
		maxFreq = Integer.parseInt(args.get(5));
		
		if (args.size() >= 7)
			growParents = Boolean.parseBoolean(args.get(6));
		else
			growParents = false;
		
		if (args.size() >= 8)
			doublingProb = Double.parseDouble(args.get(7));
		else
			doublingProb = 0;
	}
	
	public String toString() {
		StringBuilder result = new StringBuilder();
		
		result.append("Slicer variation. major: " + mainShiftFrac + " minor: " + minorShiftFrac);
		
		return result.toString();
	}
	
	/*
	 	private Cell doubleCell(Cell c) {
		Random rand = GAParameters.getParams().getRandom();

		// make identity matrix
		List<List<Integer>> coefs = new ArrayList<List<Integer>>();
		for (int i = 0; i < Constants.numDimensions; i++) {
			List<Integer> newVect = new ArrayList<Integer>();
			newVect.add(0); newVect.add(0); newVect.add(0);
			newVect.set(i, 1);
			coefs.add(newVect);
		}
		
		// double one of the dimensions
		int doubleDim = rand.nextInt(Constants.numDimensions);
		for (int i = 0; i < Constants.numDimensions; i++)
			coefs.get(doubleDim).set(i, coefs.get(doubleDim).get(i) * 2);
				
		// make the new offspring
		return Cell.getSupercell(c, coefs);
	}
	 */
	private Cell doubleCell(Cell c) {
		Random rand = GAParameters.getParams().getRandom();
		int doubleDim = rand.nextInt(Constants.numDimensions);
		
		List<Vect> newVects = new ArrayList<Vect>(c.getLatticeVectors());
		newVects.set(doubleDim, newVects.get(doubleDim).scalarMult(2.0));
		
		List<Site> newSites = new ArrayList<Site>(c.getSites());
		
		for (Site s : c.getSites())
			newSites.add(new Site(s.getElement(), s.getCoords().plus(c.getLatticeVectors().get(doubleDim))));
		
		Cell result = new Cell(newVects, newSites, c.getLabel());
		
		if (GAParameters.getParams().usingNiggliReducedCell())
			return result.getNigliReducedCell();
		else
			return result;
	}

	public Organism doVariation(Generation parents, Generation offspring, Selection sel) {
		GAParameters params = GAParameters.getParams();
		Random rand = params.getRandom();
		
		// calculate the parameters of this slicing
		amplitude = maxAmplitude * rand.nextDouble();
		etaFreq = Math.PI * rand.nextInt(2*maxFreq+1);	
		zetaFreq = Math.PI * rand.nextInt(2*maxFreq+1);
		axis = rand.nextInt(3);
		split = rand.nextDouble();
		thickness = GAUtils.renormZeroOne(thicknessMean + thicknessSigma * rand.nextGaussian());
		
		GAOut.out().stdout("amp:"+amplitude+" etaFreq:"+etaFreq+" zetaFreq:"+zetaFreq
				+" axis:"+axis+" split:"+split+" thickness:"+thickness, GAOut.DEBUG);

		//choose two random parents' Structures
		StructureOrg[] ps = (StructureOrg[])(sel.doSelection(parents, 2));
		Cell a = ps[0].getCell();
		Cell b = ps[1].getCell();
		
		// some status info
		GAOut.out().stdout("Mating "+ps[0].getID()+" (fitness " + ps[0].getFitness() + ", "
				+ a.getNumSites()+ " atoms) with " + ps[1].getID() +" (fitness "
				+ ps[1].getFitness() + ", " + b.getNumSites() + " atoms)... ", GAOut.NOTICE);
		
		// possibly double one of the parents
		if (rand.nextDouble() < doublingProb) {
			if (rand.nextBoolean()) {
				GAOut.out().stdout("(actually doubling first parent).", GAOut.NOTICE);
				a = doubleCell(a);
			} else {
				GAOut.out().stdout("(actually doubling second parent).", GAOut.NOTICE);
				b = doubleCell(b);
			}
		}
		
		// possibly take a supercell of one of the parents
		if (growParents) {
			if (a.getBasisSize() > b.getBasisSize()) {
				int sizeMultiple = a.getBasisSize() / b.getBasisSize(); // integer division
				// TODO: make be input options and similarly below
				if (sizeMultiple > 1) {
					GAOut.out().stdout("(actually using supercell of second parent with " + b.getNumSites() + "atoms).", GAOut.NOTICE);
					try {
						Cell newCell = Cell.getSupercell(b, SupercellOptimizer.getOptimalSupercell(b, false, 8 , sizeMultiple, params.getMaxLatticeLength(), false));
						b = newCell;
					} catch (Exception x) {
						GAOut.out().stdout("(actually supercelling of second parent failed).", GAOut.NOTICE);
					}
				}
			} else {
				int sizeMultiple = b.getBasisSize() / a.getBasisSize();
				if (sizeMultiple > 1) {
					GAOut.out().stdout("(actually using supercell of first parent with " + a.getNumSites() + "atoms).", GAOut.NOTICE);
					try {
						Cell newCell = Cell.getSupercell(a, SupercellOptimizer.getOptimalSupercell(a, false, 8 , sizeMultiple, params.getMaxLatticeLength(), false));
						a = newCell;
					} catch (Exception x) {
						GAOut.out().stdout("(actually supercelling of first parent failed).", GAOut.NOTICE);
					}
				}
			}
		}
				
		//average parents' bases to create a new basis
		double[] aAngles = a.getCellAnglesDegrees();
		double[] bAngles = b.getCellAnglesDegrees();
		double[] aLengths = a.getCellLengths();
		double[] bLengths = b.getCellLengths();
		double[] newAngles = new double[3];
		double[] newLengths = new double[3];
		for (int j = 0; j < 3; j++) {
			newAngles[j] = (aAngles[j] + bAngles[j]) / 2.0;
			newLengths[j] = (aLengths[j] + bLengths[j]) / 2.0;
		}	
		
		// are we shifting this time?
		double[] shiftA = new double[3];
		double[] shiftB = new double[3];
		if (rand.nextDouble() < mainShiftFrac) {
			shiftA[axis] = rand.nextDouble();
			shiftB[axis] = rand.nextDouble();
		}
		if (rand.nextDouble() < minorShiftFrac) {
			shiftA[(axis+1)%3] = rand.nextDouble();
			shiftB[(axis+1)%3] = rand.nextDouble();
		}
		if (rand.nextDouble() < minorShiftFrac) {
			shiftA[(axis+2)%3] = rand.nextDouble();
			shiftB[(axis+2)%3] = rand.nextDouble();
		}
		
		//get points/species from parent a above the split and from b below
		List<Vect> basis = (new Cell(newLengths[0], newLengths[1], newLengths[2], newAngles[0], newAngles[1], newAngles[2], null, null)).getLatticeVectors();
		ArrayList<Site> newSites = new ArrayList<Site>();
		// parent a's contribution 
		for (int j = 0; j < a.getNumSites(); j++) {
			Site s = new Site(a.getSite(j).getElement(), a.getSite(j).getCoords().plusWRT(new Vect(shiftA), basis));
			if (isAbove(s)) 
				newSites.add(s);
		}
		// parent b's contribution
		for (int j = 0; j < b.getNumSites(); j++) {
			Site s = new Site(b.getSite(j).getElement(), b.getSite(j).getCoords().plusWRT(new Vect(shiftB), basis));
			if (!isAbove(s)) 
				newSites.add(s);
		}
	
		// make the new organism
		StructureOrg newOrganism = new StructureOrg(new Cell(newLengths[0], newLengths[1], newLengths[2], newAngles[0], newAngles[1], newAngles[2], newSites, null));
			
		// some status info
		GAOut.out().stdout("("+ newOrganism.getCell().getNumSites() +" atoms)", GAOut.NOTICE, newOrganism.getID());
		
		return newOrganism;
	}
	
	// if the site is within thickness/2 of the value of f(s), then we're above (or inside)
	// the cut.  otherwise not.
	private Boolean isAbove(Site s) {
		List<Double> coords = s.getCoords().getCartesianComponents();
		double fSplit = f(s);
		return Math.abs(coords.get(axis) - fSplit) < thickness/2
				|| Math.abs(coords.get(axis) - fSplit - 1) < thickness/2
				|| Math.abs(coords.get(axis) - fSplit + 1) < thickness/2;
	}
	
	// returns a double f(site) which is value of our semi-random cell-periodic function 
	// as a function of the non-axis coordinates of site
	private double f(Site site) {
		//
		List<Double> coords = site.getCoords().getCartesianComponents();
		Vect coordsVec = new Vect(coords);
		
		// let zeta and eta be the numbers of the two axes other than axis
		int zeta = 0, eta = 0;
		do {
			zeta = (zeta + 1) % 3;
		} while (zeta == axis);
		do {
			eta = (eta + 1) % 3;
		} while (eta == zeta || eta == axis);
		//make unit vectors in the off-slice directions
		double[] etaCoords = {0d,0d,0d};
		double[] zetaCoords = {0d,0d,0d};
		etaCoords[eta] = 1;
		zetaCoords[zeta] = 1;
		Vect etaVec = new Vect(etaCoords);
		Vect zetaVec = new Vect(zetaCoords);
		
		return split + amplitude/2 * (Math.sin(coordsVec.dot(etaVec)*etaFreq))*(Math.sin(coordsVec.dot(zetaVec)*zetaFreq));
	}
	
	// just for testing
	public static void main(String[] args) {
		GAParameters params = GAParameters.getParams();
		
		StructureOrg s1 = new StructureOrg(Cell.parseCif(new File("/home/wtipton/cifs/143.cif")));
		StructureOrg s2 = new StructureOrg(Cell.parseCif(new File("/home/wtipton/cifs/44.cif")));
		s1.setFitness(-1);
		s2.setFitness(-1);
		
		Generation parents = params.makeEmptyGeneration();
		parents.addOrganism(s1);
		parents.addOrganism(s2);
		
		/*
		String[] hsArgs = {"0.5", "0.0", "0", "0.00", "0", "0"};
		Variation p = new Slicer(hsArgs);
		String[] selArgs = {"2", "0"};
		Selection sel = new ProbDistSelection(selArgs);
		
		StructureOrg o = (StructureOrg)p.doVariation(parents, null, sel);
		System.out.println(o);
		//GAUtils.writeStringToFile(o.getCIF(), new File("offspring.cif"), false); 
		o.standardize();
		//GAUtils.writeStringToFile(o.getCIF(), new File("offspring_red.cif"), false); 
		System.out.println(o); */

	}
}
