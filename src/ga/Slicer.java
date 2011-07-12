/* Genetic algorithm for crystal structure prediction.  Will Tipton.  Ceder Lab at MIT. Summer 2007 */

package ga;

import java.util.List;
import java.util.Random;
import java.util.ArrayList;

import crystallography.Cell;
import crystallography.Site;
import crystallography.SupercellOptimizer;

import java.io.*;

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
	
	public Slicer(String[] args) {
		if (args == null || args.length < 6)
			GAParameters.usage("Not enough parameters given to Slicer", true);
		
		thicknessMean = Double.parseDouble(args[0]);
		thicknessSigma = Double.parseDouble(args[1]);
		
		mainShiftFrac = Double.parseDouble(args[2]);
		minorShiftFrac = Double.parseDouble(args[3]);
		
		maxAmplitude = Double.parseDouble(args[4]);
		maxFreq = Integer.parseInt(args[5]);
		
		if (args.length >= 7)
			growParents = Boolean.parseBoolean(args[6]);
		else
			growParents = false;
	}
	
	public String toString() {
		StringBuilder result = new StringBuilder();
		
		result.append("Slicer variation. major: " + mainShiftFrac + " minor: " + minorShiftFrac);
		
		return result.toString();
	}

	public Organism doVariation(Generation parents, Generation offspring, Selection sel) {
		GAParameters params = GAParameters.getParams();
		int verbosity = params.getVerbosity();
		Random rand = params.getRandom();
		
		// calculate the parameters of this slicing
		amplitude = maxAmplitude * rand.nextDouble();
		etaFreq = Math.PI * rand.nextInt(2*maxFreq+1);	
		zetaFreq = Math.PI * rand.nextInt(2*maxFreq+1);
		axis = rand.nextInt(3);
		split = rand.nextDouble();
		thickness = GAUtils.renormZeroOne(thicknessMean + thicknessSigma * rand.nextGaussian());
		
		if(verbosity >= 4)
			System.out.println("amp:"+amplitude+" etaFreq:"+etaFreq+" zetaFreq:"+zetaFreq
				+" axis:"+axis+" split:"+split+" thickness:"+thickness);
		
		//choose two random parents' Structures
		StructureOrg[] ps = (StructureOrg[])(sel.doSelection(parents, 2));
		Cell a = ps[0].getCell();
		Cell b = ps[1].getCell();
		
		// some status info
		if (verbosity >= 3)
			System.out.print("Mating "+ps[0].getID()+" (fitness " + ps[0].getFitness() + ", "
					+ a.getNumSites()+ " atoms) with " + ps[1].getID() +" (fitness "
					+ ps[1].getFitness() + ", " + b.getNumSites() + " atoms)... ");
		
		// possibly take a supercell of one of the parents
		if (growParents) {
			if (a.getBasisSize() > b.getBasisSize()) {
				int sizeMultiple = a.getBasisSize() / b.getBasisSize(); // integer division
				// TODO: make be input options and similarly below
				if (sizeMultiple > 1) {
					b = Cell.getSupercell(b, SupercellOptimizer.getOptimalSupercell(b, false, 6 , sizeMultiple, params.getMaxLatticeLength(), false));
					if (verbosity >= 3)
						System.out.println("(actually using supercell of second parent with " + b.getNumSites() + "atoms).");
				}
			} else {
				int sizeMultiple = b.getBasisSize() / a.getBasisSize();
				if (sizeMultiple > 1) {
					a = Cell.getSupercell(a, SupercellOptimizer.getOptimalSupercell(a, false, 6 , sizeMultiple, params.getMaxLatticeLength(), false));
					if (verbosity >= 3)
						System.out.println("(actually using supercell of first parent with " + a.getNumSites() + "atoms).");
				}
			}
		}
				
		//average parents' bases to create a new basis
		double[] aAngles = a.getCellAngles();
		double[] bAngles = b.getCellAngles();
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
		if (verbosity >= 3)
			System.out.println("("+ newOrganism.getCell().getNumSites() +" atoms)");
		
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
		
		String[] hsArgs = {"0.5", "0.0", "0", "0.00", "0", "0"};
		Variation p = new Slicer(hsArgs);
		String[] selArgs = {"2", "0"};
		Selection sel = new ProbDistSelection(selArgs);
		
		StructureOrg o = (StructureOrg)p.doVariation(parents, null, sel);
		System.out.println(o);
		//GAUtils.writeStringToFile(o.getCIF(), new File("offspring.cif"), false); 
		o.standardize();
		//GAUtils.writeStringToFile(o.getCIF(), new File("offspring_red.cif"), false); 
		System.out.println(o);

	}
}
