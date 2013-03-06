/* Genetic algorithm for crystal structure prediction.  Will Tipton.  Ceder Lab at MIT. Summer 2007 */

package ga;

import java.util.*;

import utility.Vect;

import crystallography.*;

import chemistry.*;


// RandomSOCreator implements StructureOrgCreator.  It creates StructureOrgs with
// random lattice parameters and atomic positions within the hard constraints and 
// with a given number of stoichiometries.  Then, it optionally scales structures 
// to a particular initial volume.

public class RandomSOCreator implements StructureOrgCreator {
	
	private String newGenType;
	private double givenVolPerAtom;
	
	public RandomSOCreator(List<String> args) {
		GAParameters params = GAParameters.getParams();
		
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

	// makes a basis which is random except that it satisfies the hard constraints
	public static List<Vect> makeRandomLattice() {
		GAParameters params = GAParameters.getParams();
		
		if (params.getMaxLatticeAngleDegrees() == -1 || params.getMinLatticeAngleDegrees() == -1
				|| params.getMaxLatticeLength() == -1 || params.getMinLatticeLength() == -1)
			GAParameters.usage("Error: Must set lattice parameter constraints to use makeRandomLattice().", true);
		
		double maxll = params.getMaxLatticeLength();
		double minll = params.getMinLatticeLength();
		double maxla = params.getMaxLatticeAngleDegrees();
		double minla = params.getMinLatticeAngleDegrees();		
		double maxch = params.getMaxCellHeight();
		
		Random rand = params.getRandom();
		
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
		double maxZll = Math.min(maxll,2*maxch);
		
		double lc = rand.nextDouble()*(maxZll-minll)+minll;
		
		return (new Cell(la,lb,lc,adeg,bdeg,gdeg,null,null)).getLatticeVectors();
	}
	
	// fills up the population with random organisms (provided they satisfy the
	// hard constraints)
	protected StructureOrg makeRandomOrg() {
		GAParameters params = GAParameters.getParams();
		Random rand = params.getRandom();
		//HashMap<String,Integer> constituents = params.getConstituents();
		Composition comp = params.getCompSpace().getRandomIntegerCompInSpace(params.getMinNumAtoms(), params.getMaxNumAtoms());

		// make random lattice parameters satisfying hard constraints
		// the Basis for our new organism
		List<Vect> latVects = makeRandomLattice();

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

		return new StructureOrg(newStructure);
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
