/* Genetic algorithm for crystal structure prediction.  Will Tipton.  Ceder Lab at MIT. Summer 2007 */

package ga;

import java.io.File;
import java.util.*;

import utility.Pair;
import utility.Utility;

// This holds the main genetic algorithm.  It contains the evolutionary strategy
// in abstract form.  For crystal structure prediction, client code should call
// CrystalGA.crystalGA() 

public class GeneticAlgorithm {

	// Performs a generic genetic optimization procedure according to the inputed parameters
	// and returns the best Organism found.
	public static Organism doGeneticAlgorithm() {
		GAParameters params = GAParameters.getParams();
		
		// Create the objects for the algorithms
		Selection sel = params.getSelection();
		Development dev = params.getDevelopment();
		Promotion pro = params.getPromotion();
		
		// Create the initial population 
		Generation parents = null;  
		Generation offspring = params.getRecord().getCurrentGen();		
		
		// the main algorithm loop:
		while (!converged(offspring)) {
			// the wheel of life turns
			parents = offspring;
			offspring = GAParameters.getParams().makeEmptyGeneration();
			
			// do promotions (except for the first generation)
			if (params.getRecord().getGenNum() != 0)
				pro.doPromotion(parents, offspring);
			
			// make the offspring generation
			while (!madeEnough(offspring)) {
				List<Pair<Organism, Integer>> organisms = new ArrayList<Pair<Organism,Integer>>();	
				List<Thread> threads = new ArrayList<Thread>();
				
				for (int i = 0; i < params.getNumCalcsInParallel(); i++) {
					// make a new organism
					Organism newOrg = getNewOrg(parents, offspring, sel);
					if (newOrg == null)
						continue;
					organisms.add(new Pair<Organism,Integer>(newOrg, 1 + ObjectiveFunction.getNumCalculations()));
					
					// some status info
					GAOut.out().stdout(newOrg.toString(), GAOut.DEBUG, newOrg.getID());
					
					// start the energy computation
					if (!newOrg.knowsValue())
						threads.add(params.getObjectiveFunctionInstance(newOrg).evaluate());
				}
				
				// wait for the energy computations
				try {
					for (Thread t : threads)
						t.join();
				} catch (InterruptedException x) {
					GAOut.out().stdout("InterruptedException in energy calc thread: " + x.getMessage(), GAOut.WARNING);
				}
				// re-develop each of the organisms and possibly add it to the offspring generation
				for (Pair<Organism,Integer> o : organisms) {
					StructureOrg s = (StructureOrg)o.getFirst();
					if (dev != null && !dev.doDevelop(offspring, s)) {
						s.setSOCreator(null);
						continue;
					}
					GAOut.out().stdout("Adding organism " + s.getID() + " to generation " + params.getRecord().getGenNum() + ".", GAOut.NOTICE, s.getID());

					offspring.addOrganism(s);
					// if it was in the first gen, mark its creator as having successfully created
					StructureOrgCreator soc;
					if ((soc = s.getSOCreator()) != null) {
						decrementSOCUsageCount(soc);
						// get rid of the references to SoC's once theyre no longer needed so we dont serialize them
						s.setSOCreator(null);
					}		
					// write Hartke output if required
					if (GAParameters.getParams().getWriteHartkeFile()) {
						String hLine = s.getID() + " " + o.getSecond() + " " + s.getTotalEnergy() + " " + s.getValue() + "\n";
						Utility.writeStringToFile(hLine, GAParameters.getParams().getHartkeOutFile(), true);
					}
				}
			}
			
			// find the organisms' fitnesses
			GAOut.out().stdout("Starting fitness evaluations...", GAOut.WARNING);

			offspring.findFitnesses();
			
			// prints and saves status/progress info and whatnot
			GAParameters.getParams().getRecord().finishGen(offspring);
		}
		GAParameters.getParams().getRecord().cleanup();
		
		return offspring.getNthBestOrganism(1);
	}
	
	private static void decrementSOCUsageCount(StructureOrgCreator soc) {
		List<Pair<StructureOrgCreator,Integer>> ioc = GAParameters.getParams().getInitialOrgCreators();
		for (int i = 0; i < ioc.size(); i++) {
			if (soc == ioc.get(i).getFirst()) {
				ioc.set(i, new Pair<StructureOrgCreator,Integer>(ioc.get(i).getFirst(), ioc.get(i).getSecond() - 1));
			}
		}
	}
	
	private static Organism getNewOrg(Generation parents, Generation offspring, Selection sel) {
		// the 0th generation is special
		if (GAParameters.getParams().getRecord().getGenNum() == 0) {
			return makeFirstGenOrg(offspring);
		} else {	
			return makeOffspringOrg(parents, offspring, sel);
		}
	}
	
	/*private static Organism makeFirstGenOrg(Generation offspring) {
		GAParameters params = GAParameters.getParams();
		Development dev = params.getDevelopment();
		Map<StructureOrgCreator,Integer> initialOrgCreators = params.getInitialOrgCreators();
				
		StructureOrgCreator soc;
		Iterator<StructureOrgCreator> socs = initialOrgCreators.keySet().iterator();
		do {
			if (!socs.hasNext())
				return null;
			soc = socs.next();
		} while (initialOrgCreators.get(soc) <= 0);

		StructureOrg newOrg;
		// make new organisms with this soc until it returns null or until one develops correctly
		do { 
			newOrg = soc.makeOrganism(offspring);
		} while (newOrg != null && dev != null && !dev.doDevelop(offspring, newOrg));
		
		// save the organisms creator
		if (newOrg != null)
			newOrg.setSOCreator(soc);
		else
			initialOrgCreators.put(soc, 0);
		
		// some status info
		if (GAParameters.getParams().getVerbosity() >= 3)
			System.out.println("Organisms left for current StructureOrgCreator: " + initialOrgCreators.get(soc));
		
		return newOrg;
	}
	 */
	
	private static Organism makeFirstGenOrg(Generation offspring) {
		GAParameters params = GAParameters.getParams();
		Development dev = params.getDevelopment();
		List<Pair<StructureOrgCreator,Integer>> initialOrgCreators = params.getInitialOrgCreators();
				
		StructureOrg newOrg = null;

		while (newOrg == null) {
			// find the next non-empty SOC, soc
			int socIndx = 0;
			while (initialOrgCreators.size() > socIndx && initialOrgCreators.get(socIndx).getSecond() <= 0) 
				socIndx++;
			if (initialOrgCreators.size() <= socIndx)
				return null;
			StructureOrgCreator soc = initialOrgCreators.get(socIndx).getFirst();
	
			// make new organisms with this soc until it returns null or until one develops correctly
			do { 
				newOrg = soc.makeOrganism(offspring);
			} while (newOrg != null && dev != null && !dev.doDevelop(offspring, newOrg));
			
			// save the organisms creator
			if (newOrg != null)
				newOrg.setSOCreator(soc);
			else
				initialOrgCreators.set(socIndx, new Pair<StructureOrgCreator,Integer>(soc, 0));
			
			// some status info
			GAOut.out().stdout("Organisms left for current StructureOrgCreator: " + initialOrgCreators.get(socIndx).getSecond(), GAOut.NOTICE);
		}
				
		return newOrg;
	}
	
	private static Organism makeOffspringOrg(Generation parents, Generation offspring, Selection sel) {
		GAParameters params = GAParameters.getParams();
		Vector<Variation> vars = params.getVariations();
		Development dev = params.getDevelopment();
		Random rand = params.getRandom();
		
		// pick a variation, v, according to inputted probabilities
		int varNum;
		do {
			varNum = rand.nextInt(vars.size());
		} while (rand.nextDouble() > params.getVarProb(varNum));
		Variation v = vars.get(varNum);
		// apply the variation v:
		// - we have to make the prospective organism, develop, find fitness, and develop again
		// in case the evaluation makes an organism no longer satisfy hard constraints.
		// in particular, crystal relaxation will do this.  we develop once here.
		// - once we choose a variation, we make sure we consider an organism created by that
		// variation which passes the development (at least before the evaluation).  this way,
		// we won't short variations which tend to make offspring which don't fit the constraints.
		Organism newOrg;
		do { // some output
			GAOut.out().stdout("Starting variation " + v + " (done " + offspring.getNumOrganisms()
					+ " of " + params.getPopSize() + ").", GAOut.WARNING);
			newOrg = v.doVariation(parents, offspring, sel);
			if (newOrg == null) // if the chosen variation can't make a child, try another variation
				return makeOffspringOrg(parents, offspring, sel);
		} while (dev != null && !dev.doDevelop(offspring, newOrg));
		
		return newOrg;
	}
	
	private static Boolean madeEnough(Generation offspring) {
		GAParameters params = GAParameters.getParams();
		// the 0th generation is special
		if (GAParameters.getParams().getRecord().getGenNum() == 0) {
			List<Pair<StructureOrgCreator,Integer>> initialOrgCreators = params.getInitialOrgCreators();
			for (Pair<StructureOrgCreator, Integer> i : initialOrgCreators) {
				if (i.getSecond() > 0) 
					return false;
			}
			GAOut.out().stdout("Finished creating initial population.", GAOut.WARNING);
			return true;
		} else {	
			// only use the minimum population size option if we're doing parallel energy calculations
			if (params.getNumCalcsInParallel() != 1)
				return offspring.getNumOrganisms() >= Math.min(params.getMinPopSize(), params.getPopSize());
			else
				return offspring.getNumOrganisms() >= params.getPopSize();
		}
	}
	
	public static Boolean converged(Generation currentGen) {
		if(currentGen == null)
			return false;
		
		// check each of the convergence criteria
		Vector<ConvergenceCriterion> ccs = GAParameters.getParams().getCCs();
		Iterator<ConvergenceCriterion> i = ccs.iterator();
		while (i.hasNext()) {
			ConvergenceCriterion c = i.next();
			if (c.converged(currentGen))
				return true;
		}

		return false;
	}

}
