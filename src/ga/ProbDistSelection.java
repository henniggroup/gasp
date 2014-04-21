<<<<<<< HEAD
/* Genetic algorithm for crystal structure prediction.  Will Tipton.  Ceder Lab at MIT. Summer 2007 */
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

import java.io.Serializable;
import java.util.*;

// Contains the algorithm for random selection of Organisms.
// Random selection is done from a probability distribution determined by
// numSurvivors and power.  Only the numSurvivors best Organisms have a
// non-zero selection probability.  The selection probabilities of these
// are related by a power law.
// To get standard Simple selection, set numSurvivors to some fraction of
// the population and power to 0.  To get standard Roulette selection, set
// numSurvivors to the whole population and power to 1.
// Fitnesses and energies should have already been calculated.

public final class ProbDistSelection implements Selection, Serializable {
	static final long serialVersionUID = 1;
	
	private int numSurvivors;
	private double power;

	public ProbDistSelection(List<String> args) {
		if (args.size() < 2)
			GAParameters.usage("Not enough parameters given to ProbDistSelection", true);
		numSurvivors = Integer.parseInt(args.get(0));
		power = Double.parseDouble(args.get(1));	
	}
	
	public String toString() {
		StringBuilder result = new StringBuilder();
		
		result.append("ProbDistSelection selection");
		result.append("(numParents = " + numSurvivors + ", power = " + power + ")");
		
		return result.toString();
	}
	
	/*
	// set the selection probability of the worst individuals to 0, 
	// the worst being all but the top numSurvivors.  renormalizes the
	// fitnesses of the remaining organisms
	private HashMap<Organism, Double> neuterWorst(Generation g, HashMap<Organism, Double> probMap) {
		HashMap<Organism, Double> renormFitnesses = new HashMap<Organism, Double>(); 
		double bestFitness = g.getNthBestOrganism(1).getFitness();
		double fn; 
		
		// deal with the weird situations where we don't want to do any renormalization
		if (numSurvivors >= g.getNumOrganisms()
				|| (fn = g.getNthBestOrganism(numSurvivors+1).getFitness()) == bestFitness) {
			Iterator<Organism> i = g.iterator();
			while (i.hasNext()) {
				Organism o = i.next();
				renormFitnesses.put(o, o.getFitness());
			}
		} else {
			Iterator<Organism> i = g.iterator();
			while (i.hasNext()) {
				Organism o = i.next();
				if (o.getFitness() <= fn)
					probMap.put(o, 0.0);
				else
					renormFitnesses.put(o, (o.getFitness()-fn)/(bestFitness-fn));
			}
		}
		
		return renormFitnesses;
	}
	 
	
	// we let the selection probability have a nth power dependence on the fitness
	// by setting an organism's selection probability to
	// (fitness^n)/(sum of all organisms' (fitnesses^n))
	private HashMap<Organism, Double> findProbabilities(Generation g) {	
		HashMap<Organism, Double> probMap = new HashMap<Organism, Double>();
		
		// set the selection probability of all but numSurvivors organisms to 0
		HashMap<Organism, Double> renormFitnesses = neuterWorst(g, probMap);

		Iterator<Organism> i;
		// find the sum, sumF, of the fitnesses of the remainder
		Double sumF = 0.0;
		i = g.iterator();
		while (i.hasNext()) {
			Organism o = i.next();
			if (!probMap.containsKey(o))
				sumF += Math.pow(renormFitnesses.get(o), power);
		}	
		
		// set the selection probabilities of the remainder to the fitness^n/sumF
		i = g.iterator();
		while (i.hasNext()) {
			Organism o = i.next();
			if (!probMap.containsKey(o))
				probMap.put(o, Math.pow(renormFitnesses.get(o), power)/sumF);
		}	
		
		return probMap;
	} */
	
	private HashMap<Organism, Double> findProbabilities(Generation g) {	
		HashMap<Organism, Double> result = new HashMap<Organism, Double>();
		
		// make a list of all the organisms
		int numOrgs = g.getNumOrganisms();
		List<Organism> orgs = new ArrayList<Organism>();
		for (int i = 0; i < numOrgs; i++)
			orgs.add(g.getOrganism(i));
		
		// get unnormalized non-zero probabilities for top numSurvivor organisms
		// and 0 probabilities for the rest
		for (int i = 0; i < numOrgs; i++) {
			// get best organism in list by fitness
			Organism bestOrg = orgs.get(0);
			for (int j = 1; j < orgs.size(); j++) 
				if (orgs.get(j).getFitness() > bestOrg.getFitness()) 
					bestOrg = orgs.get(j);
			// remove it from the list and add it to the result
			orgs.remove(bestOrg);
			double unnormalizedProb = Math.max(numSurvivors - i, 0);
			result.put(bestOrg, unnormalizedProb);
		}
		
		// normalize stuffs 
		double sum = 0.0;
		for (Organism o : result.keySet())
			sum += Math.pow(result.get(o), power);
		
		for (Organism o : result.keySet())
			result.put(o, Math.pow(result.get(o), power) / sum);
		
		return result;
	}
	
	public Organism[] doSelection(Generation g, int n) {
		Random rand = GAParameters.getParams().getRandom();
		
		ArrayList<StructureOrg> resultList = new ArrayList<StructureOrg>();
		
		// some output
		GAOut.out().stdout("Selecting " + n + " of the top " + numSurvivors + " of " + g.getNumOrganisms() + " organisms", GAOut.NOTICE);

		
		// get the mapping from organisms to selection probabilities
		Map<Organism, Double> probMap = findProbabilities(g);
		GAOut.out().stdout(getProbMapString(probMap), GAOut.DEBUG);
		
		// select n random, distinct organisms according to the calculated selection probabilities
		for (int i = 0; i < n; i++) {
			StructureOrg o;
			do
				o = (StructureOrg)g.getOrganism(rand.nextInt(g.getNumOrganisms()));
			while (rand.nextDouble() > probMap.get(o) || resultList.contains(o));
			resultList.add(o);
		
			// some status info
			GAOut.out().stdout("Selected organism " + o.getID() + " (fitness " 
					+ o.getFitness() + ", probability " + probMap.get(o) + ")", GAOut.DEBUG, o.getID());
		}
		
		StructureOrg[] result = new StructureOrg[resultList.size()];
		result = resultList.toArray(result);
		return result;
	}
	
	private String getProbMapString(Map<Organism, Double> probMap) {
		StringBuilder result = new StringBuilder();
		
		Set<Organism> keySet = probMap.keySet();
		Iterator<Organism> i = keySet.iterator();
		while (i.hasNext()) {
			Organism o = i.next();
			result.append("Organism " + o.getID() + " has fitness " + o.getFitness()
					+" and selection probability " + probMap.get(o) + "\n");
		}
		
		return result.toString();
	}
	
	public static void main(String[] args) {
		List<String> elitistArgs = new ArrayList<String>();
		elitistArgs.add("5");
		elitistArgs.add("1");
		ProbDistSelection bob = new ProbDistSelection(elitistArgs);
		
		Structures g = new Structures();
		for (double x = -0.2; x <= 1.0; x += 0.05) {
			StructureOrg s = new StructureOrg(null);
			s.setFitness(0.7);
			g.addOrganism(s);
		}
		StructureOrg s = new StructureOrg(null);
		s.setFitness(1.0);
		g.addOrganism(s);
		
		for (int i = 0; i < 10; i++)
			System.out.println(bob.doSelection(g, 1)[0].getID());
	}
}
