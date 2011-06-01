/* Genetic algorithm for crystal structure prediction.  Will Tipton.  Ceder Lab at MIT. Summer 2007 */

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

	public ProbDistSelection(String[] args) {
		if (args.length < 2)
			GAParameters.usage("Not enough parameters given to Tournament", true);
		numSurvivors = Integer.parseInt(args[0]);
		power = Double.parseDouble(args[1]);	
	}
	
	public String toString() {
		StringBuilder result = new StringBuilder();
		
		result.append("Elitist selection");
		result.append("(numParents = " + numSurvivors + ", power = " + power + ")");
		
		return result.toString();
	}
	
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
	}
	
	public Organism[] doSelection(Generation g, int n) {
		int verbosity = GAParameters.getParams().getVerbosity();
		Random rand = GAParameters.getParams().getRandom();
		
		ArrayList<StructureOrg> resultList = new ArrayList<StructureOrg>();
		
		// some output
		if (verbosity >= 3)
			System.out.println("Selecting " + n + " of the top " + numSurvivors + " of " + g.getNumOrganisms() + " organisms");
		
		// get the mapping from organisms to selection probabilities
		Map<Organism, Double> probMap = findProbabilities(g);
		if (verbosity >= 5)
			printProbMap(probMap);
		
		// select n random, distinct organisms according to the calculated selection probabilities
		for (int i = 0; i < n; i++) {
			StructureOrg o;
			do
				o = (StructureOrg)g.getOrganism(rand.nextInt(g.getNumOrganisms()));
			while (rand.nextDouble() > probMap.get(o) || resultList.contains(o));
			resultList.add(o);
		
			// some status info
			if (verbosity >= 4)
				System.out.println("Selected organism " + o.getID() + " (fitness " 
						+ o.getFitness() + ", probability " + probMap.get(o) + ")");
		}
		
		StructureOrg[] result = new StructureOrg[resultList.size()];
		result = resultList.toArray(result);
		return result;
	}
	
	private void printProbMap(Map<Organism, Double> probMap) {
		Set<Organism> keySet = probMap.keySet();
		Iterator<Organism> i = keySet.iterator();
		while (i.hasNext()) {
			Organism o = i.next();
			System.out.println("Organism " + o.getID() + " has fitness " + o.getFitness()
					+" and selection probability " + probMap.get(o));
		}
	}
	
	public static void main(String[] args) {
		String[] elitistArgs = {"3", "1"};
		ProbDistSelection bob = new ProbDistSelection(elitistArgs);
		
		Structures g = new Structures();
		for (double x = 0.0; x <= 1.0; x += 0.2) {
			StructureOrg s = new StructureOrg(null);
			s.setFitness(x);
			g.addOrganism(s);
		}
		
		for (int i = 0; i < 1000; i++)
			System.out.println(bob.doSelection(g, 1)[0].getID());
	}
}
