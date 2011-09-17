/* Genetic algorithm for crystal structure prediction.  Will Tipton.  Ceder Lab at MIT. Summer 2007 */

package ga;

import java.io.Serializable;
import java.util.*;

// Generation contains a Vector of Organisms which are treated
// logically as a generation.  It also implements a number of operations on that group
// such as calculation of fitnesses from values. A generation is not an
// evolving population: after its creation, it rarely gains or loses members.
// This is an abstract class which is implemented by e.g. Structures.

public abstract class Generation implements Iterable<Organism>, Serializable {
	static final long serialVersionUID = 1;

	// all the organisms
	protected Vector<Organism> organisms;
	
	protected Generation() {		
		organisms = new Vector<Organism>();
	}
	
	public Iterator<Organism> iterator() {
		return organisms.iterator();
	}
	
	public int getNumOrganisms() {
		return organisms.size();
	}
	
	public void addOrganism(Organism o) {
		organisms.add(o);
	}
	
	public void removeOrganism(Organism o) {
		organisms.remove(o);
	}
	
	public Boolean contains(Organism o) {
		return organisms.contains(o);
	}
	
	public Organism getOrganism(int i) {
		return organisms.get(i);
	}
	
	public List<Organism> getOrganismsSorted() {
		
		class StructureEnergyComparator implements Comparator<Organism> {
			public int compare(Organism a, Organism b) {
				return a.getFitness() > b.getFitness() ? -1 : 1;
			}
		}
		
		List<Organism> sortedOrgs = new ArrayList<Organism>(organisms);
		Collections.sort(organisms, new StructureEnergyComparator());
		
		return sortedOrgs;
	}
	
	public Organism getOrgByID(int id) {
		for (Organism o: organisms)
			if (o.getID() == id)
				return o;
		return null;
	}
	
	// returns the number of organisms with energy within dValue of value
	public Organism[] getOrganismsOfValue(double value, double dValue) {
		Iterator<Organism> i = organisms.iterator();
		ArrayList<Organism> answerList = new ArrayList<Organism>();
		
		while (i.hasNext()) {
			Organism o = i.next();
			double v = o.getValue();
			if (v <= value+dValue && v >= value-dValue)
				answerList.add(o);
		}
		
		Organism[] answer = new Organism[answerList.size()];
		answer = answerList.toArray(answer);
		
		return answer;
	}
	
	// returns a double[2] containing {bestValue, worstValue}
	public double[] getExtremeValues() {
		Iterator<Organism> i = organisms.iterator();
		
		// do the first organism by itself to initialize best and worst
		Organism o = i.next();
		double bestValue = o.getValue();
		double worstValue = o.getValue();
		// do the rest of the organisms
		while (i.hasNext()) {
			o = i.next();
			// update best and worst
			if (o.getValue() > worstValue)
				worstValue = o.getValue();
			if (o.getValue() < bestValue)
				bestValue = o.getValue();
		}
		double[] ans = {bestValue, worstValue};
		return ans;
	}
	
	// returns the organism with the Nth highest fitness.
	// make sure to return distinct organisms for distinct n even in the case
	// that multiple structures have the same energy
	// btw, the best org is at n==1, not n==0
	public Organism getNthBestOrganism(int n) {
		return getOrganismsSorted().get(n-1);
	}
	/*}
	public Organism getNthBestOrganism(int n) {
		ArrayList<Double> fitnesses = new ArrayList<Double>();
		Iterator<Organism> i = organisms.iterator();
		while(i.hasNext()) {
			Organism c = i.next();
			fitnesses.add(new Double(c.getFitness()));
		}
		// don't try to select from more than we have
		n = Math.min(n, fitnesses.size());
		
		// find the nth best fitness
		double fn = GAUtils.doubleSelect(fitnesses, n);
		
		// find the corresponding organism
		i = organisms.iterator();
		while(i.hasNext()) {
			Organism c = i.next();
			if (fn == c.getFitness())
				return c;
		}	
		
		return null;
	} */
	
	// assume values have already been calculated
	public void findFitnesses() {
		Organism o;
		int verbosity = GAParameters.getParams().getVerbosity();
		
		// make sure we have something in our population
		if (organisms.size() < 1) {
			System.out.println("ObjectiveFunction: no organisms in population??");
			System.exit(1);
		}
		
		//Calculate the organisms' fitnesses by normalization relative to the best (1) and worst (0)
		double[] extremes = getExtremeValues();
		double bestValue = extremes[0];
		double worstValue = extremes[1];
		Iterator<Organism> j = iterator();
		while (j.hasNext()) {
			o = j.next();
			// just in case:
			if (bestValue == worstValue)
				o.setFitness(1);
			else
				o.setFitness((o.getValue() - worstValue)/(bestValue - worstValue));
		}
		
		// some output
		if (verbosity >= 3)
			for (Organism i : organisms)
				System.out.println("Org " + i.getID() + "; value: " + i.getValue() + "; fitness: "+ i.getFitness());
	}
}
