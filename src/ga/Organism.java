/* Genetic algorithm for crystal structure prediction.  Will Tipton.  Ceder Lab at MIT. Summer 2007 */

package ga;

import java.io.Serializable;

// Organism is the abstract class which represents members of a Generation
// in the genetic algorithm.  It is implemented by e.g. StructureOrg.

public abstract class Organism implements Serializable {
	static final long serialversionUID = 1l;
	
	// fitness is normalized to fall between 0 and 1, where 1 is the best.
	// we initialize to null to know we havent' computed the fitness yet. similarly
	// for the energy
	protected Double fitness = null;
	
	// value is generally non-normalized fitness, e.g. energy per atom
	protected Double value = null;
	
	private int id;
	
	public Organism() {
		id = GAParameters.getParams().getNewOrgID();
	}
	
	// organisms have a unique ID
	public int getID() {
		return id;
	}
	
	public void setID(int id_) {
		id = id_;
	}
	
	public double getFitness() {
		if (!knowsFitness())
			System.out.println("Warning: Org " + id + " using fitness w/o calculating it.");
		return fitness;
	}
	
	public double getValue() {
		if (!knowsValue())
			System.out.println("Warning: Org " + id + " using value w/o calculating it.");
		return value;
	}
	
	// since the fitness is sometimes expensive to calculate, it's useful to be
	// able to check whether we've already calculated it:
	public Boolean knowsValue() {
		return value != null;
	}
	
	public Boolean knowsFitness() {
		return fitness != null;
	}
	
	public void setFitness(double f) {
		fitness = f;
	}
	
	public void setValue(double e) {
		value = e;
	}
}
