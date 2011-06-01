/* Genetic algorithm for crystal structure prediction.  Will Tipton.  Ceder Lab at MIT. Summer 2007 */

package ga;

// Classes which implement the Energy interface can compute the total
// lattice energy of a StructureOrg.  They're generally wrappers which
// call an external code.  They're generally used by classes which 
// implement ObjectiveFunction (e.g. EnergyPerAtom).

//NB: Convention: when runs fail, energy should evaluate to Double.POSITIVE_INFINITY


public interface Energy {
	public abstract double getEnergy(StructureOrg o);
}
