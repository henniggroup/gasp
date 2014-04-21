/* Genetic algorithm for crystal structure prediction.  Will Tipton.  Ceder Lab at MIT. Summer 2007 */

package ga;

// Interface to objects containing Selection algorithms.  A Selection
// algorithm takes a generation and an integer n and returns n organisms
// from the generation.  It assumes that fitnesses have already been calculated.

public interface Selection {
	
	public Organism[] doSelection(Generation g, int n);
	
	// should also overload toString();
	
}
