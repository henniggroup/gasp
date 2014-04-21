/* Genetic algorithm for crystal structure prediction.  Will Tipton.  Ceder Lab at MIT. Summer 2007 */

package ga;

import java.io.Serializable;

// A Variation operation takes a pair of generations and uses them to create an organism to
// be added to the offspring generation.  Neither should be modified.  Application of Variations,
// along with Promotion, are how a new Generation is created from a parent Generation.

// When implementing, keep in mind that Structure, Organism, etc. mostly don't have copy constructors.

public interface Variation extends Serializable {
	
	public Organism doVariation(Generation parents, Generation offspring, Selection sel);
	
}
