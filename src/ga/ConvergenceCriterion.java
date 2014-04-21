/* Genetic algorithm for crystal structure prediction.  Will Tipton.  Ceder Lab at MIT. Summer 2007 */

package ga;

import java.io.Serializable;

// Each ConvergenceCriterion in use should be called exactly once per generation.
// It will return true if the population is determined to have "converged."
// By convention, all classes implementing ConvergenceCriterion end in "CC".
public interface ConvergenceCriterion extends Serializable {

	public Boolean converged(Generation currentGen);
	
}
