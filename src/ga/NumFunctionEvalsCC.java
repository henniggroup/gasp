/* Genetic algorithm for crystal structure prediction.  Will Tipton.  Ceder Lab at MIT. Summer 2007 */

package ga;

// NumFunctionEvalsCC implements ConvergenceCriterion.  It indicates that the algorithm
// has converged after a certain number of ObjectiveFunction evaluations have been done.
// Notice that keeping track of the number of evaluations is done by the ObjectiveFunction
// itself.

public class NumFunctionEvalsCC implements ConvergenceCriterion{
	
	int maxNumFunctionEvals;
	
	public NumFunctionEvalsCC(String[] args) {
		if (args == null || args.length < 1)
			GAParameters.usage("Not enough parameters given to NumFunctionEvalsCC", true);
		
		maxNumFunctionEvals = Integer.parseInt(args[0]);
	}
	
	public String toString() {
		return "NumFunctionEvalsCC. maxNumFunctionEvals = " + maxNumFunctionEvals;
	}
	
	public Boolean converged(Generation currentGen) {
		// get the number of function evaluations
		int numFunctionEvals = ObjectiveFunction.getNumCalculations();
		
		// check convergence
		return (numFunctionEvals >= maxNumFunctionEvals);
	}
}
