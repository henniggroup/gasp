/* Genetic algorithm for crystal structure prediction.  Will Tipton.  Ceder Lab at MIT. Summer 2007 */

package ga;

// NumGensWOImprCC implements ConvergenceCriterion.  It indicates that the algorithm has
// converged when the algorithm has run for a given number of generations without improving
// its best solution.

public class NumGensWOImprCC implements ConvergenceCriterion {
	
	int numGensWithoutProgress = 0;
	double bestValue = 0;
	
	int maxNumGensWOProgress;
	double improvementTolerance;
	
	public NumGensWOImprCC(String[] args) {
		if (args == null || args.length < 2)
			GAParameters.usage("Not enough parameters given to NumGensWOImprCC", true);
		
		maxNumGensWOProgress = Integer.parseInt(args[0]);
		improvementTolerance = Double.parseDouble(args[1]);
	}
	
	public String toString() {
		return "NumGensWOImprCC. n = " + numGensWithoutProgress;
	}
	
	public Boolean converged(Generation currentGen) {
		int currentGenNum = GAParameters.getParams().getRecord().getGenNum();
		
		// update the best value info
		double best = currentGen.getExtremeValues()[0];
		if (currentGenNum == 0 || best+improvementTolerance < bestValue) {
			bestValue = best;
			numGensWithoutProgress = 0;
		} else {
			numGensWithoutProgress++;
		}
		
		// check convergence
		return (numGensWithoutProgress >= maxNumGensWOProgress);
	}

}
