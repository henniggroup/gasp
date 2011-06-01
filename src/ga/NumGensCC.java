/* Genetic algorithm for crystal structure prediction.  Will Tipton.  Ceder Lab at MIT. Summer 2007 */

package ga;

// NumGensCC implements ConvergenceCriterion.  It indicates that the algorithm has
// converged when the algorithm has run for a given number of generations.

public class NumGensCC implements ConvergenceCriterion {

	int maxNumGens;
	
	public NumGensCC(String[] args) {
		if (args == null || args.length < 1)
			GAParameters.usage("Not enough parameters given to NumGensCC", true);
		
		maxNumGens = Integer.parseInt(args[0]);
	}
	
	public String toString() {
		return "NumGensCC. maxNumGens = " + maxNumGens;
	}
	
	public Boolean converged(Generation currentGen) {
		return GAParameters.getParams().getRecord().getGenNum() > maxNumGens;
	}

}
