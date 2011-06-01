/* Genetic algorithm for crystal structure prediction.  Will Tipton.  Ceder Lab at MIT. Summer 2007 */

package ga;

// ValueAchievedCC is a ConvergenceCriterion which indicates convergence when
// there is at least one member of the current generation with a value less
// than or equal to a given target value.

public class ValueAchievedCC implements ConvergenceCriterion {

	double valueTarget;
	
	public ValueAchievedCC(String[] args) {
		if (args == null || args.length < 1)
			GAParameters.usage("Not enough parameters given to EnergyAchievedCC", true);
		
		valueTarget = Double.parseDouble(args[0]);
	}
	
	public String toString() {
		return "EnergyAchievedCC. energyTarget = " + valueTarget;
	}
	
	public Boolean converged(Generation currentGen) {		
		// get the best value info
		double best = currentGen.getExtremeValues()[0];
		
		// check convergence
		return (best <= valueTarget);
	}

}
