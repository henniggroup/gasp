/* Genetic algorithm for crystal structure prediction.  Will Tipton.  Ceder Lab at MIT. Summer 2007 */

package ga;

// ObjectiveFunction is an abstract class which specifies the interface to
// function we want to minimize.

public abstract class ObjectiveFunction implements Runnable {
	
	protected static int numCalculations = 0;
	
	// remember to update numCalculations when implementing this
	public abstract Thread evaluate();
	
	public static int getNumCalculations() {
		return numCalculations;
	}
	
	// an ObjectiveFunction should also overload toString();
	public abstract String toString();
	
}
