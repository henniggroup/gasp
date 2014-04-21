/* Genetic algorithm for crystal structure prediction.  Will Tipton.  Ceder Lab at MIT. Summer 2007 */

package ga;

// Structures extends Generation and has list of StructureOrg. It is the type of
// Generation generally used by the genetic algorithm for the crystal structure
// prediction problem.

public final class Structures extends Generation {
	
	private RedundancyGuard rGuard;
	
	public Structures() {
		// initialize rGuard
		String rgType = GAParameters.getParams().getRedundancyGuardType();
		if (rgType.equalsIgnoreCase("both")||rgType.equalsIgnoreCase("perGeneration"))
			rGuard = new RedundancyGuard(GAParameters.getParams().getRedundancyGuardArgs());
	}
	
	public String toString() {
		StringBuilder result = new StringBuilder();
		
		result.append("Structures generation. rGuard: " + rGuard);
		
		return result.toString();
	}
	
	public void addOrganism(Organism s) {
		if (rGuard != null)
			rGuard.addStructureOrg(s);
		super.addOrganism(s);
	}
	
	public RedundancyGuard getRedundancyGuard() {
		return rGuard;
	}
}
