/* Genetic algorithm for crystal structure prediction.  Will Tipton.  Ceder Lab at MIT. Summer 2007 */

package ga;

import java.io.Serializable;
import java.util.List;

import pdvisual.IComputedEntry;

// Promotion "promotes" some number of the best organisms in one generation directly
// to the next. 

public class Promotion implements Serializable {
	static final long serialVersionUID = 1;
	
	private int n;
	
	public Promotion(List<String> args) {
		if (args.size() < 1)
			GAParameters.usage("Not enough parameters given to Promotion", true);
		
		n = Integer.parseInt(args.get(0));
	}
	
	public String toString() {
		StringBuilder result = new StringBuilder();
		
		result.append("Promotion. n = " + n);
		
		return result.toString();
	}

	// adds the n highest fitness Organism from parents to offspring
	public void doPromotion(Generation parents, Generation offspring) {
		GAParameters params = GAParameters.getParams();
		
		int verbosity = params.getVerbosity();
		
		// If we're doing a phase diagram run, then an n of 0 means to do no promotion and 
		// any nonzero n means to promote all organisms on the convex hull.
		// Otherwise, n is the number of organisms to promote.
		if (params.doingPDRun()) {
			if (n > 0)
				for (IComputedEntry e : params.getPDBuilder().getPDData().getStableEntries()) {
					StructureOrg o = (StructureOrg)e;
					if (verbosity >= 3)
						System.out.println("Promoting organism " + o.getID() + " with fitness "
								+ o.getFitness());
					o.setValue(0);
					offspring.addOrganism(o);
				}
		} else {
			for (int i = 1; i <= n; i++) {
				Organism o = parents.getNthBestOrganism(i);
				// some status info
				if (verbosity >= 3)
					System.out.println("Promoting organism " + o.getID() + " with fitness "
							+ o.getFitness());
				offspring.addOrganism(o);
			}
		}
	}

}
