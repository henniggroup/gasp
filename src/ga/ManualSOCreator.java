/* Genetic algorithm for crystal structure prediction.  Will Tipton.  Ceder Lab at MIT. Summer 2007 */

package ga;

import crystallography.Cell;


// ManualSOCreator allows for structures to be added by code which uses
// the GeneticAlgorithm.  That code should call GAParameters.setSeedStructures
// before beginning the algorithm. 

public class ManualSOCreator implements StructureOrgCreator {
	
	Cell[] seeds;

	public ManualSOCreator(String[] args) {
		// nothing
	}
	
	public String toString() {
		return "ManualSOCreator";
	}
	
	public StructureOrg makeOrganism(Generation g) {
		GAParameters params = GAParameters.getParams();
		
		// initialize seeds
		if (seeds == null)
			seeds = params.getSeedStructures();
		
		// set each seed to null as we use it
		for (int i = 0; i < seeds.length; i++)
			if (seeds[i] != null) {
				StructureOrg result = new StructureOrg(seeds[i]);
				seeds[i] = null;
				return result;
			}
		
		// if we get here, we've used all the seeds
		return null;
	}
}
