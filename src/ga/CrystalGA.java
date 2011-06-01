/* Genetic algorithm for crystal structure prediction.  Will Tipton.  Ceder Lab at MIT. Summer 2007 */

package ga;

import utility.ArgumentParser;
import utility.Utility;
import crystallography.Cell;


// Contains the top-level call for the crystal structure prediction program
// and a lot of the chemistry knowledge and logic.

// TODO: take a stoichiometry and return a Structure (long term goal.)

public class CrystalGA {

	// To seed the initial population, pass it in and have the input file use the
	// "manual" population option (e.g. population 20 structures manual).
	// When using a random initial population, the value of initialPop doesn't matter.
	public static Cell crystalGA(String inputFileName, Cell[] initialPop) {
		GAParameters params = GAParameters.getParams();
		String[] args = {"--f ", inputFileName};
		params.setArgs(args);

		params.setSeedGeneration(initialPop);
		
		StructureOrg s = (StructureOrg)GeneticAlgorithm.doGeneticAlgorithm();
		
		return s.getCell();
	}
	
	public static void main(String[] args) {
		// Parse command-line arguments

		ArgumentParser aParser = new ArgumentParser(args);
		if (aParser.hasArguments("r")) {
			System.out.println("Resuming from" +  aParser.getArgument("r"));
			GAParameters.setParams((GAParameters)(Utility.readSerializable(aParser.getArgument("r"))));
		} else {
			GAParameters.getParams().setArgs(args);
		}
		
		GeneticAlgorithm.doGeneticAlgorithm();

	}
}
