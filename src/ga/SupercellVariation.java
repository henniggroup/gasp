package ga;

import java.util.ArrayList;
import java.util.List;
import java.util.Random;
import utility.Constants;
import crystallography.Cell;

public class SupercellVariation  implements Variation {
	
	static final long serialVersionUID = 1l;

	public SupercellVariation(String[] args) {
		if (args == null || args.length != 0)
			GAParameters.usage("Wrong parameters given to SupercellVariation", true);

	}
	
	public String toString() {
		StringBuilder result = new StringBuilder();
		
		result.append("SupercellVariation variation");
		
		return result.toString();
	}
	
	public Organism doVariation(Generation parents, Generation offspring, Selection sel) {
		GAParameters params = GAParameters.getParams();
		Random rand = params.getRandom();
		int verbosity = params.getVerbosity();
		
		// pick a parent randomly 
		StructureOrg p = (StructureOrg)(sel.doSelection(parents, 1)[0]);
		Cell pStruct = p.getCell();
		
		// make identity matrix
		List<List<Integer>> coefs = new ArrayList<List<Integer>>();
		for (int i = 0; i < Constants.numDimensions; i++) {
			List<Integer> newVect = new ArrayList<Integer>();
			newVect.add(0); newVect.add(0); newVect.add(0);
			newVect.set(i, 1);
			coefs.add(newVect);
		}
		
		// double one of the dimensions
		final int scaleFactor = 2;
		int doubleDim = (int)(Math.random() * Constants.numDimensions);
		for (int i = 0; i < Constants.numDimensions; i++)
			coefs.get(doubleDim).set(i, coefs.get(doubleDim).get(i) * scaleFactor);
				
		// make the new offspring
		StructureOrg result = new StructureOrg(Cell.getSupercell(pStruct, coefs));
		
		// we dont need to recalculate energy or value
		// note that im gonna assume energy is extensive and value is intensive here!!
		result.setTotalEnergy(p.getTotalEnergy() * scaleFactor);
		result.setValue(p.getValue());
		
		if (verbosity >= 5) {
			System.out.println("SupercellVariation created new StructureOrg:");
			System.out.println(result);
		}			
		
		return result;
	}
	
	// just for testing
	/*
	public static void main(String[] args) {
		StructureOrg s1 = new StructureOrg(Cell.parseCif(new File("/home/wtipton/cifs/143.cif")));
		s1.setFitness(-1);
		
		String[] smArgs = {"0.5", "0.5", "0.5"};
		StructureMut s = new StructureMut(smArgs);
		String[] selArgs = {"1", "0"};
		Selection sel = new ProbDistSelection(selArgs);

		Generation parents = new Structures();
		parents.addOrganism(s1);
		
		System.out.println(s1);
		StructureOrg o = (StructureOrg)s.doVariation(parents, null, sel);
	//	GAUtils.writeStringToFile(o.getCIF(), new File("offspring.cif"), false);
		System.out.println(o);
	} */
}