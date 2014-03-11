/*
 * Copyright 2011-2014 Will Tipton, Richard Hennig, Ben Revard, Stewart Wenner

This file is part of the Genetic Algorithm for Structure and Phase Prediction (GASP).

    GASP is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    GASP is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with GASP.  If not, see <http://www.gnu.org/licenses/>.
    
    
    */

package ga;

import java.util.ArrayList;
import java.util.List;
import java.util.Random;
import utility.Constants;
import crystallography.Cell;

public class SupercellVariation  implements Variation {
	
	static final long serialVersionUID = 1l;
	
	final int scaleFactor = 2;
	final int maxAttempts = 20; 
	
	private boolean relaxChildren;
	
	public SupercellVariation(List<String> args) {
		if (args == null || args.size() != 1)
			GAParameters.usage("Wrong parameters given to SupercellVariation", true);
		
		relaxChildren = Boolean.parseBoolean(args.get(0));
	}
	
	public String toString() {
		StringBuilder result = new StringBuilder();
		
		result.append("SupercellVariation");
		
		return result.toString();
	}
	
	public Organism doVariation(Generation parents, Generation offspring, Selection sel) {
		GAParameters params = GAParameters.getParams();
		Random rand = params.getRandom();
		
		// pick a parent randomly 
		// but make sure it's going to work; i.e. it has less than maxnumatoms/2 atoms
		StructureOrg p = null;
		Cell pStruct = null;
		for (int i = 0; i < maxAttempts; i++) { 
			p = (StructureOrg)(sel.doSelection(parents, 1)[0]);
			if (p.getCell().getNumSites() * scaleFactor >= GAParameters.getParams().getMaxNumAtoms())
				continue;
			
			pStruct = p.getCell();
			break;
		}
		if (pStruct == null)
			return null;
		
		// make identity matrix
		List<List<Integer>> coefs = new ArrayList<List<Integer>>();
		for (int i = 0; i < Constants.numDimensions; i++) {
			List<Integer> newVect = new ArrayList<Integer>();
			newVect.add(0); newVect.add(0); newVect.add(0);
			newVect.set(i, 1);
			coefs.add(newVect);
		}
		
		// double one of the dimensions
		int doubleDim = (int)(Math.random() * Constants.numDimensions);
		for (int i = 0; i < Constants.numDimensions; i++)
			coefs.get(doubleDim).set(i, coefs.get(doubleDim).get(i) * scaleFactor);
				
		// make the new offspring
		StructureOrg result = new StructureOrg(Cell.getSupercell(pStruct, coefs));
		
		// we dont need to recalculate energy or value if relaxChildren is false
		// note that im gonna assume energy is extensive and value is intensive here!!
		if (! relaxChildren) {
			result.setTotalEnergy(p.getTotalEnergy() * scaleFactor);
			result.setValue(p.getValue());
		}
		
		GAOut.out().stdout("SupercellVariation created new StructureOrg:", GAOut.DEBUG, result.getID());
		GAOut.out().stdout(result.toString(), GAOut.DEBUG, result.getID());			
		
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
