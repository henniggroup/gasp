/* Genetic algorithm for crystal structure prediction.  Will Tipton.  Ceder Lab at MIT. Summer 2007 */
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

import java.io.*;
import java.util.List;

import utility.Utility;

import crystallography.Cell;

// FoundStructureCC is a ConvergenceCriterion which indicates convergence when
// there is a member of the current generation which is the same as a given structure.
// The target structure is given as a Cif file and the matching is done using a
// RedundancyGuard which calls a StructureFitter internally.

public class FoundStructureCC implements ConvergenceCriterion {
	
	RedundancyGuard target;
	
	public FoundStructureCC(List<String> args) {
		if (args == null || args.size() < 1)
			GAParameters.usage("Not enough parameters given to FoundStructureCC", true);
		
		String targetFilename = args.get(0);
		
		StructureOrg targetSOrg = new StructureOrg(Cell.parseCif(new File(targetFilename)));
		
		// initialize the target
		//String[] rgArgs = {"1.0", "0.001"};
		
		target = new RedundancyGuard(Utility.subList(args, 1));
		
		target.addStructureOrg(targetSOrg);
	}

	public Boolean converged(Generation currentGen) {
		for (Organism o : currentGen)
			if (target.checkStructureOrg((StructureOrg)o) != null)
				return true;
		return false;
	}

}
