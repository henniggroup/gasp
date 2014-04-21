<<<<<<< HEAD
/* Genetic algorithm for crystal structure prediction.  Will Tipton.  Ceder Lab at MIT. Summer 2007 */
=======
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
>>>>>>> 0e3189c40547bbd59ea42c4f91890d7511fb7797

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
