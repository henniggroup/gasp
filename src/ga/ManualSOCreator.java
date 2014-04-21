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

import java.util.List;

import crystallography.Cell;


// ManualSOCreator allows for structures to be added by code which uses
// the GeneticAlgorithm.  That code should call GAParameters.setSeedStructures
// before beginning the algorithm. 

public class ManualSOCreator implements StructureOrgCreator {
	
	Cell[] seeds;

	public ManualSOCreator(List<String> args) {
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
