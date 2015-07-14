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

import java.util.List;

// NumGensCC implements ConvergenceCriterion.  It indicates that the algorithm has
// converged when the algorithm has run for a given number of generations.

public class NumGensCC implements ConvergenceCriterion {

	int maxNumGens;
	
	public NumGensCC(List<String> args) {
		if (args == null || args.size() < 1)
			GAParameters.usage("Not enough parameters given to NumGensCC", true);
		
		maxNumGens = Integer.parseInt(args.get(0));
	}
	
	public String toString() {
		return "NumGensCC. maxNumGens = " + maxNumGens;
	}
	
	public Boolean converged(Generation currentGen) {
		return GAParameters.getParams().getRecord().getGenNum() > maxNumGens;
	}

}
