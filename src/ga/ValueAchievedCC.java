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

// ValueAchievedCC is a ConvergenceCriterion which indicates convergence when
// there is at least one member of the current generation with a value less
// than or equal to a given target value.

public class ValueAchievedCC implements ConvergenceCriterion {

	double valueTarget;
	
	public ValueAchievedCC(List<String> args) {
		if (args == null || args.size() < 1)
			GAParameters.usage("Not enough parameters given to EnergyAchievedCC", true);
		
		valueTarget = Double.parseDouble(args.get(0));
	}
	
	public String toString() {
		return "EnergyAchievedCC. energyTarget = " + valueTarget;
	}
	
	public Boolean converged(Generation currentGen) {		
		// get the best value info
		double best = currentGen.getExtremeValues()[0];
		
		// check convergence
		return (best <= valueTarget);
	}

}
