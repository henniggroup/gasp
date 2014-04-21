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

// NumFunctionEvalsCC implements ConvergenceCriterion.  It indicates that the algorithm
// has converged after a certain number of ObjectiveFunction evaluations have been done.
// Notice that keeping track of the number of evaluations is done by the ObjectiveFunction
// itself.

public class NumFunctionEvalsCC implements ConvergenceCriterion{
	
	int maxNumFunctionEvals;
	
	public NumFunctionEvalsCC(List<String> args) {
		if (args == null || args.size() < 1)
			GAParameters.usage("Not enough parameters given to NumFunctionEvalsCC", true);
		
		maxNumFunctionEvals = Integer.parseInt(args.get(0));
	}
	
	public String toString() {
		return "NumFunctionEvalsCC. maxNumFunctionEvals = " + maxNumFunctionEvals;
	}
	
	public Boolean converged(Generation currentGen) {
		// get the number of function evaluations
		int numFunctionEvals = ObjectiveFunction.getNumCalculations();
		
		// check convergence
		return (numFunctionEvals >= maxNumFunctionEvals);
	}
}
