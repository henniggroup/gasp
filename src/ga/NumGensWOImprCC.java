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

// NumGensWOImprCC implements ConvergenceCriterion.  It indicates that the algorithm has
// converged when the algorithm has run for a given number of generations without improving
// its best solution.

public class NumGensWOImprCC implements ConvergenceCriterion {
	
	int numGensWithoutProgress = 0;
	double bestValue = 0;
	
	int maxNumGensWOProgress;
	double improvementTolerance;
	
	public NumGensWOImprCC(List<String> args) {
		if (args == null || args.size() < 2)
			GAParameters.usage("Not enough parameters given to NumGensWOImprCC", true);
		
		maxNumGensWOProgress = Integer.parseInt(args.get(0));
		improvementTolerance = Double.parseDouble(args.get(1));
	}
	
	public String toString() {
		return "NumGensWOImprCC. n = " + numGensWithoutProgress;
	}
	
	public Boolean converged(Generation currentGen) {
		int currentGenNum = GAParameters.getParams().getRecord().getGenNum();
		
		// update the best value info
		double best = currentGen.getExtremeValues()[0];
		if (currentGenNum == 0 || best+improvementTolerance < bestValue) {
			bestValue = best;
			numGensWithoutProgress = 0;
		} else {
			numGensWithoutProgress++;
		}
		
		// check convergence
		return (numGensWithoutProgress >= maxNumGensWOProgress);
	}

}
