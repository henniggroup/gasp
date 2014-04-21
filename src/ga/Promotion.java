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

import java.io.Serializable;
import java.util.List;

import pdvisual.IComputedEntry;

// Promotion "promotes" some number of the best organisms in one generation directly
// to the next. 

public class Promotion implements Serializable {
	static final long serialVersionUID = 1;
	
	private int n;
	
	public Promotion(List<String> args) {
		if (args.size() < 1)
			GAParameters.usage("Not enough parameters given to Promotion", true);
		
		n = Integer.parseInt(args.get(0));
	}
	
	public String toString() {
		StringBuilder result = new StringBuilder();
		
		result.append("Promotion. n = " + n);
		
		return result.toString();
	}

	// adds the n highest fitness Organism from parents to offspring
	public void doPromotion(Generation parents, Generation offspring) {
		GAParameters params = GAParameters.getParams();
				
		// If we're doing a phase diagram run, then an n of 0 means to do no promotion and 
		// any nonzero n means to promote all organisms on the convex hull.
		// Otherwise, n is the number of organisms to promote.
		if (params.doingPDRun()) {
			if (n > 0)
				for (IComputedEntry e : params.getPDBuilder().getPDData().getStableEntries()) {
					StructureOrg o = (StructureOrg)e;
					GAOut.out().stdout("Promoting organism " + o.getID() + " with fitness "
							+ o.getFitness(), GAOut.NOTICE, o.getID());

					o.setValue(0);
					offspring.addOrganism(o);
				}
		} else {
			for (int i = 1; i <= n; i++) {
				Organism o = parents.getNthBestOrganism(i);
				// some status info
				GAOut.out().stdout("Promoting organism " + o.getID() + " with fitness "
						+ o.getFitness(), GAOut.NOTICE, o.getID());
				offspring.addOrganism(o);
			}
		}
	}

}
