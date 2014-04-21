<<<<<<< HEAD
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
package pdvisual;

import java.util.Iterator;
import java.util.List;

public class PDConvexable implements Convexable {
	
	int dim;
	double[][] allVertices;
	IPDPlottable pdd;
	
	public PDConvexable(IPDPlottable _pdd) {
		// make the allVertices list
		// we use the elements order in this!
		
		pdd = _pdd;
		
		List<Double[]> coords = pdd.getCompositionFractions();
		//List<Double> energies = pdd.getFormEnergiesPerAtom();
		List<Double> energies = pdd.getEnergiesPerAtom(); 

             
		/**
		 * Thu Jul 24 14:36:26 EDT 2008 Chris BUGFIX
		 * chem potential information was not used ! now updated in PDData
		 */
		allVertices = new double[coords.size()][pdd.getDimension()];
		Iterator<Double[]> iterv = coords.iterator();
		for(int i = 0; i < allVertices.length ; i++) {
			Double[] v = iterv.next();
			for (int j = 0; j < v.length; j++)
				allVertices[i][j] = v[j].doubleValue();
		}
		for(int i = 0; i < energies.size(); i++)
			allVertices[i][pdd.getDimension() - 1] = energies.get(i);
	}
	
	// getAllverticesforConvexHull should return an array of double
	// each row is a vertex and each column the coord of the vertex.
	public double[][] getAllVerticesforConvexHull() {
		return allVertices;
	}
	
	// returns the dimension in which you're building your convex hull
	public int getdimension() {
		return pdd.getDimension();
	}
}