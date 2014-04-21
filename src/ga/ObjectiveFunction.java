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

// ObjectiveFunction is an abstract class which specifies the interface to
// function we want to minimize.

public abstract class ObjectiveFunction implements Runnable {
	
	protected static int numCalculations = 0;
	
	// remember to update numCalculations when implementing this
	public abstract Thread evaluate();
	
	public static int getNumCalculations() {
		return numCalculations;
	}
	
	// an ObjectiveFunction should also overload toString();
	public abstract String toString();
	
}
