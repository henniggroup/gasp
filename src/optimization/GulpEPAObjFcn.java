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
package optimization;

import crystallography.Cell;
import ga.ObjectiveFunction;
import gulp.*;
import utility.Pair;

public class GulpEPAObjFcn extends ObjectiveFunction {
	
	String potentialFile;
	OptiSystem sys;
	int timeLimit;
	boolean cautious;
	
	public GulpEPAObjFcn(String pot, int t, OptiSystem s, boolean _cautious) {
		potentialFile = pot;
		sys = s;
		timeLimit = t;
		cautious = _cautious;
	}

	public Pair<Cell,Double> getCellAndValue(Cell c) {		
		GulpEnergy ge = new GulpEnergy(c, potentialFile, timeLimit, true, sys, cautious);
		
		return new Pair<Cell,Double>(ge.getOptimizedCell(), ge.getOptimizedEnergy() / c.getBasisSize());
	}

	@Override
	public void run() {
		// TODO Auto-generated method stub
		
	}

	@Override
	public Thread evaluate() {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public String toString() {
		// TODO Auto-generated method stub
		return null;
	}

}
