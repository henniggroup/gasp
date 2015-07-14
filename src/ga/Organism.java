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

import java.io.Serializable;

// Organism is the abstract class which represents members of a Generation
// in the genetic algorithm.  It is implemented by e.g. StructureOrg.

public abstract class Organism implements Serializable {
	static final long serialversionUID = 1l;
	
	// fitness is normalized to fall between 0 and 1, where 1 is the best.
	// we initialize to null to know we haven't computed the fitness yet. similarly
	// for the energy
	protected Double fitness = null;
	
	// value is generally non-normalized fitness, e.g. energy per atom
	protected Double value = null;
	
	private int id;
	
	public Organism() {
		id = GAParameters.getParams().getNewOrgID();
	}
	
	// organisms have a unique ID
	public int getID() {
		return id;
	}
	
	public void setID(int id_) {
		id = id_;
	}
	
	public double getFitness() {
		if (!knowsFitness())
			System.out.println("Warning: Org " + id + " using fitness w/o calculating it.");
		return fitness;
	}
	
	public double getValue() {
		if (!knowsValue())
			System.out.println("Warning: Org " + id + " using value w/o calculating it.");
		return value;
	}
	
	// since the fitness is sometimes expensive to calculate, it's useful to be
	// able to check whether we've already calculated it:
	public Boolean knowsValue() {
		return value != null;
	}
	
	public Boolean knowsFitness() {
		return fitness != null;
	}
	
	public void setFitness(double f) {
		fitness = f;
	}
	
	public void setValue(double e) {
		value = e;
	}
}
