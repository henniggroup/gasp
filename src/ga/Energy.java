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

// Classes which implement the Energy interface can compute the total
// lattice energy of a StructureOrg.  They're generally wrappers which
// call an external code.  They're generally used by classes which 
// implement ObjectiveFunction (e.g. EnergyPerAtom).

//NB: Convention: when runs fail, energy should evaluate to Double.POSITIVE_INFINITY


public interface Energy {
	public abstract double getEnergy(StructureOrg o);
	public abstract boolean cannotCompute(StructureOrg o);
}
