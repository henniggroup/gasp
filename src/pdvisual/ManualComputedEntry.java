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

import ga.GAUtils;

import java.io.Serializable;

import chemistry.Composition;
import crystallography.Cell;

public class ManualComputedEntry implements IComputedEntry,Serializable {
	static final long serialVersionUID = 1;
	
	double totalEnergy;
	Cell cell;
	
	public ManualComputedEntry(Cell c, double e) {
		cell = c;
		totalEnergy = e;
	}

	public double getEnergyPerAtom() {
		return totalEnergy / cell.getNumSites();
	}
	
	public Composition getComposition() {
		return cell.getComposition();
	}
	
	public Cell getCell() {
		return cell;
	}
	
	public double getTotalEnergy() {
		return totalEnergy;
	}
	
	public String getLabel() {
		return cell.getLabel();
	}
	
	public String toString() {
		StringBuilder result = new StringBuilder();
		String newline = GAUtils.newline();
		
		result.append("ManualComputedEntry: " + getLabel() + ":" + newline);
		result.append("Total Energy: " + totalEnergy +  newline);
		result.append(cell.toString());
		
		return result.toString();
	}
	
}
