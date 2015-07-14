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

package vasp;

import java.io.Serializable;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import pdvisual.IComputedEntry;

import chemistry.Composition;
import chemistry.Element;
import utility.Vect;
import crystallography.Cell;

public class VaspConfig implements IComputedEntry, Serializable {

	private static final long serialVersionUID = 1L;
	
	private Cell cell;
	private double totalEnergy;
	private List<Vect> forces;
	private double stress[];
	
	static final int STRESS_LEN = 6;
	
	public VaspConfig() {
		cell = null;
		totalEnergy = Double.NaN;
		forces = new ArrayList<Vect>();
		stress = new double[STRESS_LEN];
	}
	
	public void setTotalEnergy(double e) {
		totalEnergy = e;
	}
	public double getTotalEnergy() {
		return totalEnergy;
	}
	public boolean energyHasBeenSet() {
		return ! Double.isNaN(totalEnergy);
	}
	public void setStress(double s[]) {
		if (s.length != STRESS_LEN)
			throw new IllegalArgumentException("bad stress vector passed to VaspConfig.setStress()");
		for (int i = 0; i < STRESS_LEN; i++)
			stress[i] = s[i];
	}
	public double[] getStress() {
		double result[] = new double[STRESS_LEN];
		for (int i = 0; i < STRESS_LEN; i++)
			result[i] = stress[i];
		return result;
	}
	public void setForces(List<Vect> f) {
		forces = new ArrayList<Vect>(f);
	}
	public List<Vect> getForces() {
		return new ArrayList<Vect>(forces);
	}
	public void setCell(Cell c) {
		cell = c;
	}
	public Cell getCell() {
		return cell;
	}
	
	static Map<Element,Integer> typeMap = new HashMap<Element,Integer>();

	static int elemToType(Element e) {
		if (!typeMap.containsKey(e))
			typeMap.put(e, typeMap.keySet().size());
		return typeMap.get(e);
	}
	
	static void parseTypeMap(List<String> typeMapArgs) {
		if (typeMapArgs.size() % 2 != 0) {
			throw new IllegalArgumentException("ERROR: uneven number of args passed for type mapping");
		}
		
		for (int i = 0; i < typeMapArgs.size(); i = i + 2) {
			typeMap.put(Element.getElemFromSymbol(typeMapArgs.get(i)), Integer.parseInt(typeMapArgs.get(i+1)));
		}
	}
	
	public String toPotfitString(String comment) {
		StringBuilder result = new StringBuilder();
		
		String nl = System.getProperty("line.separator");
		// print header and always use forces
		result.append("#N " + getCell().getNumSites() + " 1" + nl);
		
		// print a couple comment lines
		result.append(comment+nl);
		
		// lattice vectors
		double [] lVects =  getCell().getLatticeVectorsArray();
		result.append("#X " + lVects[0] + " " + lVects[1] + " " + lVects[2] + nl);
		result.append("#Y " + lVects[3] + " " + lVects[4] + " " + lVects[5] + nl);
		result.append("#Z " + lVects[6] + " " + lVects[7] + " " + lVects[8] + nl);
		
		// energy is an energy per atom
		result.append("#E " + (getTotalEnergy() / getCell().getNumSites()) + nl);
		
		// stress tensor
		result.append("#S");
		double stress[] = getStress();
		for (int j = 0; j < VaspConfig.STRESS_LEN; j++)
			result.append(" " + stress[j]);
		result.append(nl);
		
		// weight
		result.append("#W 1.0" + nl);
		
		// ends the header
		result.append("#F" + nl);
		
		for (int j = 0; j < getCell().getNumSites(); j++) {
			result.append(elemToType(getCell().getSite(j).getElement()));
			for (double d : getCell().getSite(j).getCoords().getCartesianComponents())
				result.append(" " + d);
			for (double d : getForces().get(j).getCartesianComponents())
				result.append(" " + d);
			result.append(nl);

		}
		return result.toString();
	}

	public double getEnergyPerAtom() {
		return getTotalEnergy() / cell.getNumSites();
	}

	public Composition getComposition() {
		return cell.getComposition();
	}

	public String getLabel() {
		return cell.getLabel();
	}
}