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

import java.io.*;
import java.util.*;

import crystallography.Cell;

// RedundancyGuard is used by the algorithm to avoid considering identical
// StructureOrgs more than once.  It stores a Map of all structures "seen"
// to their Organism IDs.  Then, other algorithms (e.g. StructureDev) may
// check Organisms against the list of those already seen. 

public class RedundancyGuard implements Serializable {
	static final long serialVersionUID = 1;

	// holds all the structures the algorithm has "seen"
	// and maps them to their Organism IDs
	Map<Cell,Integer> structures;
	
	private double atomicMisfit;
	private double latticeMisfit;
	private double angleMisfit;
	private boolean usePBCs;
	
	public RedundancyGuard(List<String> args) {
		// initialize structures
		structures = new HashMap<Cell,Integer>();
		
		// parse args
		if (args == null || args.size() < 3)
			GAParameters.usage("Not enough parameters given to RedundancyGuard", true);

		atomicMisfit = Double.parseDouble(args.get(0));
		latticeMisfit = Double.parseDouble(args.get(1));
		angleMisfit = Double.parseDouble(args.get(2));
		
		if (args.size() > 3)
			usePBCs = Boolean.parseBoolean(args.get(3));
		else
			usePBCs = true;
	}
	
	public String toString() {
		return "RedundancyGuard";//: lAngleTol="+lAngleTol+", lLengthTol="+lLengthTol+", interval="+interval;
	}
	
	public void addStructureOrg(Organism o) {
		// the Organism better be a StructureOrg
		StructureOrg s = (StructureOrg)o;
		// ignore structures w/ no atoms. should maybe stop this from happening elsewhere?
		if (s.getCell().getBasisSize() < 1) {
			GAOut.out().stdout("Warning: RedundancyGuard got passed structure with no sites. Ignoring it...", GAOut.NOTICE, o.getID());
		} else {
			structures.put(s.getCell(), new Integer(s.getID()));
		}
	}
	
	/*
	public void removeStructureOrg(Organism o) {
		// the Organism better be a StructureOrg
		StructureOrg s = (StructureOrg)o;
		
		structures.remove(s);
	}*/
	
	// returns an ID of a matching StructureOrg if we've seen it before,
	// null otherwise
	public Integer checkStructureOrg(StructureOrg s1) {
		Cell s = s1.getCell();

		for (Cell t : structures.keySet()) {
			if (usePBCs) {
				if (s.matchesCell(t, atomicMisfit, latticeMisfit, angleMisfit))
					return structures.get(t);
			} else {
				if (s.matchesCellNoPBCs(t, atomicMisfit))
					return structures.get(t);
			}
		} 

		return null;
	}
	
	/*
	private Boolean sOrgEquals(Structure s, Structure t) {
		// check the lattice angles
		double[] sla = s.getCellAngles();
		double[] tla = t.getCellAngles();
		for (int i = 0; i < sla.length; i++)
			if (Math.abs(sla[i] - tla[i]) > lAngleTol)
				return false;
		
		// check the lattice lengths
		double[] sll = s.getCellLengths();
		double[] tll = t.getCellLengths();
		for (int i = 0; i < sll.length; i++)
			if (Math.abs(sll[i] - tll[i]) > lLengthTol)
				return false;
		
		// check the sites and species
		Site[] ss = s.getAllSites();
		Site[] ts = t.getAllSites();
		if (ss.length != ts.length)
			return false;
		for (int i = 0; i < ss.length; i++) {
			// check the species
			if (!ss[i].m_sp.m_symbol.equals(ts[i].m_sp.m_symbol))
				return false;
			// check atomic locations
			double[] slocs = ss[i].getCoords().getCoords();
			double[] tlocs = ts[i].getCoords().getCoords();
			for (int j = 0; j < slocs.length; j++)
				if (Math.abs(slocs[j] - tlocs[j]) > interval)
					return false;
		}
		
		return true;
	}
	*/
	
	// just for testing
	public static void main(String[] args) {
		/*
		String[] rdArgs = {"0.5", "0.5", "0.05", "false"};
		
		RedundancyGuard bob = new RedundancyGuard(rdArgs);
		
		StructureOrg temp1 = new StructureOrg(Cell.parseCif(new File("/home/wtipton/POSCAR4.cif")));
		StructureOrg temp2 = new StructureOrg(Cell.parseCif(new File("/home/wtipton/POSCAR3.cif")));	
		temp1.setTotalEnergy(1 * temp1.getCell().getBasisSize());
		temp2.setTotalEnergy(1 * temp2.getCell().getBasisSize());
//		System.out.println(bob.sOrgToString(temp1));
//		System.out.println(bob.sOrgToString(temp2));
		bob.addStructureOrg(temp1);
		System.out.print("They match: ");
		System.out.println((bob.checkStructureOrg(temp2) != null));
		*/

	}
	
}
