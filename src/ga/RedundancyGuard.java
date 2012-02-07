/* Genetic algorithm for crystal structure prediction.  Will Tipton.  Ceder Lab at MIT. Summer 2007 */

package ga;

import java.io.*;
import java.util.*;

import crystallography.Cell;

// RedundancyGuard is used by the algorithm to avoid considering identical
// StructureOrgs more than once.  It stores a Map of all structures "seen"
// to their Organism IDs.  Then, other algorithms (e.g. StructureDev) may
// check Organisms against the list of those already seen.  Structure comparisons
// are done first by a direct comparison of lattice parameters and atomic
// positions and second, if necessary, by the slower but more accurate
// StructureFitter.

public class RedundancyGuard implements Serializable {
	static final long serialVersionUID = 1;

	// holds all the structures the algorithm has "seen"
	// and maps them to their Organism IDs
	Map<Cell,Integer> structures;

	// 
	//private double lAngleTol = 0.5;
	//private double lLengthTol = 0.5;
	//private double interval = 0.05;
	
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
		
		//lAngleTol = Double.parseDouble(args[0]);
		//lLengthTol = Double.parseDouble(args[1]);
		//interval = Double.parseDouble(args[2]);

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
			if (GAParameters.getParams().getVerbosity() >= 3)
				System.out.println("Warning: RedundancyGuard got passed structure with no sites. Ignoring...");
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
	
	// returns the ID of the matching StructureOrg if we've seen it before,
	// null otherwise
	public Integer checkStructureOrg(StructureOrg s1) {
		Cell s = s1.getCell();
		Iterator<Cell> i;
		
		// loop through all the Structures we've seen and compare
/*		if (useBothComparators) {
			i = structures.keySet().iterator();
			while (i.hasNext()) {
				Structure t = i.next();
				if (sOrgEquals(s,t))
					return structures.get(t);
			}
		}
*/	
		// if all those seem ok, do the more thorough checks:
		i = structures.keySet().iterator();
		while (i.hasNext()) {
			Cell t = i.next();
			// don't fail an organism for looking like itself
		//	if (structures.get(t) != s1.getID() && s.matchesCell(t, atomicMisfit, latticeMisfit))
		//		return structures.get(t);
			
			// temporary comparison
	/*		double ENERGY_TOL = 0.01;
			double VOL_TOL = 0.1;
			Cell c1 = SpgLib.getInstance().getPrimitiveCell(s);
			Cell c2 = SpgLib.getInstance().getPrimitiveCell(t.getCell());
		*/	
		//	System.out.println(SpgLib.getInstance().getSpaceGroup(c2));
	/*		System.out.println(c1);
			System.out.println(s);
			
			System.out.println(SpgLib.getInstance().getSpaceGroup(s));
			System.out.println(SpgLib.getInstance().getSpaceGroup(t.getCell()));
		*/	
		/*	if (structures.get(t) != s1.getID() && 
					(SpgLib.getInstance().getSpaceGroup(c1) == SpgLib.getInstance().getSpaceGroup(c2)) &&
					(Math.abs(c1.getVolume() - c2.getVolume()) < VOL_TOL) &&
					(Math.abs(t.getEnergyPerAtom() - s1.getEnergyPerAtom()) < ENERGY_TOL) &&
					(c1.getBasisSize() == c2.getBasisSize()))
					*/
			if (usePBCs && s.matchesCell(t, atomicMisfit, latticeMisfit, angleMisfit))
				return structures.get(t);
			if (!usePBCs && s.matchesCellNoPBCs(t, atomicMisfit))
				return structures.get(t);
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
