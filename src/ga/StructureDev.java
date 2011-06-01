/* Genetic algorithm for crystal structure prediction.  Will Tipton.  Ceder Lab at MIT. Summer 2007 */

package ga;

import gulp.GulpSurrogate;

import java.io.Serializable;
import java.util.*;

import utility.Utility;

import crystallography.Cell;
import crystallography.Site;

// StructureDev implements Development and is the Development operation
// generally used with StructureOrgs.  It checks various hard constraints. 
// It modifies the Organism only in the case where we are using the Niggli
// reduced cell. It modifies the generation only during the dValue and 
// RedundancyGuard checks.

public final class StructureDev implements Development, Serializable {
	static final long serialVersionUID = 1;
	
	private Boolean useNiggliReducedCell;
	
	private RedundancyGuard rGuard;
	private Boolean useWholePopRG;
	private Boolean usePerGenRG;
	
	private GulpSurrogate surrogate;

	public StructureDev() {
		GAParameters params = GAParameters.getParams();
		
		useNiggliReducedCell = params.getUseNiggliReducedCell();
		
		String rgType = params.getRedundancyGuardType();
		useWholePopRG = (rgType.equalsIgnoreCase("both")||rgType.equalsIgnoreCase("wholePopulation"));
		usePerGenRG = (rgType.equalsIgnoreCase("both")||rgType.equalsIgnoreCase("perGeneration"));
		String useStructureFitter = "true";
		if (useWholePopRG)
			rGuard = new RedundancyGuard((String[])Utility.appendArray(params.getRedundancyGuardArgs(), useStructureFitter));
		
		if (params.usingSurrogateModel()) 
			surrogate = new GulpSurrogate(params.getSurrogateArgs());
	}
	
	public String toString() {
		StringBuilder result = new StringBuilder();
		
		result.append("StructureDev development. useNiggliReducedCell: "
				+ useNiggliReducedCell + ", rGuard:" + rGuard);
		
		return result.toString();
	}
	
	// remove unphysical candidate crystals based on hard constraints
	public Boolean doDevelop(Generation gen, Organism o) {
		StructureOrg s = (StructureOrg)o;
		Structures g = (Structures)gen;
		
		GAParameters params = GAParameters.getParams();
		int verbosity = params.getVerbosity();
		// double maxid = params.getMaxInteratomicDistance();
		double minid = params.getMinInteratomicDistance();
		double maxll = params.getMaxLatticeLength();
		double minll = params.getMinLatticeLength();
		double maxla = params.getMaxLatticeAngle();
		double minla = params.getMinLatticeAngle();
		int maxatoms = params.getMaxNumAtoms();
		int minatoms = params.getMinNumAtoms();
		
		// fail if cell is null.
		// this can happen if, e.g., the VASP calculation on a structure fails
		//  and so, after failing to parse the VASP output, VaspOut sets
		//  the cell to null.
		// if we added the structure to the wholepop rGuard before doing the vasp
		//  calc and then the vasp calc sets its cell to null, we'll have problems.
		//  so, in this case, we remove the structure from the rGuard.
		if (s.getCell() == null) {
			if (verbosity >= 3)
				System.out.println("Organism " + s.getID() + " failed: had null cell.");
			if (rGuard != null)
				rGuard.removeStructureOrg(s);
			return false;
		}
		
		// use the Niggli reduced cell 
		if (useNiggliReducedCell ) {
			s.standardize();
			if (!s.isReduced()) {
				if (verbosity >= 3)
					System.out.println("Organism " + s.getID() + " failed Niggli reduction.");
				return false;
			}
		}
		
		Cell structure = s.getCell();
		
		// check the number of atoms
		if (minatoms != -1 && structure.getNumSites() < minatoms) {
			if (verbosity >= 3)
				System.out.println("Organism " + s.getID() + " failed min number of atoms constraint: (numAtoms = "
						+ structure.getNumSites() + ").");
			return false;
		}
		if (maxatoms != -1 && structure.getNumSites() > maxatoms) {
			if (verbosity >= 3)
				System.out.println("Organism " + s.getID() + " failed max number of atoms constraint: (numAtoms = "
						+ structure.getNumSites() + ").");
			return false;
		}
		
		// sanity check on the energy
		if (Double.isInfinite(s.getEnergyPerAtom())) {
			if (verbosity >= 3)
				System.out.println("Organism " + s.getID() + " discarded: energy is inf.");
			return false;
		}
		if (params.getDoNonnegativityConstraint()) {
			if (s.knowsValue() && s.getValue() >= 0) {
				if (verbosity >= 3)
					System.out.println("Organism " + s.getID() + " failed value nonnegativity constraint: (value = "
							+ s.getValue() + ").");
				return false;
			}
		} 
		
		// scale the volume if we haven't yet done the energy calculation
		// and we're beyond the 0th generation
		if (params.getOptimizeDensity() && !s.knowsValue() && GAParameters.getParams().getRecord().getGenNum() != 0) {
			double volume = GAParameters.getParams().getRecord().getBestDensityEstimate() * s.getCell().getNumSites();
			s.setCell(s.getCell().scaleTo(volume));
		}

		// check the lattice lengths
		double[] lLengths = structure.getCellLengths();
		for(int i = 0; i <= 2; i++) {
			if(minll != -1 && lLengths[i] < minll) {
				if (verbosity >= 3)
					System.out.println("Organism " + s.getID() + " failed min lattice length constraint: (lLengths["
							+ i + "] == " + lLengths[i] + ").");
				return false;
			}
			if(maxll != -1 && lLengths[i] > maxll) {
				if (verbosity >= 3)
					System.out.println("Organism " + s.getID() + " failed max lattice length constraint: (lLengths["
							+ i + "] == " + lLengths[i] + ").");
				return false;
			}
		}
		
		// check the lattice angles
		double[] lAngles = GAUtils.anglesToDegrees(structure.getCellAngles());
		for(int i = 0; i <= 2 ; i++) {
			if(minla != -1 && lAngles[i] < minla) {
				if (verbosity >= 3)
					System.out.println("Organism " + s.getID() + " failed min lattice angle constraint (lAngles["
							+ i + "] == " + lAngles[i] + ").");
				return false;
			}
			if(maxla != -1 && lAngles[i] > maxla) {
				if (verbosity >= 3)
					System.out.println("Organism " + s.getID() + " failed max lattice angle constraint (lAngles["
							+ i + "] == " + lAngles[i] + ").");
				return false;
			}
		}

		// Interatomic Distances
		if (minid != -1)
			for (int i = 0; i < structure.getNumSites(); i++) {
				// it's no good if there are any other atoms in the minimum radius sphere
				if (structure.getAtomsInSphereSorted(structure.getSite(i).getCoords(), minid).size() > 1) {
					if (verbosity >= 3)
						System.out.println("Organism " + s.getID() + " failed minimum interatomic distance constraint.");
					return false;
				}
			}	
		
		// check the stoichiometry
		/*
		if (fixStoich) {
			HashMap<String, Integer> constituents = params.getConstituents();
			// choose the first constituent in the crystal
			Set<String> symbols = constituents.keySet();
			Iterator<String> i = symbols.iterator();
			String firstSymbol = i.next();
			// find the number of stoichiometries worth of this element in the crystal
			Element spec = Element.getElemFromSymbol(firstSymbol);
			double numSymbol = s.getCell().getNumSitesWithElement(spec);
			double numStoichs = numSymbol / constituents.get(firstSymbol);
			// make sure it equals that of all the other elements
			while (i.hasNext()) {
				String currentSymbol = i.next();
				Element currentSpec = Element.getElemFromSymbol(currentSymbol);
				numSymbol = s.getCell().getNumSitesWithElement(currentSpec);
				if (numStoichs != (numSymbol/constituents.get(currentSymbol))) {
					if (verbosity >= 3)
						System.out.println("Organism " + s.getID() + " failed stoichiometry constraint.");
					return false;}
			}
		}*/
		// make sure composition is in compSpace
		if (!params.getCompSpace().contains(s.getCell().getComposition())) {
			if (verbosity >= 3)
				System.out.println("Organism " + s.getID() + " failed stoichiometry constraint.");
			return false;
		}
		
		// check nearest neighbors
		if (params.getMaxNearestNeighborLength() != 0) {
			List<Set<String>> excludedPairs = params.getNotNearestNeighbors();
			if (excludedPairs != null) {
				// for all the atoms in the structure
				List<Site> sites = structure.getSites();
				for (int i = 0; i < sites.size(); i++) {
					String symbol = sites.get(i).getElement().getSymbol();
					// check if there is a constraint on this species.  set otherSymb to the
					// symbol of the species the current one is not supposed to be next to
					String otherSymb = null;
					for (Set<String> pair : excludedPairs)
					{
						Iterator <String> k = pair.iterator();
						String first = k.next();
						String second = k.next();
						if (symbol == first)
							otherSymb = second;
						else if (symbol == second)
							otherSymb = first;
						else
							continue;
					}
					if (otherSymb == null)
						continue;
					// the current atom shouldnt have nearest neighbor with symbol otherSymb
					List<Site> atomsInSphere = structure.getAtomsInSphereSorted(sites.get(i).getCoords(), params.getMaxNearestNeighborLength());
					if (atomsInSphere.size() < 2)
						continue;
					String nearestNeighbor = atomsInSphere.get(1).getElement().getSymbol();
					if (otherSymb.equalsIgnoreCase(nearestNeighbor)) {
						if (verbosity >= 3)
							System.out.println("Organism " + s.getID() + " failed nearest neighbor constraint.");
						return false;
					}
				}
			}
		}
		
		// check redundancy guards ;
		// check against the perGeneration RG
		if (usePerGenRG) {
			Integer orgID = g.getRedundancyGuard().checkStructureOrg(s);
			if (orgID != null) {
				if (verbosity >= 3)
					System.out.println("Organism " + s.getID() + " failed perGeneration redundancy guard (looks like org "
							 + orgID + ").");
				// keep it anyway if it's better than the old one, but remove the old one
				Organism org = g.getOrgByID(orgID);
				// org can equal null if it has already been replaced in this generation
				if (org != null && s.knowsValue() && s.getValue() < org.getValue()) {
					g.removeOrganism(org);
					g.addOrganism(s);
					if (verbosity >= 4)
						System.out.println("perGeneration RedundancyGuard removing " + orgID
								+ " and replacing with " + s.getID());
				}
				return false;
			} 
		}
		// check against the wholePopulation RG only if we haven't done the energy calculation
		if (useWholePopRG && !s.knowsValue()) {
			Integer orgID = rGuard.checkStructureOrg(s);
			if (orgID != null) {
				if (verbosity >= 3)
					System.out.println("Organism " + s.getID() + " failed wholePopulation redundancy guard (looks like org "
							 + orgID + ").");
				return false;
			} 
		}

		// the dValue rule
		double dValue = params.getDValue();
		if (dValue != 0)
			if (o.knowsValue()) {
				Organism[] orgs = g.getOrganismsOfValue(o.getValue(), dValue);
				if (orgs.length >= 1) {
					if (verbosity >= 3)
						System.out.println("Organism " + s.getID() + " failed dValue constraint.");
					
					// keep it anyway if it's better than the old one, but remove the old one
					if (o.getValue() < orgs[0].getValue()) {
						g.removeOrganism(orgs[0]);
						g.addOrganism(o);
						if (verbosity >= 4)
							System.out.println("dValue rule removing " + orgs[0].getID()
									+ " and replacing with " + o.getID());
					}
					return false;
				}
			}
		
		// check w/ surrogate model
		if (s.knowsValue() && params.usingSurrogateModel()) {
			surrogate.addEntry(s);
			if (surrogate.fails(s))
				return false;
		}
		
		// if we're using the RedundancyGuards, add o to the list
		// (orgs added to the perGen RG in Structures.addOrganism())
		if (useWholePopRG)
			rGuard.addStructureOrg(s);
		
		return true;
	}
}