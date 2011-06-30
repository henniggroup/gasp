package ga;

import java.io.File;
import java.util.ArrayList;
import java.util.LinkedList;
import java.util.List;
import java.util.Random;

import Jama.Matrix;

import chemistry.Composition;
import chemistry.Element;

import utility.Constants;
import utility.RandomNumbers;
import utility.Utility;
import utility.Vect;
import crystallography.Cell;
import crystallography.Site;

public class UnitsSOCreator implements StructureOrgCreator {
	
	private String[] atoms;
	private double[][] coords;
	private int[] numAtoms;
	private int[] numUnits;
	private int difUnits;
	private static boolean unitsOnly;
	private List<Site> sites;
	private List<Vect> refLoc;
	private int m = 1; //For testing only
	private int[] units;
	private double targetDensity;
	private double densityTol;
	private int numSites;
	
// ------------------------------------------------------------
	
	public UnitsSOCreator(String[] args) {
		
		// Length of arguments
		int length = args.length;
		
		// Parse unitsOnly
		unitsOnly = Boolean.parseBoolean(args[length-1]);
		
		// Parse density parameters
		densityTol = Double.parseDouble(args[length-2]);
		targetDensity = Double.parseDouble(args[length-3]);
		
		// Parse number of different molecules
		difUnits = Integer.parseInt(args[0]);
		
		// Marks location of coordinates
		int cStart = difUnits + 1;
		
		numAtoms = new int[difUnits]; numUnits = new int[difUnits];
		
		// Parse in number of atoms per molecule
		numSites = 0;
		for (int i=0; i<difUnits; i++) {
			int l = Integer.parseInt(args[1+i]);
			numAtoms[i] = l;
			numSites = numSites + l;
		}
		
		// Marks location after coordinates end
		int eStart = cStart + 4*numSites;
		
		// Parse number of each molecule to be added
		for (int i=0; i<difUnits; i++) {
			numUnits[i] = Integer.parseInt(args[eStart+i]);
		}
		
		// String array of all atoms
		atoms = new String[numSites];
		for (int i=0; i<numSites; i++) {
			atoms[i] = args[cStart + i*4];
		}
		
		// Double array of atom coordinates
		coords = new double[numSites][Constants.numDimensions];
		int pos = 1;
		for (int r=0; r<numSites; r++) {
			for (int c=0; c<Constants.numDimensions; c++) {
				if (pos % 4 != 0) {
					coords[r][c] = Double.parseDouble(args[cStart + pos]);
					pos++;
				}
				else {
					coords[r][c] = Double.parseDouble(args[cStart + pos + 1]);
					pos = pos + 2;
				}
			}
		}
		
		// Creates list of sites and list of (unique) elements involved
		sites = new LinkedList<Site>();
		for (int i=0; i<numSites; i++) {
			Vect v = new Vect(coords[i][0],coords[i][1],coords[i][2]);
			Element e = Element.getElemFromSymbol(atoms[i]);
			Site s = new Site(e,v);
			sites.add(s);
		}
		
	}

/*
	// TODO: make this work for multiple units (giant pain in the ass)
	public String toString() {
		boolean u = unitsOnly;
		
		String a = ""; String c = ""; String n = ""; String d = ""; String f = "";
		
		for (int i=0; i<numAtoms; i++) {
			if (i<numAtoms-1) {
				a = a + atoms[i] + ", ";
			}
			else {
				a = a + atoms[i];
			}
			if (i<numAtoms-1) {
				c = c + "(" + coords[i][0] + "," + coords[i][1] + "," + coords[i][2] + "), ";
			}
			else {
				c = c + "(" + coords[i][0] + "," + coords[i][1] + "," + coords[i][2] + ")";
			}
		}
		
		if (numUnits == 0) {
			n = n + "random";
		}
		else {
			n = n + numUnits;
		}
		
		if (targetDensity == 0) {
			d = d + "generated";
		}
		else {
			d = d + targetDensity + " g/cm^3";
		}
		
		if (densityTol == 0) {
			f = f + "20%";
		}
		else {
			double perc = densityTol*100;
			f = f + perc + "%";
		}
		
		return "UnitsSOCreator: atoms " + a + " with coordinates " + c + " and unitsOnly " + u + ", # molecules = " + n + ", target density " + d + " and tolerance " + f;
	}

*/
	
	// Fills up the population with unit organisms (provided they satisfy the
	// hard constraints)
	protected StructureOrg makeUnitOrg() {
		GAParameters params = GAParameters.getParams();
		Random rand = params.getRandom();
		
/*		// Checks interatomic distances within molecule
		for (int n=0; n<numAtoms; n++) {
			for (int m=0; m<numAtoms; m++) {
				if (n!=m) {
					Vect nvect = sites.get(n).getCoords(); Vect mvect = sites.get(m).getCoords();
					Vect dif = nvect.subtract(mvect);
					if (dif.length() < params.getMinInteratomicDistance()) {
						GAParameters.usage("Minimum interatomic distance is greater than distances within molecule", true);
					}
				}
			}
		}
*/	
		
		// make random lattice parameters satisfying hard constraints
		// the Basis for our new organism
		List<Vect> latVects = RandomSOCreator.makeRandomLattice();

		// make lists to hold the points and species which will make our new Structure
		ArrayList<Site> sitesList = new ArrayList<Site>();

		final int maxFails = 100;
		int failCount = 0;
		

		int totNum = 0;
		for (int r=0; r<difUnits; r++) {
			totNum = totNum + numUnits[r]*numAtoms[r];
		}
		if (totNum > params.getMaxNumAtoms()) {
			GAParameters.usage("Number of units exceeds max number of atoms constraint", true);
		}	
		
		// Parse or generate target number of molecules of each type, ensuring total # atoms is within constraint
		int targetAtoms = RandomNumbers.getUniformIntBetween(numSites,params.getMaxNumAtoms()); int totAtoms = 0;
		double[] fracValues = new double[difUnits]; double fracTotal = 0.0; int[] target = new int[difUnits]; units = new int[difUnits];
		if (numUnits[0] != 0) {
			for (int t=0; t<difUnits; t++) {
				target[t] = numUnits[t];
				units[t] = target[t];
			}
		}
		else {
			for (int y=0; y<difUnits; y++) {
				fracValues[y] = RandomNumbers.getUniformDoubleBetween(0.0, 1.0);
				fracTotal = fracTotal + fracValues[y];
			}
			for (int y=0; y<difUnits; y++) {
				fracValues[y] = fracValues[y]/fracTotal;
				target[y] = (int) Math.round(fracValues[y]*targetAtoms);
				target[y] = target[y]/numAtoms[y];
				units[y] = target[y];
				totAtoms = totAtoms + target[y]*numAtoms[y];
			}
		}
		
/*		units = new int[difUnits]; int[] target = new int[difUnits]; int totAtoms = 0;
		while (totAtoms==0 || totAtoms > params.getMaxNumAtoms()) {
			for (int y=0; y<difUnits; y++) {
				if (numUnits[y] != 0) {
					target[y] = numUnits[y];
					units[y] = target[y];
					totAtoms = totAtoms + numUnits[y]*numAtoms[y];
				}
				else {
				target[y] = RandomNumbers.getUniformIntBetween(1, params.getMaxNumAtoms());
				units[y] = target[y];
				totAtoms = totAtoms + target[y]*numAtoms[y];
				}
			}
		}
*/
		
/*		int totAtoms = 0; units = new int[difUnits]; int[] target = new int[difUnits];
		for (int i=0; i<difUnits; i++) {
			if (numUnits[i] != 0) {
				target[i] = numUnits[i];
				units[i] = target[i];
				totAtoms = totAtoms + numUnits[i]*numAtoms[i];
			}
			else {
				int maxAtoms = params.getMaxNumAtoms(); int minAtoms = params.getMinNumAtoms();
				int randInt = RandomNumbers.getUniformIntBetween(numAtoms[i], (maxAtoms - totAtoms));
				target[i] = randInt/numAtoms[i];
				units[i] = target[i];
				totAtoms = totAtoms + target[i]*numAtoms[i];
			}
			if (target[i] == 0 && (params.getMaxNumAtoms() - totAtoms) >= numAtoms[i]) {
				int randInt = RandomNumbers.getUniformIntBetween(numAtoms[i], params.getMaxNumAtoms() - totAtoms);
				target[i] = randInt/numAtoms[i];
				units[i] = target[i];
				totAtoms = totAtoms + target[i]*numAtoms[i];				
			}			
		}
*/


		// For testing only: creates directory for pre-relaxation structures
		File outDir = new File("/home/skw57/structures");
		outDir.mkdir();
		
		// List of locations of previously placed units
		refLoc = new LinkedList<Vect>();
		
		// Loop through the different types of molecules
		for (int i = 0; i < difUnits; i++) {
			
			// Loop cycles through unit additions until it hits target number
			for (int k = 0; k < target[i]; k++) {
				
				int currentSite = 0;
				for (int m=0; m<i; m++) {
					currentSite = currentSite + numAtoms[m];
				}
				
				// Initializes and creates list of vectors from input file (original unit geometry)
				List<Vect> basis = new LinkedList<Vect>();
				for (int n=0; n<numAtoms[i]; n++) {
					Site s = sites.get(n + currentSite);
					basis.add(s.getCoords());
				}
				
				// Returns list of randomly rotated vectors
				List<Vect> rotatedBasis = getRandomRotation3D(basis,i);
				
				// Creates list of sites from new vectors
				List<Site> newSites = new LinkedList<Site>();
				for (int p=0; p<numAtoms[i]; p++) {
					Site s = sites.get(p + currentSite);
					newSites.add(new Site(s.getElement(), rotatedBasis.get(p)));
				}
				
				// Generates potential location for unit
				Vect potentialLocation = new Vect(rand.nextDouble(),rand.nextDouble(),rand.nextDouble(), latVects);
				
				// TODO: This code seems ineffective, probably just a waste of resources
				// If this is not the first unit placement
				if (refLoc.size() > 0) {
					boolean regen = true;
					while (regen) {
						regen = false;
						// Calculate distances between potential location and previous placements
						for (Vect v : refLoc) {
							double dist = potentialLocation.subtract(v).length();
							if (dist < params.getMinInteratomicDistance()) {
								regen = true;
							}
							// If distance is less than minInteratomicDistance, regenerate potential location
							potentialLocation = new Vect(rand.nextDouble(),rand.nextDouble(),rand.nextDouble(), latVects);
						}
					}
				}
				
				
				if ((new Cell(latVects,sitesList)).getAtomsInSphereSorted(potentialLocation, params.getMinInteratomicDistance()).size() == 0) {
					// Actual "placement" of atoms occurs here
					for (Site r: newSites) {
						Element e = r.getElement();
						Vect relativeLoc = r.getCoords();
						// To preserve unit geometry, adds (rotated) input vector to potentialLocation
						Vect v = potentialLocation.plus(relativeLoc);
						sitesList.add(new Site(e,v)); // Adds atom to list of sites
						refLoc.add(potentialLocation);
						
						// for debugging purposes
//						System.out.println("potentialLocation " + potentialLocation + " and relativeLoc " + relativeLoc + " and finalLoc " + v);
					}	
				}
				else if (failCount < maxFails) {
					failCount++;
					k--;
				}
			}
		}
			

			// Add in atoms to match composition space if applicable
			if (!unitsOnly) {
				Composition comp = params.getCompSpace().getRandomIntegerCompInSpace(params.getMinNumAtoms(), (params.getMaxNumAtoms() - totAtoms));
				for (Element e : comp.getElements())
					// add a stoichiometry's worth of that species
					for (int k = 0; k < comp.getOrigAmount(e); k++) {
						Vect potentialLocation = new Vect(rand.nextDouble(),rand.nextDouble(),rand.nextDouble(), latVects);
						if ((new Cell(latVects,sitesList)).getAtomsInSphereSorted(potentialLocation, params.getMinInteratomicDistance()).size() == 0)
							sitesList.add(new Site(e,potentialLocation));
						else if (failCount < maxFails) {
							failCount++;
							k--;
						}
					}
			}
			
		// Creates a cell based on the generated sites and vectors
		Cell newStructure = new Cell(latVects, sitesList);
		
		// Optimizes cell density
		Cell optimizedStructure = optimizeDensity(newStructure);
		
		// For testing purposes only: writes pre-relaxation structure to file
		Utility.writeStringToFile(optimizedStructure.getCIF(), "/home/skw57/structures/struct" + m +".cif");
		m++;

		return new StructureOrg(optimizedStructure);
	}
	
	public StructureOrg makeOrganism(Generation g) {
		
		StructureOrg o = makeUnitOrg();
		
		return o;
	}
	
	public List<Vect> getRandomRotation3D(List<Vect> basis, int q) {
		
		int dim = Constants.numDimensions;
		
		int currentSite = 0;
		for (int m=0; m<q; m++) {
			currentSite = currentSite + numAtoms[m];
		}
		
		// Generates random angles (degrees)
		double randx = RandomNumbers.getUniformDoubleBetween(0, 360);
		double randy = RandomNumbers.getUniformDoubleBetween(0, 360);
		double randz = RandomNumbers.getUniformDoubleBetween(0, 360);
		
		// Change from degrees to radians
		randx = randx*(Math.PI/180); randy = randy*(Math.PI/180); randz = randz*(Math.PI/180);
		
		// Initialize X,Y,Z rotation matrices
		Matrix X = new Matrix(dim,dim); Matrix Y = new Matrix(dim,dim); Matrix Z = new Matrix(dim,dim);
		X.set(0,0,1); X.set(0,1,0); X.set(0,2,0); X.set(1,0,0); X.set(1,1,Math.cos(randx)); X.set(1,2,-Math.sin(randx)); X.set(2,0,0); X.set(2,1,Math.sin(randx)); X.set(2,2,Math.cos(randx));
		Y.set(0,0,Math.cos(randy)); Y.set(0,1,0); Y.set(0,2,Math.sin(randy)); Y.set(1,0,0); Y.set(1,1,1); Y.set(1,2,0); Y.set(2,0,-Math.sin(randy)); Y.set(2,1,0); Y.set(2,2,Math.cos(randy));
		Z.set(0,0,Math.cos(randz)); Z.set(0,1,-Math.sin(randz)); Z.set(0,2,0); Z.set(1,0,Math.sin(randz)); Z.set(1,1,Math.cos(randz)); Z.set(1,2,0); Z.set(2,0,0); Z.set(2,1,0); Z.set(2,2,1);
		
		// Creates full 3D rotation matrix
		Matrix R = X.times(Y).times(Z);
		
		// Form matrix from list of vectors
		Matrix V = new Matrix(dim,numAtoms[q]);
			for (int r=0; r<dim; r++) {
				for (int c=0; c<numAtoms[q]; c++) {
					V.set(r, c, coords[c+currentSite][r]);
				}
			}
		
		// Generated rotated matrix
		Matrix P = R.times(V);
		
		// Extract transformed list of vectors from matrix
		List<Vect> newBasis = new LinkedList<Vect>();
		for (int i=0; i<numAtoms[q]; i++) {
			Vect v = new Vect(P.get(0,i),P.get(1,i),P.get(2,i));
			newBasis.add(v);
		}
		
		return newBasis;

	}	
		
	public Cell optimizeDensity(Cell structure) {
		// Obtain list of sites and lattice vectors from cell
		List<Site> sites = structure.getSites();
		List<Vect> vects = structure.getLatticeVectors();
		
		// Conversion from amu/angstrom^3 to g/cm^3
		double amuang_to_gmcm = 1.6606;
		
		// Tolerance level (i.e. +/- percentage of target density allowable)
		double tol = 0;
		if (densityTol == 0) {
			tol = 0.2;
		}
		else {
			tol = densityTol;
		}
		
		// Calculate initial volume
		Vect va = vects.get(0); Vect vb = vects.get(1); Vect vc = vects.get(2);
		double initialVol = Math.abs(va.dot(vb.cross(vc)));
		
		// Calculate mass of atoms currently placed
		double mass = 0.0;
		for (Site s : sites) {
			double m = s.getElement().getAtomicMass();
			mass = mass + m; // amu
		}
		
		// Calculate current density
		double density = mass/initialVol; // amu/angstrom^3
		density = density*amuang_to_gmcm; // g/cm^3
		
		// Calculate target density
		double target = 0;
		if (targetDensity == 0) {
			double total = 0.0;
			int count = 0;
			for (Site s : sites) {
				double d = s.getElement().getDensity();
				total = total + d;
				count++;
			}
			target = total/count; // g/cm^3
		}
		else {
			target = targetDensity;
		}
		
		double minTarget = (1-tol)*target; double maxTarget = (1+tol)*target;
		
		GAParameters params = GAParameters.getParams();
		Random rand = params.getRandom();
		
		List<Vect> newVects = new LinkedList<Vect>();
		
		// Create new cell and shift unit locations, recalculate density until it is within tolerance
		int c = 0;
//		int maxFails = 500;
		
		while (density < minTarget || density > maxTarget) {
			
//			if (c > maxFails) {
//				break;
//			}
			
			if (density == 0) {
				break;
			}
			
			// Generate new lattice parameters
			newVects = RandomSOCreator.makeRandomLattice();
			
		// Generate random vector, shift current molecules to these locations
			int currentSite = 0;
			for (int u=0; u<difUnits; u++) {	
				for (int i=0; i<units[u]; i++) {
				
					Vect mainLoc = sites.get(currentSite).getCoords();

					// Generates potential location for unit
					Vect potentialLoc = new Vect(rand.nextDouble(),rand.nextDouble(),rand.nextDouble(), newVects);
				
					// Shift is difference between site location and potential location
					Vect locShift = mainLoc.subtract(potentialLoc);
				
					// Cycles through individual atoms in unit (adding the same shift to each)
					for (int k=0; k<numAtoms[u]; k++) {
						Vect v = sites.get(currentSite+k).getCoords();
						Element e = sites.get(currentSite+k).getElement();
						// Add the position shift vector
						Vect vNew = v.plus(locShift);
						Site newSite = new Site(e,vNew);
						sites.set(currentSite+k, newSite);
					}
					currentSite = currentSite + numAtoms[u];
				}
			}
			
			// Generate random vector, shift current individual atoms to these locations (if !unitsOnly)
			if (!unitsOnly) {
				for (int k=currentSite; k<sites.size(); k++) {
					Vect mainLoc = sites.get(k).getCoords();
					Vect potentialLoc = new Vect(rand.nextDouble(),rand.nextDouble(),rand.nextDouble(), newVects);
					Vect locShift = mainLoc.subtract(potentialLoc);
					Element e = sites.get(k).getElement();
					Vect vNew = mainLoc.plus(locShift);
					Site newSite = new Site(e,vNew);
					sites.set(k, newSite);
				}
			}
			
			// Recalculate cell volume
			Vect va2 = newVects.get(0); Vect vb2 = newVects.get(1); Vect vc2 = newVects.get(2);
			double currentVol = Math.abs(va2.dot(vb2.cross(vc2)));
			
			// Recalculate cell density
			density = mass/currentVol; // amu/angstrom^3
			density = density*amuang_to_gmcm; // g/cm^3
			
			c++;
			
			System.out.println("Density attempt: " + density + ". Target density: " + target);
		}
		
		if (c == 0) {
			newVects = vects;
		}
		
//		System.out.println("Final density: " + density + ". Target density: " + target);
		return new Cell(newVects, sites);
		
	}
	
	
/*	public static void main(String[] args) {
		String[] params = { "--minInteratomicDistance", "0.8", "--maxLatticeLength", "35", "--minLatticeLength", "1", "--maxLatticeAngle", "140", "--minLatticeAngle", "40", "--maxNumAtoms", "20", "minNumAtoms", "1", "doNonnegativityConstraint", "false", "--compositionSpace", "2", "H", "O", "0", "1", "1", "0", "--ObjectiveFunction", "epa", "gulp", ""};
		GAParameters.getParams().setArgs(params); 
		String[] t = { "O", "0.0", "0.0", "0.0", "H", "0.94", "0.0", "0.0", "H", "0.0", "0.94", "0.0", "true", "3"};
		UnitsSOCreator will = new UnitsSOCreator(t);
//		System.out.println(will.makeUnitOrg().getCell().getCIF());
//		Utility.writeStringToFile(will.makeUnitOrg().getCell().getCIF(), "/home/skw57/stew.cif");
	}
*/	
	
// -----------------------------GETTERS------------------------------
	
	public String[] getAtoms() {
		return atoms;
	}
	
	public double[][] getCoords() {
		return coords;
	}
	
	public static boolean getUnitsOnly() {
		return unitsOnly;
	}
	
}