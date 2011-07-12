/* Genetic algorithm for crystal structure prediction.  Will Tipton.  2010. */

package ga;

import java.io.*;
import java.text.DecimalFormat;
import java.util.*;

import chemistry.Element;

import utility.*;

import crystallography.Cell;
import crystallography.Site;

// LammpsEnergy computes the total energy of a StructureOrg using Lammps and the given potential.
// It contains all of the methods and utilities that are specific to Lammps.

public class LammpsEnergy implements Energy {
	
	private String potlStr;
	
	private static final String dataFileName = "data.in";
	private static final String dumpFileName = "dump.atom";
	private static final String inFileName = "in.min";

	public LammpsEnergy(String[] args)
	{
		if (args == null || args.length < 1)
			GAParameters.usage("Not enough parameters given to LammpsEnergy", true);

		// read in the LAMMPS potential to use
		File potlFile = new File(args[0]);
		potlStr = GAUtils.readStringFromFile(potlFile);

	}

	public String toString() {
		StringBuilder result = new StringBuilder();
		
		result.append("Lammps total energy:" + GAUtils.newline());
		result.append(potlStr + GAUtils.newline());		
		
		return result.toString();
	}
	
	public double getEnergy(StructureOrg c) {
		return lammpsRun(c);
	}

	// returns a structure representation in format parse-able by Lammps
	private static String getLammpsDataFile(StructureOrg c) {
		
		Cell cell = c.getCell().getCellRotatedIntoPrincDirs();
		//example:
		/* 
header

3     atoms
1       atom types
-1.26815 1.24047      xlo xhi
-1.26815 1.24047      ylo yhi
-1.26815 1.24047      zlo zhi
0.1 0.1 0.1 xy xz yz

Masses

1 12 

Atoms

1       1       -2.9378454      -4.4592615      -4.8109196
2       1        5.6222143      -2.7335026      -1.7157569
3       1       -2.6614623      -5.5431059       1.6353686

*/
		
		StringBuilder result = new StringBuilder();
		String newline = GAUtils.newline();

		DecimalFormat df = new DecimalFormat();
		df.setMaximumFractionDigits(6);

		result.append("Written by LammpsEnergy" + newline + newline);
		result.append(cell.getBasisSize() + " atoms" + newline);
		result.append(cell.getComposition().getNumElements() + " atom types" + newline);
		
		// The parallelepiped has its "origin" at (xlo,ylo,zlo) and is defined by 3 edge vectors 
		// starting from the origin given by A = (xhi-xlo,0,0); B = (xy,yhi-ylo,0); C = (xz,yz,zhi-zlo).
		double xlo = 0, xhi = cell.getLatticeVectors().get(0).getCartesianComponents().get(0);
		double ylo = 0, yhi = cell.getLatticeVectors().get(1).getCartesianComponents().get(1);
		double zlo = 0, zhi = cell.getLatticeVectors().get(2).getCartesianComponents().get(2);
		double xy = cell.getLatticeVectors().get(1).getCartesianComponents().get(0);
		double xz = cell.getLatticeVectors().get(2).getCartesianComponents().get(0);
		double yz = cell.getLatticeVectors().get(2).getCartesianComponents().get(1);
		
		result.append(xlo + " " + xhi + " xlo xhi" + newline);
		result.append(ylo + " " + yhi + " ylo yhi" + newline);
		result.append(zlo + " " + zhi + " zlo zhi" + newline);
		result.append(xy + " " + xz + " " + yz + " xy xz yz" + newline + newline);
		
		// Masses
		result.append("Masses" + newline + newline);
		List<Element> elems = cell.getComposition().getElements();
		for (int i = 0; i < elems.size(); i++)
			result.append(i+1 + " " + elems.get(i).getAtomicMass() + newline);
		
		result.append(newline + "Atoms" + newline + newline);
		
        for (int i = 0; i < cell.getSites().size(); i++) {
        	Site s = cell.getSites().get(i);
        	
        	result.append(i+1 + " " + (1 + cell.getComposition().getElements().indexOf(cell.getSite(i).getElement())) + " ");
        	//for (double d : s.getCoords().getComponentsWRTBasis(c.getCell().getLatticeVectors()))
        	for (double d : s.getCoords().getCartesianComponents())
        		result.append(df.format(d) + " ");
        	result.append(newline);
        }
        
		return result.toString();
	}

	private static String getLammpsInputFile(StructureOrg c, String potlStr) {
		StringBuilder ans = new StringBuilder();
		String newline = GAUtils.newline();

		ans.append("units		metal" + newline);
		ans.append("dimension	3" + newline);
		ans.append("atom_style	atomic" + newline);
		ans.append("boundary	p p p" + newline);
		ans.append("read_data	" + dataFileName + newline);

		/*
		# GA user supplies this part
		pair_style	lj/cut 2.5
		pair_coeff	1 1 1.0 1.0 2.5
		pair_modify	shift yes
		*/
		ans.append(newline + potlStr + newline);

		ans.append("minimize 0.0 1.0e-8 1 1" + newline);
		ans.append("dump myDump all atom 100000000000000 " + dumpFileName + newline);
		ans.append("dump_modify myDump sort 1 scale no" + newline);
		ans.append("fix 1 all box/relax iso 1e4 vmax 0.001" + newline);
		ans.append("minimize 0.0 1.0e-8 10000 100000 " + newline);
		
		return ans.toString();
	}

	// runs Lammps on the input file given and returns the results in a String
	private static String runLAMMPS(String runDir) {
		int verbosity = GAParameters.getParams().getVerbosity();
		StringBuilder lammpsOutput = new StringBuilder();

		String s = null;
		try {
			// run the gulp command. GULP will only take input from it's stdin,
			// not from a file. in order to avoid having
			// to open inputFile and feed it to GULP, we use a wrapper script,
			// callgulp, which basically is just gulp <$1 .
			Process p = Runtime.getRuntime().exec("calllammps " + runDir);

			BufferedReader stdInput = new BufferedReader(new InputStreamReader(
					p.getInputStream()));
			BufferedReader stdError = new BufferedReader(new InputStreamReader(
					p.getErrorStream()));

			try {
				// read the output
				while ((s = stdInput.readLine()) != null) {
					lammpsOutput.append(s + GAUtils.newline());
					if (verbosity >= 5)
						System.out.println(s);
				}
	
				// print out any errors
				while ((s = stdError.readLine()) != null) {
					System.out.println(s);
				}
			} finally {
				stdInput.close();
				stdError.close();
			}

		} catch (IOException e) {
			System.out.println("IOException in LammpsEnergy.runLammps: " + e.getMessage());
			System.exit(-1);
		} 

		return lammpsOutput.toString();
	}

	// Does a Lammps run on the StructureOrg c. 
	// And we update c to be the local minimum found.
	private double lammpsRun(StructureOrg c) {
		GAParameters params = GAParameters.getParams();
		double finalEnergy = Double.POSITIVE_INFINITY;
		int verbosity = GAParameters.getParams().getVerbosity();
		
		// make temp directory
		String outDirPath = params.getTempDirName() + "/" + params.getRunTitle() + "." + c.getID();
		File outDir = new File(outDirPath);
		outDir.mkdir();
		
		// write unrelaxed cell to disk
		c.getCell().writeCIF(outDirPath + "/" + c.getID() + ".unrelaxed.cif");
		
		utility.Utility.writeStringToFile(getLammpsInputFile(c, potlStr), outDirPath + "/" + inFileName);
		utility.Utility.writeStringToFile(getLammpsDataFile(c), outDirPath + "/" + dataFileName);
		if (verbosity >= 3)
			System.out.println("Starting Lammps computation on organism " + c.getID());
		
		String lammpsOutput = runLAMMPS(outDir.getAbsolutePath());

		// store the results of the run
		// update c to be the structure in Lammps's output
		Cell a = parseOutputStructure(c.getCell(), outDirPath + "/" + dumpFileName);
		if (a == null) {
			if (verbosity >= 3)
				System.out.println("Warning: bad Lammps output.  Not updating structure.");
		} else {
			c.setCell(a);
		}
		
		// write cell to disk
		//Utility.writeStringToFile(c.getCIF(), outDirPath + "/" + c.getID() + ".relaxed.cif");

		finalEnergy = parseFinalEnergy(lammpsOutput);
		
		if (verbosity >= 3)
			System.out.println("Energy of org " + c.getID() + ": " + finalEnergy + " ");

		return finalEnergy;
	}
	
	private Cell parseOutputStructure(Cell origCell, String outFile) {
		// make sure out file was created successfully
		
		if (! (new File(outFile)).exists()) {
			if (GAParameters.getParams().getVerbosity() >= 4)
				System.out.println("In LammpsEnergy.parseOutputStructure: dump file doesn't exist.");
			return null;
		}
		
		String output = Utility.readStringFromFile(outFile);

		List<Vect> basis = new LinkedList<Vect>();
		List<Site> sites = new LinkedList<Site>();
		
		String lines[] = output.split("\n");
		
		/* should look like: (w/o the numbers)
0: ITEM: TIMESTEP
1: 466
2: ITEM: NUMBER OF ATOMS
3: 1
4: ITEM: BOX BOUNDS xy xz yz
5: 0.463004 1.537 9.54527e-17
6: 0.463004 1.537 1.0302e-16
7: 0.463004 1.537 9.09363e-17
8: ITEM: ATOMS id type xs ys zs
9: 1 1 0 0 0
10-etc: more atoms like line 9
		 */
		
		int numAtoms = Integer.parseInt(lines[3]);
		
		String xs[] = lines[5].split(" ");
		double xlo = Double.parseDouble(xs[0]);
		double xhi = Double.parseDouble(xs[1]);
		double xy = Double.parseDouble(xs[2]);
		String ys[] = lines[6].split(" ");
		double ylo = Double.parseDouble(ys[0]);
		double yhi = Double.parseDouble(ys[1]);
		double xz = Double.parseDouble(ys[2]);
		String zs[] = lines[7].split(" ");
		double zlo = Double.parseDouble(zs[0]);
		double zhi = Double.parseDouble(zs[1]);
		double yz = Double.parseDouble(zs[2]);
		
		basis.add(new Vect(xhi-xlo,0.0,0.0));
		basis.add(new Vect(xy, yhi-ylo, 0.0));
		basis.add(new Vect(xz, yz, zhi-zlo));
		
		for (int i = 0; i < numAtoms; i++) {
			String as[] = lines[9+i].split(" ");
			int type = Integer.parseInt(as[1]) - 1;
			double x = Double.parseDouble(as[2]);
			double y = Double.parseDouble(as[3]);
			double z = Double.parseDouble(as[4]);
			sites.add(new Site(origCell.getComposition().getElements().get(type),
					   new Vect(x - xlo,y - ylo,z - zlo)));
		}
		
		if (origCell.getBasisSize() != numAtoms) {
			System.out.println("ERROR: Lammps calculation lost atoms");
			System.out.println(outFile);
		}
		
		return new Cell(basis, sites, origCell.getLabel());
	}
	
	private double parseFinalEnergy(String lammpsOutput) {
		Double finalEnergy = Double.POSITIVE_INFINITY;
		int verbosity = GAParameters.getParams().getVerbosity();
		
		String lines[] = lammpsOutput.split("\n");
		int i;
		for (i = 0; i < lines.length; i++) {
			if (lines[i].matches(".*Energy .* final.*")) {
				String energies[] = lines[i+1].split("  *");
				finalEnergy = Double.parseDouble(energies[2]);
			}
		}

		// toss sketchy structures:
		/*if (cautious 
				&& (lammpsOutput.contains("warning") || lammpsOutput.contains("failed") || lammpsOutput.contains("caution")
						|| lammpsOutput.contains("Maximum number of function calls"))) {
			if (verbosity >= 4)
				System.out.println("Lammps run failed.");
			return 0;
		}*/
				
		return finalEnergy;
	}
	
	// just for testing:
	public static void main(String[] args) {
		String[] geArgs = {"/home/wtipton/temp/lammpspot"};
		LammpsEnergy bob = new LammpsEnergy(geArgs);
		StructureOrg c = new StructureOrg(Cell.parseCif(new File("/home/wtipton/POSCAR4.cif")));
//		Utility.writeStringToFile(c.getCIF(), "bob.cif");
		double energy = bob.getEnergy(c);
		System.out.println(energy);
		System.out.println(c.getCell());
		
		
	}
}
