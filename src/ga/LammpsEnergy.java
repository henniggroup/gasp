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

import java.io.*;
import java.text.DecimalFormat;
import java.util.*;

import chemistry.CompositionSpace;
import chemistry.Element;

import utility.*;
import vasp.VaspIn;

import crystallography.Cell;
import crystallography.Site;

// LammpsEnergy computes the total energy of a StructureOrg using Lammps and the given potential.
// It contains all of the methods and utilities that are specific to Lammps.

public class LammpsEnergy implements Energy {
	
	private String potlStr;
	private String unitsStr;
	private boolean relax_box;
	
	private static final String dataFileName = "data.in";
	private static final String dumpFileName = "dump.atom";
	private static final String inFileName = "in.min"; 

	public LammpsEnergy(List<String> args)
	{
		if (args == null || args.size() < 2)
			GAParameters.usage("Not enough parameters given to LammpsEnergy", true);

		// read in the LAMMPS potential to use
		File potlFile = new File(args.get(0));
		potlStr = GAUtils.readStringFromFile(potlFile);
		unitsStr = args.get(1);
		
		if (args.size() > 2)
			relax_box = Boolean.parseBoolean(args.get(1));
		else
			relax_box = true;

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
	public static String getLammpsDataFile(Cell cell) {
		
		cell.rotatedIntoPrincDirs();
		CompositionSpace compSpace = GAParameters.getParams().getCompSpace();
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
		// turns out we need to tell it there's an atom type for each entry in the potential file
	//	result.append(cell.getComposition().getNumElements() + " atom types" + newline);
		result.append(compSpace.getElements().size() + " atom types" + newline);
		
		// The parallelepiped has its "origin" at (xlo,ylo,zlo) and is defined by 3 edge vectors 
		// starting from the origin given by A = (xhi-xlo,0,0); B = (xy,yhi-ylo,0); C = (xz,yz,zhi-zlo).
		//
		// Basically, if the 3 lattice vectors are (a1,a2,a3),(b1,b2,b3),(c1,c2,c3), then we're now
		// writing them as (a1,0,0),(b1,b2,0),(c1,c2,c3) which is not a problem since we called
		// cell.getCellRotatedIntoPrincDirections() so that a2=a3=b3=0.
		
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
		List<Element> elems = compSpace.getElements();
		for (int i = 0; i < elems.size(); i++)
			result.append(i+1 + " " + elems.get(i).getAtomicMass() + newline);
		
		result.append(newline + "Atoms" + newline + newline);
		
        for (int i = 0; i < cell.getSites().size(); i++) {
        	Site s = cell.getSites().get(i);
        	
        	int charge = 0; //TODO: change?
        	result.append((i+1) + " " + (1 + compSpace.getElements().indexOf(cell.getSite(i).getElement())) + " " + charge + " ");
        	//for (double d : s.getCoords().getComponentsWRTBasis(c.getCell().getLatticeVectors()))
        	for (double d : s.getCoords().getCartesianComponents())
        		result.append(df.format(d) + " ");
        	result.append(newline);
        }
        
		return result.toString();
	}

	private static String getLammpsInputFile(StructureOrg c, String potlStr, String units, boolean relax) {
		StringBuilder ans = new StringBuilder();
		String newline = GAUtils.newline();

		//ans.append("units		metal" + newline);
		ans.append("units		" + units + newline);
		ans.append("dimension	3" + newline);
		//ans.append("atom_style	atomic" + newline);
		ans.append("atom_style	charge" + newline);
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
		if (relax)
			ans.append("fix 1 all box/relax tri 1e4 vmax 0.001" + newline);
		ans.append("minimize 0.0 1.0e-8 10000 100000 " + newline);
		ans.append("dump myDump all atom 100000000000000 " + dumpFileName + newline);
		ans.append("dump_modify myDump sort 1 scale no" + newline);
		if (relax)
			ans.append("fix 1 all box/relax tri 0 vmax 0.001" + newline);
		ans.append("#min_modify line quadratic" + newline);
		ans.append("minimize 0.0 1.0e-8 10000 100000 " + newline);
		
		return ans.toString();
	}

	// runs Lammps on the input file given and returns the results in a String
	private static String runLAMMPS(String runDir) {
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
					GAOut.out().stdout(s, GAOut.DEBUG);
				}
	
				// print out any errors
				while ((s = stdError.readLine()) != null) {
					System.out.println(s);
				}
			} finally {
				stdInput.close();
				stdError.close();
			    p.getOutputStream().close();
			    p.getInputStream().close();
			    p.getErrorStream().close();
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
		
		// make temp directory
		String outDirPath = params.getTempDirName() + "/" + params.getRunTitle() + "." + c.getID();
		File outDir = new File(outDirPath);
		outDir.mkdir();
		
		// write unrelaxed cell to disk
		VaspIn.writePoscar(c.getCell(), outDirPath + "/" + c.getID() + ".unrelaxed.POSCAR", false);
		
		utility.Utility.writeStringToFile(getLammpsInputFile(c, potlStr, unitsStr, relax_box), outDirPath + "/" + inFileName);
		utility.Utility.writeStringToFile(getLammpsDataFile(c.getCell()), outDirPath + "/" + dataFileName);
		GAOut.out().stdout("Starting Lammps computation on organism " + c.getID(), GAOut.NOTICE, c.getID());
		
		String lammpsOutput = runLAMMPS(outDir.getAbsolutePath());

		// store the results of the run
		// update c to be the structure in Lammps's output
		Cell a = parseOutputStructure(c.getCell(), outDirPath + "/" + dumpFileName);
		if (a == null) {
			GAOut.out().stdout("Warning: bad Lammps output.  Not updating structure.", GAOut.NOTICE, c.getID());
		} else {
			c.setCell(a);
			finalEnergy = parseFinalEnergy(lammpsOutput);
		}
		
		// write cell to disk
		//Utility.writeStringToFile(c.getCIF(), outDirPath + "/" + c.getID() + ".relaxed.cif");
		
		GAOut.out().stdout("Energy of org " + c.getID() + ": " + finalEnergy + " ", GAOut.NOTICE, c.getID());

		return finalEnergy;
	}
	
	public static Cell parseOutputStructure(Cell origCell, String outFile) {
		// make sure out file was created successfully
		
		if (! (new File(outFile)).exists() || Utility.readStringFromFile(outFile).isEmpty()) {
			GAOut.out().stdout("In LammpsEnergy.parseOutputStructure: dump file doesn't exist or is empty.", GAOut.INFO);
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
		double xlo_bound = Double.parseDouble(xs[0]);
		double xhi_bound = Double.parseDouble(xs[1]);
		double xy = Double.parseDouble(xs[2]);
		String ys[] = lines[6].split(" ");
		double ylo_bound = Double.parseDouble(ys[0]);
		double yhi_bound = Double.parseDouble(ys[1]);
		double xz = Double.parseDouble(ys[2]);
		String zs[] = lines[7].split(" ");
		double zlo_bound = Double.parseDouble(zs[0]);
		double zhi_bound = Double.parseDouble(zs[1]);
		double yz = Double.parseDouble(zs[2]);
		
		// invert according to http://lammps.sandia.gov/doc/Section_howto.html#howto_12
		double xlo = xlo_bound - Math.min(Math.min(0.0, xy), Math.min(xz, xy+xz));
		double xhi = xhi_bound - Math.max(Math.max(0.0, xy), Math.max(xz, xy+xz));
		double ylo = ylo_bound - Math.min(0.0, yz);
		double yhi = yhi_bound - Math.max(0.0, yz);
		double zlo = zlo_bound;
		double zhi = zhi_bound;
		
		basis.add(new Vect(xhi-xlo,0.0,0.0));
		basis.add(new Vect(xy, yhi-ylo, 0.0));
		basis.add(new Vect(xz, yz, zhi-zlo));
		
		for (int i = 0; i < numAtoms; i++) {
			String as[] = lines[9+i].split(" ");
			int type = Integer.parseInt(as[1]) - 1;
			double x = Double.parseDouble(as[2]);
			double y = Double.parseDouble(as[3]);
			double z = Double.parseDouble(as[4]);
			sites.add(new Site(GAParameters.getParams().getCompSpace().getElements().get(type),
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
		
		String lines[] = lammpsOutput.split("\n");
		int i;
		/*
		for (i = 0; i < lines.length; i++) {
			if (lines[i].matches(".*Energy .* final.*")) {
				String energies[] = lines[i+1].split("  *");
				finalEnergy = Double.parseDouble(energies[2]);
			}
		} */
		for (i = 0; i < lines.length; i++) {
			if (lines[i].matches(".*Step *Temp *E_pair *E_mol *TotEng *Press.*")) {
				String energies[] = lines[i+2].trim().split("  *");
				try {
					finalEnergy = Double.parseDouble(energies[4]);
				}	catch (NumberFormatException x) {
					System.out.println("Warning: malformed LAMMPS output when trying to parse energy.");
				}
			//	break;
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
	
	
	public boolean cannotCompute(StructureOrg o) {
		// make sure we don't expect the error:
		// Triclinic box skew is too large
	    //     The displacement in a skewed direction must be less than half the box length in that dimension. 
		//     E.g. the xy tilt must be between -half and +half of the x box length. 
		
		Cell cell = o.getCell();
		cell.rotatedIntoPrincDirs();
		
		double xlo = 0, xhi = cell.getLatticeVectors().get(0).getCartesianComponents().get(0);
		double ylo = 0, yhi = cell.getLatticeVectors().get(1).getCartesianComponents().get(1);
		double zlo = 0, zhi = cell.getLatticeVectors().get(2).getCartesianComponents().get(2);
		double xy = cell.getLatticeVectors().get(1).getCartesianComponents().get(0);
		double xz = cell.getLatticeVectors().get(2).getCartesianComponents().get(0);
		double yz = cell.getLatticeVectors().get(2).getCartesianComponents().get(1);
		
		if ((xhi < Math.abs(2*xy)) || (xhi < Math.abs(2*xz)) || (yhi < Math.abs(2*yz)))
			return true;
		
		return false;
	}
	
	// just for testing:
	public static void main(String[] args) {
		/*
		String[] geArgs = {"/home/wtipton/temp/lammpspot"};
		LammpsEnergy bob = new LammpsEnergy(geArgs);
		StructureOrg c = new StructureOrg(Cell.parseCif(new File("/home/wtipton/POSCAR4.cif")));
//		Utility.writeStringToFile(c.getCIF(), "bob.cif");
		double energy = bob.getEnergy(c);
		System.out.println(energy);
		System.out.println(c.getCell());
		
		*/
	}
}
