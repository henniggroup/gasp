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

import utility.*;
import vasp.VaspIn;

import crystallography.Cell;
import crystallography.Site;

// OhmmsEnergy computes the total energy of a StructureOrg using OHMMS and the given potential.
// It contains all of the methods and utilities that are specific to OHMMS.

public class OhmmsEnergy implements Energy {
	
	private String footerStr;
	private String headerStr;
	private Boolean cautious;
	
	private final String outStructureFileName = "restart.xml";

	public OhmmsEnergy(List<String> args)
	{
		if (args == null || args.size() < 3)
			GAParameters.usage("Not enough parameters given to OhmmsEnergy", true);
		
		// read in the GULP header to use
		File headerFile = new File(args.get(0));
		headerStr = GAUtils.readStringFromFile(headerFile);

		// read in the GULP potential to use
		File footerFile = new File(args.get(1));
		footerStr = GAUtils.readStringFromFile(footerFile);
		
		cautious = Boolean.parseBoolean(args.get(2));
	}

	public String toString() {
		StringBuilder result = new StringBuilder();
		
		result.append("OHMMS total energy:" + GAUtils.newline());
		result.append(headerStr + GAUtils.newline());
		result.append(footerStr + GAUtils.newline());
		
		
		return result.toString();
	}
	
	public double getEnergy(StructureOrg c) {
		return ohmmsRun(c);
	}

	// returns a structure representation in format parse-able by OHMMS
	private static String structureToString(StructureOrg c) {
		//example:
		/*    <particlesets>
 <particleset name="BulkSi" size="8">
        <!-- Description of the cubic unit cell -->
        <unitcell>
          <parameter name="scale" Unit="AA">5.432 </parameter>
          <parameter name="lattice">  1.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 1.0 </parameter>
          <parameter name="bconds"> p p p </parameter>
        </unitcell>
        <!-- creating particles using generic attrib 
        name = the name of a particleattrib in an applicatoin 
        id   = special keyword to determine the type of data
        inunit = 1 in the unit, 0 for cartesian (depends on id)
        -->
        <attrib name="ionid" datatype="stringArray">
          Si Si Si Si Si Si Si Si
        </attrib>
        <attrib name="position" datatype="posArray" condition="1">
          0.0 0.0 0.0 
          0.20 0.25 0.25 
          0.5 0.5 0.0 
          0.75 0.75 0.25 
          0.5 0.0 0.5 
          0.75 0.25 0.75 
          0.0 0.5 0.5 
          0.25 0.75 0.75 
        </attrib>
      </particleset>

    </particlesets>
*/
		
		
		StringBuilder result = new StringBuilder();
		String newline = GAUtils.newline();

		DecimalFormat df = new DecimalFormat();
		df.setMaximumFractionDigits(6);
		// cell parameters
		result.append("<particlesets>" + newline);
		result.append(newline);
		result.append("<particleset name=\"OhmssEnergy\" size=\"" + c.getCell().getNumSites() + "\">" + newline);
		result.append("<unitcell>" + newline);
		result.append("<parameter name=\"scale\" Unit=\"AA\"> 1.0 </parameter>" + newline);
		result.append("<parameter name=\"lattice\"> ");
		for (double d : c.getCell().getLatticeVectorsArray())
			result.append(df.format(d) + " ");
		result.append(" </parameter>");
        result.append("<parameter name=\"bconds\"> p p p </parameter>" + newline);
        result.append("</unitcell>" + newline);

		// atoms
        result.append("<attrib name=\"ionid\" datatype=\"stringArray\">" + newline);
        for (Site s : c.getCell().getSites())
        	result.append(s.getElement().getSymbol() + " ");
        result.append(newline +  "</attrib>" + newline);
        result.append("<attrib name=\"position\" datatype=\"posArray\" condition=\"1\">" + newline);
        for (Site s : c.getCell().getSites()) {
        	for (double d : s.getCoords().getComponentsWRTBasis(c.getCell().getLatticeVectors()))
        	//for (double d : s.getCoords().getCartesianComponents())
        		result.append(df.format(d) + " ");
        	result.append(newline);
        }
        result.append("</attrib>" + newline + "</particleset>" + newline + "</particlesets>");
        
		return result.toString();
	}

	// Returns the name of a OHMMS input file for the crystal c
	// optimize the structure or not.
	private static String ohmmsInputFile(StructureOrg c, String headerStr, String footerStr, String path) {
		GAParameters params = GAParameters.getParams();
		String ans = new String("Maybe file creation failed.");
		String newline = GAUtils.newline();
		
		try {
			// create the OHMMS input file
			File f = new File(path,params.getRunTitle() + "." + c.getID() + ".xml");
			//File f = new File(c.getID() + ".xml");
			// store the filename
			ans = f.getPath();

			// write the file
			BufferedWriter out = new BufferedWriter(new FileWriter(f));

			out.write(headerStr + newline);		
			out.write(newline);
			out.write(structureToString(c));
			out.write(newline);
			out.write(footerStr);
			out.close();
		} catch (IOException e) {
			GAOut.out().stdout("OhmmsEnergy IOException: " + e.getMessage(), GAOut.WARNING, c.getID());
		}

		return ans;
	}

	// runs OHMMS on the input file given and returns the results in a String
	private static String runOHMMS(String runDir) {
		StringBuilder ohmmsOutput = new StringBuilder();

		String s = null;
		try {
			// run the gulp command. GULP will only take input from it's stdin,
			// not from a file. in order to avoid having
			// to open inputFile and feed it to GULP, we use a wrapper script,
			// callgulp, which basically is just gulp <$1 .
			Process p = Runtime.getRuntime().exec("callohmms " + runDir);

			BufferedReader stdInput = new BufferedReader(new InputStreamReader(
					p.getInputStream()));
			BufferedReader stdError = new BufferedReader(new InputStreamReader(
					p.getErrorStream()));

			// read the output
			while ((s = stdInput.readLine()) != null) {
				ohmmsOutput.append(s + GAUtils.newline());
				GAOut.out().stdout(s, GAOut.DEBUG);
			}

			// print out any errors
			while ((s = stdError.readLine()) != null) {
				System.out.println(s);
			}

		} catch (IOException e) {
			System.out.println("IOException in OhmmsEnergy.runOhmms: " + e.getMessage());
			System.exit(-1);
		}

		return ohmmsOutput.toString();
	}

	// Does a OHMMS run on the StructureOrg c. 
	// And we update c to be the local minimum found.
	private double ohmmsRun(StructureOrg c) {
		GAParameters params = GAParameters.getParams();
		double finalEnergy = 0;
		
		// make temp directory
		File outDir = new File(params.getTempDirName() + "/" + params.getRunTitle() + "." + c.getID());
		outDir.mkdir();
		
		ohmmsInputFile(c, headerStr, footerStr, outDir.getAbsolutePath() + "/");
		GAOut.out().stdout("Starting OHMMS computation on organism " + c.getID(), GAOut.NOTICE, c.getID());
		
		String ohmmsOutput = runOHMMS(outDir.getAbsolutePath());

		// store the results of the run
		// update c to be the structure in OHMMS's output
		Cell a = parseOutputStructure(c.getCell(), outDir.getAbsolutePath() + "/" + outStructureFileName);
		if (a == null) {
			GAOut.out().stdout("Warning: bad OHMMS output.  Not updating structure.", GAOut.NOTICE, c.getID());
		} else {
			c.setCell(a);
		}
		
		// write cell to disk
		VaspIn.writePoscar(c.getCell(), outDir.getAbsolutePath() + "/" + c.getID() + ".relaxed.POSCAR", false);

		finalEnergy = parseFinalEnergy(ohmmsOutput);
		
		GAOut.out().stdout("Energy of org " + c.getID() + ": " + finalEnergy + " ", GAOut.NOTICE, c.getID());

		return finalEnergy;
	}
	
	private Cell parseOutputStructure(Cell origCell, String xmlOutputFile) {
		String output = Utility.readStringFromFile(xmlOutputFile);
		List<Vect> basis = new LinkedList<Vect>();
		List<Site> sites = new LinkedList<Site>();
		
		String lines[] = output.split("\n");
		int lineNum = 0;
		for ( ; lineNum < lines.length; lineNum++)
			if (lines[lineNum].matches(".*<parameter name=\"lattice\".*"))
				break;
		
		for (int i = 1; i <= 3; i++) {
			StringTokenizer vectTok = new StringTokenizer(lines[lineNum + i]);
			basis.add(new Vect(Double.parseDouble(vectTok.nextToken())
							, Double.parseDouble(vectTok.nextToken())
							, Double.parseDouble(vectTok.nextToken())));
		}
		
		for ( ; lineNum < lines.length; lineNum++)
			if (lines[lineNum].matches(".*<attrib name=\"position\".*"))
				break;
		
		List<Site> origSites = origCell.getSites();
		for (int i = 0; i < origSites.size(); i++) {
			StringTokenizer siteTok = new StringTokenizer(lines[lineNum + i + 1]);
			sites.add(new Site(origSites.get(i).getElement(),
					   new Vect(Double.parseDouble(siteTok.nextToken())
							, Double.parseDouble(siteTok.nextToken())
							, Double.parseDouble(siteTok.nextToken()), basis)));
		}
		
		return new Cell(basis, sites, origCell.getLabel());
	}
	
	private double parseFinalEnergy(String ohmmsOutput) {
		Double finalEnergy = 0.0;

		// parse the ohmmsOutput to find the final energy 
		// (looking for the 2nd number in the last line)
		String lastLine = null;
		String line = null;
		BufferedReader r = new BufferedReader(new StringReader(ohmmsOutput));

		try {
			while ((line = r.readLine()) != null)
				lastLine = line;
		} catch (IOException x) {
			GAOut.out().stdout("OhmmsEnergy.ohmmsRun: IOException: " + x.getLocalizedMessage(), GAOut.CRITICAL);
		}

		String[] tokens = lastLine.trim().split("[ 	][ 	]*");

		try {
			finalEnergy = Double.parseDouble(tokens[1]);
		} catch(NumberFormatException x) {
			System.out.println("NumberFormatException in OhmmsEnergy.parseFinalEnergy: " + x.getMessage());
			return 0;
		}


		// toss sketchy structures:
		/*if (cautious 
				&& (ohmmsOutput.contains("warning") || ohmmsOutput.contains("failed") || ohmmsOutput.contains("caution")
						|| ohmmsOutput.contains("Maximum number of function calls"))) {
			if (verbosity >= 4)
				System.out.println("OHMMS run failed.");
			return 0;
		}*/
				
		return finalEnergy;
	}
	
	
	public boolean cannotCompute(StructureOrg o) {
		return false;
	}
	
	// just for testing:
	public static void main(String[] args) {
		/*
		String[] geArgs = {"header", "footer", "true"};
		OhmmsEnergy bob = new OhmmsEnergy(geArgs);
		StructureOrg c = new StructureOrg(Cell.parseCif(new File("in.cif")));
		c.getCell().writeCIF("bob.cif");
	//	Utility.writeStringToFile(c.getCell().getCIF(), "bob.cif");
		double energy = bob.getEnergy(c);
		System.out.println(energy);
		System.out.println(c.getCell());
		*/
		
	}
}
