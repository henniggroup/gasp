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

import ga.*;
import chemistry.*;
import java.io.*;

import org.xml.sax.SAXException;

import crystallography.Cell;

import java.util.*;

import utility.Utility;

import javax.xml.parsers.ParserConfigurationException;

// VaspEnergy implements Energy.  It computes the total energy of a 
// StructureOrg using Vasp and the given directory of pseudopotentials.
// It contains all of the methods and utilities in ga that 
// are specific to VASP.

public class VaspEnergy implements Energy {
	
	private String kpointsFile;
	private Map<Element,String> potcarFileMap;
	private String incarFile;
	private boolean cautious;
		
	public VaspEnergy(List<String> args) {
		if (args == null || args.size() < 5 || args.size() % 2 != 1)
			GAParameters.usage("Not enough or malformed parameters given to VaspEnergy", true);
		
		cautious = Boolean.parseBoolean(args.get(0));
		kpointsFile = args.get(1);
		incarFile = args.get(2);
		potcarFileMap = new HashMap<Element,String>();
		for (int i = 3; i < args.size(); i = i+2) {
			Element e = Element.getElemFromSymbol(args.get(i));
			String potcarFile = args.get(i+1);
			potcarFileMap.put(e, potcarFile);
		}

	}
	
	// runs VASP on the input file given and returns the results in a String
	private String runVasp(String inputDir) {
		StringBuilder vaspOutput = new StringBuilder();

		String s = null;
		Process p = null;
		try {
			// run the vasp command. in order to avoid hard-coding things which are
			// probably system-dependent, we call a wrapper script which is probably
			// just "cd $1; vasp" for simple set ups
			p = Runtime.getRuntime().exec("callvasp " + inputDir);

			BufferedReader stdInput = new BufferedReader(new InputStreamReader(
					p.getInputStream()));
			BufferedReader stdError = new BufferedReader(new InputStreamReader(
					p.getErrorStream()));

			try {
				// read the output
				while ((s = stdInput.readLine()) != null) {
					vaspOutput.append(s + GAUtils.newline());
					GAOut.out().stdout(s, GAOut.DEBUG);
				}
	
				// print out any errors
				while ((s = stdError.readLine()) != null) {
					System.out.println(s);
				}
			} finally {
				stdInput.close();
				stdError.close();
		    	  try {
		    		  p.getErrorStream().close();
		    	  } catch (IOException e) {
		    		  e.printStackTrace();
		    	  }
		    	  try {
		    		  p.getOutputStream().close();
		    	  } catch (IOException e) {
		    		  e.printStackTrace();
		    	  }
		    	  try {
		    		  p.getInputStream().close();
		    	  } catch (IOException e) {
		    		  e.printStackTrace();
		    	  }
			}

		} catch (IOException e) {
			System.out.println("IOException in VaspEnergy.runVasp: " + e.getMessage());
			System.exit(-1);
		}

		return vaspOutput.toString();
	}
	
	private double vaspRun(StructureOrg o) {
		GAParameters params = GAParameters.getParams();
		
		// some output
		GAOut.out().stdout("Starting VASP computation on organism "
				+ o.getID() + "... ", GAOut.NOTICE, o.getID());
		
		// make temp directory
		File outDir = new File(params.getTempDirName() + "/" + params.getRunTitle() + "." + o.getID());
		outDir.mkdir();
		
		// make the vasp files
		VaspIn vaspin = new VaspIn(o.getCell(), kpointsFile, incarFile, potcarFileMap);
		vaspin.makeINCAR(outDir.getAbsolutePath() + "/");
		vaspin.makePOSCAR(outDir.getAbsolutePath() + "/");
		vaspin.makeKPOINTS(outDir.getAbsolutePath() + "/");
		vaspin.makePOTCAR(outDir.getAbsolutePath() + "/");
	
		// run vasp
		String vaspOutput = runVasp(outDir.getAbsolutePath());
		
		GAOut.out().stdout(vaspOutput, GAOut.DEBUG, o.getID());

		// store the relaxed structure back into o
		Cell newCell = VaspOut.getPOSCAR(outDir.getAbsolutePath() + "/CONTCAR");
		if (newCell != null)
			o.setCell(newCell);
		
		// return the energy
		double finalEnergy = VaspOut.getFinalEnergy(outDir.getAbsolutePath() + "/OUTCAR", cautious);
		GAOut.out().stdout("Energy of org " + o.getID() + ": " + finalEnergy + " ", GAOut.NOTICE, o.getID());
		
		return finalEnergy; 
	}

	public double getEnergy(StructureOrg o) {
		return vaspRun(o);
	}
	
	
	public boolean cannotCompute(StructureOrg o) {
		return false;
	}

	// just for testing:
	public static void main(String[] args) {
	/*	try {
			VaspXMLReader vaspReader = new VaspXMLReader("vasprun.xml", System.err);
			Cell bob = vaspReader.getFinalStructure();
			vaspReader.getFinalEnergy();
			System.out.println(bob);
		} catch (SAXException x) {
			System.out.println(x.getMessage());
		} catch (ParserConfigurationException x) {
			System.out.println(x.getMessage());
		} catch (IOException x) {
			System.out.println(x.getMessage());
		} */
		
	}
}
