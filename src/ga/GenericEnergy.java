/* Genetic algorithm for crystal structure prediction.  */
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

import crystallography.Cell;
import java.util.*;
import vasp.*;

public class GenericEnergy implements Energy {
		
	public GenericEnergy(List<String> args) {
		if (args == null || args.size() > 0 )
			GAParameters.usage("Parameters given to GenericEnergy?", true);
	}
	
	// runs VASP on the input file given and returns the results in a String
	private String runGeneric(String inputDir) {
		StringBuilder output = new StringBuilder();

		String s = null;
		try {
			//  we call a wrapper script which is probably
			// just, e.g.,  "cd $1; vasp" for simple set ups
			Process p = Runtime.getRuntime().exec("callgeneric " + inputDir);

			BufferedReader stdInput = new BufferedReader(new InputStreamReader(
					p.getInputStream()));
			BufferedReader stdError = new BufferedReader(new InputStreamReader(
					p.getErrorStream()));

			try {
				// read the output
				while ((s = stdInput.readLine()) != null) {
					output.append(s + GAUtils.newline());
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
			System.out.println("IOException in GenericEnergy.runGeneric: " + e.getMessage());
			System.exit(-1);
		}

		return output.toString();
	}
	
	private double genericRun(StructureOrg o) {
		GAParameters params = GAParameters.getParams();
		
		// some output
		GAOut.out().stdout("Starting Generic computation on organism "
				+ o.getID() + "... ", GAOut.NOTICE, o.getID());
		
		// make temp directory
		File outDir = new File(params.getTempDirName() + "/" + params.getRunTitle() + "." + o.getID());
		outDir.mkdir();
		
		// make the poscar file
		VaspIn vaspin = new VaspIn(o.getCell(), null, null, null);
		vaspin.makePOSCAR(outDir.getAbsolutePath() + "/");
	
		// run 
		String output = runGeneric(outDir.getAbsolutePath());
		
		GAOut.out().stdout(output, GAOut.DEBUG, o.getID());

		// store the relaxed structure back into o
		Cell newCell = VaspOut.getPOSCAR(outDir.getAbsolutePath() + "/CONTCAR");
		if (newCell != null)
			o.setCell(newCell);
		
		// return the energy
		double finalEnergy = VaspOut.getFinalEnergy(outDir.getAbsolutePath() + "/OUTCAR", false);
		GAOut.out().stdout("Energy of org " + o.getID() + ": " + finalEnergy + " ", GAOut.NOTICE, o.getID());
		
		return finalEnergy; 
	}

	public double getEnergy(StructureOrg o) {
		return genericRun(o);
	}
	
	
	public boolean cannotCompute(StructureOrg o) {
		return false;
	}

	// just for testing:
	public static void main(String[] args) {

		
	}
}
