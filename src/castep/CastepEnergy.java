/* Genetic algorithm for crystal structure prediction.  */

package castep;

import ga.*;
import chemistry.*;
import java.io.*;
import org.xml.sax.SAXException;

import crystallography.Cell;

import java.util.*;

import utility.Utility;

import javax.xml.parsers.ParserConfigurationException;

// CastepEnergy implements Energy.  It computes the total energy of a 
// StructureOrg using Castep and the given directory of pseudopotentials.
// It contains all of the methods and utilities in ga that 
// are specific to CASTEP.

public class CastepEnergy implements Energy {
	
	private double kpointSpacing;
	private double pressure;
	private Map<Element,String> ppFileMap;
	private String paramFile;
	private boolean cautious;
	
	public static String castepPrefix = "in";
		
	public CastepEnergy(List<String> args) {
		if (args == null || args.size() < 6 || args.size() % 2 != 0)
			GAParameters.usage("Not enough or malformed parameters given to CastepEnergy", true);
		
		cautious = Boolean.parseBoolean(args.get(0));
		kpointSpacing = Double.parseDouble(args.get(1));
		pressure = Double.parseDouble(args.get(2));
		paramFile = args.get(3);
		ppFileMap = new HashMap<Element,String>();
		for (int i = 4; i < args.size(); i = i+2) {
			Element e = Element.getElemFromSymbol(args.get(i));
			String potcarFile = args.get(i+1);
			ppFileMap.put(e, potcarFile);
		}

	}
	
	// runs CASTEP on the input file given and returns the results in a String
	private String runCastep(String inputDir) {
		StringBuilder castepOutput = new StringBuilder();

		String s = null;
		try {
			// run the castep command. in order to avoid hardcoding thing which are
			// probably system-dependent, we call a wrapper script which is probably
			// just "cd $1; castep"
			Process p = Runtime.getRuntime().exec("callcastep " + inputDir);

			BufferedReader stdInput = new BufferedReader(new InputStreamReader(
					p.getInputStream()));
			BufferedReader stdError = new BufferedReader(new InputStreamReader(
					p.getErrorStream()));

			try {
				// read the output
				while ((s = stdInput.readLine()) != null) {
					castepOutput.append(s + GAUtils.newline());
					GAOut.out().stdout(s, GAOut.DEBUG);
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
			System.out.println("IOException in CastepEnergy.runCastep: " + e.getMessage());
			System.exit(-1);
		}

		return castepOutput.toString();
	}
	
	private double castepRun(StructureOrg o) {
		GAParameters params = GAParameters.getParams();
		
		// some output
		GAOut.out().stdout("Starting CASTEP computation on organism "
				+ o.getID() + "... ", GAOut.NOTICE, o.getID());
		
		// make temp directory
		File outDir = new File(params.getTempDirName() + "/" + params.getRunTitle() + "." + o.getID());
		outDir.mkdir();
		
		// make the castep files
		CastepIn castepin = new CastepIn(o.getCell(), kpointSpacing, pressure, paramFile, ppFileMap);

		castepin.makeParam(outDir.getAbsolutePath() + "/");
		castepin.makeCell(outDir.getAbsolutePath() + "/");

		// run castep
		String castepOutput = runCastep(outDir.getAbsolutePath());
		
		GAOut.out().stdout(castepOutput, GAOut.DEBUG, o.getID());

		// store the relaxed structure back into o
		Cell newCell = CastepOut.getCell(outDir.getAbsolutePath() + "/" + castepPrefix + ".castep");
		if (newCell != null)
			o.setCell(newCell);
		
		// return the energy
		double finalEnergy = CastepOut.getFinalEnergy(outDir.getAbsolutePath() + "/" + castepPrefix + ".castep", cautious);
		
		GAOut.out().stdout("Energy of org " + o.getID() + ": " + finalEnergy + " ", GAOut.NOTICE, o.getID());
		
		return finalEnergy; 
	}

	public double getEnergy(StructureOrg o) {
		return castepRun(o);
	}
	

	// just for testing:
	public static void main(String[] args) {
	/*	try {
			CastepXMLReader castepReader = new CastepXMLReader("casteprun.xml", System.err);
			Cell bob = castepReader.getFinalStructure();
			castepReader.getFinalEnergy();
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
