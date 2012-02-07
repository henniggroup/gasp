package ga;

import java.io.BufferedReader;
import java.io.File;
import java.io.FilenameFilter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.HashMap;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;

import chemistry.Element;

import utility.Utility;
import utility.Vect;

import crystallography.Cell;
import crystallography.Site;

public class DFTPPEnergy implements Energy {
	
	private String inputStr;
	private boolean cautious;
	private Map<Element,String> ppsMap;
	
	static final String inFileName = "dft.in";
	static final String optionsFileName = "options.in";
	static final String outFileName = "out";
	static final double bohrPerAngstrom = 1.889725989;
	static final String ionPosFileString = "ionpos";

	public DFTPPEnergy(List<String> args)
	{
		if (args == null || args.size() < 4 || (args.size()) % 2 != 0)
			GAParameters.usage("Wrong parameters given to DFTPPEnergy", true);

		inputStr = GAUtils.readStringFromFile(new File(args.get(0)));
		cautious = Boolean.parseBoolean(args.get(1));

		ppsMap = new HashMap<Element,String>();
		for (int i = 2; i < args.size() - 1; i = i+2) {
			Element e = Element.getElemFromSymbol(args.get(i));
			ppsMap.put(e, args.get(i+1));
		}
		

	}
	
	public double getEnergy(StructureOrg c) {
		return dftppRun(c);
	}
	
	private String getInputFile(StructureOrg c) {
		StringBuilder ans = new StringBuilder();
		String newline = GAUtils.newline();

		/*
		# The supercell vectors for our system:
			# the units are Bohr radii, and the supercell vectors
			# are the columns of the matrix specified.

			lattice  0.                5                 5  \
			         5                 0.                5  \
			         5                 5                 0.

			# We're doing a molecule, so its a sensible to use
			# absolute, cartesian coords

			coords-type cartesian

			# Ionic species needed for Water

			ion-species Hydrogen 1.0 1.00797  H/H.pot none
			ion-species Oxygen   6.0 15.9994  O/O.pot none

			# These lines say where are atoms are located.
			# The "1" flag at the end says
			# the ions are movable.

			ion Oxygen    0.000000  0.000000  0.000000  0
			ion Hydrogen  1.470365  1.173033  0.000000  1
			ion Hydrogen -1.470365  1.173033  0.000000  1
		*/
		
		ans.append("dump End IonicPositions" + GAUtils.newline() + GAUtils.newline());

		// removed in new version of dft++:
		//ans.append("relax-ions" + GAUtils.newline() + GAUtils.newline());
		//ans.append("calculate-forces" + GAUtils.newline() + GAUtils.newline());
		//ans.append("electronic-minimization" + GAUtils.newline() + GAUtils.newline());

		ans.append("latt-scale 1 1 1" + GAUtils.newline());
		ans.append("lattice ");
		for (double d : c.getCell().getLatticeVectorsArray())
			ans.append(d * bohrPerAngstrom + " ");
		ans.append(GAUtils.newline() + GAUtils.newline());
		
		ans.append("coords-type cartesian" + GAUtils.newline() + GAUtils.newline());
		
		for (Element e : c.getCell().getComposition().getElements()) {
			ans.append("ion-species " + ppsMap.get(e) + GAUtils.newline());
			// removed in newer version of dftpp
			// ans.append("ion-relax-param " + e.getName() + " 0 0.5" + GAUtils.newline());
		}
		ans.append(GAUtils.newline());
		
		// fix the position of the first ion
		boolean writingFirstIon = true;
		for (Site s : c.getCell().getSites()) {
			ans.append("ion " + s.getElement().getSymbol() + " ");
			for (Double d : s.getCoords().getCartesianComponents())
				ans.append(d * bohrPerAngstrom + " ");
			if (writingFirstIon) {
				ans.append("0" + GAUtils.newline());
				writingFirstIon = false;
			} else {
				ans.append("1" + GAUtils.newline());
			}
		}
		ans.append(GAUtils.newline());
		
		ans.append("include options.in" + GAUtils.newline());
		
		ans.append(GAUtils.newline());

		return ans.toString();
	}

	private double dftppRun(StructureOrg c) {
		GAParameters params = GAParameters.getParams();
		double finalEnergy = Double.POSITIVE_INFINITY;
		int verbosity = GAParameters.getParams().getVerbosity();
		
		// make temp directory
		String outDirPath = params.getTempDirName() + "/" + params.getRunTitle() + "." + c.getID();
		File outDir = new File(outDirPath);
		outDir.mkdir();
		
		// write unrelaxed cell to disk
		c.getCell().writeCIF(outDirPath + "/" + c.getID() + ".unrelaxed.cif");
		
		utility.Utility.writeStringToFile(getInputFile(c), outDirPath + "/" + inFileName);
		utility.Utility.writeStringToFile(inputStr, outDirPath + "/" + optionsFileName);
		if (verbosity >= 3)
			System.out.println("Starting DFT++ computation on organism " + c.getID());
		
		String dftppOutput = runDFTPP(outDir.getAbsolutePath());

		// store the results of the run
		// update c to be the structure in DFT++ output
		Cell a = parseOutputStructure(c.getCell(), outDirPath);
		if (a == null) {
			if (verbosity >= 3)
				System.out.println("Warning: bad DFT++ output.  Not updating structure.");
		} else {
			c.setCell(a);
			finalEnergy = parseFinalEnergy(Utility.readStringFromFile(outDirPath + "/" + outFileName));
		}
		
		// write cell to disk
		//Utility.writeStringToFile(c.getCIF(), outDirPath + "/" + c.getID() + ".relaxed.cif");
		
		if (verbosity >= 3)
			System.out.println("Energy of org " + c.getID() + ": " + finalEnergy + " ");

		return finalEnergy;
	}
	
	// runs DFT++ on the input file given and returns the results in a String
	private static String runDFTPP(String runDir) {
		int verbosity = GAParameters.getParams().getVerbosity();
		StringBuilder output = new StringBuilder();

		String s = null;
		try {
			// run the gulp command. GULP will only take input from it's stdin,
			// not from a file. in order to avoid having
			// to open inputFile and feed it to GULP, we use a wrapper script,
			// callgulp, which basically is just gulp <$1 .
			Process p = Runtime.getRuntime().exec("calldftpp " + runDir);

			BufferedReader stdInput = new BufferedReader(new InputStreamReader(
					p.getInputStream()));
			BufferedReader stdError = new BufferedReader(new InputStreamReader(
					p.getErrorStream()));

			try {
				// read the output
				while ((s = stdInput.readLine()) != null) {
					output.append(s + GAUtils.newline());
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
			System.out.println("IOException in DFTPPEnergy.runDFTPP: " + e.getMessage());
			System.exit(-1);
		} 

		return output.toString();
	}
	
	public static Cell parseOutputStructure(Cell origCell, String outDir) {
		// make sure out file was created successfully
		
		if (! (new File(outDir)).exists()) {
			if (GAParameters.getParams().getVerbosity() >= 4)
				System.out.println("In DFTPPEnergy: output dir doesn't exist?");
			return null;
		}
				
		List<Site> sites = new LinkedList<Site>();
		
		// get outpos file 
		File dir = new File(outDir);
		File outPosFiles[] = dir.listFiles(new OutPosFileFilter());
		
		if ((outPosFiles == null || outPosFiles.length == 0) && GAParameters.getParams().getVerbosity() >= 2) {
			System.out.println("Error: dftpp output positions not found in directory " + outDir);
			return null;
		}
		if (outPosFiles != null && outPosFiles.length > 1) 
			System.out.println("Warning: dftpp found >1 output positions file in directory " + outDir);
		
		String newPositions = Utility.readStringFromFile(outPosFiles[0].getAbsolutePath());
	
		String lines[] = newPositions.split("\n");
		
		for (int i = 0; i < lines.length; i++) {
			String tokens[] = lines[i].trim().split("  *");
			if (tokens.length != 6 || tokens[0].startsWith("#"))
				continue;
			String elementSym = tokens[1];
			double x = Double.parseDouble(tokens[2]) / bohrPerAngstrom;
			double y = Double.parseDouble(tokens[3]) / bohrPerAngstrom;
			double z = Double.parseDouble(tokens[4]) / bohrPerAngstrom;

			sites.add(new Site(Element.getElemFromSymbol(elementSym), new Vect(x,y,z)));
		} 
		
		return new Cell(origCell.getLatticeVectors(), sites, origCell.getLabel());
	}
	
	private static class OutPosFileFilter implements FilenameFilter {
		public boolean accept(File f, String name) {
			return name.contains(ionPosFileString);
		}
	}
	
	private static double parseFinalEnergy(String output) {
		Double finalEnergy = Double.POSITIVE_INFINITY;
		
		/*
		 * from Kendra:
		 * 
		 *  - apparently dft++ prints out free energy "F " only if its doing a calc w/ fillings
		 *    and prints out Etot otherwise.
		 *  - so we look for "F " first and if it's not found, look for Etot
		 
	       to get energy for non-metallic
	
	       grep Etot *.out|tail -1|awk '{print $3}'
	
	       to get energy for metallic (with fillings)
	
	       grep "F " *.out|tail -1|awk '{print $3}'
		 */
		
		String lines[] = output.split("\n");
		
		/*
		for (int i = 0; i < lines.length; i++) {
			if (lines[i].matches("F .*")) {
				String energies[] = lines[i].trim().split("  *");
				try {
					finalEnergy = Double.parseDouble(energies[4]);
				}	catch (NumberFormatException x) {
					System.out.println("Warning: malformed DFT++ output when trying to parse F.");
				}
			}
		} */

		// IonicMinimize: Iter:   4  Etot: -1.726088327116451e+01  |grad|_K:  4.057e-04  alpha:  2.761e+00, linmin =  6.784e-01, cgtest = -8.019e-01
		if (finalEnergy == Double.POSITIVE_INFINITY) {
			for (int i = 0; i < lines.length; i++) {
				if (lines[i].matches("IonicMinimize: Iter: *[0-9][0-9]* *Etot:.*")) {
					String energies[] = lines[i].trim().split("  *");
					try {
						finalEnergy = Double.parseDouble(energies[4]);
					}	catch (NumberFormatException x) {
						System.out.println("Warning: malformed DFT++ output when trying to parse Etot.");
					}
				}
			} 
		}
				
		return finalEnergy;
	}
	
	// just for testing
	public static void main(String args[]) {
		String output = Utility.readStringFromFile("/home/wtipton/out");
		
		System.out.println(parseFinalEnergy(output));
		
	//	Cell origCell = Cell.parseCif(new File("/home/wtipton/gpu_test/4.unrelaxed.cif"));
	//	System.out.println(parseOutputStructure(origCell, "/home/wtipton/gpu_test/"));
	}
	
}
