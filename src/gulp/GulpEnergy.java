/* Genetic algorithm for crystal structure prediction.  Will Tipton.  Ceder Lab at MIT. Summer 2007 */

package gulp;

import ga.Energy;
import ga.GAOut;
import ga.GAParameters;
import ga.GAUtils;
import ga.StructureOrg;
import ga.UnitsSOCreator;

import java.io.*;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.List;
import java.util.StringTokenizer;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import utility.Utility;
import utility.Vect;

import crystallography.Cell;
import crystallography.Site;


// GulpEnergy computes the total energy of a StructureOrg using GULP and the given potential.
// It contains all of the methods and utilities that are specific to GULP.

//NB: Convention: when runs fail, energy should evaluate to Double.POSITIVE_INFINITY

public class GulpEnergy implements Energy {
	
	private String potlStr;
	private String headerStr;
	private ArrayList<String> speciesWithShell;
	private Boolean cautious;

	public GulpEnergy(List<String> args)
	{
		if (args == null || args.size() < 3)
			GAParameters.usage("Not enough parameters given to GulpEnergy", true);
		
		// read in the GULP header to use
		File headerFile = new File(args.get(0));
		headerStr = GAUtils.readStringFromFile(headerFile);

		// read in the GULP potential to use
		File potlFile = new File(args.get(1));
		potlStr = GAUtils.readStringFromFile(potlFile);
		
		cautious = Boolean.parseBoolean(args.get(2));
		
		speciesWithShell = new ArrayList<String>();
		for (int i = 4; i < args.size(); i++)
			speciesWithShell.add(args.get(i));
	}

	public String toString() {
		StringBuilder result = new StringBuilder();
		
		result.append("GULP total energy:" + GAUtils.newline());
		result.append(headerStr + GAUtils.newline());
		result.append(potlStr + GAUtils.newline());		
		
		return result.toString();
	}
	
	public double getEnergy(StructureOrg c) {
		Boolean relax = true;
		return gulpRun(c, potlStr, relax, speciesWithShell);
	}

	// returns a structure representation in format parse-able by GULP
	public static String structureToString(Cell c) {
		return structureToString(c, new ArrayList<String>());
	}
	public static String structureToString(Cell c, List<String> speciesWithShell) {
		StringBuilder result = new StringBuilder();
		String newline = GAUtils.newline();
		
		DecimalFormat df = new DecimalFormat();
		df.setMaximumFractionDigits(6);
		// cell parameters
		result.append("cell");
		result.append(newline);
		// lengths
		double[] lLengths = c.getCellLengths();
		for (int i = 0; i < 3; i++)
			result.append(df.format(lLengths[i]) + " ");
		// angles in degrees
		double[] lAngles = c.getCellAnglesDegrees();
		for (int i = 0; i < 3; i++)
			result.append(df.format(lAngles[i]) + " ");
		result.append(newline);

		result.append("cart");
		result.append(newline);
		// atoms
		// if UFF is being used (assuming potl file will include "uff" in file name or path)
		// TODO: this is not a good assumption --Will
	/*	if (potentialName.toLowerCase().contains("uff")) {
			String[] newLabels = getGulpFormat(c);

			for (int i = 0; i < c.getNumSites(); i++) {
				Site s = c.getSite(i);
				String symbol = newLabels[i];
				List<Double> coords = s.getCoords().getCartesianComponents();
				result.append(symbol + " ");
				result.append(coords.get(0) + " " + coords.get(1) + " " + coords.get(2));
				result.append(newline);
			}
		}
		else { */
			for (int i = 0; i < c.getNumSites(); i++) {
				Site s = c.getSite(i);
				String symbol = s.getElement().getSymbol();
				List<Double> coords = s.getCoords().getCartesianComponents();
				// the core
				result.append(symbol + " core ");
				result.append(coords.get(0) + " " + coords.get(1) + " " + coords.get(2));
				result.append(newline);
				// the shell
				if (GAUtils.listHasString(speciesWithShell,symbol)) {
					result.append(symbol + " shell ");
					result.append(coords.get(0) + " " + coords.get(1) + " " + coords.get(2));
					result.append(newline);
				}
			}
	//	}

		return result.toString();
	}

	// Returns the name of a GULP input file for the crystal c, where we can
	// either
	// optimize the structure or not.
	private String gulpInputFile(StructureOrg c, String potlStr) {
		GAParameters params = GAParameters.getParams();
		String ans = new String("Maybe file creation failed.");
		String newline = GAUtils.newline();
		
		StringBuilder out = new StringBuilder();
		
		// create the GULP input file
		//File f = File.createTempFile(params.getRunTitle() + "." + c.getID() + ".", ".gin",new File(params.getTempDirName()));
		File f = new File(params.getTempDirName(),params.getRunTitle() + "." + c.getID() + ".gin");
		// store the filename
		ans = f.getPath();

		// only output a cif if we're optimizing the structure
		out.append(headerStr + newline);		
		out.append(newline);
		out.append("output cif " + ans + ".cif");
		out.append(newline);
		out.append(structureToString(c.getCell()));
		out.append(newline);
		out.append(potlStr);

		// write the file
		Utility.writeStringToFile(out.toString(), f.getPath());

		return ans;
	}

	// runs GULP on the input file given and returns the results in a String
	public static String runGULP(String inputFile) {
		StringBuilder gulpOutput = new StringBuilder();

		String s = null;
		BufferedReader stdInput = null;
		BufferedReader stdError = null;
		try {
			// run the gulp command. GULP will only take input from it's stdin,
			// not from a file. in order to avoid having
			// to open inputFile and feed it to GULP, we use a wrapper script,
			// callgulp, which basically is just gulp <$1 .
			Process p = Runtime.getRuntime().exec("callgulp " + inputFile);

			stdInput = new BufferedReader(new InputStreamReader(
					p.getInputStream()));
			stdError = new BufferedReader(new InputStreamReader(
					p.getErrorStream()));

			// read the output
			while ((s = stdInput.readLine()) != null) {
				gulpOutput.append(s + GAUtils.newline());
				GAOut.out().stdout(s, GAOut.DEBUG);
			}

			// print out any errors
			while ((s = stdError.readLine()) != null) {
				System.out.println(s);
			}

		} catch (IOException e) {
			System.out.println("IOException in GulpEnergy.runGulp: " + e.getMessage());
			System.exit(-1);
		} finally {
			if (stdInput != null) 
				try{ stdInput.close(); } catch (Exception x) { } //ignore
			if (stdError != null) 
				try{ stdError.close(); } catch (Exception x) { } //ignore
		}

		return gulpOutput.toString();
	}

	// Does a Lennard-Jones GULP run on the StructureOrg c. If optimize is true,
	// we update c to be the local minimum found by GULP.
	public double gulpRun(StructureOrg c, String potlStr,
			Boolean relax, ArrayList<String> swShell) {
		double finalEnergy = Double.POSITIVE_INFINITY;
		
		speciesWithShell = swShell;

		String inputFile = gulpInputFile(c, potlStr);
		GAOut.out().stdout("Starting GULP computation on organism " + c.getID(), GAOut.NOTICE, c.getID());
	
		String gulpOutput = runGULP(inputFile);

		// if we're optimizing the structure, store the results of the run
		if (relax) {
			// update c to be the structure in GULP's output
			String cifFileName = inputFile + ".cif";
			Cell a = Cell.parseCif(new File(cifFileName));
			if (a == null) {
				GAOut.out().stdout("Warning: bad GULP CIF.  Not updating structure.", GAOut.NOTICE, c.getID());
			} else {
				c.setCell(a);
			}

		}

		finalEnergy = parseFinalEnergy(gulpOutput, cautious);
		
		GAOut.out().stdout("Energy of org " + c.getID() + ": " + finalEnergy + " ", GAOut.NOTICE, c.getID());

		return finalEnergy;
	}
	
	public static double parseFinalEnergy(String gulpOutput, boolean cautious) {
		Double finalEnergy = Double.POSITIVE_INFINITY;
		
		// parse the gulpOutput to find the final energy
		String line = null;
		Pattern energyPattern = Pattern.compile("Final en[a-z]* = .*");
		Matcher energyMatcher = energyPattern.matcher(gulpOutput);
		BufferedReader r = new BufferedReader(new StringReader(gulpOutput));
		try {
			while ((line = r.readLine()) != null) {
				energyMatcher.reset(line);
				// when we find the final energy line, store the final energy
				// and exit the loop
				if (energyMatcher.find()) {
					// the energy should be the 4th string on the line which
					// looks like, e.g.:
					// Final energy = -9571.00691982 eV
					StringTokenizer t = new StringTokenizer(energyMatcher.group());
					t.nextToken();
					t.nextToken();
					t.nextToken();
					try {
						finalEnergy = Double.parseDouble(t.nextToken());
					} catch (NumberFormatException x) {
						GAOut.out().stdout("GulpEnergy.parseFinalEnergy: " + x.getMessage(), GAOut.NOTICE);
						GAOut.out().stdout("GULP output follows:", GAOut.DEBUG);
						GAOut.out().stdout(gulpOutput, GAOut.DEBUG);
					}
					break;
				}
			}
			line = r.readLine();
		} catch (IOException x) {
			GAOut.out().stdout("GulpEnergy.gulpRun: IOException: " + x.getLocalizedMessage(), GAOut.CRITICAL);
		}
		
		if (line == null)
			return Double.POSITIVE_INFINITY;

		// the final gNorm should be on the next line:
		StringTokenizer t = new StringTokenizer(line);
		t.nextToken();
		t.nextToken();
		t.nextToken();
		double finalGNorm = 0;
		try {
			finalGNorm = Double.parseDouble(t.nextToken());
		} catch (NumberFormatException x) {
			GAOut.out().stdout("GulpEnergy.parseFinalEnergy: " + x.getMessage(), GAOut.NOTICE);
			if (cautious)
				return Double.POSITIVE_INFINITY;
		}
		// toss sketchy structures:
		//  - the GNorm check only lets us off the hook for newt-raphs optimizer.  
		//	  however, if we have a gnorm < 0.1, then we've hopefully switched to
		//	  the newt-raph optimizer.  The number 0.1 comes from
		//    a suggestion in GULP's output.
		if (cautious 
				&& !(finalGNorm < 0.1 && gulpOutput.contains("BFGS"))
				&& (gulpOutput.contains("warning") || gulpOutput.contains("failed") || gulpOutput.contains("caution")
						|| gulpOutput.contains("Maximum number of function calls"))) {
			GAOut.out().stdout("WARNING: Being cautious - GULP run failed.", GAOut.INFO);
			return Double.POSITIVE_INFINITY;
		}
				
		return finalEnergy;
	}
	
	//TODO: currently ignores square planar -- how to identify? also when # bonds > 6 just leaves blank -- change to octahedral?
	//TODO: doesn't work when both given range and running in parallel, due to UnitsSOCreator.getTargets() being overwritten... how to fix?
	public static String[] getUFFSiteStrings(Cell c) {
		List<Site> sites = c.getSites();
		String[] newLabels = new String[sites.size()];
		
		int difUnits = UnitsSOCreator.getDifUnits();
		int[] numUnits = UnitsSOCreator.getTargets();
		int[] numAtoms = UnitsSOCreator.getNumAtoms();
		
		int loc = 0;
		for (int r=0; r<difUnits; r++) {
			for (int s=0; s<numUnits[r]; s++) {
				// TODO: this makes some kind of sketchy assumptions about the order of the sites list --Will
				List<Site> subsites = sites.subList(loc, loc + numAtoms[r]);

				for (Site k: subsites) {
					Vect v1 = k.getCoords();
					int counter = 0;
					for (Site l: subsites) {
						if (k != l) {
							// TODO: this ignores periodic boundary conditions
							//		 consider using Cell.getAtomsInSphereSorted() 
							//		 and then crossreferencing those results w/ those known to be in the same molecular unit --Will
							Vect v2 = l.getCoords();
							if (v1.getCartDistanceTo(v2) <= 1.85) 	// TODO: this magic number needs to be configurable or documented
								counter++;		//		 or something --Will
						}
					}
					if (counter == 0 || counter > 6) {
						if (k.getElement().getSymbol().length() == 2)
							newLabels[loc] = k.getElement().getSymbol();
						else
							newLabels[loc] = k.getElement().getSymbol() + "_";
					}
					else if (counter == 1 || counter == 2) {
						newLabels[loc] = k.getElement().getSymbol() + "_1";
					}
					else if (counter == 3) {
						newLabels[loc] = k.getElement().getSymbol() + "_2";
					}
					else if (counter == 4) {
						newLabels[loc] = k.getElement().getSymbol() + "_3";
					}
					else if (counter == 5) {
						newLabels[loc] = k.getElement().getSymbol() + "_5";
					}
					else if (counter == 6) {
						newLabels[loc] = k.getElement().getSymbol() + "_6";
					}
					if (k.getElement().getSymbol().equals("H")) {
						newLabels[loc] = k.getElement().getSymbol() + "_";
					}
					loc = loc + 1;
				}
			}
		}
		return newLabels;
		
				
				
				
				
/*				
			for (Site i: sites) {
			Vect v1 = i.getCoords();
			int counter = 0;
			for (Site j: sites) {
				double dist = 0.0;
				if (i != j) {
					Vect v2 = j.getCoords();
					dist = v1.getCartDistanceTo(v2);
					if (dist <= 1.85) {
						counter = counter + 1;
					}
				}
			}
			if (counter == 0 || counter > 6) {
				if (i.getElement().getSymbol().length() == 2)
					newLabels[loc] = i.getElement().getSymbol();
				else
					newLabels[loc] = i.getElement().getSymbol() + "_";
			}
			else if (counter == 1 || counter == 2) {
				newLabels[loc] = i.getElement().getSymbol() + "_1";
			}
			else if (counter == 3) {
				newLabels[loc] = i.getElement().getSymbol() + "_2";
			}
			else if (counter == 4) {
				newLabels[loc] = i.getElement().getSymbol() + "_3";
			}
			else if (counter == 5) {
				newLabels[loc] = i.getElement().getSymbol() + "_5";
			}
			else if (counter == 6) {
				newLabels[loc] = i.getElement().getSymbol() + "_6";
			}
			if (i.getElement().getSymbol().equals("H")) {
				newLabels[loc] = i.getElement().getSymbol() + "_";
			}
			loc = loc + 1;
		}
		return newLabels;
*/		
	}
	
	// just for testing:
	public static void main(String[] args) {
		/*
		String[] geArgs = {"/home/wtipton/projects/ga_for_crystals/gulp_header", "/home/wtipton/projects/ga_for_crystals/gulppotls/gulppotl_alcu", "true"};
		GulpEnergy bob = new GulpEnergy(geArgs);
		Cell c = Cell.parseCif(new File("/home/wtipton/projects/ga_for_crystals/test_runs/alcu_compspace/refstates/8.cif"));
		
		System.out.println(bob.getEnergy(new StructureOrg(c)));
		*/
	//	String output = GulpEnergy.runGULP("mno2_poly.gin");
	//	System.out.println(output);
	//	System.out.println(GulpEnergy.parseFinalEnergy(output, bob.cautious));
	//	System.out.println(output.contains("failed"));
	}
}