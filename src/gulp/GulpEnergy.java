/* Genetic algorithm for crystal structure prediction.  Will Tipton.  Ceder Lab at MIT. Summer 2007 */

package gulp;

import ga.Energy;
import ga.GAParameters;
import ga.GAUtils;
import ga.StructureOrg;

import java.io.*;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.List;
import java.util.StringTokenizer;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import utility.Utility;

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

	public GulpEnergy(String[] args)
	{
		if (args == null || args.length < 3)
			GAParameters.usage("Not enough parameters given to GulpEnergy", true);
		
		// read in the GULP header to use
		File headerFile = new File(args[0]);
		headerStr = GAUtils.readStringFromFile(headerFile);

		// read in the GULP potential to use
		File potlFile = new File(args[1]);
		potlStr = GAUtils.readStringFromFile(potlFile);
		
		cautious = Boolean.parseBoolean(args[2]);
		
		speciesWithShell = new ArrayList<String>();
		for (int i = 4; i < args.length; i++)
			speciesWithShell.add(args[i]);
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
		double[] lAngles = GAUtils.anglesToDegrees(c.getCellAngles());
		for (int i = 0; i < 3; i++)
			result.append(df.format(lAngles[i]) + " ");
		result.append(newline);

		// atoms
		result.append("cart");
		result.append(newline);
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
		int verbosity = GAParameters.getParams().getVerbosity();
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
				if (verbosity >= 5)
					System.out.println(s);
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
		int verbosity = GAParameters.getParams().getVerbosity();
		
		speciesWithShell = swShell;

		String inputFile = gulpInputFile(c, potlStr);
		if (verbosity >= 3)
			System.out.println("Starting GULP computation on organism "
					+ c.getID());
		String gulpOutput = runGULP(inputFile);

		// if we're optimizing the structure, store the results of the run
		if (relax) {
			// update c to be the structure in GULP's output
			String cifFileName = inputFile + ".cif";
			Cell a = Cell.parseCif(new File(cifFileName));
			if (a == null) {
				if (verbosity >= 3)
					System.out.println("Warning: bad GULP CIF.  Not updating structure.");
			} else {
				c.setCell(a);
			}

		}

		finalEnergy = parseFinalEnergy(gulpOutput, cautious);
		
		if (verbosity >= 3)
			System.out.println("Energy of org " + c.getID() + ": " + finalEnergy + " ");

		return finalEnergy;
	}
	
	public static double parseFinalEnergy(String gulpOutput, boolean cautious) {
		Double finalEnergy = Double.POSITIVE_INFINITY;
		int verbosity = GAParameters.getParams().getVerbosity();
		
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
						if (verbosity >= 3) 
							System.out.println("GulpEnergy.parseFinalEnergy: " + x.getMessage());
						if (verbosity >= 5) {
							System.out.println("GULP output follows:");
							System.out.println(gulpOutput);
						}
					}
					break;
				}
			}
			line = r.readLine();
		} catch (IOException x) {
			if (verbosity >= 1)
				System.out.println("GulpEnergy.gulpRun: IOException.");
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
			if (verbosity >= 3) 
				System.out.println("GulpEnergy.parseFinalEnergy: " + x.getMessage());
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
			if (verbosity >= 4)
				System.out.println("WARNING: Being cautious - GULP run failed.");
			return Double.POSITIVE_INFINITY;
		}
				
		return finalEnergy;
	}
	
	// just for testing:
	public static void main(String[] args) {
		String[] geArgs = {"/home/wtipton/projects/ga_for_crystals/gulp_header", "/home/wtipton/projects/ga_for_crystals/gulppotls/gulppotl_alcu", "true"};
		GulpEnergy bob = new GulpEnergy(geArgs);
		Cell c = Cell.parseCif(new File("/home/wtipton/projects/ga_for_crystals/test_runs/alcu_compspace/refstates/8.cif"));
		
		System.out.println(bob.getEnergy(new StructureOrg(c)));
		
	//	String output = GulpEnergy.runGULP("mno2_poly.gin");
	//	System.out.println(output);
	//	System.out.println(GulpEnergy.parseFinalEnergy(output, bob.cautious));
	//	System.out.println(output.contains("failed"));
	}
}

/*
package gulp;

import crystallography.*;
import java.io.*;
import java.util.*;
import java.util.regex.*;
import optimization.*;

import chemistry.Element;

import utility.*;

public class GulpEnergy {
	
	private Cell origCell;
	private String potential;
	private int timeLimit;
	private Cell optimizedCell;
	private Double optimizedEnergy;
	private OptiSystem sys;
	
	// TODO: configurable?
	private String gulpExe = "gulp";
	private Boolean debug = true;
	private Boolean cautious = true;
	
	public GulpEnergy(Cell c, String pot, int tl, Boolean potStrIsFilename, OptiSystem _sys, boolean _cautious) {
		origCell = c;
		cautious = _cautious;

		if (potStrIsFilename) {
			StringBuilder potBuilder = new StringBuilder();
			// read in potential string
			FileInputStream fis = null;
		    BufferedInputStream bis = null;
		    BufferedReader dis = null;

		    try {
		      fis = new FileInputStream(pot);

		      // Here BufferedInputStream is added for fast reading.
		      bis = new BufferedInputStream(fis);
		      dis = new BufferedReader(new InputStreamReader(bis));

	          String s;
	          while ((s = dis.readLine()) != null)  {
	              potBuilder.append(s + "\n");
	          }
	          potential = potBuilder.toString();

		      // dispose all the resources after using them.
		      fis.close();
		      bis.close();
		      dis.close();

		    } catch (FileNotFoundException e) {
		      e.printStackTrace();
		    } catch (IOException e) {
		      e.printStackTrace();
		    }
		} else {
			potential = pot;
		}
		
		optimizedCell = null;
		optimizedEnergy = Double.NaN;
		sys = _sys;
		timeLimit = tl;
	}
	
	private boolean elementNeedsShell(Element e) {
		//search the String potential for "e.getSymbol  .* shel"
		String energyRegexp = e.getSymbol() + " * shel";
		Pattern energyPattern = Pattern.compile(energyRegexp);
		Matcher energyMatcher = energyPattern.matcher(potential);
		return energyMatcher.find();
	}
	
	private String getGulpInput() {
		StringBuilder answer = new StringBuilder();
		
		// the header
		if (sys.getOptimizeCell() && sys.getOptimizeSites())
			answer.append("opti conp conj \n");
		else if (sys.getOptimizeCell())
			answer.append("opti conp cellonly rfo\n");
		else if (sys.getOptimizeSites())
			answer.append("opti conv\n");
		else
			answer.append("conv\n");
		
		answer.append("time " + timeLimit + "\n");
//		answer.append("switch_minimiser rfo gnorm 0.3 \n\n");
		
		
		
		// lattice
		answer.append("vectors\n");
		List<Vect> lVectors = origCell.getLatticeVectors();
		for (int i = 0; i < 3; i++) {
			Vect v = lVectors.get(i);
			List<Double> components = v.getCartesianComponents();
			// TODO: number formatting necessary here?
			for (int j = 0; j < 3; j++) 
				answer.append(components.get(j) + " ");
			answer.append("\n");
		}
		
		// basis
		answer.append("\ncart\n");
		for (Site s : origCell.getSites()) {
			List<Double> coord = s.getCoords().getCartesianComponents();
			answer.append(s.getElement().getSymbol() + " core " + coord.get(0)
					+ " " + coord.get(1)
					+ " " + coord.get(2) + "\n");
			if (elementNeedsShell(s.getElement()))
				answer.append(s.getElement().getSymbol() + " shel " + coord.get(0)
						+ " " + coord.get(1)
						+ " " + coord.get(2) + "\n");
		}
		
		// potential
		answer.append(potential);
		
		return answer.toString();
	}
	
	// runs gulp on input file and returns output
	private String runGulp(String gulpInput) {
		// thanks to: http://www.devdaily.com/java/edu/pj/pj010016
        StringBuilder output = new StringBuilder();

        try {
        	// run the gulp command
            // using the Runtime exec method:
            Process p = Runtime.getRuntime().exec(gulpExe);
            
            OutputStreamWriter gulpStdIn = new OutputStreamWriter(p.getOutputStream());
            gulpStdIn.write(gulpInput);
            gulpStdIn.close();
            
            BufferedReader stdInput = new BufferedReader(new 
                 InputStreamReader(p.getInputStream()));

            // read the output from the command
            String s;
            while ((s = stdInput.readLine()) != null)  {
                output.append(s + "\n");
                if (debug)
                	System.out.println(s);
            }
            stdInput.close();

            //BufferedReader stdError = new BufferedReader(new 
            //        InputStreamReader(p.getErrorStream()));
            // read any errors from the attempted command
            //while ((s = stdError.readLine()) != null) {
            //    System.out.println(s);
            //}
            
        }
        catch (IOException e) {
            System.out.println("ERROR: runGulp(): " + e.getMessage());
            e.printStackTrace();
        }
        return output.toString();
	}
	
	// return energy
	private double parseOptiRunOutputEnergy(String output) {
		
		// get energy: (kinda a hack.)
		double energy = 0;
		try {
			String energyRegexp = "Final en[a-z]* *= *[-\\.0-9]* eV";
			Pattern energyPattern = Pattern.compile(energyRegexp);
			Matcher energyMatcher = energyPattern.matcher(output);
			energyMatcher.find();
			StringTokenizer energyLineTok = new StringTokenizer(energyMatcher.group(energyMatcher.groupCount()));
			energyLineTok.nextToken();energyLineTok.nextToken();energyLineTok.nextToken();
			energy = Double.parseDouble(energyLineTok.nextToken());
		} catch (IllegalStateException x) {
			System.out.println("ERROR: parseOptiRunOutputEnergy(): " + x.getMessage());
		}
		
		return energy;
	}
	
	private List<Vect> parseOptiRunOutputCell(String output) {
		
		// Get lattice vectors
		List<Vect> latticeVectors = new LinkedList<Vect>();
		BufferedReader lvecsReader = new BufferedReader(new StringReader(output));
		try {
			// read until "Final Cartesian lattice vectors"
			while (! lvecsReader.readLine().contains("Final Cartesian lattice vectors"))
				;
			// read a blank line
			lvecsReader.readLine();
			// read in three lattice vectors
			for (int i = 0; i < 3; i++) {
				StringTokenizer tokens = new StringTokenizer(lvecsReader.readLine());
				double x1 = Double.parseDouble(tokens.nextToken());
				double x2 = Double.parseDouble(tokens.nextToken());
				double x3 = Double.parseDouble(tokens.nextToken());
				latticeVectors.add(new Vect(x1,x2,x3));
			}		
			lvecsReader.close();
		} catch (IOException x) {
			System.out.println("ERROR: parseOptiRunOutputCell (lattice vectors part): " + x.getMessage());
		}
		
		return latticeVectors;
	}
	
	private boolean gulpRunFailed(String output) {
		//strings to search for:
		List<String> errorStrs = new LinkedList<String>();
		errorStrs.add("failed"); errorStrs.add("Failed");
		errorStrs.add("Too many reciprocal lattice vectors needed");
		errorStrs.add("WARNING"); errorStrs.add("Warning");
		errorStrs.add("Error"); errorStrs.add("ERROR"); errorStrs.add("error");
		errorStrs.add("Final energy = \\*");
		
		//search the String potential for "e.getSymbol  .* shel"
		for (String s : errorStrs) {
			if (Utility.stringContains(output, s))
				return true;
		}
		
		// finally, possibly check for the gradient norm issue specifically:
		if (Utility.stringContains(output, "Conditions for a minimum have not been satisfied")) {
			// get gnorm: 
			double maxAcceptableGnorm = 10;
			double gnorm = 0;
			try {
				String energyRegexp = "Final Gnorm *= *[-\\.0-9]*";
				Pattern energyPattern = Pattern.compile(energyRegexp);
				Matcher energyMatcher = energyPattern.matcher(output);
				energyMatcher.find();
				StringTokenizer energyLineTok = new StringTokenizer(energyMatcher.group(energyMatcher.groupCount()));
				energyLineTok.nextToken();energyLineTok.nextToken();energyLineTok.nextToken();
				gnorm = Double.parseDouble(energyLineTok.nextToken());
			} catch (IllegalStateException x) {
				System.out.println("ERROR: gulpRunFailed(): " + x.getMessage());
			} catch (NoSuchElementException x) { //Final Gnorm = # didn't match anything 
				gnorm = 9999999999999.0;
			}
			if (gnorm > maxAcceptableGnorm)
				return true;
		}
		
		return false;
	}
	
	private List<Site> parseOptiRunOutputSites(String output, List<Vect> latticeVectors) {
		// Get sites
		List<Site> sites = new LinkedList<Site>();
		BufferedReader sitesReader = new BufferedReader(new StringReader(output));
		try {
			// read lines until we see "Final fractional coordinates of atoms"
			while (! sitesReader.readLine().contains("Final fractional coordinates of atoms"))
				;
			// read 5 informationless lines
			sitesReader.readLine();sitesReader.readLine();sitesReader.readLine();sitesReader.readLine();sitesReader.readLine();
			// lines should match the pattern:
			//  <int> <element symbol> <char> <xfrac> <yfrac> <zfrac> <double>
			for (int i = 0; i < origCell.getBasisSize(); i++) {
				StringTokenizer tokens = new StringTokenizer(sitesReader.readLine());
				tokens.nextToken();
				String symbol = tokens.nextToken();
				Element elem = Element.getElemFromSymbol(symbol);
				tokens.nextToken();
				double x1 = Double.parseDouble(tokens.nextToken());
				double x2 = Double.parseDouble(tokens.nextToken());
				double x3 = Double.parseDouble(tokens.nextToken());
				sites.add(new Site(elem, new Vect(x1,x2,x3,latticeVectors)));
			}
			sitesReader.close();
		} catch (IOException x) {
			System.out.println("ERROR: parseOptiRunOutputCell (sites part): " + x.getMessage());
		}

		return sites;
	}
	
	private void doRun() {
		// get the input file
		String input = getGulpInput();
		
		if (debug)
			System.out.println(input);
		
	    // write some output files maybe
	    if (sys.getWriteTempFiles()) 
			Utility.writeStringToFile(input ,sys.getOutDir() + "/" + origCell.getLabel() + ".gin");
			

		// run GULP
	    String gulpOutput = runGulp(input);
	    
	    // write some output files maybe
	    if (sys.getWriteTempFiles()) 
			Utility.writeStringToFile(gulpOutput ,sys.getOutDir() + "/" + origCell.getLabel() + ".gout");
		
		// parse out final energy and optimizedCell
	    if (cautious && gulpRunFailed(gulpOutput)) {
	    	optimizedEnergy = 123456.0 * origCell.getBasisSize();
	    	optimizedCell = origCell;
	    } else {
		    optimizedEnergy = parseOptiRunOutputEnergy(gulpOutput);
		    
		    List<Vect> latticeVectors = null;
		    if (sys.getOptimizeCell())
		    	latticeVectors = parseOptiRunOutputCell(gulpOutput);
		    else
		    	latticeVectors = origCell.getLatticeVectors();
		    
		    List<Site> sites = null;
		    if (sys.getOptimizeSites())
		    	sites = parseOptiRunOutputSites(gulpOutput, latticeVectors);
		    else
		    	sites = origCell.getSites();
		    
		    optimizedCell = new Cell(latticeVectors, sites, origCell.getLabel());
	    }
	}
	
	public double getOptimizedEnergy() {
		if (optimizedEnergy.isNaN())
			doRun();
	    
	    if (debug)
	    	System.out.println("Energy: " + optimizedEnergy);
		
		return optimizedEnergy;
	}
	
	public Cell getOptimizedCell() {
		if (optimizedCell == null)
			doRun();
	    
	    if (debug)
	    	System.out.println("Cell: " + optimizedCell);
		
		return optimizedCell;
	}

}

*/