/* Genetic algorithm for crystal structure prediction.  Will Tipton.  Ceder Lab at MIT. Summer 2007 */

package avogadro;

import ga.Energy;
import ga.GAParameters;
import ga.GAUtils;
import ga.StructureOrg;

import java.io.*;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.List;
import java.util.Scanner;
import java.util.StringTokenizer;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import utility.Utility;

import crystallography.Cell;
import crystallography.Site;


// AvogadroEnergy computes the total energy of a StructureOrg using Avogadro and the given potential.
// It contains all of the methods and utilities that are specific to Avogadro.

//NB: Convention: when runs fail, energy should evaluate to Double.POSITIVE_INFINITY

public class AvogadroEnergy implements Energy {
	
	private String headerStr;
	private Boolean cautious;

	public AvogadroEnergy(String[] args)
	{
		if (args == null || args.length < 2)
			GAParameters.usage("Not enough parameters given to AvogadroEnergy", true);
		
		// read in the Avogadro header to use
		File headerFile = new File(args[0]);
		headerStr = GAUtils.readStringFromFile(headerFile);
		
		cautious = Boolean.parseBoolean(args[1]);
	}

	public String toString() {
		StringBuilder result = new StringBuilder();
		
		result.append("Avogadro total energy:" + GAUtils.newline());
		result.append(headerStr + GAUtils.newline());		
		
		return result.toString();
	}
	
	public double getEnergy(StructureOrg c) {
		Boolean relax = true;

		return avogadroRun(c, relax);
	}

	// Returns the name of a Avogadro input file for the crystal c, where we can
	// either optimize the structure or not.
	private String avogadroInputFile(StructureOrg c) {
		GAParameters params = GAParameters.getParams();
		int verbosity = params.getVerbosity();
		String ans = new String("Maybe file creation failed.");
		String newline = GAUtils.newline();
		
		try {
			// create the Avogadro input file
			//File f = File.createTempFile(params.getRunTitle() + "." + c.getID() + ".", ".gin",new File(params.getTempDirName()));
			File f = new File(params.getTempDirName(),params.getRunTitle() + "." + c.getID() + ".py");
			// store the filename
			ans = f.getPath();

			// write the file
			BufferedWriter out = new BufferedWriter(new FileWriter(f));
			out.write(headerStr + newline);		
			out.write("conv.ReadFile(mol,\"" + params.getTempDirName() + "/orig." + c.getID() + ".cif\")");
			out.write(newline);
			out.write("log = open('" + params.getTempDirName() + "/log." + c.getID() + "','w')");
			out.write(newline);
			out.write(newline);
			out.write("pFF = openbabel.OBForceField.FindForceField(\"UFF\")");
			out.write(newline);
			out.write("pFF.Setup(mol)");
			out.write(newline);
			out.write("e = pFF.Energy()");
			out.write(newline);
			out.write(newline);
			out.write("pFF.ConjugateGradients(1000)");
			out.write(newline);
			out.write("pFF.GetCoordinates(mol)");
			out.write(newline);
			out.write("e = pFF.Energy()");
			out.write(newline);
			out.write("log.write('Final energy = ' + str(e))");
			out.write(newline);
			out.write(newline);
			out.write("conv.WriteFile(mol,\"" + params.getTempDirName() + "/" + params.getRunTitle() + "." + c.getID() + ".py.cif\")");
			out.close();
		} catch (IOException e) {
			if (verbosity >= 2)
				System.out.println("AvogadroEnergy IOException: " + e.getMessage());
		}

		return ans;
	}

	// runs Avogadro on the input file given and prints out any errors
	public static void runAvogadro(String inputFile) {
		int verbosity = GAParameters.getParams().getVerbosity();
		StringBuilder avogadroOutput = new StringBuilder();

		String s = null;
		try {
			// run the avogadro commands
			Process p = Runtime.getRuntime().exec("python " + inputFile);

			BufferedReader stdInput = new BufferedReader(new InputStreamReader(
					p.getInputStream()));
			BufferedReader stdError = new BufferedReader(new InputStreamReader(
					p.getErrorStream()));

			// read the output
			while ((s = stdInput.readLine()) != null) {
				avogadroOutput.append(s + GAUtils.newline());
//				if (verbosity >= 5)
//					System.out.println(s);
			}

			// print out any errors
			while ((s = stdError.readLine()) != null) {
				System.out.println(s);
			}

		} catch (IOException e) {
			System.out.println("IOException in AvogadroEnergy.runAvogadro: " + e.getMessage());
			System.exit(-1);
		}
	}

	// Does an Avogadro run on the StructureOrg c. If optimize is true,
	// we update c to be the local minimum found by Avogadro.
	public double avogadroRun(StructureOrg c, Boolean relax) {
		double finalEnergy = Double.POSITIVE_INFINITY;
		int verbosity = GAParameters.getParams().getVerbosity();

		GAParameters params = GAParameters.getParams();
		
		String inputFile = avogadroInputFile(c);
		Utility.writeStringToFile(c.getCell().getCIF(), params.getTempDirName() + "/orig." + c.getID() +".cif");
		if (verbosity >= 3)
			System.out.println("Starting Avogadro computation on organism "
					+ c.getID());
		
		runAvogadro(inputFile);

//		System.out.println("this is the org number: " + c.getID());
	// cifs	
		// if we're optimizing the structure, store the results of the run
		if (relax) {
			// update c to be the structure in Avogadro's output
			String cifFileName = inputFile + ".cif";
			Cell a = Cell.parseAvogCif(new File(cifFileName));
			a.writeCIF(params.getTempDirName() + "/orig." + c.getID() + ".cif");
//			a.writeCIF("/home/stewart/testdirectory/" + c.getID() + ".cif");
			if (a == null) {
				if (verbosity >= 3)
					System.out.println("Warning: bad Avogadro CIF.  Not updating structure.");
			} else {
				c.setCell(a);
			}
		}

		File file = new File(params.getTempDirName() + "/log." + c.getID());
		String line = null;
        try {
            Scanner scan = new Scanner(file);
            while (scan.hasNextLine()) {
                line = scan.nextLine();
                System.out.println(line);
            }
        } catch (FileNotFoundException e) {
            e.printStackTrace();
        }
		
//		System.out.println(line);
        
		finalEnergy = parseFinalEnergy(line, cautious);
		
		if (verbosity >= 3)
			System.out.println("Energy of org " + c.getID() + ": " + finalEnergy + " ");

		return finalEnergy;
	}
	
	public static double parseFinalEnergy(String avogadroOutput, boolean cautious) {
		Double finalEnergy = Double.POSITIVE_INFINITY;
		int verbosity = GAParameters.getParams().getVerbosity();
		
		// parse the avogadroOutput to find the final energy
		String line = null;
		Pattern energyPattern = Pattern.compile("Final en[a-z]* = .*");
		Matcher energyMatcher = energyPattern.matcher(avogadroOutput);
		BufferedReader r = new BufferedReader(new StringReader(avogadroOutput));
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
							System.out.println("AvogadroEnergy.parseFinalEnergy: " + x.getMessage());
						if (verbosity >= 5) {
							System.out.println("Avogadro output follows:");
							System.out.println(avogadroOutput);
						}
					}
					break;
				}
			}
//			line = r.readLine();
		} catch (IOException x) {
			if (verbosity >= 1)
				System.out.println("AvogadroEnergy.avogadroRun: IOException.");
		}
		
		if (line == null)
			return Double.POSITIVE_INFINITY;
		
/*		// the final gNorm should be on the next line:
		StringTokenizer t = new StringTokenizer(line);
		t.nextToken();
		t.nextToken();
		t.nextToken();
		double finalGNorm = 0;
		try {
			finalGNorm = Double.parseDouble(t.nextToken());
		} catch (NumberFormatException x) {
			if (verbosity >= 3) 
				System.out.println("AvogadroEnergy.parseFinalEnergy: " + x.getMessage());
			if (cautious)
				return Double.POSITIVE_INFINITY;
		}

		if (cautious 
				&& !(finalGNorm < 0.1 && avogadroOutput.contains("BFGS"))
				&& (avogadroOutput.contains("warning") || avogadroOutput.contains("failed") || avogadroOutput.contains("caution")
						|| avogadroOutput.contains("Maximum number of function calls"))) {
			if (verbosity >= 4)
				System.out.println("WARNING: Being cautious - Avogadro run failed.");
			return Double.POSITIVE_INFINITY;
		}
*/				
		return finalEnergy;
	}
	
/*	// just for testing:
	public static void main(String[] args) {
		String[] geArgs = {"gulppotl_mno2", "true"};
		AvogadroEnergy bob = new AvogadroEnergy(geArgs);
		String output = AvogadroEnergy.runAvogadro("mno2_poly.gin");
		System.out.println(output);
		System.out.println(AvogadroEnergy.parseFinalEnergy(output, bob.cautious));
		System.out.println(output.contains("failed"));
	}
*/	
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