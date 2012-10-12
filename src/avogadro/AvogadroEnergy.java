/* Genetic algorithm for crystal structure prediction.  Will Tipton.  Ceder Lab at MIT. Summer 2007 */

package avogadro;

import ga.Energy;
import ga.GAOut;
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

import org.openbabel.*;


// AvogadroEnergy computes the total energy of a StructureOrg using Avogadro and the given potential.
// It contains all of the methods and utilities that are specific to Avogadro.

// NB: Convention: when runs fail, energy should evaluate to Double.POSITIVE_INFINITY

public class AvogadroEnergy implements Energy {
	
	private String headerStr;
	private static String unrelaxed_cifpath;
	private String logpath;
	private static String errorpath;

	public AvogadroEnergy(List<String> args) {		
		
		if (args == null || args.size() < 1)
			GAParameters.usage("Not enough parameters given to AvogadroEnergy", true);
		
		// read in the Avogadro header to use
		File headerFile = new File(args.get(0));
		headerStr = GAUtils.readStringFromFile(headerFile);		

	}

	public String toString() {
		StringBuilder result = new StringBuilder();
		
		result.append("Avogadro total energy:" + GAUtils.newline());
		result.append(headerStr + GAUtils.newline());		
		
		return result.toString();
	}
	
	public double getEnergy(StructureOrg c) {
		GAParameters params = GAParameters.getParams();
		
		// creates file paths for cleaner coding
		unrelaxed_cifpath = params.getTempDirName() + "/orig." + c.getID() + ".cif";
		logpath = params.getTempDirName() + "/log." + c.getID();
		errorpath = params.getTempDirName() + "/errors/error." + c.getID() + ".cif";
		// creates directory for "bad" CIFs
		File errorDir = new File(params.getTempDirName() + "/errors");
		errorDir.mkdir();

		return avogadroRun(c);
	}

	// Returns the name of a Avogadro input file for the crystal c, where we can
	// either optimize the structure or not.
	private String avogadroInputFile(StructureOrg c) {
		GAParameters params = GAParameters.getParams();
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
			out.write("log = open('" + logpath + "','w')");
			out.write(newline);
			out.write(newline);
			out.write("uff = openbabel.OBForceField.FindForceField(\"UFF\")");
			out.write(newline);
			out.write("uff.Setup(mol)");
			out.write(newline);
			out.write("e = uff.Energy()");
			out.write(newline);
			out.write(newline);
			out.write("uff.ConjugateGradients(1000)");
			out.write(newline);
			out.write("uff.GetCoordinates(mol)");
			out.write(newline);
			out.write("e = uff.Energy()");
			out.write(newline);
			out.write("log.write('Final energy = ' + str(e))");
			out.write(newline);
			out.write(newline);
			out.write("conv.WriteFile(mol,\"" + params.getTempDirName() + "/" + params.getRunTitle() + "." + c.getID() + ".py.cif\")");
			out.close();
		} catch (IOException e) {
			GAOut.out().stdout("AvogadroEnergy IOException: " + e.getMessage(), GAOut.WARNING, c.getID());
		}

		return ans;
	}

	// runs Avogadro on the input file given and prints out any errors
	public static Boolean runAvogadro(String inputFile) {		
		GAParameters params = GAParameters.getParams();
		
		Boolean status = true;
		
		StringBuilder avogadroOutput = new StringBuilder();

		String s = null;
		
		try {
			// run the avogadro commands
			Process p = Runtime.getRuntime().exec("python " + inputFile);

			BufferedReader stdInput = new BufferedReader(new InputStreamReader(
					p.getInputStream()));
			BufferedReader stdError = new BufferedReader(new InputStreamReader(
					p.getErrorStream()));

			// print out any errors
			// if error, kill process, copy bad CIF to errors directory
			while ((s = stdError.readLine()) != null) {
//			if (verbosity >= 5)	
				if (s.contains("glibc")) {
					System.out.println("python crashed. skipping calculation");
					Process e = Runtime.getRuntime().exec("cp " + unrelaxed_cifpath + " " + errorpath);
					p.destroy();
					status = false;
					continue;
				}
				System.out.println(s);
			}
			//TODO: change this once we find a fix
			// read the output
			while ((s = stdInput.readLine()) != null) {
				avogadroOutput.append(s + GAUtils.newline());
				System.out.println("hi: "+s);
			}

		} catch (IOException e) {
			System.out.println("IOException in AvogadroEnergy.runAvogadro: " + e.getMessage());
//			System.exit(-1);
		}
		
		return status;
		
	}

	// Does an Avogadro run on the StructureOrg c. If optimize is true,
	// we update c to be the local minimum found by Avogadro.
	public double avogadroRun(StructureOrg c) {
		double finalEnergy = Double.POSITIVE_INFINITY;
		
		String inputFile = avogadroInputFile(c);
		Utility.writeStringToFile(c.getCell().getCIF(), unrelaxed_cifpath);
		
		GAOut.out().stdout("Starting Avogadro computation on organism " + c.getID(), GAOut.NOTICE, c.getID());

		// Execute the python script
		if (!runAvogadro(inputFile))
			return Double.POSITIVE_INFINITY;

		// update c to be the structure in Avogadro's output
		String cifFileName = inputFile + ".cif";
		Cell a = Cell.parseAvogCif(new File(cifFileName));
		if (a == null) {
			GAOut.out().stdout("Bad Avogadro CIF.  Not updating structure.", GAOut.NOTICE, c.getID());
		} else {
			c.setCell(a);
		}

		// reads string from log
		String line = Utility.readStringFromFile(logpath);

		// parse energy from log
		finalEnergy = parseFinalEnergy(line);
		
		GAOut.out().stdout("Energy of org " + c.getID() + ": " + finalEnergy + " ", GAOut.NOTICE, c.getID());
	
		return finalEnergy;			
	}
	
	public static double parseFinalEnergy(String avogadroOutput) {
		Double finalEnergy = Double.POSITIVE_INFINITY;
		
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
					// the energy should be the only line, looks like:
					// Final energy = -9571.00691982
					StringTokenizer t = new StringTokenizer(energyMatcher.group());
					t.nextToken();
					t.nextToken();
					t.nextToken();
					try {
						finalEnergy = Double.parseDouble(t.nextToken());
					} catch (NumberFormatException x) {
						GAOut.out().stdout("AvogadroEnergy.parseFinalEnergy: " + x.getMessage(), GAOut.NOTICE);
						GAOut.out().stdout("Avogadro output follows: \n" + avogadroOutput, GAOut.DEBUG);
					}
					break;
				}
			}
//			line = r.readLine();
		} catch (IOException x) {
			GAOut.out().stdout("AvogadroEnergy.avogadroRun: IOException.", GAOut.CRITICAL);
		}
		
		if (line == null)
			return Double.POSITIVE_INFINITY;
		
		return finalEnergy;
	}
	
	public boolean cannotCompute(StructureOrg o) {
		return false;
	}

}