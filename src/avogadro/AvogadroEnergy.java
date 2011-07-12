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

// NB: Convention: when runs fail, energy should evaluate to Double.POSITIVE_INFINITY

public class AvogadroEnergy implements Energy {
	
	private String headerStr;
	private Boolean cautious;
	private String unrelaxed_cifpath;
	private String logpath;

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
		GAParameters params = GAParameters.getParams();
		
		unrelaxed_cifpath = params.getTempDirName() + "/orig." + c.getID() + ".cif";
		logpath = params.getTempDirName() + "/log." + c.getID();

		return avogadroRun(c);
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
			out.write("log = open('" + logpath + "','w')");
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

			// print out any errors
			while ((s = stdError.readLine()) != null) {
//			if (verbosity >= 5)	
				System.out.println(s);
			}
			
			// read the output
			while ((s = stdInput.readLine()) != null) {
				avogadroOutput.append(s + GAUtils.newline());
				System.out.println("hi: "+s);
			}

		} catch (IOException e) {
			System.out.println("IOException in AvogadroEnergy.runAvogadro: " + e.getMessage());
//			System.exit(-1);
		}
		
	}

	// Does an Avogadro run on the StructureOrg c. If optimize is true,
	// we update c to be the local minimum found by Avogadro.
	public double avogadroRun(StructureOrg c) {
		double finalEnergy = Double.POSITIVE_INFINITY;
		int verbosity = GAParameters.getParams().getVerbosity();

		GAParameters params = GAParameters.getParams();
		
		String inputFile = avogadroInputFile(c);
		Utility.writeStringToFile(c.getCell().getCIF(), unrelaxed_cifpath);
		if (verbosity >= 3)
			System.out.println("Starting Avogadro computation on organism "
					+ c.getID());
		
		runAvogadro(inputFile);

		// update c to be the structure in Avogadro's output
		String cifFileName = inputFile + ".cif";
		Cell a = Cell.parseAvogCif(new File(cifFileName));
//		a.writeCIF(unrelaxed_cifpath);
		if (a == null) {
			if (verbosity >= 3) {
				System.out.println("Warning: bad Avogadro CIF.  Not updating structure.");
			}
		} else {
			c.setCell(a);
		}

		String line = Utility.readStringFromFile(logpath);
	/*	File file = new File(logpath);
		String line = null;
        try {
            Scanner scan = new Scanner(file);
            while (scan.hasNextLine()) {
                line = scan.nextLine();
//              System.out.println(line);
            }
        } catch (FileNotFoundException e) {
            e.printStackTrace();
        }
		
//		System.out.println(line);
*/        
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
		
		return finalEnergy;
	}
}
