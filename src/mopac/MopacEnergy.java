package mopac;

import ga.Energy;
import ga.GAParameters;
import ga.GAUtils;
import ga.StructureOrg;
import ga.UnitsSOCreator;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.LinkedList;
import java.util.List;
import java.util.StringTokenizer;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import utility.Constants;
import utility.Utility;
import utility.Vect;
import crystallography.Cell;
import crystallography.Site;

public class MopacEnergy implements Energy {

	private static String execpath;
	
	public MopacEnergy(String[] args) {		
		
		if (args == null || args.length < 1)
			GAParameters.usage("Not enough parameters given to MopacEnergy", true);
		
		// read path to executable, ensure it points directly to the executable
		execpath = args[0];
/*		int length = execpath.length();
		if (!execpath.substring(length-4).equals("exe")) {
			if (execpath.substring(length-1).equals("/")) {
				execpath = execpath + "MOPAC2009.exe";
			}
			else {
				execpath = execpath + "/MOPAC2009.exe";
			}
		}
*/		
	}
	
	public double getEnergy(StructureOrg c) {
		GAParameters params = GAParameters.getParams();	
		
		// Copy original structure (for testing)
		c.getCell().writeCIF(params.getTempDirName() + "/orig" + c.getID() + ".cif");
		
		// Run MOPAC
		runMopac(c);
		
		// Parse final energy
		Double energy = parseFinalEnergy(params.getTempDirName() + "/" + c.getID() + ".out");
		
		return energy;
	}
	
	private static String writeInput(StructureOrg c) {
		GAParameters params = GAParameters.getParams();
		
		String outdir = params.getTempDirName() + "/" + c.getID() + ".mop";
		//TODO: add "LET DDMIN=0.0"? will help with "NUMERICAL PROBLEMS IN BRACKETING LAMBDA" error
		// uses PM6, overrides interatomic distance check, uses all cartesian coordinates
		String keywds = "PM6 GEO-OK XYZ \n";
		String title = "Structure " + c.getID() + "\n\n";
		
		List<Vect> latVects = c.getCell().getLatticeVectors();
		List<Site> sites = c.getCell().getSites();
		
		// Creates list of atomic sites
		String atoms = "";
		for (Site s: sites) {
			List<Double> coords = s.getCoords().getCartesianComponents();
			atoms = atoms + s.getElement().getSymbol() + "  ";
			for (int k=0; k<Constants.numDimensions; k++) {
				atoms = atoms + coords.get(k) + "      ";
			}
			atoms = atoms + "\n";
		}
		
		// Creates list of lattice vectors
		String lattice = "";
		for (Vect v: latVects) {
			List<Double> xyz = v.getCartesianComponents();
			lattice = lattice + "Tv   ";
			for (int i=0; i<Constants.numDimensions; i++) {
				lattice = lattice + xyz.get(i) + "     ";
			}
			lattice = lattice + "\n";
		}
		
		String total = keywds + title + atoms + lattice;
		Utility.writeStringToFile(total, outdir);
		
		return outdir;
	}
	
	public static void runMopac(StructureOrg c) {
		GAParameters params = GAParameters.getParams();
		int verbosity = params.getVerbosity();
		
		// Write input files
		String input = writeInput(c);		
		
		// Execute MOPAC
		String s = null;
		BufferedReader stdInput = null;
		BufferedReader stdError = null;
		try {
			Process p = Runtime.getRuntime().exec(execpath + " " + input);

			stdInput = new BufferedReader(new InputStreamReader(
					p.getInputStream()));
			stdError = new BufferedReader(new InputStreamReader(
					p.getErrorStream()));
			
			// read the output
			while ((s = stdInput.readLine()) != null) {
				if (verbosity >= 5)
					System.out.println(s);
			}

			// print out any errors
			while ((s = stdError.readLine()) != null) {
				System.out.println(s);
			}		
		
		} catch (IOException e) {
			System.out.println("IOException in MopacEnergy.runMopac: " + e.getMessage());
			System.exit(-1);
		} finally {
			if (stdInput != null) 
				try{ stdInput.close(); } catch (Exception x) { } //ignore
			if (stdError != null) 
				try{ stdError.close(); } catch (Exception x) { } //ignore
		}
			
		// Parse final structure, set as return structure
		if (parseStructure(c, params.getTempDirName() + "/" + c.getID() + ".out") == null) {
			if (verbosity >= 3) {
				System.out.println("Warning: bad MOPAC CIF.  Not updating structure.");
			}
		} else {
			Cell p = parseStructure(c, params.getTempDirName() + "/" + c.getID() + ".out");
			c.setCell(p);
		}
		
	}
	
	public static Cell parseStructure(StructureOrg c, String output) {
		GAParameters params = GAParameters.getParams();
		int verbosity = params.getVerbosity();
		
		List<Site> sites = c.getCell().getSites();
		List<Vect> newVects = new LinkedList<Vect>();
		List<Site> newSites = new LinkedList<Site>();
		
		// check for case where structure is already optimized
		String line = null;
		Pattern goodPattern = Pattern.compile("     GRADIENTS WERE INITIALLY ACCEPTABLY SMALL");
		Matcher goodMatcher = goodPattern.matcher(output);
		try {
			BufferedReader t = new BufferedReader(new FileReader(output));
			try {
				while ((line = t.readLine()) != null) {
					goodMatcher.reset(line);
					if (goodMatcher.find()) {
						System.out.println("Structure already optimized, returning original");
						return c.getCell();
					}
				}
			} catch (IOException x) {
				if (verbosity >= 1)
					System.out.println("MopacEnergy: IOException.");
			}
		} catch (FileNotFoundException e) {
			System.out.println("MopacEnergy.parseStructure: .out not found");
			return c.getCell();
		}
		
		// parse the output to return a structure
		line = null;
		Pattern coordsPattern = Pattern.compile("       FINAL  POINT  AND  DERIVATIVES");
		Matcher coordsMatcher = coordsPattern.matcher(output);
		try {
			BufferedReader r = new BufferedReader(new FileReader(output));			
			try {
				while ((line = r.readLine()) != null) {
					coordsMatcher.reset(line);
//					System.out.print("line: " + line + "\n");

					if (coordsMatcher.find()) {
//						System.out.println("here's the line: " + line);
						r.readLine(); r.readLine();
						//TODO: is this worth improving (i.e. not random-looking numbers)?
						for (int e=0; e<(3*sites.size()+16); e++) {
							r.readLine();
						}
						line = r.readLine();
						coordsMatcher.reset(line);
						try {
							// read in atomic locations
							for (Site s: sites) {
								StringTokenizer t = new StringTokenizer(line);
//								System.out.println("this is the token zone: " + line);
								t.nextToken(); t.nextToken();
								Double x = Double.parseDouble(t.nextToken()); t.nextToken();
								Double y = Double.parseDouble(t.nextToken()); t.nextToken();
								Double z = Double.parseDouble(t.nextToken());
								Vect v = new Vect(x,y,z);
								newSites.add(new Site(s.getElement(),v));
								line = r.readLine();
							}
							
							// read in lattice vectors
							for (int k=0; k<Constants.numDimensions; k++) {
								StringTokenizer m = new StringTokenizer(line);
//								System.out.println("this is the lattice token zone: " + line);
								m.nextToken(); m.nextToken();
								Double x = Double.parseDouble(m.nextToken()); m.nextToken();
								Double y = Double.parseDouble(m.nextToken()); m.nextToken();
								Double z = Double.parseDouble(m.nextToken());
								Vect v = new Vect(x,y,z);
								newVects.add(v);
								line = r.readLine();
							}
							
						} catch (NumberFormatException x) {
							if (verbosity >= 3) 
								System.out.println("MopacEnergy.parseStructure: " + x.getMessage());
							if (verbosity >= 5) {
								System.out.println("MOPAC output follows: ");
								System.out.println(output);
							}
						}
						break;
					}
				}
			} catch (IOException x) {
				if (verbosity >= 1)
					System.out.println("MopacEnergy: IOException.");
			}
		} catch (FileNotFoundException e) {
			System.out.println("MopacEnergy.parseStructure: .out not found");
			return c.getCell();
		}

		
		Cell p = new Cell(newVects, newSites);
		
		//TODO: remove this line, is just for testing
//		p.writeCIF(params.getTempDirName() + "/" + c.getID() + "/parsed" + c.getID() + ".cif");
				
		return p;
	}
	
	//TODO: parsing from the .arc file rather than .out might be cleaner, though not necessary
	public static Double parseFinalEnergy(String output) {
		GAParameters params = GAParameters.getParams();
		int verbosity = params.getVerbosity();
		
		Double finalEnergy = Double.POSITIVE_INFINITY;
		
		// parse the output to return the final energy
		String line = null;
		Pattern coordsPattern = Pattern.compile("TOTAL ENERGY");
		Matcher coordsMatcher = coordsPattern.matcher(output);
		try {
			BufferedReader r = new BufferedReader(new FileReader(output));			
			try {
				while ((line = r.readLine()) != null) {
					coordsMatcher.reset(line);

					if (coordsMatcher.find()) {
//						System.out.println("indicator line: " + line);
						try {
							StringTokenizer m = new StringTokenizer(line);
//							System.out.println("energy line: " + line);
							m.nextToken(); m.nextToken(); m.nextToken();
							finalEnergy = Double.parseDouble(m.nextToken());
						} catch (NumberFormatException x) {
							if (verbosity >= 3) 
								System.out.println("MopacEnergy.parseFinalEnergy: " + x.getMessage());
							if (verbosity >= 5) {
								System.out.println("MOPAC output follows: ");
								System.out.println(output);
							}
						}
						break;
					}
				}
			} catch (IOException x) {
				if (verbosity >= 1)
					System.out.println("MopacEnergy: IOException.");
			}
		} catch (FileNotFoundException e) {
			System.out.println("MopacEnergy.parseFinalEnergy: .out not found");
		}
		
		return finalEnergy;
	}


	//testing
/*	public static void main(String[] args) {
		String[] argsIn = {"/home/skw57/Downloads/dl_poly_4.02/execute/"};
		DLPolyEnergy bob = new DLPolyEnergy(argsIn);
		Cell c = Cell.parseCif(new File("/home/skw57/2.cif"));
		
		bob.getEnergy(new StructureOrg(c));
		
	//	String output = GulpEnergy.runGULP("mno2_poly.gin");
	//	System.out.println(output);
	//	System.out.println(GulpEnergy.parseFinalEnergy(output, bob.cautious));
	//	System.out.println(output.contains("failed"));
	}
*/
	
	
}
