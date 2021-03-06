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

package mopac;

import ga.Energy;
import ga.GAOut;
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
import java.util.ArrayList;
import java.util.LinkedList;
import java.util.List;
import java.util.StringTokenizer;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import utility.Constants;
import utility.Utility;
import utility.Vect;
import vasp.VaspOut;
import crystallography.Cell;
import crystallography.Site;

public class MopacEnergy implements Energy {

	private static String execpath;

	public MopacEnergy(List<String> args) {		

		if (args == null || args.size() < 1)
			GAParameters.usage("Not enough parameters given to MopacEnergy", true);

		// read path to executable
		execpath = args.get(0);
	}

	public double getEnergy(StructureOrg c) {
		GAParameters params = GAParameters.getParams();	

		// Run MOPAC
		runMopac(c);

		// Parse final energy
		Double energy = parseFinalEnergy(params.getTempDirName() + "/" + c.getID() + ".out");

		return energy;
	}

	private static String writeInput(StructureOrg c) {
		GAParameters params = GAParameters.getParams();

		String outdir = params.getTempDirName() + "/" + c.getID() + ".mop";
		// uses PM7, overrides interatomic distance check, uses all cartesian coordinates
		String keywds = "PM7 GEO-OK XYZ T=500.00M RELSCF=0.01 RMIN=-10\n"; 
		String title = "Structure " + c.getID() + "\n\n";

		List<Vect> latVects = c.getCell().getLatticeVectors();
		List<Site> sites = c.getCell().getSites();

		// Creates list of atomic sites
		String atoms = "";
		for (Site s: sites) {
			List<Double> coords = s.getCoords().getCartesianComponents();
			atoms = atoms + s.getElement().getSymbol() + "  ";
			for (int k=0; k<Constants.numDimensions; k++) {
				atoms = atoms + coords.get(k) + "   1   "; // 0 means dont optimize coordinate; 1 means do opt
			}
			atoms = atoms + "\n";
		}

		// Creates list of lattice vectors
		String lattice = "";
		for (Vect v: latVects) {
			List<Double> xyz = v.getCartesianComponents();
			lattice = lattice + "Tv   ";
			for (int i=0; i<Constants.numDimensions; i++) {
				lattice = lattice + xyz.get(i) + "   1   ";
			}
			lattice = lattice + "\n";
		}

		String total = keywds + title + atoms + lattice;
		Utility.writeStringToFile(total, outdir);

		return outdir;
	}

	public static void runMopac(StructureOrg c) {
		GAParameters params = GAParameters.getParams();

		// Write input files
		String input = writeInput(c);		

		// Execute MOPAC
		String s = null;
		BufferedReader stdInput = null;
		BufferedReader stdError = null;
		Process p = null;
		try {
			p = Runtime.getRuntime().exec(execpath + " " + input);

			stdInput = new BufferedReader(new InputStreamReader(
					p.getInputStream()));
			stdError = new BufferedReader(new InputStreamReader(
					p.getErrorStream()));

			// read the output
			while ((s = stdInput.readLine()) != null) {
				GAOut.out().stdout(s, GAOut.DEBUG);
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

		// Parse final structure, set as return structure
		if (parseStructure(c, params.getTempDirName() + "/" + c.getID() + ".out") == null) {
			GAOut.out().stdout("Warning: bad MOPAC CIF.  Not updating structure.", GAOut.WARNING, c.getID());
		} else {
			Cell finalStructure = parseStructure(c, params.getTempDirName() + "/" + c.getID() + ".out");
			c.setCell(finalStructure);
		}

	}

	public static Cell parseStructure(StructureOrg c, String output) {
		GAParameters params = GAParameters.getParams();

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
				GAOut.out().stdout("MopacEnergy: IOException: " + x.getLocalizedMessage(), GAOut.CRITICAL, c.getID());
			}
		} catch (FileNotFoundException e) {
			System.out.println("MopacEnergy.parseStructure: .out not found");
			return c.getCell();
		}

		// parse the output to return a structure
		line = null;
		boolean finalGeometryFound = false;
		Pattern coordsPattern = Pattern.compile(" *ATOM *CHEMICAL *X *Y *Z *");
		Matcher coordsMatcher = coordsPattern.matcher(output);
		try {
			BufferedReader r = new BufferedReader(new FileReader(output));			
			try {
				while ((line = r.readLine()) != null) {
					coordsMatcher.reset(line);
					//					System.out.print("line: " + line + "\n");

					if (coordsMatcher.find()) {
						while ((line = r.readLine()) != null) {
							coordsMatcher.reset(line);
							if (coordsMatcher.find()) {
								finalGeometryFound = true;
								//						System.out.println("here's the line: " + line);
								r.readLine(); r.readLine();
								//TODO: is this worth improving (i.e. not random-looking numbers)?
								//						for (int e=0; e<16; e++) {
								//							r.readLine();
								//						}
								line = r.readLine().replace("*", ""); // sometimes the asterisks are included in the output, sometimes not

								//System.out.println(line);
								coordsMatcher.reset(line);
								try {
									// read in atomic locations
									for (Site s: sites) {
										StringTokenizer t = new StringTokenizer(line);
										//								System.out.println("this is the token zone: " + line);
										t.nextToken(); t.nextToken();
										Double x = Double.parseDouble(t.nextToken());
										Double y = Double.parseDouble(t.nextToken());
										Double z = Double.parseDouble(t.nextToken());
										Vect v = new Vect(x,y,z);
										newSites.add(new Site(s.getElement(),v));
										line = r.readLine().replace("*", "");
										//System.out.println(line);
									}

									// read in lattice vectors
									for (int k=0; k<Constants.numDimensions; k++) {
										StringTokenizer m = new StringTokenizer(line);
										//								System.out.println("this is the lattice token zone: " + line);
										m.nextToken(); m.nextToken();
										Double x = Double.parseDouble(m.nextToken()); 
										Double y = Double.parseDouble(m.nextToken()); 
										Double z = Double.parseDouble(m.nextToken());
										Vect v = new Vect(x,y,z);
										newVects.add(v);
										line = r.readLine().replace("*", "");
										//System.out.println(line);
									}

								} catch (NumberFormatException x) {
									GAOut.out().stdout("NumberFormatException in MopacEnergy.parseStructure: " + x.getMessage(), GAOut.NOTICE, c.getID());
									GAOut.out().stdout("MOPAC output follows: ", GAOut.DEBUG, c.getID());
									GAOut.out().stdout(output, GAOut.DEBUG, c.getID());
								}
								break;
							}
						}
					}
				}
			} catch (IOException x) {
				GAOut.out().stdout("MopacEnergy: IOException: " + x.getLocalizedMessage(), GAOut.CRITICAL, c.getID());
			}
		} catch (FileNotFoundException e) {
			System.out.println("MopacEnergy.parseStructure: .out not found");
			return c.getCell();
		}

		if (finalGeometryFound)
			return new Cell(newVects, newSites);
		else
			return null;
	}

	//TODO: parsing from the .arc file rather than .out might be cleaner, though not necessary
	public static Double parseFinalEnergy(String output) {
		GAParameters params = GAParameters.getParams();

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
							GAOut.out().stdout("MopacEnergy.parseFinalEnergy: " + x.getMessage(), GAOut.NOTICE);
							GAOut.out().stdout("MOPAC output follows: ", GAOut.DEBUG);
							GAOut.out().stdout(output, GAOut.DEBUG);
						}
						break;
					}
				}
			} catch (IOException x) {
				GAOut.out().stdout("MopacEnergy: IOException: " + x.getLocalizedMessage(), GAOut.CRITICAL);
			}
		} catch (FileNotFoundException e) {
			System.out.println("MopacEnergy.parseFinalEnergy: .out not found");
		}

		return finalEnergy;
	}

	public boolean cannotCompute(StructureOrg o) {
		
		double minSphereRadius = 4.5; // necessary to get accurate results, according to jimmy stewart
		
		List<Vect> vects = o.getCell().getLatticeVectors();
		Vect a = vects.get(0);
		Vect b = vects.get(1);
		Vect c = vects.get(2);

		Vect n1 = a.cross(b);
		Vect n2 = c.cross(b);
		Vect n3 = a.cross(c);
		
		Vect center = (a.plus(b).plus(c)).scalarMult(0.5);
		
		double d1 = Math.abs(n1.dot(center)) / n1.length();
		double d2 = Math.abs(n1.dot(center)) / n2.length();
		double d3 = Math.abs(n1.dot(center)) / n3.length();
	
		if (d1 < minSphereRadius || d2 < minSphereRadius || d3 < minSphereRadius)
			return true;
		else		
			return false;
	}
	
	// just for testing
	public static void main(String args[]) {
	//	StructureOrg a = new StructureOrg(VaspOut.getPOSCAR("/home/bcr48/GA/2D/garun_Si7/926059.POSCAR"));
	//	StructureOrg b = new StructureOrg(VaspOut.getPOSCAR("/home/bcr48/GA/2D/garun_Si7/1678920.POSCAR"));
	//	runMopac(a);
	//	runMopac(b);
	//	System.out.println(a);
	//	System.out.println(b);
		//	String output = "/home/wtipton/59.out";
	//	List<String> as = new ArrayList<String>();
	//	as.add("");
	//	StructureOrg c = new StructureOrg(VaspOut.getPOSCAR("/home/wtipton/1.POSCAR"));
	//	System.out.println((new MopacEnergy(as)).cannotCompute(c));
	}

}
