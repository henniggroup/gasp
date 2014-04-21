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


package dlpoly;

import ga.Energy;
import ga.GAOut;
import ga.GAParameters;
import ga.GAUtils;
import ga.StructureOrg;
import ga.UnitsSOCreator;
import gulp.GulpEnergy;

import java.io.*;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.LinkedList;
import java.util.List;
import java.util.Scanner;
import java.util.StringTokenizer;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import utility.Constants;
import utility.Utility;
import utility.Vect;

import crystallography.Cell;
import crystallography.Site;

import org.openbabel.*;

public class DLPolyEnergy implements Energy {

	public DLPolyEnergy(List<String> args) {		
		
		if (args == null || args.size() < 2)
			GAParameters.usage("Not enough parameters given to DLPolyEnergy", true);
		
		// read in the location of DLPoly files
		String loc = args.get(0);
		
		// read in the location of the potl file
		String potl = args.get(1);
		
/*		// create a directory for calldlpoly files
		File calldl = new File("/home/skw57/bin/calldl");
		calldl.mkdir();
*/		
		// create file containing location
		if (loc.charAt(loc.length()-1) == '/') {
			loc = loc + "execute/DLPOLY.Z";
		}
		else
			loc = loc + "/execute/DLPOLY.Z";
		
		Utility.writeStringToFile(loc, GAParameters.getParams().getOutDirName() + "/loc");
	}
	
	public double getEnergy(StructureOrg c) {
		GAParameters params = GAParameters.getParams();		
		
		// Prepare directory for this structure
		String subPath = params.getTempDirName() + "/" + c.getID();
		File newDir = new File(subPath);
		newDir.mkdir();
		
		// Copy original structure (for testing)
		c.getCell().writeCIF(params.getTempDirName() + "/" + c.getID() + "/orig" + c.getID() + ".cif");
		
		runDLPoly(c);
		
		//Parse final energy
		Double energy = parseFinalEnergy(params.getTempDirName() + "/" + c.getID() + "/OUTPUT");
		
		return energy;
	}
	
	private static void writeConfig(StructureOrg c) {
		GAParameters params = GAParameters.getParams();
		
		String outdir = params.getTempDirName() + "/" + c.getID() + "/CONFIG";
//		String outdir = "/home/skw57/polyfiles/CONFIG";
		String title = "Structure " + c.getID() + "\n";
		// '0' -- coordinates only, '3' -- parallelepiped boundary conditions
		String record = "0         3" + "\n";
		
		
		List<Vect> latVects = c.getCell().getLatticeVectors();
		List<Site> sites = c.getCell().getSites();
		
		String lattice = "";
		for (Vect v: latVects) {
			List<Double> xyz = v.getCartesianComponents();
			for (int i=0; i<Constants.numDimensions; i++) {
				lattice = lattice + xyz.get(i) + "     ";
			}
			lattice = lattice + "\n";
		}
		
		String atoms = "";
		int n = 1;
		for (Site s: sites) {
			List<Double> coords = s.getCoords().getCartesianComponents();
			atoms = atoms + s.getElement().getSymbol() + "         " + n + "\n";
			for (int k=0; k<Constants.numDimensions; k++) {
				atoms = atoms + coords.get(k) + "      ";
			}
			atoms = atoms + "\n";
			n++;
		}
		String total = title + record + lattice + atoms;
		Utility.writeStringToFile(total, outdir);
	}
	
	public static void writeControl(StructureOrg c) {
		GAParameters params = GAParameters.getParams();
		String newline = GAUtils.newline();
		
		// Calculate maximum cutoff value (minimum half-cell width)
		Vect x = c.getCell().getLatticeVectors().get(0);
		Vect y = c.getCell().getLatticeVectors().get(1);
		Vect z = c.getCell().getLatticeVectors().get(2);
		Vect xy = x.cross(y); Vect xz = x.cross(z); Vect yz = y.cross(z);
		Double one = yz.magnitudeOfProjection(x); Double two = xz.magnitudeOfProjection(y); Double three = xy.magnitudeOfProjection(z);
//		System.out.println("lengths: " + 0.5*one + " " + 0.5*two + " " + 0.5*three);
		Double l = 0.5*Math.min(Math.min(one, two), Math.min(two,three)); l = 0.99*l; String m = l.toString().substring(0, 6); l = Double.parseDouble(m);
//		System.out.println("cutoff: " + l);
		
		try {
			File f = new File(params.getTempDirName() + "/" + c.getID(),"CONTROL");
//			File f = new File("/home/skw57/polyfiles/","CONTROL");

			// write the file
			BufferedWriter out = new BufferedWriter(new FileWriter(f));
			out.write(params.getRunTitle() + newline + newline);
			out.write("minimise        energy 1");
			out.write(newline);
			out.write("ensemble nvt evans");
			out.write(newline);
//			out.write("ewald precision");
//			out.write(newline);
			out.write("equilibration        5000");
			out.write(newline);
			out.write("temperature            1");
			out.write(newline);
			out.write("steps               10000");
			out.write(newline);
			out.write("print                 50");
			out.write(newline);
			out.write("stack                 10");
			out.write(newline);
			out.write("stats                 50");
			out.write(newline);
			out.write("timestep          5.0E-7");
			out.write(newline);
			out.write("cutoff             " + l);
			out.write(newline + newline);
			out.write("job time       1000000.0");
			out.write(newline);
			out.write("close time         500.0");
			out.write(newline + newline);
			out.write("finish");
			out.close();
		} catch (IOException e) {
				System.out.println("DLPolyEnergy IOException: " + e.getMessage());
		}
	}
	
	//TODO: might be a problem running in parallel (again with the #units deal)
		// point code towards a potential file a la GULP
	private static void writeField(StructureOrg c) {
		GAParameters params = GAParameters.getParams();
		int difUnits = UnitsSOCreator.getDifUnits();
		int[] numUnits = UnitsSOCreator.getTargets();
		int[] numAtoms = UnitsSOCreator.getNumAtoms();
		List<Site> sites = c.getCell().getSites();
		
		String outdir = params.getTempDirName() + "/" + c.getID() + "/FIELD";
//		String outdir = "/home/skw57/polyfiles/FIELD";
		String title = "Structure " + c.getID() + "\n";
		String units = "UNITS eV \n";
		
		String num = "MOLECULES " + difUnits + "\n";
		
		String molecules = "";
		int loc = 0;
		for (int k=0; k<difUnits; k++) {
			molecules = molecules + "mol" + (k+1) + "\n";
			molecules = molecules + "NUMMOLS " + numUnits[k] + "\n";
			molecules = molecules + "ATOMS " + numAtoms[k] + "\n";
			for (int l=0; l<numAtoms[k]; l++) {
				//TODO: get rid of this junk(?)
				Double n = sites.get(loc+l).getElement().getAtomicMass(); int y = n.toString().length(); int dif = 8 - y; String mass = "";
				for (int m=0; m<dif; m++) {
					mass = mass + " ";
				}
				mass = mass + n.toString();
				
				// TODO: partial charge is necessary, and should probably be inputted
				n = sites.get(loc+l).getElement().getCharge(); y = n.toString().length(); dif = 6 - y; String charge = "";
				for (int m=0; m<dif; m++) {
					charge = charge + " ";
				}
				charge = charge + n.toString();
				molecules = molecules + "    " + sites.get(loc+l).getElement().getSymbol() + "     " + mass + "     " + charge + "\n";
			}
			molecules = molecules + "FINISH \n";
			loc = loc + numAtoms[k]*numUnits[k];
		}
		
//		String potl = "VDW 1 \nO     O     LJ    0.16    3.196 \n";
		String potl = "VDW 2 \nO     O     buck  15764   0.149  27.88 \nH     O     buck  311.97  0.25  0.0 \n";
//		String potl = "VDW 3 \nAl    Al    LJ   0.392   1.4472 \nCu    Cu    LJ   0.415   2.277 \nAl    Cu    LJ   0.4   1.76 \n";
		
		
		String total = title + units + num + molecules + potl + "CLOSE";
		Utility.writeStringToFile(total, outdir);
	}
	
	public static void runDLPoly(StructureOrg c) {
		GAParameters params = GAParameters.getParams();
		
		// Write input files
		writeConfig(c);
		writeControl(c);
		writeField(c);
		
		//TODO: why is this necessary? seems like a problem with the code, very weird
		// Modify and save
		String control = Utility.readStringFromFile(params.getTempDirName() + "/" + c.getID() + "/CONTROL");
		control.concat("\n ");
		Utility.writeStringToFile(control, params.getTempDirName() + "/" + c.getID() + "/CONTROL");
		
		String field = Utility.readStringFromFile(params.getTempDirName() + "/" + c.getID() + "/FIELD");
		field.concat("\n ");
		Utility.writeStringToFile(field, params.getTempDirName() + "/" + c.getID() + "/FIELD");
		
		//TODO: need to change this to write multiple files so running in parallel is possible
		// Write calldlpoly script
		String poly = "cd " + params.getTempDirName() + "/" + c.getID() + "\n \n" + Utility.readStringFromFile(params.getOutDirName() + "/loc");
//		Utility.writeStringToFile(poly, "/home/skw57/bin/calldl/calldlpoly" + c.getID());
		Utility.writeStringToFile(poly, "/home/skw57/bin/calldlpoly");
		
		
		// Execute DL_Poly
		String s = null;
		BufferedReader stdInput = null;
		BufferedReader stdError = null;
		try {
//			Process p = Runtime.getRuntime().exec("sh calldlpoly" + c.getID());
			Process p = Runtime.getRuntime().exec("calldlpoly");

			stdInput = new BufferedReader(new InputStreamReader(
					p.getInputStream()));
			stdError = new BufferedReader(new InputStreamReader(
					p.getErrorStream()));
			
			// read the output
			while ((s = stdInput.readLine()) != null) {
				GAOut.out().stdout(s, GAOut.DEBUG, c.getID());
			}

			// print out any errors
			while ((s = stdError.readLine()) != null) {
				System.out.println(s);
			}			
			
		} catch (IOException e) {
			System.out.println("IOException in DLPolyEnergy.runDLPoly: " + e.getMessage());
			System.exit(-1);
		} finally {
			if (stdInput != null) 
				try{ stdInput.close(); } catch (Exception x) { } //ignore
			if (stdError != null) 
				try{ stdError.close(); } catch (Exception x) { } //ignore
		}
		
//		System.out.println("here's the path: " + params.getTempDirName() + "/REVCON");
		
/*		//TODO: code is trying to access files before they've been written. how to wait more effectively?
		// Wait to allow REVCON to be generated
		try {
		Thread.sleep(10000);
		}
		catch(InterruptedException e) {
		e.printStackTrace();
		}
*/
		
		// Parse final structure, set as return structure
		if (parseStructure(c, params.getTempDirName() + "/" + c.getID() + "/REVCON") == null) {
			GAOut.out().stdout("Warning: bad DLPoly CIF.  Not updating structure.", GAOut.NOTICE, c.getID());
		} else {
			Cell p = parseStructure(c, params.getTempDirName() + "/" + c.getID() + "/REVCON");
			c.setCell(p);
		}
		
	}
	
	public static Cell parseStructure(StructureOrg c, String revcon) {
		GAParameters params = GAParameters.getParams();
		
		List<Site> sites = c.getCell().getSites();
		List<Vect> newVects = new LinkedList<Vect>();
		List<Site> newSites = new LinkedList<Site>();
		
		// parse the output to return a structure
		String line = null;
		Pattern coordsPattern = Pattern.compile("Structure");
		Matcher coordsMatcher = coordsPattern.matcher(revcon);
		try {
			BufferedReader r = new BufferedReader(new FileReader(revcon));			
			try {
				while ((line = r.readLine()) != null) {
					coordsMatcher.reset(line);

					if (coordsMatcher.find()) {
//						System.out.println("here's the line: " + line);
						r.readLine(); line = r.readLine();
						coordsMatcher.reset(line);
						try {
							for (int k=0; k<Constants.numDimensions; k++) {
								StringTokenizer m = new StringTokenizer(line);
//								System.out.println("this is the lattice token zone: " + line);
								Vect v = new Vect(Double.parseDouble(m.nextToken()), Double.parseDouble(m.nextToken()), Double.parseDouble(m.nextToken()));
								newVects.add(v);
								line = r.readLine();
								if (k==2)
									line = r.readLine();
							}
							
							int n = 0;
							for (Site s: sites) {
								StringTokenizer t = new StringTokenizer(line);
//								System.out.println("this is the token zone: " + line);
								Vect v = new Vect(Double.parseDouble(t.nextToken()), Double.parseDouble(t.nextToken()), Double.parseDouble(t.nextToken()));
								newSites.add(new Site(s.getElement(),v));
								if (n < sites.size())
									r.readLine(); r.readLine(); r.readLine(); line = r.readLine();
							}
						} catch (NumberFormatException x) {
							GAOut.out().stdout("DLPolyEnergy.parseStructure: " + x.getMessage(), GAOut.NOTICE, c.getID());
							GAOut.out().stdout("DLPoly output follows: ", GAOut.DEBUG, c.getID());
							GAOut.out().stdout(revcon, GAOut.DEBUG, c.getID());
						}
						break;
					}
				}
			} catch (IOException x) {
				GAOut.out().stdout("DLPolyEnergy: IOException: " + x.getLocalizedMessage(), GAOut.CRITICAL, c.getID());
			}
		} catch (FileNotFoundException e) {
			System.out.println("DLPolyEnergy.parseStructure: REVCON not found");
			return c.getCell();
		}

		
		Cell p = new Cell(newVects, newSites);
		
		//TODO: remove this line, is just for testing
//		p.writeCIF(params.getTempDirName() + "/" + c.getID() + "/parsed" + c.getID() + ".cif");
				
		return p;
	}
	
	public static Double parseFinalEnergy(String output) {
		GAParameters params = GAParameters.getParams();
		
		Double finalEnergy = Double.POSITIVE_INFINITY;
		
		// parse the output to return the final energy
		String line = null;
		Pattern coordsPattern = Pattern.compile("run terminated after");
		Matcher coordsMatcher = coordsPattern.matcher(output);
		try {
			BufferedReader r = new BufferedReader(new FileReader(output));			
			try {
				while ((line = r.readLine()) != null) {
					coordsMatcher.reset(line);

					if (coordsMatcher.find()) {
//						System.out.println("indicator line: " + line);
						r.readLine(); r.readLine(); r.readLine(); r.readLine(); r.readLine(); r.readLine(); r.readLine(); r.readLine(); r.readLine(); line = r.readLine();
						coordsMatcher.reset(line);
						try {
							StringTokenizer m = new StringTokenizer(line);
//							System.out.println("energy line: " + line);
							m.nextToken();
							finalEnergy = Double.parseDouble(m.nextToken());
						} catch (NumberFormatException x) {
							GAOut.out().stdout("DLPolyEnergy.parseFinalEnergy: " + x.getMessage(), GAOut.NOTICE);
							GAOut.out().stdout("DLPoly output follows: ", GAOut.DEBUG);
							GAOut.out().stdout(output, GAOut.DEBUG);
						}
						break;
					}
				}
			} catch (IOException x) {
				GAOut.out().stdout("DLPolyEnergy: IOException: " + x.getLocalizedMessage(), GAOut.CRITICAL);
			}
		} catch (FileNotFoundException e) {
			GAOut.out().stdout("DLPolyEnergy.parseFinalEnergy: OUTPUT not found", GAOut.NOTICE);
		}
		
		return finalEnergy;
	}
	
	
	public boolean cannotCompute(StructureOrg o) {
		return false;
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