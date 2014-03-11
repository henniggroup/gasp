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

package vasp;

import crystallography.Cell;
import ga.GAOut;
import ga.GAParameters;

import java.io.*;
import java.util.*;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import javax.xml.parsers.ParserConfigurationException;

import org.xml.sax.SAXException;

import utility.*;
import crystallography.Site;
import chemistry.Element;

/* Reads files in a vasp directory into appropriate data structures
 * 
 * Will Tipton
 * 
 * TODO:
 *  - make use of other vasp input datas?
 *  - In the future, this should really just do the parsing and then
 *    create a VaspData object which holds all the parameters, etc. that
 *    go into a vasp run.  Then VaspOut can take a VaspData file and 
 *    write vasp input files to disk.
 */ 

//NB: Convention: when runs fail, energy should evaluate to Double.POSITIVE_INFINITY

public class VaspOut {
	
	static final String vaspSuccessString = "reached required accuracy";

	public static Cell getPOSCAR(String poscarInFile) {
		String description = null;	/* Will hold first line of the POSCAR */
		
		/**************** Parse me a POSCAR *******************/
		List<Vect> latticeVectors = new LinkedList<Vect>();
		List<Site> basis = new LinkedList<Site>();
		
		StringTokenizer tokens = null;
		BufferedReader poscarReader = null;
		try {
			poscarReader = new BufferedReader(new FileReader(poscarInFile));
			/* Read the comment line */
			description = poscarReader.readLine() + "; Read from " + poscarInFile; 
			/* Get the scaling factor */
			double scalingFactor = Double.parseDouble(poscarReader.readLine()); 
			/* Get lattice vectors */
			for (int i = 0; i < Constants.numDimensions; i++) {
				tokens = new StringTokenizer(poscarReader.readLine());
				Double x = Double.parseDouble(tokens.nextToken()) * scalingFactor;
				Double y = Double.parseDouble(tokens.nextToken()) * scalingFactor;
				Double z = Double.parseDouble(tokens.nextToken()) * scalingFactor;
				latticeVectors.add(new Vect(x,y,z));
			}
			/* Get components */
			tokens = new StringTokenizer(poscarReader.readLine());
			List<Element> elements = new LinkedList<Element>();
			while (tokens.hasMoreTokens()) 
				elements.add(Element.getElemFromSymbol(tokens.nextToken()));
			
			/* Get numbers of components */
			tokens = new StringTokenizer(poscarReader.readLine());
			List<Integer> numsOfComponents = new LinkedList<Integer>();
			
			while (tokens.hasMoreTokens()) 
				numsOfComponents.add(Integer.parseInt(tokens.nextToken()));
			
			/* Skip Selective Dynamics and get Direct or Cartesian */
			String line = poscarReader.readLine();
			if (line.startsWith("S") || line.startsWith("s"))
				line = poscarReader.readLine();
			Boolean usesCartesianCoords = line.startsWith("C") || line.startsWith("c") 
											|| line.startsWith("k") || line.startsWith("K");
			/* Get site positions */
			for (int elemNum = 0; elemNum < numsOfComponents.size(); elemNum++) {
				for (int i = 0; i < numsOfComponents.get(elemNum); i++) {
					tokens = new StringTokenizer(poscarReader.readLine());
					Double x = Double.parseDouble(tokens.nextToken());
					Double y = Double.parseDouble(tokens.nextToken());
					Double z = Double.parseDouble(tokens.nextToken());
					Vect siteLoc = null;
					if (usesCartesianCoords)
						siteLoc = new Vect(x,y,z);
					else
						siteLoc = new Vect(x,y,z,latticeVectors);
					basis.add(new Site(elements.get(elemNum), siteLoc));
				}
			}
		} catch (Exception x) {
			GAOut.out().stdout("Warning: VaspOut.getCell() failed: " + x.getMessage(), GAOut.NOTICE);

			return null;
		} finally {
			if (poscarReader != null)
				try {
					poscarReader.close();
				} catch (IOException x) {
					GAOut.out().stdout("Warning: VaspOut.getCell() failed to close poscarReader: " + x.getMessage(), GAOut.NOTICE);
				}
		}
		
		/* Make the VaspData object which for now just has a Cell */
		return new Cell(latticeVectors, basis, description);

	}
	
	public static List<VaspConfig> getConfigs(String outcarFileName, Cell origCell) {
		List<VaspConfig> result = new ArrayList<VaspConfig>();
		
		BufferedReader outcarReader = null;
		try {
			outcarReader = new BufferedReader(new FileReader(outcarFileName));
			String line;
			
			int currentConfigNum = 0;
			
			Pattern newConfigPattern = Pattern.compile("--* *Iteration *"+(currentConfigNum+1));
			Pattern stressPattern = Pattern.compile("in kB");
			Pattern energyPattern = Pattern.compile(" energy  without entropy"); // two spaces between energy and without :p
			Pattern atomsPattern = Pattern.compile("POSITION *TOTAL-FORCE");
			Pattern lvectsPattern = Pattern.compile("direct lattice vectors");
			
			VaspConfig currentConf = null;
			List<Vect> lVects = new ArrayList<Vect>();
			List<Site> sites = new ArrayList<Site>();
			while ((line = outcarReader.readLine()) != null) {
				if (newConfigPattern.matcher(line).find()) {
					if (currentConfigNum > 0) {
						Cell newCell = new Cell(lVects, sites, "Configuration " + currentConfigNum + " from " + outcarFileName);
						currentConf.setCell(newCell);
						result.add(currentConf);
					}
					currentConf = new VaspConfig();
					lVects = new ArrayList<Vect>();
					sites = new ArrayList<Site>();
					currentConfigNum++;
					newConfigPattern = Pattern.compile("--* *Iteration *"+(currentConfigNum+1));
				} else
				
				if (atomsPattern.matcher(line).find()) {
					outcarReader.readLine(); // skip a line of hyphens
					List<Vect> forces = new ArrayList<Vect>();
					for (int i = 0; i < origCell.getNumSites(); i++) {	
						StringTokenizer atomLineTok = new StringTokenizer(outcarReader.readLine());
						double x = Double.parseDouble(atomLineTok.nextToken());
						double y = Double.parseDouble(atomLineTok.nextToken());
						double z = Double.parseDouble(atomLineTok.nextToken());
						double fx = Double.parseDouble(atomLineTok.nextToken());
						double fy = Double.parseDouble(atomLineTok.nextToken());
						double fz = Double.parseDouble(atomLineTok.nextToken());
						sites.add(new Site(origCell.getSite(i).getElement(),new Vect(x,y,z)));
						forces.add(new Vect(fx,fy,fz));
					}
					currentConf.setForces(forces);
				} else
				
				if (lvectsPattern.matcher(line).find()) {
					for (int i = 0; i < 3; i++) {
						StringTokenizer vectLineTok = new StringTokenizer(outcarReader.readLine());
						double x = Double.parseDouble(vectLineTok.nextToken());
						double y = Double.parseDouble(vectLineTok.nextToken());
						double z = Double.parseDouble(vectLineTok.nextToken());
						lVects.add(new Vect(x,y,z));
					}
				} else
				
				if (stressPattern.matcher(line).find()) {
					StringTokenizer energyLineTok = new StringTokenizer(line);
					energyLineTok.nextToken();energyLineTok.nextToken(); // want 6 tokens starting with the 3rd
					double stress[] = new double[VaspConfig.STRESS_LEN];
					for (int i = 0; i < VaspConfig.STRESS_LEN; i++)
						stress[i] = Double.parseDouble(energyLineTok.nextToken()) / 1602.0;
					// TODO: the 1602 is copied from the potfit parser. no idea where it came from.
					currentConf.setStress(stress);
				} else
				
				if (energyPattern.matcher(line).find()) {
					StringTokenizer energyLineTok = new StringTokenizer(line);
					energyLineTok.nextToken();energyLineTok.nextToken(); // we want the 7th token
					energyLineTok.nextToken();energyLineTok.nextToken();
					energyLineTok.nextToken();energyLineTok.nextToken();
					currentConf.setTotalEnergy(Double.parseDouble(energyLineTok.nextToken()));
				}	
			}
			if (currentConf != null && currentConf.energyHasBeenSet()) {
				Cell newCell = new Cell(lVects, sites, "Configuration " + currentConfigNum + " from " + outcarFileName);
				currentConf.setCell(newCell);
				result.add(currentConf);
			}
			
		} catch (Exception x) {
			GAOut.out().stdout("Warning: VaspOut.getConfigs() failed: " + x.getMessage(), GAOut.NOTICE);
		} finally {
			if (outcarReader != null)
				try {
					outcarReader.close();
				} catch (IOException x) {
					GAOut.out().stdout("Warning: VaspOut.getConfigs() failed to close outcarReader: " + x.getMessage(), GAOut.NOTICE);
				}
		}
		
		return result;
	}

	
	public static double getFinalEnergy(String outcarFileName, boolean cautious) {
		
		double energy = Double.POSITIVE_INFINITY;
		boolean runConvergedSuccessfully = false;
		
		BufferedReader outcarReader = null;
		try {
			outcarReader = new BufferedReader(new FileReader(outcarFileName));
			String line;
			String energyRegexp = "TOTEN *= *[-\\.0-9]* eV";
			Pattern energyPattern = Pattern.compile(energyRegexp);
			while ((line = outcarReader.readLine()) != null) {
				Matcher energyMatcher = energyPattern.matcher(line);
				if (energyMatcher.find()) {
					StringTokenizer energyLineTok = new StringTokenizer(energyMatcher.group(energyMatcher.groupCount()));
					energyLineTok.nextToken();energyLineTok.nextToken();
					energy = Double.parseDouble(energyLineTok.nextToken());
				}
				if (line.contains(vaspSuccessString))
					runConvergedSuccessfully = true;
			}
		} catch (Exception x) {
			GAOut.out().stdout("Warning: VaspOut.getCell() failed: " + x.getMessage(), GAOut.NOTICE);
		} finally {
			if (outcarReader != null)
				try {
					outcarReader.close();
				} catch (IOException x) {
					GAOut.out().stdout("Warning: VaspOut.getCell() failed to close outcarReader: " + x.getMessage(), GAOut.NOTICE);
				}
		}		
		
		
		if (cautious && !runConvergedSuccessfully)
			return Double.POSITIVE_INFINITY;
		
		return energy;
	}
	/*
	public static double getFinalEnergy(String outcarFile, boolean cautious) {
		String outcar = Utility.readStringFromFile(outcarFile);
		
		double energy = Double.POSITIVE_INFINITY;
		if (cautious && !vaspRunConverged(outcar))
			return energy;
		
		// get energy: (kinda a hack.)
		try {
			String energyRegexp = "TOTEN *= *[-\\.0-9]* eV";
			Pattern energyPattern = Pattern.compile(energyRegexp);
			Matcher energyMatcher = energyPattern.matcher(outcar);
			StringTokenizer energyLineTok = null;
			while (energyMatcher.find())
				energyLineTok = new StringTokenizer(energyMatcher.group(energyMatcher.groupCount()));
			energyLineTok.nextToken();energyLineTok.nextToken();
			energy = Double.parseDouble(energyLineTok.nextToken());
		} catch (Exception x) {
			GAOut.out().stdout("Warning: VaspOut.getFinalEnergy() failed: " + x.getMessage(), GAOut.NOTICE);
		}
		
		return energy;
	}*/
	
	/*
	public static boolean vaspRunConverged(String outcar) {
		return outcar.contains(vaspSuccessString);
	} */

	public static void main(String args[]) {
		for (VaspConfig i : getConfigs("/home/wtipton/projects/ga_for_crystals/oldruns/garun_vasp1/temp/testrun.17/OUTCAR", 
				getPOSCAR("/home/wtipton/projects/ga_for_crystals/oldruns/garun_vasp1/temp/testrun.17/POSCAR"))) {
			System.out.println(i.getTotalEnergy());
			for (double d : i.getStress())
				System.out.print(d + " ");
			System.out.println("");
			System.out.println(i.getCell());
		}
	}
}
