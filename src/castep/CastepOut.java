<<<<<<< HEAD
=======
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

>>>>>>> 0e3189c40547bbd59ea42c4f91890d7511fb7797
package castep;

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

/* Reads files in a castep directory into appropriate data structures
 * 
 * Will Tipton
 * 
 * TODO:
 *  - make use of other castep input datas?
 *  - In the future, this should really just do the parsing and then
 *    create a CastepData object which holds all the parameters, etc. that
 *    go into a castep run.  Then CastepOut can take a CastepData file and 
 *    write castep input files to disk.
 */

//NB: Convention: when runs fail, energy should evaluate to Double.POSITIVE_INFINITY

public class CastepOut {
	
	static final String castepSuccessString = "Geometry optimization completed successfully";

	public static Cell getCell(String castepOutFile) {
		String description = null;	/* Will hold first line of the POSCAR */
		
		/**************** Parse me a POSCAR *******************/
		List<Vect> latticeVectors = new LinkedList<Vect>();
		List<Site> basis = new LinkedList<Site>();
		
		StringTokenizer tokens = null;
		try {
			BufferedReader poscarReader = new BufferedReader(new FileReader(castepOutFile));
			/* Read the comment line */
			description = poscarReader.readLine(); 
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
			GAOut.out().stdout("Warning: CastepOut.getCell() failed: " + x.getMessage(), GAOut.NOTICE);

			return null;
		} 
		
		/* Make the CastepData object which for now just has a Cell */
		return new Cell(latticeVectors, basis, null);

	}
	
	public static double getFinalEnergy(String castepOutFile, boolean cautious) {
		String out = Utility.readStringFromFile(castepOutFile);
		
		double energy = Double.POSITIVE_INFINITY;
		if (cautious && !castepRunConverged(castepOutFile))
			return energy;
		
		// get energy: (kinda a hack.)
		try {
			String energyRegexp = "Final Enthalpy *= *[-\\.0-9]* eV";
			Pattern energyPattern = Pattern.compile(energyRegexp);
			Matcher energyMatcher = energyPattern.matcher(castepOutFile);
			StringTokenizer energyLineTok = null;
			while (energyMatcher.find())
				energyLineTok = new StringTokenizer(energyMatcher.group(energyMatcher.groupCount()));
			energyLineTok.nextToken();energyLineTok.nextToken();
			energy = Double.parseDouble(energyLineTok.nextToken());
		} catch (Exception x) {
			GAOut.out().stdout("Warning: CastepOut.getFinalEnergy() failed: " + x.getMessage(), GAOut.NOTICE);
		}
		
		return energy;
	}
	
	public static boolean castepRunConverged(String castepOutFile) {
		return castepOutFile.contains(castepSuccessString);
	}

}
