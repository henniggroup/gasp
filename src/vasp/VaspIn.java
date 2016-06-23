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

import crystallography.*;
import java.io.*;
import java.text.NumberFormat;
import java.util.*;
import utility.*;

import chemistry.Element;
// writes internal data structures into VASP input files

public class VaspIn {
		
	private Cell cell;
	private List<Element> elements; /* Holds the elements involved in this
	 								 * VASP calculation in the order they come
	 								 * in the POTCAR file. */
	String description;	/* Holds the comment in the first line of the POSCAR file. */
	
	String kpointsFile;
	String incarFile;
	Map<Element,String> potcarFileMap;
	
	public VaspIn(Cell _cell,  String kpointsFile, String incarFile, Map<Element,String> potcarFileMap) {
		this(_cell, "Created by VaspData", kpointsFile, incarFile, potcarFileMap);
	}
	
	public VaspIn(Cell _cell, String _description, String _kpointsFile, String _incarFile, Map<Element,String> _potcarFileMap) {
		cell = _cell;
		elements = cell.getComposition().getElements();
		description = _description;
		kpointsFile = _kpointsFile;
		potcarFileMap = _potcarFileMap;
		incarFile = _incarFile;
	}
	
	public Cell getCell() {
		return cell;
	}
	
	public List<Element> getElements() {
		List<Element> result = new LinkedList<Element>();
		result.addAll(elements);
		return result;
	}
	
	public String getDescription()	{
		return description;
	}
		
	public static void writePoscar(Cell c, String outPoscar, Boolean useCartesianCoords) {
		(new VaspIn(c, null, null, null)).writePoscar(outPoscar, useCartesianCoords);
	}
		
	public void writePoscar(String outPoscar, Boolean useCartesianCoords) {
		/* Note that this will overwrite outPoscar. If this is not desired behavior,
		 * the calling routine should check if the output file already exists */
		
		/* Write the output */
		try {
			NumberFormat nf = NumberFormat.getInstance();
			nf.setMaximumFractionDigits(8);
			nf.setMinimumFractionDigits(8);
			
			BufferedWriter writer = new BufferedWriter(new FileWriter(outPoscar));
			writer.write(description + "\n");
			writer.write("1.0 \n");
			List<Vect> latticeVectors = cell.getLatticeVectors();
			for (Vect vec : latticeVectors) {
				List<Double> comps = vec.getCartesianComponents();
				for (int i = 0; i < Constants.numDimensions; i++) {
					String coordStr = nf.format(comps.get(i));
					writer.write(coordStr + " ");
				}
				writer.write("\n");
			}
			/* Write elements */
			for (Element e : elements)
				writer.write(e.getSymbol() + " ");
			writer.write("\n");
			/* Write number of each type of site */
			for (Element e : elements)
				writer.write(cell.getNumSitesWithElement(e) + " ");
			writer.write("\n");
			writer.write("Selective Dynamics\n"); // For specifying which atoms to allow to move
			if (useCartesianCoords)
				writer.write("Cartesian\n");
			else
				writer.write("Direct\n");
			/* Make sure we're printing these out in the right order */
			List<Site> basis = cell.getSites();	
			for (Element e : elements)
				for (Site s : basis) {
					if (s.getElement().equals(e)) {
						List<Double> coords = null;
						if (useCartesianCoords)
							coords = s.getCoords().getCartesianComponents();
						else
							coords = s.getCoords().getComponentsWRTBasis(cell.getLatticeVectors());
						for (int i = 0; i < Constants.numDimensions; i++) {
							String coordStr = nf.format(coords.get(i));
							writer.write(coordStr + " ");
						}
						
						// check to see if we should let this atom move or not
						if (s.getRelax())  
							writer.write("T T T");
						else
							writer.write("F F F");
						
						writer.write("\n");
					}
				}
			
			/* Write it all to disk */
			writer.flush();
		} catch (IOException x) {
			System.out.println("IOException in VaspOut.writePoscar(): " + x.getMessage());
		}
	}

	public void makeINCAR(String directory) {
		String incarStr = Utility.readStringFromFile(incarFile);
		Utility.writeStringToFile(incarStr, directory + "/INCAR");
		
	}

	public void makePOSCAR(String directory) {
		writePoscar(directory + "/POSCAR", false);
	}

	public void makeKPOINTS(String directory) {
		String kpointsStr = Utility.readStringFromFile(kpointsFile);
		Utility.writeStringToFile(kpointsStr, directory + "/KPOINTS");
	}

	public void makePOTCAR(String directory) {
		StringBuilder potcar = new StringBuilder();
		
		for (Element e : elements) {
//	System.out.println("potcar for " + e.getSymbol() + " : " + potcarFileMap.get(e));
			String potcarStr = Utility.readStringFromFile(potcarFileMap.get(e));
			potcar.append(potcarStr);
		}
		
		Utility.writeStringToFile(potcar.toString(), directory + "/POTCAR");
	}

}
