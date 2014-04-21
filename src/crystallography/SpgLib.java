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

package crystallography;

import java.io.File;
import java.util.*;

import chemistry.Element;

import spglib.*;
import utility.Constants;
import utility.Vect;
import ga.StructureOrg;

// singleton
public class SpgLib {
	
	private static final double symprec = 0.001; // TODO: dunno if this is a reasonable number
	
	static
	{
	    System.loadLibrary("spglib");
	    System.loadLibrary("symspg");
	}
	
	private static SpgLib instance = null;
	
	private SpgLib() {
		//System.loadLibrary("spglib");
	}
	
	public static SpgLib getInstance() {
		if (instance == null)
			instance = new SpgLib();
		return instance;
	}

//	private native int spg_find_primitive(double lattice[][], double position[][],
//            int types[], int num_atom, double symprec);
	
	public int getSpaceGroup(Cell c) {
		int num_atoms = c.getBasisSize();
		
		SWIGTYPE_p_a_3__double lvects = spglib.getVectsArray();
		for (int i= 0; i < Constants.numDimensions; i++)
			for (int j= 0; j < Constants.numDimensions; j++)
				spglib.putInDoubleNby3Array(lvects, i, j, c.getLatticeVectors().get(i).getCartesianComponents().get(j));

		SWIGTYPE_p_a_3__double positions = spglib.getPositionsArray(num_atoms);
		for (int i= 0; i < num_atoms; i++)
			for (int j= 0; j < Constants.numDimensions; j++)
				spglib.putInDoubleNby3Array(positions, i, j, c.getSite(i).getCoords().getComponentsWRTBasis(c.getLatticeVectors()).get(j));
		
		SWIGTYPE_p_int types = spglib.getTypesArray(num_atoms);
		for (int i = 0; i < num_atoms; i++)
			spglib.putInIntArray(types, i, c.getComposition().getElements().indexOf(c.getSite(i).getElement()));
		
		String symbol = "12345678910";
		int result = spglib.spg_get_international(symbol, lvects, positions, types, num_atoms, symprec);
		
		spglib.deleteDoubleNby3Array(lvects);
		spglib.deleteIntArray(types);
		spglib.deleteDoubleNby3Array(positions);
		
		return result;
	}
	
	public Cell getPrimitiveCell(Cell c) {
		int num_atoms = c.getBasisSize();
		List<Element> elements = c.getComposition().getElements();
				
		SWIGTYPE_p_a_3__double lvects = spglib.getVectsArray();
		for (int i= 0; i < Constants.numDimensions; i++)
			for (int j= 0; j < Constants.numDimensions; j++)
				spglib.putInDoubleNby3Array(lvects, i, j, c.getLatticeVectors().get(i).getCartesianComponents().get(j));

		SWIGTYPE_p_a_3__double positions = spglib.getPositionsArray(num_atoms);
		for (int i= 0; i < num_atoms; i++)
			for (int j= 0; j < Constants.numDimensions; j++)
				spglib.putInDoubleNby3Array(positions, i, j, c.getSite(i).getCoords().getComponentsWRTBasis(c.getLatticeVectors()).get(j));
		
		SWIGTYPE_p_int types = spglib.getTypesArray(num_atoms);
		for (int i = 0; i < num_atoms; i++)
			spglib.putInIntArray(types, i, elements.indexOf(c.getSite(i).getElement()));
		
		int num_prim_atoms = spglib.spg_find_primitive(lvects, positions, types, num_atoms, symprec);
		
/*		System.out.println(num_prim_atoms);
		
		for (int i = 0; i < 3; i++) {
			for (int j = 0; j < 3; j++)
				System.out.print(spglib.getFromDoubleNby3Array(lvects, i, j) + "	");
			System.out.println("");
		}
		System.out.println("");
		for (int i = 0; i < num_prim_atoms; i++) {
			System.out.print(spglib.getFromIntArray(types, i) + " ");
			for (int j = 0; j < 3; j++)
				System.out.print(spglib.getFromDoubleNby3Array(positions, i, j) + "	");
			System.out.println(" ");
		}
		System.out.println("");
		
		spglib.spg_show_symmetry(lvects, positions, types, num_atoms, symprec);
*/
		
		// if the function returned 0, it didn't find a primitive cell, or the current cell is one already
		if (num_prim_atoms == 0)
			return c;
		
		List<Vect> vects = new LinkedList<Vect>();
		vects.add(new Vect(	spglib.getFromDoubleNby3Array(lvects, 0, 0), 
							spglib.getFromDoubleNby3Array(lvects, 0, 1), 
							spglib.getFromDoubleNby3Array(lvects, 0, 2)));
		vects.add(new Vect(	spglib.getFromDoubleNby3Array(lvects, 1, 0), 
							spglib.getFromDoubleNby3Array(lvects, 1, 1), 
							spglib.getFromDoubleNby3Array(lvects, 1, 2)));
		vects.add(new Vect(	spglib.getFromDoubleNby3Array(lvects, 2, 0), 
							spglib.getFromDoubleNby3Array(lvects, 2, 1), 
							spglib.getFromDoubleNby3Array(lvects, 2, 2)));
		
		List<Site> sites = new LinkedList<Site>();
		for (int i = 0; i < num_prim_atoms; i++) {
			int type = spglib.getFromIntArray(types, i);
			sites.add(new Site(elements.get(type),
					  new Vect(	spglib.getFromDoubleNby3Array(positions, i, 1),
							  	spglib.getFromDoubleNby3Array(positions, i, 2),
								spglib.getFromDoubleNby3Array(positions, i, 3))));
		}
		
		spglib.deleteDoubleNby3Array(lvects);
		spglib.deleteIntArray(types);
		spglib.deleteDoubleNby3Array(positions);
		
		return new Cell(vects, sites, c.getLabel());
	}

	// for testing
	public static void main(String [] args) {
		SpgLib s = SpgLib.getInstance();
		
		//Cell c = StructureOrg.parseCif(new File("/home/wtipton/projects/ga_for_crystals/test_runs/alcu_compspace/refstates/77.cif"));
		Cell c = Cell.parseCif(new File("/home/wtipton/POSCAR4.cif"));
		
		System.out.println(c);
		
		System.out.println(s.getPrimitiveCell(c));
	}
	
}
