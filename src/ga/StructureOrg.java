<<<<<<< HEAD
/* Genetic algorithm for crystal structure prediction.  Will Tipton.  Ceder Lab at MIT. Summer 2007 */
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

package ga;

import java.io.*;
import java.util.*;

import pdvisual.IComputedEntry;

import chemistry.Composition;
import chemistry.Element;

import utility.Vect;

import crystallography.Cell;
import crystallography.Site;

// StructureOrg extends Organism and has a Structure. It is the type of
// Organism generally used by the genetic algorithm for the crystal structure
// prediction problem.

public final class StructureOrg extends Organism implements IComputedEntry, Serializable {
	
	static final long serialVersionUID = 1l;
	
	private Cell structure;
	private Boolean reduced = false;
	private StructureOrgCreator SOCreator = null;
	
	private double totalEnergy = Double.NaN;
	
	public StructureOrg(Cell s) {
		structure = s;
	}
	
	public Cell getCell() {
		return structure;
	}
	
	
	public void setCell(Cell s) {
		setCell(s, true);
	}
		
	public void setCell(Cell s, boolean invalidateEnergies) {
		structure = s;
		// changing the structure invalidates the energy and fitness and satisfiesConstraints
		if (invalidateEnergies) {
			value = null;
			fitness = null;
		}
	}
	
	public void standardize() {
		structure = structure.getCellWithAllAtomsInCell().getNigliReducedCell();
		reduced = true;
	}
	
<<<<<<< HEAD
	public void standardize2D() {
		structure = structure.getCellWithAllAtomsInCell().getNigliReduced2DCell();
		reduced = true;
	}
	
=======
>>>>>>> 0e3189c40547bbd59ea42c4f91890d7511fb7797
	public double getTotalEnergy() {
		return totalEnergy;
	}
	
	public void setTotalEnergy(double t) {
		totalEnergy = t;
	}
	
	public String toString() {
		StringBuilder result = new StringBuilder();
		String newline = GAUtils.newline();
		
		result.append("StructureOrg " + getID() + ":" + newline);
		result.append("Total Energy: " + totalEnergy + " Value: " + value + newline);
		result.append(structure.toString());
		
		return result.toString();
	}
	
	public Boolean isReduced() {
		return reduced;
	}
	
	// to test niggliReduceStructure
	public static void main(String[] args) throws IOException {
/*		Cell bobStruct = StructureOrg.parseCif(new File("/home/wtipton/projects/ga_for_crystals/testing/ljc/seeds18/105.cif"));
		String[] rgArgs = {"0.5", "0.5", "0.05", "true"};
		RedundancyGuard rg = new RedundancyGuard(rgArgs);
		
		double[][] ba = {{1,0,0},{0,2,0},{0,0,3}};
		double[][] ca = {{1,0,0},{0,1,0},{0,0,1}};
		Basis b = new Basis(ba);
		Basis cart = new Basis(ca);
		double[] po = {1,2,3};
		Point p = new Point(po, b);
		
		Structure s = new ExplicitStructure(b);
		s.addSiteChecked(s.new Site(p, new AtomicSpecies(ElementCollection.ELEMENTS.getElement("AA"),0.0)));
		
		
		System.out.println(s);
		*/
	
		//bobStruct = new ExplicitStructure(com.cmcweb.db.postgresql.StructureExtractor.getGivenCell(3,5.196,2,(103.0+55.0/60.0)*Math.PI/180.0,(109.0+28.0/60.0)*Math.PI/180.0,(134.0+53.0/60.0)*Math.PI/180.0));
/*		Vaspin vaspin = new Vaspin(bobStruct, "/home/wtipton");
		vaspin.makePOSCAR("/home/wtipton/");
		
		
		StructureOrg bob = new StructureOrg(bobStruct);
		rg.addStructureOrg(bob);
		bob = new StructureOrg(bobStruct.getNigliReducedCell());
		//bob.niggliReduceStructure2();
		GAUtils.writeStringToFile(bob.getCIF(), new File("bob.cif"), false);

		System.out.println(rg.checkStructureOrg(bob)); */
	}
	
	public StructureOrgCreator getSOCreator() {
		return SOCreator;
	}
	
	public void setSOCreator(StructureOrgCreator s) {
		SOCreator = s;
	}

	@Override
	public double getEnergyPerAtom() {
		return getTotalEnergy() / structure.getNumSites();
	}

	@Override
	public Composition getComposition() {
		return structure.getComposition();
	}

	@Override
	public String getLabel() {
	//	return toString();
		return Integer.toString(getID());
	}
}
