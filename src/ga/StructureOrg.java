/* Genetic algorithm for crystal structure prediction.  Will Tipton.  Ceder Lab at MIT. Summer 2007 */

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
		structure = s;
		// changing the structure invalidates the energy and fitness and satisfiesConstraints
		value = null;
		fitness = null;
	}
	
	public void standardize() {
		structure = structure.getNigliReducedCell();
		reduced = true;
	}
	
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
