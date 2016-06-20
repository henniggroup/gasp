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

package ga;

import java.io.File;
import java.util.*;

import utility.RandomNumbers;
import utility.Vect;

import chemistry.*;
import crystallography.*;

// NumStoichsMut implements Variation.  It creates an offspring StructureOrg from
// a parent by adding or removing a randomly chosen number of stoichiometries worth
// of species.

public class NumStoichsMut implements Variation {
	
	double meanNum;
	double sigmaNum;
	
	public NumStoichsMut(List<String> args) {
		if (args == null || args.size() < 2)
			GAParameters.usage("Not enough parameters given to NumStoichsMut", true);
		
		meanNum = Double.parseDouble(args.get(0));
		sigmaNum = Double.parseDouble(args.get(1));
	}
	
	private void addAtoms(int numAtoms, List<Site> newSites, List<Vect> newVects) {
		GAParameters params = GAParameters.getParams();
		Composition randComp = params.getCompSpace().getRandomIntegerCompInSpace(numAtoms,numAtoms);
		Random rand = params.getRandom();
		// add n atoms
		for (Element e : randComp.getElements()) {
			for (int i = 0; i < randComp.getOrigAmount(e); i++)
				// make the new point randomly and add them to the list
				newSites.add(new Site(e, new Vect(rand.nextDouble(),rand.nextDouble(),rand.nextDouble(), newVects)));
		}
	}

	private void removeAtoms(int numAtoms, List<Site> newSites) {
		GAParameters params = GAParameters.getParams();
		Random rand = params.getRandom();
		// remove atoms n times
		for (int i = 0; i < numAtoms; i++) {
			newSites.remove(RandomNumbers.getUniformIntBetweenInclusive(0, newSites.size() - 1));
		}
	}
	
	/*
	private Boolean hasAStoich(List<Site> speciesList, Map<Element,Double> constituents) {
		Set<Element> elements = constituents.keySet();
		Iterator<Element> i = elements.iterator();
		while (i.hasNext()) {
			Element el = i.next();
			double numToRemove = constituents.get(el);
			Iterator<Site> j = speciesList.iterator();
			while(j.hasNext() && numToRemove > 0) {
				if (!j.next().getElement().getSymbol().toLowerCase().startsWith(el.getSymbol().toLowerCase()))
					continue;
				numToRemove--;
			}
			if (!j.hasNext() && numToRemove > 0)
				return false;
		}
		return true;
	}*/

	public Organism doVariation(Generation parents, Generation offspring, Selection sel) {
		GAParameters params = GAParameters.getParams();
		Random rand = params.getRandom();
		
		// pick a parent randomly 
		StructureOrg p = (StructureOrg)(sel.doSelection(parents, 1)[0]);
		Cell pStruct = p.getCell();
		// copy the parent's vectors and sites
	//	Basis newBasis = new Basis(pStruct.getBasis());
		List<Site> newSites = new ArrayList<Site>();
		List<Vect> newVects= new ArrayList<Vect>(pStruct.getLatticeVectors());
		for (Site s : pStruct.getSites()) 
			newSites.add(s);
		
		// we'll add (or remove) abs(n) atoms to the structure.  force (n != 0)
		int n;
		do {
			n = Math.round(Math.round(meanNum + rand.nextGaussian()*sigmaNum));
		} while (n == 0 || n + pStruct.getBasisSize() > params.getMaxNumAtoms()
						|| n + pStruct.getBasisSize() < params.getMinNumAtoms());
		
		// some output:
		GAOut.out().stdout("Creating new StructureOrg by adding " + n + 
				" atoms to StructureOrg " + p.getID(), GAOut.NOTICE, p.getID());

		// add or remove species and their points, depending on the sign of n
		if (n > 0)
			addAtoms(n, newSites, newVects);
		else
			removeAtoms(Math.abs(n), newSites);	
		
		StructureOrg newOrg = new StructureOrg(new Cell(newVects, newSites));
		
		// if we're using the island objective function, preserve the parent's location and interlayer distance
		if (params.getObjFcnArgs().get(0) == "island") {
			newOrg.setLocation(p.getLocation());
			newOrg.setInterlayerDist(p.getInterlayerDist());
		}
		
		return newOrg;
	}

	// just for testing
	public static void main(String[] args) {
		/*
		StructureOrg s1 = new StructureOrg(Cell.parseCif(new File("/home/wtipton/cifs/143.cif")));
		GAParameters params = GAParameters.getParams();
		s1.setFitness(1);
		
		String[] gaArgs = {"--f", "/home/wtipton/test"};
		params.setArgs(gaArgs);
		
		String[] smArgs = {"0", "100"};
		NumStoichsMut s = new NumStoichsMut(smArgs);
		String[] selArgs = {"1", "0"};
		Selection sel = new ProbDistSelection(selArgs);

		Generation parents = new Structures();
		parents.addOrganism(s1);
		
		System.out.println(s1);
		
		StructureOrg o = (StructureOrg)s.doVariation(parents, null, sel);
	//	System.out.println(o.getCell().getNumSites());
	//	GAUtils.writeStringToFile(o.getCIF(), new File("offspring.cif"), false);
		System.out.println(o); */
	}
}
