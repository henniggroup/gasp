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

import java.util.*;

import chemistry.Element;

import utility.Utility;
import utility.Vect;

import crystallography.Cell;
import crystallography.Site;

// Permutation is a Variation operation which creates an offspring from
// one parent by randomly swapping some of its atoms.  The number of swaps
// is normally distributed with mean and standard deviation taken by the
// function as arguments.

public class Permutation implements Variation {
	
	private double meanExchanges;
	private double sigmaExchanges;
	List<String> pairStrings;
	private List<Set<String>> pairs = null;
	
	public Permutation(List<String> args) {
		if (args.size() < 3)
			GAParameters.usage("Not enough parameters given to Permutation", true);
		
		meanExchanges = Double.parseDouble(args.get(0));
		sigmaExchanges = Double.parseDouble(args.get(1));
		
		if (meanExchanges == 0 && sigmaExchanges == 0)
			GAParameters.usage("Bad mean/sigma given to Permutation", true);
		
		pairStrings = Utility.subList(args, 2);
	}
	
	public String toString() {
		StringBuilder result = new StringBuilder();
		
		result.append("Permutation variation (" + meanExchanges + ", " + sigmaExchanges + ") swaps: ");
		for (String p : pairStrings) 
			result.append(p + " ");
		
		return result.toString();
	}

	// selects and copies an organism from parents, permutes numExchanges of its
	// sites, and returns it
	public Organism doVariation(Generation parents, Generation offspring, Selection sel) {
		GAParameters params = GAParameters.getParams();
		Random rand = params.getRandom();
		
		if (pairs == null)
			pairs = GAUtils.parsePairs(pairStrings);
		
		// get a parent and its basis
		StructureOrg p = (StructureOrg)(sel.doSelection(parents, 1)[0]);
		Cell pStruct = p.getCell();
		
		// copy the parent's vectors and sites
		List<Site> newSites = new ArrayList<Site>();
		List<Vect> newVects= new ArrayList<Vect>(pStruct.getLatticeVectors());
		for (Site s : pStruct.getSites()) 
			newSites.add(s);
		
		// find a (nonzero) number of exchanges to do (mean meanExchanges, standard dev. sigmaExchanges)
		int numExchanges;
		do {
			numExchanges = Math.round(Math.round(rand.nextGaussian()*sigmaExchanges + meanExchanges));
		} while (numExchanges == 0);
		
		// permute numExchanges of them at random:
		for (int i = 0; i < numExchanges; i++) {
			// select a pair of species
			Set<String> pair = pairs.get(rand.nextInt(pairs.size()));
			Iterator<String> j = pair.iterator();
			String firstSymbol = j.next();
			String secondSymbol = j.next();
			
			// find one site of each symbol at random. make sure we can first
			if (pStruct.getNumSitesWithElement(Element.getElemFromSymbol(firstSymbol)) == 0 
					|| pStruct.getNumSitesWithElement(Element.getElemFromSymbol(secondSymbol)) == 0)
				continue;
			int indexA, indexB;
			do {
				indexA = rand.nextInt(newSites.size());
			} while (!newSites.get(indexA).getElement().getSymbol().startsWith(firstSymbol));
			do {
				indexB = rand.nextInt(newSites.size());
			} while (!newSites.get(indexB).getElement().getSymbol().startsWith(secondSymbol));
			
			// swap the species (but not the points)
			Element elemA = newSites.get(indexA).getElement();
			Element elemB = newSites.get(indexB).getElement();
			newSites.set(indexA, new Site(elemB, newSites.get(indexA).getCoords()));
			newSites.set(indexB, new Site(elemA, newSites.get(indexB).getCoords()));

			// some output
			GAOut.out().stdout("Permuting " + newSites.get(indexA).getElement() + " and " + newSites.get(indexB).getElement(), GAOut.INFO);

		}
		
		StructureOrg newOrg = new StructureOrg(new Cell(newVects, newSites));
		
		// if we're using the island objective function, preserve the parent's location and interlayer distance
		if (params.getObjFcnArgs().get(0) == "island") {
			newOrg.setLocation(p.getLocation());
			newOrg.setInterlayerDist(p.getInterlayerDist());
		}
		
		return newOrg;	
	}
}
