/* Genetic algorithm for crystal structure prediction.  Will Tipton.  Ceder Lab at MIT. Summer 2007 */

package ga;

import java.util.*;

import chemistry.Element;

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
	String[] pairStrings;
	private List<Set<String>> pairs = null;;
	
	public Permutation(String[] args) {
		if (args.length < 3)
			GAParameters.usage("Not enough parameters given to Permutation", true);
		
		meanExchanges = Double.parseDouble(args[0]);
		sigmaExchanges = Double.parseDouble(args[1]);
		
		pairStrings = GAUtils.subArray(args, 2);
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
		int verbosity = params.getVerbosity();
		Random rand = params.getRandom();
		
		if (pairs == null)
			pairs = GAUtils.parsePairs(pairStrings);
		
		// get a parent and its basis
		Cell pStruct = ((StructureOrg)(sel.doSelection(parents, 1)[0])).getCell();
		
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
			// find one site of each symbol at random
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
			if (verbosity >= 4)
				System.out.println("Permuting " + newSites.get(indexA).getElement() + " and " + newSites.get(indexB).getElement());
		}
		return new StructureOrg(new Cell(newVects, newSites));	
	}
}
