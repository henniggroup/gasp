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
package crystallography;

import java.util.*;

import Jama.Matrix;

import chemistry.*;
import pdvisual.*;
import utility.Utility;


/**
 * Thu Jul 24 18:11:10 EDT 2008 Chris BUGFIX
 * Class needed update to fix 'general' problem with current "Entry" code.
 * Entry objects should be extensive (i.e., the energy of Ti2O4 should be larger
 * than TiO2, but somewhere along the way we converted to "per atom" quantities).
 * 
 * a mixed phase could be a decomposition such as the following 
 * 
 * ABO2 = AO + 1/3 B2O + 1/3 BO2
 * 
 * or it could be
 * 
 * 1/4 ABO2 = 1/4 ( AO + 1/3 B2O + 1/3 BO2 )
 * 
 * but we getEnergy() for these two to be different. Also, getNumAtoms() should
 * be different for the two. 
 * 
 * 
 * 
 * 
 * @author wtipton, Will Tipton 
 * @author ccfisch, Chris Fischer <ccfisch@gmail.com>
 * Date: Jul 24, 2008
 *
 */
public class MixedPhase {
	
	private Double energy = Double.NaN;
	private Double energyPerAtom = Double.NaN; 
	
	private Double volume = Double.NaN;
	private Double volumePerAtom = Double.NaN;
	
	private static double EPSILON = 0.000001;
	
	// we represent the mixed phase primarily as a mapping from the entries
	// which constitute it to the fractions of the mixed phase which they constitute
	/**
	 * Chris: 
	 * the fractions now does NOT have to sum to one.
	 */
	private Map<IComputedEntry, Double> phases;
	
	private Composition composition;
	
	public MixedPhase(Map<IComputedEntry, Double> _phases) {
		phases = _phases;
		
		// make the composition
		Map<Element, Double> comps = new TreeMap<Element, Double>();
		for (IComputedEntry entry : phases.keySet()) {
			Composition eComp = entry.getComposition();
			for (Element e : eComp.getElements())
				if (comps.containsKey(e)) {
					Double currentAmount = comps.get(e);
					Double amount = currentAmount + phases.get(entry) * eComp.getOrigAmount(e); 
					comps.put(e, amount);
				} else {
					comps.put(e, phases.get(entry) * eComp.getOrigAmount(e));
				}
		}
		Element elements[] = new Element[comps.keySet().size()];
		Double amounts[] = new Double[comps.keySet().size()];
		elements = comps.keySet().toArray(elements);
		for (int i = 0; i < elements.length; i++)
			amounts[i] = comps.get(elements[i]);
		composition = new Composition(elements, amounts);
	}
	
	/**
	 *  Constructor: makes a mixed phase with composition c out of 
	 *  a list of given entries
	 *  Takes a composition C and a list of compositions of length N: {C1,C2,...,CN}.
	 *  Computes a list of (positive) coefficients a1,a2,...,aN such that 
	 *     C = a1*C1 + a2*c2 + ... + aN*cN
	 *  We then normalize so that sum(a'_i) = 1.0 (Chris: WE DON'T DO THIS ANY LONGER)
	 *  note that a_i >= 0 
	 *  
	 */
	public MixedPhase(Composition c, List<IComputedEntry> origBasis) {
		composition = c;
		/**
		 * Sun Jul 20 14:02:26 EDT 2008 Chris BUGFIX
		 * You have to be careful here, because dim(c) might be less than
		 * max_{i=1,N} dim(c_i). In that case it might be possible to write
		 * C = a1*C1 + ... + aN*CN in an imbalanced way if you don't embed C in 
		 * the higher dimensional composition space. For example, if you are
		 * writing the composition Ti1 O2 as a linear combination of Ti3 O5, Ti1 O2,
		 * and Li1 Ti1 O2.
		 * 
		 * so, make sure the given C (and all C's for that matter) are embedded
		 * in the highest dimensional composition space in the origBasis
		 */
		// get Compositions
		List<MixedPhase> cBasis = new LinkedList<MixedPhase>();
		Set <Element> elemSet = new HashSet<Element>();
		
		for (IComputedEntry e : origBasis) {
			Map<IComputedEntry,Double> phaseMap = new HashMap<IComputedEntry,Double>();
			phaseMap.put(e, 1.0);
			cBasis.add(new MixedPhase(phaseMap));
			for (Element elem : e.getComposition().getElements())
				elemSet.add(elem);
	/*		if(e.getComposition().getElements().length > 
			   cMaxDim.getElements().length ){
			    cMaxDim = e.getComposition(); 
		}*/	
		}
		// now embed given c into possibly higher dim space
		Double[] amts = new Double[elemSet.size()];
		Element[] elements = new Element[elemSet.size()];
		int idx = 0; 
		for(Element e : elemSet) {
			amts[idx] = c.getOrigAmount(e);
			elements[idx++] = e;
		}
		composition = new Composition(elements, amts); 
		
		// make phases map, omitting phases that don't contribute
		phases = new HashMap<IComputedEntry, Double>();
		Double[] compFracs = makeCompoFromBasis(composition, cBasis, elements);
		for (int i = 0; i < origBasis.size(); i++)
//			if (compFracs[i] > PDAnalyzer.COMP_TOL)
				phases.put(origBasis.get(i), compFracs[i]);

		//normalizeFractions();
	}
	
	
	/**
	 * returns a list of coefficients (that sum to 1) indicating the linear combination of the given
	 * basis that makes up the given composition.  order to return value corresponds to order of given basis
	 * And not only that, but if there is more than one way to make up the given composition
	 * from the given basis, we need to return the one that gives the lowest energy.
	 * 
	 * Sun Jul 20 14:14:42 EDT 2008 Chris COMMENT
	 * 
	 * An example or two of a case where it is necessary to iterate over 'more than
	 * one way to make up a given composition' would be nice -- I'm not sure when
	 * one would encounter such a situation, but I suppose the 'best' among them
	 * is the combination with the lowest energy
	 *  Will: it happens whenever the input basis is degenerate.
	 * 
	 * CF: note this method signature seems to depend on @see PseudoPDData -- 
	 *     I'm not sure why you have to supply a List<MixedPhase> ? 
	 *   Will: We can't simply take a list a Compositions because of the degenerate
	 *      case where there's more than one linear combination of the basis compounds
	 *      that can be used to make up the desired composition.  In this case, we
	 *      choose from all the possibilities that of the lowest energy.  To have an
	 *      energy we need a MixedPhase, as opposed to just a Composition.
	 */
	public static Double[] makeCompoFromBasis(Composition c, List<MixedPhase> origBasis, Element units[]) {
        // make a column vector representing the amount of elements in
        // our target composition    
	    double[] Cvals = new double[units.length];
        for (int i = 0; i < units.length; i++)
            Cvals[i] = c.getOrigAmount(units[i]); 

		List<List<MixedPhase>> basisPowerSet = Utility.getPowerSet(origBasis);
		
		Double result[] = null;
		
		List<Double[]> possibleResults = new LinkedList<Double[]>();

		for (List<MixedPhase> basis : basisPowerSet) {
			try {		
				result = new Double[origBasis.size()];
				int nCols = basis.size();
				int nRows = units.length;

				double[][] Cnvals = new double[nRows][nCols];
				for (int i = 0; i < nRows; i++)
					for (int j = 0; j < nCols; j++)
						// fill in the matrix with the composition fraction of each element
						// in each basis entry
						Cnvals[i][j] = basis.get(j).getComposition().getOrigAmount(units[i]);
				// make the matrix m a Jama.Matrix Cn
				Matrix Cn = new Matrix(Cnvals);
				
				// and make it a matrix, C
				Matrix C = new Matrix(Cvals, nRows);
				
				// check if the matrix is singular
				if ((Cn.transpose().times(Cn)).det() < EPSILON)
					throw new RuntimeException("Singular matrix.");
				/**
				 * Sun Jul 20 15:58:19 EDT 2008: Chris (Update)
				 * Cn.solve(C) will default to a least squares solution when 
				 * rank(C) < nRows; then check residual
				 */
				// then, our solution a satisfies Cn*a = C.  so, a = (Cn^T Cn)^-1 Cn^T * C
				//Matrix a = Cn.transpose().times(Cn).inverse().times(Cn.transpose()).times(C);
				Matrix a = Cn.solve(C);
				double resid = Cn.times(a).minus(C).normInf(); //picks off maximum row
				if(resid > EPSILON){
				    throw new RuntimeException("basis insufficient to express composition"); 
				}
				
				// zero out the result array, then fill it in
				for (int i = 0; i < result.length; i++)
					result[i] = 0.0;
				for (int i = 0; i < nCols; i++) {
					double d = a.get(i, 0);
					if (Math.abs(d) < EPSILON)
						d = 0;
					else if (Math.abs(d-1) < EPSILON)
						d = 1;
					else if (d < 0)
						throw new RuntimeException("Negative composition coefficient");
					for (int j = 0; j < origBasis.size(); j++)
						if ((Object)(basis.get(i)) == ((Object)origBasis.get(j)))
							result[j] = d;
				}
				
				possibleResults.add(result);
			} catch(RuntimeException x) {
				//
			    //System.out.println(x); 
			}			
		}
		
		if (possibleResults.isEmpty())
			throw new RuntimeException("Can't make composition out of given basis");
		// return the possible result with the lowest energy
		Double lowestEnergy = Double.POSITIVE_INFINITY;
		for (Double possibleResult[] : possibleResults) {
			Double resultEnergy = 0.0;
			for (int i = 0; i < origBasis.size(); i++)
				resultEnergy += possibleResult[i] * origBasis.get(i).getEnergy();
			if (resultEnergy < lowestEnergy) {
				lowestEnergy = resultEnergy;
				result = possibleResult;
			}
		}
		if (lowestEnergy == Double.POSITIVE_INFINITY)
			System.out.println("ERROR: lowestEnergy == Double.POSITIVE_INFINITY in MixedPhase.makeCompFromBasis indicates something is wrong with structures' energies.");
		return result;
	}
	
	// TODO: do we really want to (compute?) return the amount of each phase in terms
	//		 of number of compositions worth of stuff? by mass?
	public Map<IComputedEntry, Double> getPhases() {
		return phases;
	}
	
	public Double getEnergy() {
		// if we've computed it before, don't do it again
		if (energy.isNaN()) {	
			// the energy is a linear combination of energies of its constituents:
			energy = 0.0;
			for (IComputedEntry e : phases.keySet())
				energy += e.getTotalEnergy() * phases.get(e);
		}
		
		return energy;
	}
	
	public Double getVolume() {
	    if(volume.isNaN()){
	        volume = 0.0;
	        for(IComputedEntry e : phases.keySet()){
	            volume += e.getCell().getVolume() * phases.get(e); 
	        }
	    }
	    return volume; 
	}
	
	public Double getEnergyPerAtom(){
		// if we've computed it before, don't do it again
	/*	if (energyPerAtom.isNaN()) {	
			// the energy is a linear combination of energies of its constituents:
			energyPerAtom = 0.0;
			for (ComputedEntry e : phases.keySet())
				energyPerAtom += (e.getEnergyPerAtom()) * (phases.get(e)/getOrigAmountSum());
		}
		
		return energyPerAtom; 
		*/
		
		
		double numAtoms = 0.0;
		for (IComputedEntry e : phases.keySet())
			numAtoms += e.getComposition().getNumComponents() * phases.get(e);
			
		return getEnergy() / numAtoms;
	
/*		double epa2 = 0.0;
		double sum = 0.0;
		for (ComputedEntry e : phases.keySet()) {
			epa2 += e.getEnergyPerAtom()* phases.get(e)*e.getComposition().getOrigAmountSum();
			sum += phases.get(e)*e.getComposition().getOrigAmountSum();
		}
		epa2 /= sum;
*/		
	//	System.out.println("should have " + getEnergy()/numAtoms + " = " + epa2);
		
	//	return epa2;
	}

	public Double getVolumePerAtom(){
	    if(volumePerAtom.isNaN()){
            double numAtoms = 0.0;
            for (IComputedEntry e : phases.keySet()){
                numAtoms += e.getComposition().getNumComponents() * phases.get(e);
            }
            volumePerAtom = getVolume() / numAtoms;
	    }
	    return volumePerAtom; 
	}
	
	// make a string representation of our mixed phase by printing entries
	// and the stoichiometry
	public String toString() {
		StringBuilder result = new StringBuilder();
		
		result.append("MixedPhase: Energy = ");
		result.append(getEnergy().toString());
		result.append("\n");
		for (IComputedEntry e : phases.keySet()) {
			result.append(phases.get(e));
			result.append(": ");
			result.append(e.getComposition());
			//result.append(e.toString());
			result.append("\n");
		}
		
		return result.toString();
	}
	
	public Composition getComposition() {
		return composition;
	}
	
}
