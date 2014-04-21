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

package pdvisual;

import java.util.*;

import chemistry.*;
import crystallography.*;
import com.sun.org.apache.bcel.internal.classfile.Utility;

public class PseudoPDData implements IPDPlottable {
	
	// all the phases in the pseudopddata
	private List<MixedPhase> allPhases;
	private List<List<Integer>> facets;
	private List<Integer> basisIndxs;
	private List<Integer> indxUnstableEntries;
	
	PDData m_pdd;
	
	public PseudoPDData(List<MixedPhase> _allPhases, List<Integer> b, PDData p) {	
		allPhases = _allPhases;
		basisIndxs = b;
		m_pdd = p;
		
	// do we really want to do this?  i'd rather filter them out in the visualization. --will
	//	removeEntriesWPositiveFE();
	}
	
	public PDData getPDD() {
		return m_pdd;
	}
	
	public List<Integer> getPhasesWSameCompAs(int indxEntry) {
		List<Integer> result = new LinkedList<Integer>();
		
		for (int i = 0; i < allPhases.size(); i++)
			if (allPhases.get(i).getComposition().equals(allPhases.get(indxEntry).getComposition()))
				result.add(i);
//		for now we'll just sort the results by energy per atom
		List<Integer> sorted=new LinkedList<Integer>();
		
		Integer chosenOne=null;
		while(result.size()!=0)
		{
			double energymin=Double.MAX_VALUE;
			for(Integer i:result)
			{
				if(energymin>allPhases.get(i).getEnergyPerAtom())
				{
					energymin=allPhases.get(i).getEnergyPerAtom();
					chosenOne=i;
					System.out.println(chosenOne);
				}
			}
			result.remove(chosenOne);
			sorted.add(chosenOne);
		}
	
		return sorted;

	}

    public List<Element> getElements(){
        throw new RuntimeException("not yet implemented!");
    }
	
	public int getNumEntries() {
		return allPhases.size();
	}
	
	public List<Composition> getCompositions() {
		List<Composition> result = new LinkedList<Composition>();
		
		for (MixedPhase m : allPhases)
			result.add(m.getComposition());
		
		return result;
	}
	
	private void removeEntriesWPositiveFE() {
		// remove entries with positive formation energies
		List<Double> fEnergies = getFormEnergiesPerAtom();
		for (int i = fEnergies.size() - 1; i >= 0; i--)
			if (fEnergies.get(i) > 0)
				allPhases.remove(i);
	}
	
	// returns the dimension of the plottable pd object
	public int getDimension() {
		return basisIndxs.size();
	}
	
	public List<MixedPhase> getAllPhases() {
		// return a copy of the phases list
		return new LinkedList<MixedPhase>(allPhases);
	}
	
	public List<List<Integer>> getIndxFacets() {
		return facets;
	}
	
	public void setFacets(List<List<Integer>> fs) {
		facets = fs;
	}
	
	// composition fractions in terms of the basis IComputedEntries' compositions (NORMALIZED)
	public List<Double[]> getCompositionFractions() {
		List<Double[]> coords = new LinkedList<Double[]>();

		for (MixedPhase m : allPhases)
			coords.add(getCompositionFraction(m));

		return coords;
	}
	
	public Double[] getCompositionFraction(MixedPhase m) {
		Double[] unnormFrac = getPseudoCompositionCoefs(m);
		
		for (int i = 0; i < unnormFrac.length; i++) {
			unnormFrac[i] *= allPhases.get(basisIndxs.get(i)).getComposition().getNumComponents();
		}
		
		double sum = 0.0;
		for (Double d : unnormFrac)
			sum += d;
		for (int i = 0; i < unnormFrac.length; i++)
			unnormFrac[i] /= sum;
		
		return unnormFrac;
	}
	
	public List<Double> getEnergies() {
		List<Double> energies = new LinkedList<Double>();

		for (MixedPhase m : allPhases)
			energies.add(m.getEnergy());

		return energies;
	}
	
	public List<Double> getFormEnergiesPerAtom() {
		// NB: can't immediately shortcircuit this calculation b/c of its use in constructor
		List<Double> energies = getEnergiesPerAtom();
		List<Double[]> compFracs = getCompositionFractions();
		List<Double> fEnergies = new LinkedList<Double>();
		
		for(int i = 0; i < allPhases.size(); i++) {
			Double e = new Double(energies.get(i));
			Double[] cFrac = compFracs.get(i);
			for (int j = 0; j < cFrac.length; j++)
				e -= (cFrac[j])*allPhases.get(basisIndxs.get(j)).getEnergyPerAtom();
//			e /= allPhases.get(i).getComposition().getOrigAmountSum();
			fEnergies.add(e);
		}
		
		return fEnergies;
	}
	
	// returns UNNORMALIZED coefficients of composition of a MixedPhase in terms of the basis of
	// MixedPhases from this PseudoPDData
	public Double[] getPseudoCompositionCoefs(MixedPhase m) {
		List<MixedPhase> basis = new LinkedList<MixedPhase>();
		for (Integer i : basisIndxs)
			basis.add(allPhases.get(i));

		Double fracs[] =  MixedPhase.makeCompoFromBasis(m.getComposition(), basis, m_pdd.getElementsArray());

		return fracs;		
	}
	
	public List<String> getAxisLabels() {
		List<String> result = new LinkedList<String>();
		
		Iterator<Integer> i = basisIndxs.iterator();
		while (i.hasNext()) {
			Composition c = allPhases.get(i.next()).getComposition();
			result.add(c.toString());
		}
		
		return result;
	}
	
	public boolean isPurePhase(int i) {
		MixedPhase m = getAllPhases().get(i);
		
		Map<IComputedEntry, Double> phases = new HashMap(m.getPhases());
		Iterator<IComputedEntry> iter = phases.keySet().iterator();
		while (iter.hasNext()) {
			IComputedEntry e = iter.next();
			if (phases.get(e) == 0)
				iter.remove();
		}			
		
		// it's a pure entry if it consists of exactly one icomputedentry
		return phases.keySet().size() == 1;
	}
	
	public List<String> getVertexLabels() {
		List<String> result = new LinkedList<String>();
		
		for (MixedPhase m : allPhases) {
			//result.add(com.cmcweb.util.Utility.getNormalizedFormula(m.getComposition().getformulasum()));
			result.add(m.getComposition().toString());
		}
			//result.add(m.getComposition().toString());
		
		
		return result;
	}
	
	public List<Integer> getIndxUnstableEntries() {
		// start out with a list of all entries, and remove ones which lie on a facet
		if (indxUnstableEntries == null) {
			indxUnstableEntries = new LinkedList<Integer>();
			
			for (int i = 0; i < allPhases.size(); i++)
				indxUnstableEntries.add(i);
			for (List<Integer> f : facets)
				for (Integer i : f)
					indxUnstableEntries.remove(i);
		}
		
		return indxUnstableEntries;
	}
	
	public Composition getComposition(int i) {
		return getAllPhases().get(i).getComposition();
	}

    /**
     * @return
     *
     * @see com.cmcweb.analysis.IPDPlottable#getPerAtomEnergies()
     */
    public List<Double> getEnergiesPerAtom() {
        List<Double> ret = new ArrayList<Double>();
        for(MixedPhase m : allPhases){
            ret.add(m.getEnergyPerAtom()); 
        }
        return ret;
    }
    
    public List<Integer> getBasisEntryIndxs() {
		return basisIndxs;
    }

}
