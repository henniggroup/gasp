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

import java.util.List;

import chemistry.Composition;
import chemistry.Element;


public interface IPDPlottable {
	
	// returns the dimension of the plottable pd object
	public int getDimension();
	public int getNumEntries();
    public List<Element> getElements();
	public List<Integer> getIndxUnstableEntries();
	public List<Double> getEnergiesPerAtom(); 
	public List<Double> getFormEnergiesPerAtom();
	public List<List<Integer>> getIndxFacets();
	public List<String> getVertexLabels();
	public List<String> getAxisLabels();
	public boolean isPurePhase(int i);
	// returns the indices of the "basis" entries, e.g. the elements in a normal phase diagram
	public List<Integer> getBasisEntryIndxs();
	// get (REAL) entrie with same composition as entry i
	public List<Integer> getPhasesWSameCompAs(int i);
	public Composition getComposition(int i);
	public List<Composition> getCompositions();
	public List<Double[]> getCompositionFractions();
	

}
