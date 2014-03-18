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

package randomsearch;

import java.util.*;

import chemistry.Element;
import crystallography.*;
import optimization.OptiSystem;
import pdvisual.*;
import utility.ArgumentParser;
import utility.Pair;

public class Scratch {

	public static void main(String[] args) {

		ArgumentParser aparser = new ArgumentParser(args);
		
		OptiSystem sys = new OptiSystem(args);
		
		// make some cells
		List<Cell> cells = new LinkedList<Cell>();
		for (int i = 0; i < 150; i++)
			cells.add(RandomSearchMethod.getRandomCell(sys, i));
		
		// make some entries
		List<ComputedEntry> entries = new LinkedList<ComputedEntry>();
		for (Cell c : cells) {
			Pair<Cell,Double> result = sys.getObjFcn().getCellAndValue(c);
			entries.add(new ComputedEntry(result.getFirst(), result.getSecond()));
		}
		
		for (ComputedEntry c : entries)
			System.out.println("Entry " + c.getLabel() + " has epa " + c.getEnergyPerAtom());

		
		List<Element> elements = sys.getCompositionSpace().getElements();
		Map<Element, Double> chempots = new HashMap<Element,Double>();
		for (Element e : elements) {
			chempots.put(e, 0.0);
		}
		
		PDBuilder pdbuilder = new PDBuilder(entries, elements, chempots);
		
		PDData pddata = pdbuilder.getPDData();
		
		TernPD3D pdviz = new TernPD3D(pddata);
		
	}

}
