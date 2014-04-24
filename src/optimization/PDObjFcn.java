
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

package optimization;

import ga.*;

import java.util.*;
import pdvisual.*;
import utility.*;
import crystallography.*;
import chemistry.*;


public class PDObjFcn extends ObjectiveFunction {
	
	ObjectiveFunction energyFcn;
	
	StructureOrg org;
	
/*	public PDObjFcn(ObjectiveFunction _objFcn, List<Element> elements, PDData _pddata) {
		List<ComputedEntry> entries = new LinkedList<ComputedEntry>();
		Map<Element,Double> chempots = new HashMap<Element,Double>();
		pdbuilder = new PDBuilder(entries, elements, chempots);
	} */
	
	public PDObjFcn(List<String> subArray, Organism o) {
		// TODO Auto-generated constructor stub
		energyFcn = new EnergyPerAtom(subArray, o);
		org = (StructureOrg)o;
	}

	// TODO: remember to update numCalculations when implementing this
	public Thread evaluate() {		
		// short circuit here if we've done the calculation already
		if (org.knowsValue()) // TODO: have to update due to changing phase diagram ???
			return null;
		
		// another total energy calculation:
		// not necessary to increment here since it's done when we use EnergyPerAtom under the covers
	//	numCalculations++;
		
		// start the calculation and return the Thread
		Thread t = new Thread(this);
		t.start();
		return t;
	}
	
	public void run() {
		
		// relax the cell - have to wait for it to finish before using results
		Thread t = energyFcn.evaluate();
		try {
				t.join();
		} catch (InterruptedException x) {
			if (GAParameters.getParams().getVerbosity() >= 3)
				System.out.println("InterruptedException in energy calc thread in PDObjFcn: " + x.getMessage());
		}
		
		// updating structure w/ relaxed version is this done in EPA already
		org.setValue((new PDAnalyzer(GAParameters.getParams().getPDBuilder().getPDData())).getEnergyPerAtomAboveHull(org));
		
	}
	
	// an ObjectiveFunction should also overload toString();
	public String toString() {
		return "PDObjFcn"; // TODO: more? what does cellobjfcn do?
	}
}
