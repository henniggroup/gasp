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

import chemistry.*;
import ga.ObjectiveFunction;
import crystallography.*;
import utility.*;
import java.util.*;

public class OptiSystem {
	
	ObjectiveFunction objFcn;
	
	// parameters
	private double minLatticeLength;
	private double maxLatticeLength;
	private double minInteratomicDistance;
	private double minLatticeAngle;
	private double maxLatticeAngle;
	private int minNumAtoms;
	private int maxNumAtoms;
	private boolean optimizeCell;
	private boolean optimizeSites;
	
	private boolean writeTempFiles;
	private String outDir;
	
	// convergence criteria
	private int numStructures;
	
	
	//TODO: fix when optisystem becomes non-fixed composition
	CompositionSpace compositionSpace;
	
	public OptiSystem(String args[]) {
		// parse input args
		ArgumentParser aparser = new ArgumentParser(args);
		
		// parse out the composition: //TODO: fix when OptiSystem becomes non-fixed composition
		//	first, some input validation
	/*	if (!aparser.hasOption("composition")) 
			throw new IllegalArgumentException("No --composition option passed.");
		List<String> compTokens = aparser.getArguments("composition");
		if (compTokens.size() % 2 != 0 || compTokens.size() == 0)
			throw new IllegalArgumentException("Malformed composition passed.");
		//  now, read off element, amount pairs
		int numElements = compTokens.size() / 2;
		Map <Element, Double> components = new HashMap<Element, Double>();
		for (int i = 0; i < numElements; i++) {
			Element newElem = Element.getElemFromSymbol( compTokens.get(2 * i));
			double num = Double.parseDouble(compTokens.get(2 * i + 1));
			components.put(newElem, num);
		}
		composition = new Composition(components);
		*/
		if (!aparser.hasOption("compositionSpace")) 
			throw new IllegalArgumentException("No --compositionSpace option passed.");
		compositionSpace = new CompositionSpace(aparser.getArguments("compositionSpace"), false);
		
		// parse out some basic parameters
		// first, set some defaults
		minLatticeLength = 3;
		maxLatticeLength = 15;
		minInteratomicDistance = 0.8;
		minLatticeAngle = Math.PI / 5;
		maxLatticeAngle = 4 * Math.PI / 5;
		minNumAtoms = 10;
		maxNumAtoms = 20;
		optimizeSites = true;
		optimizeCell = true;
		if (aparser.hasArguments("minLatticeLength"))
			minLatticeLength = Double.parseDouble(aparser.getArgument("minLatticeLength"));
		if (aparser.hasArguments("maxLatticeLength"))
			maxLatticeLength = Double.parseDouble(aparser.getArgument("maxLatticeLength"));
		if (aparser.hasArguments("minInteratomicDistance"))
			minInteratomicDistance = Double.parseDouble(aparser.getArgument("minInteratomicDistance"));
		if (aparser.hasArguments("minLatticeAngle"))
			minLatticeAngle = Double.parseDouble(aparser.getArgument("minLatticeAngle")) * Math.PI/180;
		if (aparser.hasArguments("maxLatticeAngle"))
			maxLatticeAngle = Double.parseDouble(aparser.getArgument("maxLatticeAngle")) * Math.PI/180;
		if (aparser.hasArguments("minNumAtoms"))
			minNumAtoms = Integer.parseInt(aparser.getArgument("minNumAtoms"));
		if (aparser.hasArguments("maxNumAtoms"))
			maxNumAtoms = Integer.parseInt(aparser.getArgument("maxNumAtoms"));
		if (aparser.hasArguments("optimizeCell"))
			optimizeCell = Boolean.parseBoolean(aparser.getArgument("optimizeCell"));
		if (aparser.hasArguments("optimizeSites"))
			optimizeSites = Boolean.parseBoolean(aparser.getArgument("optimizeSites"));
		
		// convergence criteria
		numStructures = 100;
		if (aparser.hasArguments("numStructures"))
			numStructures = Integer.parseInt(aparser.getArgument("numStructures"));
		
		// parse out the objective function
		//		first, some input validation
		if (!aparser.hasOption("objfcn"))
			throw new IllegalArgumentException("No --ObjFcn option passed.");
		List<String> objfcnTokens = aparser.getArguments("objfcn");
		String objFcnType = objfcnTokens.get(0);
		if (objfcnTokens.size() < 4)
			throw new IllegalArgumentException("Malformed objfcn passed.");
		if (objFcnType.compareToIgnoreCase("gulpepa") == 0) {
			String potentialFile = objfcnTokens.get(1);
			int timeLimit = Integer.parseInt(objfcnTokens.get(2));
			boolean cautious = Boolean.parseBoolean(objfcnTokens.get(3));
			objFcn = new GulpEPAObjFcn(potentialFile, timeLimit, this, cautious);
		} else if (objFcnType.compareToIgnoreCase("vaspepa") == 0) {
			throw new RuntimeException("vasp not supported yet srys.");
		} else {
			throw new RuntimeException("Unknown ObjFcn passed: " + objFcnType);
		}
		
		// other stuffs
		writeTempFiles = false;
		outDir = null;
		if (aparser.hasOption("outputTempFiles"))
			writeTempFiles = true;
		if (aparser.hasOption("outDir"))
			outDir = aparser.getArgument("outDir");

	}
	
	public ObjectiveFunction getObjFcn() {
		return objFcn;
	}
	
	public double getMinLatticeLength() {
		return minLatticeLength;
	}
		
	public double getMaxLatticeLength() {
		return maxLatticeLength;
	}
	
	public double getMinInteratomicDistance() {
		return minInteratomicDistance;
	}
	
	public double getMinLatticeAngle() {
		return minLatticeAngle;
	}
	
	public double getMaxLatticeAngle() {
		return maxLatticeAngle;
	}
	
	public CompositionSpace getCompositionSpace() {
		//TODO: throw exception or fix or something when optisystem becomes non-fixed composition
		return compositionSpace;
	}
	
	public int getMinNumAtoms() {
		return minNumAtoms;
	}
	
	public int getMaxNumAtoms() {
		return maxNumAtoms;
	}
	
	public boolean getOptimizeSites() {
		return optimizeSites;
	}
	
	public boolean getOptimizeCell() {
		return optimizeCell;
	}
	
	public boolean getWriteTempFiles() {
		return writeTempFiles;
	}
	
	public String getOutDir() {
		return outDir;
	}
	
	public boolean isConverged(Map<Cell,Double> cells) {
		return cells.size() >= numStructures;
	}
}
