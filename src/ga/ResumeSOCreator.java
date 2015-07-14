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

import java.util.ArrayList;
import java.util.Collections;
import java.util.Iterator;
import java.io.*;

import crystallography.Cell;

// ResumeSOCreator implements StructureOrgCreator.  It takes as input the directory
// of a previous run of the GA.  It loads the last generation of StructureOrgs as well
// as their IDs and, optionally, their energies.  However, it does not load the other
// parameters of the algorithm and thus is not a complete "resume"-type function.

public class ResumeSOCreator implements StructureOrgCreator {
	
	String dirName;
	
	ArrayList<StructureOrg> orgList;
	
	Boolean recalcEnergies;
	
	public ResumeSOCreator(String[] args) {
		if (args == null || args.length < 2)
			GAParameters.usage("Not enough arguments passed to ResumeSOCreator.", true);
		dirName = new String(args[0]);	
		recalcEnergies = Boolean.parseBoolean(args[1]);
	}
	
	public String toString() {
		return "ResumeSOCreator: directory = " + dirName;
	}
	
	private void resume() {
		orgList = new ArrayList<StructureOrg>();
		
		// read the index file
		String index = GAUtils.readStringFromFile(GAParameters.getParams().getRecord().getOutFile());
//		String index = GAUtils.readStringFromFile(new File(dirName + "/index"));
		// find the last generation and its filenames
		ArrayList<String> filenames = new ArrayList<String>();
		ArrayList<Integer> ids = new ArrayList<Integer>();
		ArrayList<Double> values = new ArrayList<Double>();
		String[] lines = index.split("\n");
		int generation = 0;
		for (int i = 0; i < lines.length; i++) {
			if (lines[i].startsWith("generation")) {
				// store the generation
				String[] tokens = lines[i].split(" ");
				generation = Integer.parseInt(tokens[1]);
				// clear the filenames and ids arrays
				filenames = new ArrayList<String>();
				ids = new ArrayList<Integer>();
				continue;
			}
			String[] fields = lines[i].split(" ");
			filenames.add(fields[2]);
			values.add(new Double(Double.parseDouble(fields[1])));
			ids.add(new Integer(Integer.parseInt(fields[0])));
		}
		
		// make StructureOrgs from each of the filenames and set their ID's
		Iterator<String> i = filenames.iterator();
		Iterator<Integer> j = ids.iterator();
		Iterator<Double> d = values.iterator();
		while (i.hasNext()) {
			String filename = i.next();
			StructureOrg s = new StructureOrg(Cell.parseCif(new File(dirName + "/" + filename)));
			s.setID(j.next());
			if (!recalcEnergies)
				s.setValue(d.next());
			orgList.add(s);
		}
		// seed the IDs and generation
		GAParameters.getParams().seedIDs(Collections.max(ids));
		GAParameters.getParams().getRecord().setCurrentGenNum(generation);
	}
	
	// returns an array of all the structures in the most recent generation on disk.
	// sets their ID appropriately.  also sets the current generation number.
	// we assume the format of the "index" file. 
	public StructureOrg makeOrganism(Generation g) {
		if (orgList == null)
			resume();
		
		if (orgList.size() == 0)
			return null;
		
		return orgList.get(0);
	}

}
