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
package ga;

import java.io.Serializable;
import java.util.ArrayList;
import java.util.List;

public class GAOut implements Serializable {
	final static long serialVersionUID = 1l;
	
	private static GAOut instance;
	
	// output levels
	public static final int DEBUG = 5;
	public static final int INFO = 4;
	public static final int NOTICE = 3;
	public static final int WARNING = 2;
	public static final int CRITICAL = 1;
	
	private List<Integer> orgsSeen;
	
	private GAOut() {
		orgsSeen = new ArrayList<Integer>();
	}
	
	// singleton
	public static GAOut out() {
		if (instance == null)
			instance = new GAOut();

		return instance;
	}
	
	public void stdout(String message, int level) {
		stdout(message, level, -1);
	}
	
	private String possiblyColorifyStringByStructID(String s, int id) {
		
		// valid ansi color codes go from 30 to 38
		int colorNum = (id % 9) + 30;
		if (GAParameters.getParams().getColorOutput())		
			return "\033[1;" + colorNum + "m" + s + "\033[m";
		else
			return s;

	}
	
	public void stdout(String message, int level, int structureID) {
		if (level <= GAParameters.getParams().getVerbosity()) {
			if (structureID > 0) {
				if (!orgsSeen.contains(structureID)) {
					orgsSeen.add(structureID);
					System.out.println(possiblyColorifyStringByStructID("Organism " + structureID,structureID));
				}
				System.out.println(possiblyColorifyStringByStructID("   " + message,structureID));
			} else {
				System.out.println(message);
			}
		}
	}
	
}
