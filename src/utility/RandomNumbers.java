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

package utility;

import ga.GAParameters;

public class RandomNumbers {
	
	public static double getUniformDoubleBetween(double a, double b) {
		double lesser, greater;
		if (a > b) {
			greater = a;
			lesser = b;
		} else {
			greater = b;
			lesser = a;
		}
		return lesser + GAParameters.getParams().getRandom().nextDouble() * (greater - lesser);
	}
	
	public static int getUniformIntBetweenInclusive(int a, int b) {
		
		// TODO: wait does this work?
		int lesser, greater;
		if (a > b) {
			greater = a;
			lesser = b;
		} else {
			greater = b;
			lesser = a;
		}
		return lesser + GAParameters.getParams().getRandom().nextInt(greater - lesser + 1);
	//	return lesser + Math.round(Math.round(GAParameters.getParams().getRandom().nextDouble() * (greater - lesser)));
	}
	
	public static void main(String args[]) {
		for(int i = 0; i < 10; i++)
			System.out.println(getUniformIntBetweenInclusive(0,4));
	}

}
