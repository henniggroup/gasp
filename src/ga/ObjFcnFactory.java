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

import java.util.List;

import optimization.PDObjFcn;
import utility.Utility;

public class ObjFcnFactory {
	
	public static ObjectiveFunction getObjectiveFunctionInstance(Organism o, List<String> objFcnArgs) {
		if (objFcnArgs == null || objFcnArgs.size() < 1)
			GAParameters.usage("Not enough parameters passed to --ObjectiveFunction", true);
		
		ObjectiveFunction obj = null;
				
		String objFcnType = objFcnArgs.get(0);
		if (objFcnType.equalsIgnoreCase("epa")) 
			obj = new EnergyPerAtom(Utility.subList(objFcnArgs, 1), o);
		else if (objFcnType.equalsIgnoreCase("pd")) {
		//	if (pdbuilder == null)
		//		pdbuilder = makePDBuilder();
			if (GAParameters.getParams().getRecord().getGenNum() < 1) {
				obj = new EnergyPerAtom(Utility.subList(objFcnArgs, 1), o);
			} else {
				obj = new PDObjFcn(Utility.subList(objFcnArgs, 1), o);
			}
		} else if (objFcnType.equals("cluster")) {
			obj = new ClusterObjFcn(Utility.subList(objFcnArgs, 1), o);
		} else if (objFcnType.equals("surface")) {
			obj = new SurfaceObjFcn(Utility.subList(objFcnArgs, 1), o);
//		} else if (objFcnType.equals("substrate")) {
//			obj = new SubstrateObjFcn(Utility.subList(objFcnArgs, 1), o);
		} else {
			GAParameters.usage("Unknown objective function " + objFcnType, true);	
		}
		
		return obj;
	}

}
