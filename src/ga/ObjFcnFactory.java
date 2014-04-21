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
		} else {
			GAParameters.usage("Unknown objective function " + objFcnType, true);	
		}
		
		return obj;
	}

}
