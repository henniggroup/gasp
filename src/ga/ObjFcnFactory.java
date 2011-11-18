package ga;

import optimization.PDObjFcn;

public class ObjFcnFactory {
	
	public static ObjectiveFunction getObjectiveFunctionInstance(Organism o, String[] objFcnArgs) {
		if (objFcnArgs == null || objFcnArgs.length < 1)
			GAParameters.usage("Not enough parameters passed to --ObjectiveFunction", true);
		
		ObjectiveFunction obj = null;
				
		String objFcnType = objFcnArgs[0];
		if (objFcnType.equalsIgnoreCase("epa")) 
			obj = new EnergyPerAtom(GAUtils.subArray(objFcnArgs, 1), o);
		else if (objFcnType.equalsIgnoreCase("pd")) {
		//	if (pdbuilder == null)
		//		pdbuilder = makePDBuilder();
			if (GAParameters.getParams().getRecord().getGenNum() < 1) {
				obj = new EnergyPerAtom(GAUtils.subArray(objFcnArgs, 1), o);
			} else {
				obj = new PDObjFcn(GAUtils.subArray(objFcnArgs, 1), o);
			}
		} else if (objFcnType.equals("cluster")) {
			obj = new ClusterObjFcn(GAUtils.subArray(objFcnArgs, 1), o);
		} else {
			GAParameters.usage("Unknown objective function " + objFcnType, true);	
		}
		
		return obj;
	}

}
