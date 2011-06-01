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
	
	public static int getUniformIntBetween(int a, int b) {
		
		// TODO: wait does this work?
		int lesser, greater;
		if (a > b) {
			greater = a;
			lesser = b;
		} else {
			greater = b;
			lesser = a;
		}
		return lesser + Math.round(Math.round(GAParameters.getParams().getRandom().nextDouble() * (greater - lesser)));
	}

}
