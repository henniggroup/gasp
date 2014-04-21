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
