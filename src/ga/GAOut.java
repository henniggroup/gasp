package ga;

import java.io.Serializable;

public class GAOut implements Serializable {
	final static long serialVersionUID = 1l;
	
	private static GAOut instance;
	
	// output levels
	public static final int DEBUG = 5;
	public static final int INFO = 4;
	public static final int NOTICE = 3;
	public static final int WARNING = 2;
	public static final int CRITICAL = 1;
	
	// singleton
	public static GAOut out() {
		if (instance == null)
			instance = new GAOut();

		return instance;
	}
	
	public void stdout(String message, int level) {
		stdout(message, level, -1);
	}
	
	public void stdout(String message, int level, int structureID) {
		
		/*
		for (int i = 30; i < 38; i++) {
			System.out.println("\033[" + i + "mHello " +
			 "\033[1;" + i + "mWorld!\033[m");
			} */
		
		if (level <= GAParameters.getParams().getVerbosity())
			System.out.println(message);
	}
	
}
