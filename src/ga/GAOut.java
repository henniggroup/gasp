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
	
	private String colorifyStringByStructID(String s, int id) {
		
		// valid ansi color codes go from 30 to 38
		int colorNum = (id % 9) + 30;
		
		return "\033[1;" + colorNum + "m" + s + "\033[m";
	}
	
	public void stdout(String message, int level, int structureID) {
		if (level <= GAParameters.getParams().getVerbosity()) {
			if (structureID > 0) {
				if (!orgsSeen.contains(structureID)) {
					orgsSeen.add(structureID);
					System.out.println(colorifyStringByStructID("Organism " + structureID,structureID));
				}
				System.out.println(colorifyStringByStructID("   " + message,structureID));
			} else {
				System.out.println(message);
			}
		}
	}
	
}
