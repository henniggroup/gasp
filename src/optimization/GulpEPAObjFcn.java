package optimization;

import crystallography.Cell;
import ga.ObjectiveFunction;
import gulp.*;
import utility.Pair;

public class GulpEPAObjFcn extends ObjectiveFunction {
	
	String potentialFile;
	OptiSystem sys;
	int timeLimit;
	boolean cautious;
	
	public GulpEPAObjFcn(String pot, int t, OptiSystem s, boolean _cautious) {
		potentialFile = pot;
		sys = s;
		timeLimit = t;
		cautious = _cautious;
	}

	public Pair<Cell,Double> getCellAndValue(Cell c) {		
		GulpEnergy ge = new GulpEnergy(c, potentialFile, timeLimit, true, sys, cautious);
		
		return new Pair<Cell,Double>(ge.getOptimizedCell(), ge.getOptimizedEnergy() / c.getBasisSize());
	}

	@Override
	public void run() {
		// TODO Auto-generated method stub
		
	}

	@Override
	public Thread evaluate() {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public String toString() {
		// TODO Auto-generated method stub
		return null;
	}

}
