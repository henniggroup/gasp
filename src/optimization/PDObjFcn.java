package optimization;

import ga.*;

import java.util.*;
import pdvisual.*;
import utility.*;
import crystallography.*;
import chemistry.*;


public class PDObjFcn extends ObjectiveFunction {
	
	ObjectiveFunction energyFcn;
	PDBuilder pdbuilder;
	
	StructureOrg org;
	
/*	public PDObjFcn(ObjectiveFunction _objFcn, List<Element> elements, PDData _pddata) {
		List<ComputedEntry> entries = new LinkedList<ComputedEntry>();
		Map<Element,Double> chempots = new HashMap<Element,Double>();
		pdbuilder = new PDBuilder(entries, elements, chempots);
	} */
	
	public PDObjFcn(String[] subArray, Organism o, PDBuilder pdb) {
		// TODO Auto-generated constructor stub
		energyFcn = new EnergyPerAtom(subArray, o);
		org = (StructureOrg)o;
		pdbuilder = pdb;
	}

	// TODO: remember to update numCalculations when implementing this
	public Thread evaluate() {		
		// short circuit here if we've done the calculation already
		if (org.knowsValue()) // TODO: have to update due to changing phase diagram :/
			return null;
		
		// another total energy calculation:
		numCalculations++;
		
		// start the calculation and return the Thread
		Thread t = new Thread(this);
		t.start();
		return t;
	}
	
	public void run() {
		
		// relax the cell - have to wait for it to finish before using results
		Thread t = energyFcn.evaluate();
		try {
				t.join();
		} catch (InterruptedException x) {
			if (GAParameters.getParams().getVerbosity() >= 3)
				System.out.println("InterruptedException in energy calc thread in PDObjFcn: " + x.getMessage());
		}
		
		// updating structure w/ relaxed version is this done in EPA already
		org.setValue((new PDAnalyzer(pdbuilder.getPDData())).getEnergyPerAtomAboveHull(org));
		
	}
	
	// an ObjectiveFunction should also overload toString();
	public String toString() {
		return "PDObjFcn"; // TODO: more? what does cellobjfcn do?
	}
}
