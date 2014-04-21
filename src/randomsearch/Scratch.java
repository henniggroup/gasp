package randomsearch;

import java.util.*;

import chemistry.Element;
import crystallography.*;
import optimization.OptiSystem;
import pdvisual.*;
import utility.ArgumentParser;
import utility.Pair;

public class Scratch {

	public static void main(String[] args) {

		ArgumentParser aparser = new ArgumentParser(args);
		
		OptiSystem sys = new OptiSystem(args);
		
		// make some cells
		List<Cell> cells = new LinkedList<Cell>();
		for (int i = 0; i < 150; i++)
			cells.add(RandomSearchMethod.getRandomCell(sys, i));
		
		// make some entries
		List<ComputedEntry> entries = new LinkedList<ComputedEntry>();
		for (Cell c : cells) {
			Pair<Cell,Double> result = sys.getObjFcn().getCellAndValue(c);
			entries.add(new ComputedEntry(result.getFirst(), result.getSecond()));
		}
		
		for (ComputedEntry c : entries)
			System.out.println("Entry " + c.getLabel() + " has epa " + c.getEnergyPerAtom());

		
		List<Element> elements = sys.getCompositionSpace().getElements();
		Map<Element, Double> chempots = new HashMap<Element,Double>();
		for (Element e : elements) {
			chempots.put(e, 0.0);
		}
		
		PDBuilder pdbuilder = new PDBuilder(entries, elements, chempots);
		
		PDData pddata = pdbuilder.getPDData();
		
		TernPD3D pdviz = new TernPD3D(pddata);
		
	}

}
