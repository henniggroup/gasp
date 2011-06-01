package pdvisual;

import java.util.*;

import chemistry.Element;

import crystallography.Cell;
import utility.ArgumentParser;
import utility.Utility;
import vasp.VaspOut;

public class PDViewer {

	public static void usage () {
		System.out.println("Usage: pdviewer --i <input.pdb.tgz>");
		System.out.println("       --removeZeroEnergies");
		System.out.println("       --m POSCAR1 ENERGY1 POSCAR2 ENERGY2...");
	}
	
	public static void main(String[] args) {
		ArgumentParser aparser = new ArgumentParser(args);

		if (!aparser.hasArguments("i") && !aparser.hasArguments("m")) {
			usage();
			return;
		}
		
		java.util.List<IComputedEntry> entries = new LinkedList<IComputedEntry>();
		java.util.List<Element> elements = new LinkedList<Element>();
		java.util.Map<Element,Double> chempots = new HashMap<Element,Double>();
		
		if (aparser.hasArguments("i")) {
			PDData pd = (PDData)(Utility.readSerializable(aparser.getArgument("i")));
		
			for (IComputedEntry e : pd.getAllEntries())
				if (e.getTotalEnergy() != 0.0 || ! aparser.hasOption("removeZeroEnergies"))
					entries.add(e);
			
			elements.addAll(pd.getElements());
			chempots.putAll(pd.getChemPots());
		}

		if (aparser.hasArguments("m")) {
			List<String> mArgs = aparser.getArguments("m");
			if (mArgs.size() % 2 != 0) {
				usage();
				return;
			}
			for (int i = 0; i < mArgs.size(); i = i + 2) {
				Cell cell = VaspOut.getPOSCAR(mArgs.get(i));
				double energy = Double.parseDouble(mArgs.get(i+1));
				entries.add(new ManualComputedEntry(cell, energy));
				for (Element e : cell.getComposition().getElements())
					if (!elements.contains(e)) {
						elements.add(e);
						chempots.put(e, 0.0);
					}
			}
		}
				
		new TernPD3D((new PDBuilder(entries, elements, chempots)).getPDData());
	}

}
