package pdvisual;

import ga.StructureOrg;

import java.util.*;

import chemistry.Element;

import crystallography.Cell;
import crystallography.Site;
import utility.ArgumentParser;
import utility.Utility;
import utility.Vect;
import vasp.VaspOut;

public class PDViewer {

	public static void usage () {
		System.out.println("Usage: pdviewer --i <input.pdb.tgz>");
		System.out.println("       --removeZeroEnergies");
		System.out.println("       --m POSCAR1 ENERGY1 POSCAR2 ENERGY2...");
		System.out.println("       --v    - just print volume and exit");
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
		
		if (aparser.hasOption("r")) {
			List<IComputedEntry> newEnts = new ArrayList<IComputedEntry>(entries);

			for (String s : aparser.getArguments("r")) {
				int sid = Integer.parseInt(s);
				for (IComputedEntry ent : entries)
					if (((StructureOrg)ent).getID() == sid)
						newEnts.remove(ent);
			}
			
			entries = newEnts;
		}

		if (aparser.hasArguments("m")) {
			List<String> mArgs = aparser.getArguments("m");
			if (mArgs.size() % 2 != 0) {
				System.out.println(mArgs.size() +" args given to --m");
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
				System.out.println("Using ManualComputedEntry with composition " + cell.getComposition() + " and total energy " + energy);
			}
		}
		
		PDBuilder pdb = new PDBuilder(entries, elements, chempots);
		PDData pdd = pdb.getPDData();
				
		if (aparser.hasOption("v")) {
			PDHartkeMaker.getHartkeVolume(pdd);
			System.exit(0);
		}
		
		new TernPD3D(pdd);
	}

}
