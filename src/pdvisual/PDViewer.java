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
				System.out.println("Using ManualComputedEntry with composition " + cell.getComposition() + " and total energy " + energy);
			}
		}
		
		PDBuilder pdb = new PDBuilder(entries, elements, chempots);
		PDData pdd = pdb.getPDData();
				
		if (aparser.hasOption("v")) {
			// remove entries with energy > 0 and add dummy entries w/ energy = 0 at each vertex
			
			java.util.List<IComputedEntry> newEntries = new ArrayList<IComputedEntry>();
			for (IComputedEntry e : entries)
				if (e.getEnergyPerAtom() <= 0)
					newEntries.add(e);
			
			for (Element e : elements) {
				List<Site> els = new ArrayList<Site>();
				els.add(new Site(e, new Vect(0.0,0.0,0.0)));
				StructureOrg dummy = new StructureOrg(new Cell(1.,1.,1.,90.,90.,90.,els,""));
				dummy.setTotalEnergy(0);
				newEntries.add(dummy);
			}
			
			System.out.println((new PDBuilder(newEntries, elements, chempots)).getPDData().getCHullVolume());
			System.exit(0);
		}
		
		new TernPD3D(pdd);
	}

}
