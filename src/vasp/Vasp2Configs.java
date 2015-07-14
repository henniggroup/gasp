/*
 * Copyright 2011-2014 Will Tipton, Richard Hennig, Ben Revard, Stewart Wenner
 
This file is part of the Genetic Algorithm for Structure and Phase Prediction (GASP).

    GASP is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.
    
    GASP is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.
    
    You should have received a copy of the GNU General Public License
    along with GASP.  If not, see <http://www.gnu.org/licenses/>.
    
    
    */

package vasp;

import ga.RedundancyGuard;
import ga.StructureOrg;

import java.io.File;
import java.io.FilenameFilter;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;

import pdvisual.IComputedEntry;
import pdvisual.PDAnalyzer;
import pdvisual.PDBuilder;
import pdvisual.PDData;

import utility.ArgumentParser;
import utility.Utility;
import utility.Vect;

import chemistry.Element;

public class Vasp2Configs {
	
	
	public static void usage() {
		System.out.println("vasp2config --i <vasp run directories with OUTCAR and POSCARs>+");
		System.out.println("  --map (<symbol> <int>)+ : mapping of element to type id");
		System.out.println("  --mid <number> : discard configs which violate this minimum interatomic distance");
		System.out.println("  --maxeform <fepa> : generate chull from all input configs and discard ones this far above chull");
		System.out.println("  --maxsamples <num >= 2> : max configs from a single OUTCAR (always takes the ground state and some sampling of others)");
		System.out.println("  --writePDB <filename> : write configs to a .pdd file instead of printing them");
		System.out.println("  --writeConfigs : write configs to stdout in v2fout format");
		System.out.println("  --writePOSCARS : write POSCAR files to current dir");
		System.out.println("  --eliminateDups <atomicMisfit> <latticeMisfit> <angleMisfit> : explicitly check for and remove duplicate structures");
		System.out.println("  --forceSort : sort according to sum of forces");
		System.out.println("  --shuf : randomly shuffle outputs");
	}
	
	public static void applyMidConstraint(double mid, List<VaspConfig> configs) {
		Iterator<VaspConfig> it = configs.iterator();
		while (it.hasNext()) {
		   VaspConfig i = it.next();
		   if (! i.getCell().satisfiesMinInteratomicDistance(mid))
			   it.remove();
		}
	}

	public static void applyMaxsamplesConstraint(int maxsamples, List<VaspConfig> configs) {
		// let k = configs.size() - maxsamples be the number of configs we want to remove
		// the list of indices we want to remove is then
		//      floor((1:k)*(n-2)/k)

		int n = configs.size();
		int k = n - maxsamples;
		if (k <= 0)
			return;
		
		for (int i = k; i >= 1; i--)
			configs.remove((int) Math.floor(i * (n-2)/k));
	}
	
	public static PDBuilder makePDBFromConfigs(List<VaspConfig> configs) {
		// make the PDD and PDA from the list of configs
		List<IComputedEntry> entries = new LinkedList<IComputedEntry>();
		List<Element> elements = new LinkedList<Element>();
		Map<Element,Double> chempots = new HashMap<Element,Double>();
		for (VaspConfig c : configs) {
			entries.add(c);
			for (Element e : c.getComposition().getElements())
				if (!elements.contains(e)) {
					elements.add(e);
					chempots.put(e, 0.0);
				}
		}	
		
		return new PDBuilder(entries, elements, chempots);
	}
	
	public static void applyMaxeformConstraint(double maxeform, List<VaspConfig> configs) {

		PDAnalyzer pda = new PDAnalyzer(makePDBFromConfigs(configs).getPDData());
		
		// remove bad configs
		Iterator<VaspConfig> it = configs.iterator();
		while (it.hasNext()) {
		   VaspConfig i = it.next();
		   if (pda.getEnergyPerAtomAboveHull(i) > maxeform)
			   it.remove();
		}
	}

	public static void main(String args[]) {
		
		ArgumentParser aParser = new ArgumentParser(args);
		if (! aParser.hasArguments("i")) {
			usage();
			System.exit(0);
		}
		
		if (aParser.hasArguments("map")) {
			VaspConfig.parseTypeMap(aParser.getArguments("map"));
		}
		
		List<VaspConfig> configs = new ArrayList<VaspConfig>();
		
		// get all the configs from each input vasp dir
		for (String vaspdir : aParser.getArguments("i")) { 
			File dir = new File(vaspdir);
			if (!dir.exists())
				System.out.println("ERROR: directory " + dir + " does not exist.");
			File[] outcarFiles = dir.listFiles(new FilenameFilter() {
				public boolean accept(File f, String name) {
					return name.startsWith("OUTCAR");
				}
			});

			for (File oc : outcarFiles) {
			
				List<VaspConfig> iconfigs = VaspOut.getConfigs(oc.getAbsolutePath(), VaspOut.getPOSCAR(vaspdir + "/POSCAR"));
				
				// remove configs that fail mid and maxsamples constraints
				if (aParser.hasArguments("mid"))
					applyMidConstraint(Double.parseDouble(aParser.getArgument("mid")), iconfigs);
				if (aParser.hasArguments("maxsamples"))
					applyMaxsamplesConstraint(Integer.parseInt(aParser.getArgument("maxsamples")),iconfigs);
				
				configs.addAll(iconfigs);		
			}
		}
		
		// remove configs that violate maxeform constraint
		if (aParser.hasArguments("maxeform"))
			applyMaxeformConstraint(Double.parseDouble(aParser.getArgument("maxeform")), configs);
		
		// possibly check for duplicate configs
		if (aParser.hasArguments("eliminateDups")) {
			Iterator<VaspConfig> it = configs.iterator();
			RedundancyGuard rg = new RedundancyGuard(aParser.getArguments("eliminateDups"));
			while (it.hasNext()) {
			   VaspConfig i = it.next();
			   StructureOrg si = new StructureOrg(i.getCell());
			   if (rg.checkStructureOrg(si) != null) { // if i fails the check
				   it.remove();
			   } 
			   rg.addStructureOrg(si);
			}
		}
		
		// sort maybe
		if(aParser.hasOption("forceSort")) {
			Collections.sort(configs, new Comparator<VaspConfig>() {

				@Override
				public int compare(VaspConfig c1, VaspConfig c2) {
					double c1sum = 0;
					for (Vect v : c1.getForces())
						c1sum += v.length();
					double c2sum = 0;
					for (Vect v : c2.getForces())
						c2sum += v.length();
						
					return (c1sum - c2sum) > 0 ? 1 : -1;
				}
				
			});
		}

		// outputs
		if (aParser.hasArguments("writePDB")) {
			String fName = aParser.getArgument("writePDB");
			Utility.writeSerializable(makePDBFromConfigs(configs).getPDData(), fName);
			
		} 
		if (aParser.hasOption("shuf")) {
			Collections.shuffle(configs);
		}
		if(aParser.hasOption("writeConfigs")) {
			int numCells = 0;
			for (VaspConfig i : configs) {
				String comment = "## cell " + (numCells++) + System.getProperty("line.separator") 
									+ "## " + i.getCell().getLabel();
				System.out.println(i.toPotfitString(comment));
			}
		}
		if(aParser.hasOption("writePOSCARS")) {
			int numCells = 0;
			for (VaspConfig i : configs) {
				String fname = (numCells++) + ".POSCAR";
				VaspIn.writePoscar(i.getCell(), fname, true);
			}
		}
	}

}