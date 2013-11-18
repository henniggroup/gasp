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
		System.out.println("Usage: pdviewer --i <input.pdb.tgz>+");
		System.out.println("       --removeZeroEnergies");
		System.out.println("       --m POSCAR1 ENERGY1 POSCAR2 ENERGY2...");
		System.out.println("       --v    - just print volume and exit");
		System.out.println("       --c <file.pdb.tgz>	- compare other structures to input");
		System.out.println("       --minEn <num>  - minimum total energy");
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
			List<String> pddfnames = aparser.getArguments("i");
			for (String fname : pddfnames) {
				PDData pd = (PDData)(Utility.readSerializable(fname));
			
				for (IComputedEntry e : pd.getAllEntries())
					if (e.getTotalEnergy() != 0.0 || ! aparser.hasOption("removeZeroEnergies"))
						entries.add(e);
				
				for (int i = 0; i < pd.getElements().size(); i++) {
					if (!elements.contains(pd.getElements().get(i))) {
						Element e = pd.getElements().get(i);
						elements.add(e);
						chempots.put(e, pd.getChemPots().get(e));
					}
				}
			}
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
		
		if (aparser.hasArguments("minEn")) {
			double minEn = Double.parseDouble(aparser.getArgument("minEn"));
			Iterator<IComputedEntry> it = entries.iterator();
			while (it.hasNext()) {
				IComputedEntry e = it.next();
				if (e.getTotalEnergy() < minEn)
					it.remove();
			}
		}
		
		PDBuilder pdb = new PDBuilder(entries, elements, chempots);
		PDData pdd = pdb.getPDData();
		
		if (aparser.hasOption("c")) {
			PDData pd_to_compare = (PDData)(Utility.readSerializable(aparser.getArgument("c")));
			PDAnalyzer pda = new PDAnalyzer(pdd);
			for (IComputedEntry e : pd_to_compare.getAllEntries()) {
				System.out.println(e.toString());
				System.out.println("Energy above chull: " + pda.getEnergyPerAtomAboveHull(e));
			}
			System.exit(0);
		}
				
		if (aparser.hasOption("v")) {
			PDHartkeMaker.getHartkeVolume(pdd);
			System.exit(0);
		}
		
		if (aparser.hasOption("makeChullAGR")) {
			// assume ternary for now
		    Element e1 = pdd.getElements().get(0);
		    Element e2 = pdd.getElements().get(1);
		    Element e3 = pdd.getElements().get(2);
		    
		    /* xmgrace data set looks like:
@target G0.S0
@type xy
0.375 0.21650635
0.42857143 0.24743583
0.25 0.4330127
0.375 0.21650635
&
		     */
		    
		    StringBuilder output = new StringBuilder();

		    output.append("@version 50123\n");
		    output.append("@with g0\n");
		    output.append("@    world -0.05, -0.05, 1, 1\n");
		    output.append("@    view 0.300000, 0.150000, 1.000000, 0.850000\n");
		    output.append("@    xaxis tick off\n");
		    output.append("@    yaxis tick off\n");
		    output.append("@    xaxis ticklabel off\n");
		    output.append("@    yaxis ticklabel off\n");
		    output.append("@    frame color 0\n");
		    int counter = 0;
		   
		    for (int i = 0; i < 3; i++) {
		    	int indx1 = i;
		    	int indx2;
		    	int indx3;
		    	if (i == 0) {
		    		indx2 = 2;
		    		indx3 = 1;
		    	} else if (i == 1) {
		    		indx2 = 0;
		    		indx3 = 2;
		    	} else {
		    		indx2 = 1;
		    		indx3 = 0;
		    	}
		    	
		    	for (int j = 1; j <= 9; j++) {
		    		double[] cartCoords = new double[3];
		    		output.append("@    s"+counter+" line color 7\n");
			    	output.append("@target G0.S"+counter+"\n");
			    	output.append("@type xy\n");
	
		    		cartCoords[indx1] = j / 10.0;
		    		cartCoords[indx2] = 1 - j / 10.0;
		    		cartCoords[indx3] = 0;
			    	double[] ternCoords = TernPD3DProj.getTernaryXYCoord(cartCoords);
			     
			    	if (j%2 == 0) {
				    	double sqrt3o2 = 0.5 * Math.sqrt(3);
				    	double nub_scale = 0.02;
				    	if (indx1 == 1) {
				    		ternCoords[0] += nub_scale;
				    	} else if (indx1 == 0) {
				    		ternCoords[0] += - nub_scale * 0.5;
				    		ternCoords[1] += - nub_scale * sqrt3o2;
				    	} else {
				    		ternCoords[0] += - nub_scale * 0.5;
				    		ternCoords[1] += nub_scale * sqrt3o2;
				    	}
			    	}
			    	
			    	output.append(ternCoords[0] + " " + ternCoords[1] + "\n");
			    	
		    		cartCoords[indx1] = j / 10.0;
		    		cartCoords[indx2] = 0;
		    		cartCoords[indx3] = 1 - j / 10.0;
			    	ternCoords = TernPD3DProj.getTernaryXYCoord(cartCoords);
			    	output.append(ternCoords[0] + " " + ternCoords[1] + "\n");
			    	output.append("&\n");

			    	counter++;
		    	}
		    }
		    
		    for (List<Integer> facet : pdd.getIndxFacets()) {
	    		output.append("@    s"+counter+" line color 1\n");
	    		output.append("@    s"+counter+" symbol 1\n");
	    		output.append("@    s"+counter+" symbol size 0.500000\n");
	    		output.append("@    s"+counter+" symbol color 1\n");
	    		output.append("@    s"+counter+" symbol fill color 1\n");
	    		output.append("@    s"+counter+" symbol fill pattern 1\n");

		    	output.append("@target G0.S"+counter+"\n");
		    	output.append("@type xy\n");
		    	
		    	// assume ternary
		    	double[] cartCoords = new double[3];
		    	for (int i = 0; i < 3; i++) {
			    	cartCoords[0] = pdd.getAllEntries().get(facet.get(i)).getComposition().getFractionalCompo(e1);
			    	cartCoords[1] = pdd.getAllEntries().get(facet.get(i)).getComposition().getFractionalCompo(e2);
			    	cartCoords[2] = pdd.getAllEntries().get(facet.get(i)).getComposition().getFractionalCompo(e3);
			    	double[] ternCoords = TernPD3DProj.getTernaryXYCoord(cartCoords);
			    	output.append(ternCoords[0] + " " + ternCoords[1] + "\n");
		    	}
		    	cartCoords[0] = pdd.getAllEntries().get(facet.get(0)).getComposition().getFractionalCompo(e1);
		    	cartCoords[1] = pdd.getAllEntries().get(facet.get(0)).getComposition().getFractionalCompo(e2);
		    	cartCoords[2] = pdd.getAllEntries().get(facet.get(0)).getComposition().getFractionalCompo(e3);
		    	double[] ternCoords = TernPD3DProj.getTernaryXYCoord(cartCoords);
		    	output.append(ternCoords[0] + " " + ternCoords[1] + "\n");
		    	output.append("&\n");
		    	
		    	counter++;
		    }
	    	
	    	Utility.writeStringToFile(output.toString(), "output_pd.agr");
			
			System.exit(0);
		}
		
		new TernPD3D(pdd);
	}

}
