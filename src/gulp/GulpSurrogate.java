package gulp;

import ga.GAParameters;
import ga.Generation;
import ga.RandomSOCreator;
import ga.StructureOrg;
import ga.Structures;

import java.io.Serializable;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.StringTokenizer;

import com.sun.org.apache.xml.internal.serializer.utils.Utils;

import chemistry.Composition;
import chemistry.CompositionSpace;
import chemistry.Element;

import crystallography.Cell;
import crystallography.Site;

import pdvisual.IComputedEntry;
import utility.Pair;
import utility.Utility;

public class GulpSurrogate implements Serializable {
	private static final long serialVersionUID = 1l;
	
	private List<IComputedEntry> entries;
	List<GulpPotential> pots;
	private String headerStr;
	private boolean usingShift = true;
	private double shift = 0.0;
	private int refit_freq;
	
	public GulpSurrogate(String args[]) {
		
		if (args.length != 3) {
			System.out.println("ERROR: wrong args given to GulpSurrogate:");
			for (String s : args)
				System.out.println(s);
		}
		
		headerStr = Utility.readStringFromFile(args[0]);
		
		entries = new ArrayList<IComputedEntry>();
 
		parsePotentialsString(Utility.readStringFromFile(args[1]));
		
		refit_freq = Integer.parseInt(args[2]);
	}
	
	private void parsePotentialsString(String potlsStr) {
		pots = new ArrayList<GulpPotential>();
		StringTokenizer tok = new StringTokenizer(potlsStr, "\n");
		
		int numPotentials = Integer.parseInt(tok.nextToken());
		
		for (int i = 0; i < numPotentials; i++) {
			// type
			String type = tok.nextToken();
			// species
			List<Element> species = new ArrayList<Element>();
			for (String s : tok.nextToken().split("[ ][ ]*"))
				species.add(Element.getElemFromSymbol(s));
			// cutoffs
			String cutoffStrs[] = tok.nextToken().split("[ ][ ]*");
			double cutoff_min = Double.parseDouble(cutoffStrs[0]);
			double cutoff_max = Double.parseDouble(cutoffStrs[1]);
			// num params
			int numParams = Integer.parseInt(tok.nextToken());
			// param, min, max, opt_order
			List<Double> parameters = new ArrayList<Double>();
			List<Double> parameters_min = new ArrayList<Double>();
			List<Double> parameters_max = new ArrayList<Double>();
			List<Integer> optOrders = new ArrayList<Integer>();
			for (int j = 0; j < numParams; j++) {
				String paramStrs[] = tok.nextToken().split("[ ][ ]*");
				parameters.add(Double.parseDouble(paramStrs[0]));
				parameters_min.add(Double.parseDouble(paramStrs[1]));
				parameters_max.add(Double.parseDouble(paramStrs[2]));
				optOrders.add(Integer.parseInt(paramStrs[3]));
			}			
			// add the potential
			pots.add(new GulpPotential(type, species, parameters, parameters_min, parameters_max, optOrders, cutoff_min, cutoff_max));
		}
	}
	
	private String getFittingRunInput(List<List<Boolean>> fitParms) {
		StringBuilder result = new StringBuilder();
		
		DecimalFormat format = new DecimalFormat("0.000000");
		
		String nl = System.getProperty("line.separator");
		result.append("fit c6 noflag" + nl);
		result.append("title" + nl);
		result.append("fitting run by GulpSurrogate" + nl);
		result.append("end" + nl);
		if (usingShift)
			result.append("shift 1.0" + nl);
		for (IComputedEntry ent : entries) {
			result.append("# Entry" + nl);
			result.append("cell" + nl);
			List<Double> lps = ent.getCell().getLatticeParameters();
			//  need to convert radians->degrees
			result.append(format.format(lps.get(0)) + " " + format.format(lps.get(1)) + " " + format.format(lps.get(2)) + " " +
					format.format(lps.get(3)*180/Math.PI) + " " + format.format(lps.get(4)*180/Math.PI) + " " + format.format(lps.get(5)*180/Math.PI) + nl);
			result.append("fractional" + nl);
			for (Site s : ent.getCell().getSites()) {
				List<Double> fcoords = s.getCoords().getComponentsWRTBasis(ent.getCell().getLatticeVectors());
				result.append(s.getElement().getSymbol() + " core " +
						format.format(fcoords.get(0)) + " " + format.format(fcoords.get(1)) + " " + format.format(fcoords.get(2)) + nl);
			}
			result.append("observable" + nl);
			result.append("energy eV" + nl);
			result.append(format.format(ent.getTotalEnergy()) + nl);
			
		//	result.append
			
			result.append("end" + nl);
		}
		if (usingShift) {
			result.append("vary" + nl);
			result.append("shift" + nl);
			result.append("end" + nl);
		}

		for (int i = 0; i < pots.size(); i++)
			result.append(pots.get(i).toGulpPotlStr(fitParms.get(i)) + nl);
		
		return result.toString();
	}
	
	/*
	private int getTotalNumParams() {
		int result = 0;
		if (usingShift)
			result++;
		for (Pair<GulpPotential,List<Boolean>> p : pots)
			result += getNumFittedParams(p.getSecond());
		return result;
	} */
	
	private int getMaxOptIndex() {
		int result = 0;
		for (GulpPotential p : pots)
			for (Integer indx : p.getOptIndexes())
				result = Math.max(result, indx);
		return result;
	}
	
	private double evaluateTotalEnergy(Cell c) {
		StringBuilder gulpInput = new StringBuilder();
		String nl = System.getProperty("line.separator");
		gulpInput.append(headerStr + nl);
		gulpInput.append(GulpEnergy.structureToString(c) + nl);
		for (int i = 0; i < pots.size(); i++)
			gulpInput.append(pots.get(i).toGulpPotlStr(null) + nl);
			
		//TODO: better temporary file generation
		Utility.writeStringToFile(gulpInput.toString(), "temp");
		return GulpEnergy.parseFinalEnergy(GulpEnergy.runGULP("temp"), false) + shift;
	}
	
	private List<List<Integer>> getInitPinnedParms() {
		List<List<Integer>> newPinnedParams = new ArrayList<List<Integer>>();
		for (int i = 0; i < pots.size(); i++) {
			List<Integer> pp = new ArrayList<Integer>();
			for (int j = 0; j < pots.get(i).getOptIndexes().size(); j++) 
				pp.add(0);
			newPinnedParams.add(pp);
		}
		return newPinnedParams;
	}
	
	// returns bool indicating whether or not this optimization run
	// caused any parameters to be pinned
	private void updateAllParameters(List<List<Integer>> pinnedParams) {			
//printPinnedParams(pinnedParams);
		for (int i = 1; i <= getMaxOptIndex(); i++) { // loop over optimization priorities
			List<List<Boolean>> fitParms = new ArrayList<List<Boolean>>();
			for (int j = 0; j < pots.size(); j++) {
				GulpPotential p = pots.get(j);// loop over all the potentials
				// set params to optimize for all relevant parameters
				List<Boolean> fp = new ArrayList<Boolean>();
				for (int k = 0; k < p.getOptIndexes().size(); k++) {
					int optOrder = p.getOptIndexes().get(k);
					if (optOrder == i && (pinnedParams.get(j).get(k) == 0))
						fp.add(true);
					else
						fp.add(false); 
				}
				fitParms.add(fp);				
			}
			
			// TODO: unique temporary file or whatnot
			String fname = "temp" + i;
			Utility.writeStringToFile(getFittingRunInput(fitParms), fname);
			List<GulpPotential> newPots = parseFittingRunOutput(GulpEnergy.runGULP(fname), fitParms);
			
			if (newPots != null) {
				pots = newPots;
				for (int j = 0; j < pots.size(); j++) {
					for (int k = 0; k < pots.get(j).getOptIndexes().size(); k++) 
						if (pots.get(j).isPinned(k) && pots.get(j).getOptIndexes().get(k) == i)
							pinnedParams.get(j).set(k, pinnedParams.get(j).get(k) + 1);
					pots.get(j).enforceConstraints();
				}
			}
			
		}


//printPinnedParams(pinnedParams);
	}
	
	private boolean constraintsViolated() {
		for (GulpPotential p : pots) 
			if (p.violatesConstraints()) 
				return true;

		return false;
	}
	
	private boolean needsRerun(List<List<Integer>> pinnedParams) {
		for (List<Integer> l : pinnedParams)
			for (Integer b : l)
				if (b == 1)
					return true;
		return false;
	}
	
	private void printPinnedParams(List<List<Integer>> pinnedParams) {
		for (List<Integer> l : pinnedParams) {
			for (Integer b : l)
				System.out.print(b + " ");
			System.out.println("");
		}
	}
	
	private void updateSurrogateModel() {
		if (GAParameters.getParams().getVerbosity() >= 4)
			System.out.println("Updating surrogate model...");
		
		List<List<Integer>> pinnedParms = getInitPinnedParms();
		updateAllParameters(pinnedParms);
		
		if (GAParameters.getParams().getVerbosity() >= 4)
			for (GulpPotential p : pots)
				System.out.println(p.toGulpPotlStr(null));
		
		while (needsRerun(pinnedParms)) {
			if (GAParameters.getParams().getVerbosity() >= 4)
				System.out.println("Constraints violated. Reoptimizing...");
			updateAllParameters(pinnedParms);
			if (GAParameters.getParams().getVerbosity() >= 4)
				for (GulpPotential p : pots)
					System.out.println(p.toGulpPotlStr(null));
		}
	}
	
	private int getTotalNumParams(List<List<Boolean>> fitParms) {
		int result = 0;
		for (List<Boolean> l: fitParms)
			for (Boolean b : l)
				if (b)
					result++;
		return result;
	}
	
	private List<GulpPotential> parseFittingRunOutput(String output, List<List<Boolean>> fitParms) {
		
		List<GulpPotential> result = null;
		
		try {
			List<Double> finalParameters = new ArrayList<Double>();
			int indx = output.lastIndexOf("Final values of parameters :");
			if (indx < 0)
				System.out.println("ERROR: final parameters not found in gulp output");
			output = output.substring(indx);
			
			// get all the final parameters
			StringTokenizer tok = new StringTokenizer(output);
			int numParams = getTotalNumParams(fitParms);
			if (usingShift)
				numParams++;
			for (Integer paramNum = 1; paramNum <= numParams; paramNum++) {
				while (tok.nextToken().compareTo(paramNum.toString()) != 0)
					;
				tok.nextToken();
				finalParameters.add(Double.parseDouble(tok.nextToken()));
			}
			
			// update the shift and potentials
			if (usingShift) {
				shift = finalParameters.remove(0);
			}
			result = new ArrayList<GulpPotential>();
			for (int i = 0; i < pots.size(); i++) {
				GulpPotential currPot = pots.get(i);
				List<Double> newParams = new ArrayList<Double>();
				for (int j = 0; j < getNumFittedParams(fitParms.get(i)); j++)
					newParams.add(finalParameters.remove(0));
				result.add(new GulpPotential(currPot, newParams, fitParms.get(i)));
			}
		} catch (Exception x) {
			if (GAParameters.getParams().getVerbosity() >= 2)
				System.out.println("WARNING: gulp surrogate model fitting failed. Not updating potential model.");
			result = null;
		}
		
		return result;
	}
	
	public void addEntry(IComputedEntry e) {
		entries.add(e);
		
		if (entries.size() % refit_freq == 0)
			updateSurrogateModel();
	}
	
	public boolean fails(IComputedEntry e) {
		
		// TODO: something better and also maybe make sure we dont use this or w/e before fitting model
		StructureOrg s = (StructureOrg)e;
		
		System.out.println("Surrogate model energy of org " + s.getID() + ": " + evaluateTotalEnergy(e.getCell()));
		
		return false;
	}
	
	public static int getNumFittedParams(List<Boolean> fitParams) {
		int result = 0;
		
		for (Boolean b : fitParams)
			if (b)
				result++;
		
		return result;
	}
	
	// just for testing
	public static void main (String args[]) {
		// set up the stuff needed for RandomSOCreator
		Generation g = new Structures();
		GAParameters.getParams().setMaxNumAtoms(10);
		GAParameters.getParams().setMinNumAtoms(1);
		GAParameters.getParams().setMaxLatticeAngle(140);
		GAParameters.getParams().setMinLatticeAngle(40);
		GAParameters.getParams().setMaxLatticeLength(15);
		GAParameters.getParams().setMinLatticeLength(2);
		List<Composition> comps = new LinkedList<Composition>();
		Element els[] = {Element.getElemFromSymbol("Al"), Element.getElemFromSymbol("Cu")};
		Double cs1[] = {1.0, 6.0}; Double cs2[] = {0.0, 0.0};
		comps.add(new Composition(els,cs1));
		comps.add(new Composition(els,cs2));
		CompositionSpace cspace = new CompositionSpace(comps);
		GAParameters.getParams().setCompSpace(cspace);
		String arg[] = {"randomVol"};
		RandomSOCreator soc = new RandomSOCreator(arg);
	
		String geArgs[] = {"/home/wtipton/projects/ga_for_crystals/gulp_header", "/home/wtipton/projects/ga_for_crystals/gulppotls/gulppotl_alcu", "false"};
		GulpEnergy energy = new GulpEnergy(geArgs);
		String surrArgs[] = {"/home/wtipton/projects/ga_for_crystals/gulp_header", "/home/wtipton/projects/materials/potsSpec", "15"};
		GulpSurrogate s = new GulpSurrogate(surrArgs);
		
		// randomly generate structures, evaluate them w/ "real" energy model, add them to surrogate model object
		List<StructureOrg> entries = new LinkedList<StructureOrg>();
		for (int i = 0; i < 30; i++) {
			StructureOrg cell = soc.makeOrganism(g);
			cell.setTotalEnergy(energy.getEnergy(cell));
			cell.getCell().writeCIF("/home/wtipton/structs/" + i + ".cif");
			s.addEntry(cell);
			entries.add(cell);
		}
		
		// randomly generate new structures
		// evaluate their energies with "real" model and with surrogate model
		for (int i = 0; i < 10; i++) {
	//		StructureOrg cell = soc.makeOrganism(g);
	StructureOrg cell = entries.get(i);
			System.out.println("Trial structure " + i + 
					". Real energy: " + energy.getEnergy(cell) +
					". Surrogate model energy: " + s.evaluateTotalEnergy(cell.getCell()));
		} 
		
	}	
}
