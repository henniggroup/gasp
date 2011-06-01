package randomsearch;

import optimization.*;
import crystallography.*;
import utility.RandomNumbers;
import java.util.*;
import java.io.*;
import chemistry.*;
import utility.*;

public class RandomSearchMethod implements OptimizationMethod {
	
	private Map<Cell,Double> cells;
	private OptiSystem system;
	private String outDirPath;
	
	public RandomSearchMethod(OptiSystem _system) {
		this(_system, null);
	}
	
	public RandomSearchMethod(OptiSystem _system, String outDir) {
		// initialize the list of cells
		cells = new HashMap<Cell,Double>();
		outDirPath = outDir;
		system = _system;
	}
	
	// return minimum distance between Vect and the locations of a List<Site>
	private static double minDistanceHelper(Vect v, List<Site> sites, List<Vect> basis) {
		double min = Double.MAX_VALUE;
		List<List<Integer>> coefs = new LinkedList<List<Integer>>();
		List<Integer> first = new LinkedList<Integer>();
		first.add(2); first.add(0);first.add(0);
		List<Integer> second = new LinkedList<Integer>();
		second.add(0); second.add(2); second.add(0);
		List<Integer> third = new LinkedList<Integer>();
		third.add(0); third.add(0); third.add(2);
		coefs.add(first); coefs.add(second); coefs.add(third);
		Cell superCell = (new Cell(basis,sites)).getSupercell(coefs);
		
		Vect translatedVect = v;
		for (Vect vec : basis)
			translatedVect = translatedVect.plus(vec);
		
		for (Site s : superCell.getSites())
			min = Math.min(min, translatedVect.getCartDistanceTo(s.getCoords()));
		
		return min;
	}
	
	public static Cell getRandomCell(OptiSystem sys, int numCreated) {
		
		double a = RandomNumbers.getUniformDoubleBetween(sys.getMinLatticeLength(), sys.getMaxLatticeLength());
		double b = RandomNumbers.getUniformDoubleBetween(sys.getMinLatticeLength(), sys.getMaxLatticeLength());
		double c = RandomNumbers.getUniformDoubleBetween(sys.getMinLatticeLength(), sys.getMaxLatticeLength());
		double alpha = RandomNumbers.getUniformDoubleBetween(sys.getMinLatticeAngle(), sys.getMaxLatticeAngle());
		double beta = RandomNumbers.getUniformDoubleBetween(sys.getMinLatticeAngle(), sys.getMaxLatticeAngle());
		// so, um, the sum alpha+beta+gamma can't exceed 360 degrees
		// and also can't have gamma >  alpha + beta 
		double gammaMax = Math.min(sys.getMaxLatticeAngle(), 2 * Math.PI - alpha - beta);
		gammaMax = Math.min(gammaMax, alpha + beta);
		// similarly, can't have gamma + alpha < beta (that is, we need gamma > beta - alpha)
		//   				  or gamma + beta < alpha (that is, we need gamma > alpha - beta)
		double gammaMin = Math.max(sys.getMinLatticeAngle(), Math.abs(alpha-beta));
		double gamma = RandomNumbers.getUniformDoubleBetween(gammaMin, gammaMax);
		List<Vect> latticeVectors = Cell.getVectorsfromLParams(a,b,c,alpha,beta,gamma);
		
		// OptiSystem became non-fixed-composition
		List<Site> basis = new LinkedList<Site>();
		// TODO: do a better job of getting the right number of atoms
		Composition comp = null;
		do {
			comp = sys.getCompositionSpace().getRandomIntegerCompInSpace(10);
		} while (comp.isEmpty()); 
		Map<Element,Double> stoichUnit = comp.getStoichiometricUnit();
		int numAtoms = RandomNumbers.getUniformIntBetween(sys.getMinNumAtoms(), sys.getMaxNumAtoms());
		
		while (basis.size() < numAtoms) { // add a stoichiometric unit
			for (Element e : stoichUnit.keySet()) {
				// add numOfElement randomly placed sites to our list of sites
				for (int i = 0; i < stoichUnit.get(e); i++) {
					// make sure we choose a location that respects minInteratomicDistance
					Vect potentialLocation = null;
					int numTries = 0;
					do {
						potentialLocation = new Vect(Math.random(), Math.random(), Math.random(), latticeVectors);
						if (numTries++ > 50) // give up eventually if lattice vectors are just too small and try again
							return getRandomCell(sys, numCreated);
					} while(minDistanceHelper(potentialLocation,basis,latticeVectors) < sys.getMinInteratomicDistance());
					
					basis.add(new Site(e,potentialLocation));
				}
			}
		}
		
		String label = Integer.toString(numCreated);
			
		return new Cell (latticeVectors, basis, label);
	}
	
	private boolean isConverged(OptiSystem s, Map<Cell, Double> cells) {
		//TODO: do better? use criteria in the optisystem?
		if (s.isConverged(cells))
			return true;
		return false;
	}

	// doOpt: make a bunch of trial cells, evaluate them, and return the best one
	public Cell doOpt() {
		
		Cell relaxedCell = null;
		try {
			// open the output file
			String outFilePath = outDirPath + "/output";
			BufferedWriter out = new BufferedWriter(new FileWriter(outFilePath));
			
			// make a bunch of trial cells
			int numCreated = 1;
			while (!isConverged(system, cells)) {
				// get a random cell
				Cell trialCell = getRandomCell(system, numCreated);
				// relax it and add it to the list
				Pair<Cell,Double> result = system.getObjFcn().getCellAndValue(trialCell);
				relaxedCell = result.getFirst();
				Double energy = result.getSecond();
				cells.put(relaxedCell, energy);
				// output the relaxed and unrelaxed (maybe) cells
				if (system.getWriteTempFiles())
					trialCell.writeCIF(outDirPath + "/" + trialCell.getLabel() + "_unrel.cif");
				relaxedCell.writeCIF(outDirPath + "/" + relaxedCell.getLabel() + "_rel.cif");
				// append it's energy to the output file
				out.write(relaxedCell.getLabel() + " " + energy + " ");
				for (Element e : system.getCompositionSpace().getElements())
					out.write(e.getSymbol() + " " + relaxedCell.getComposition().getStoichiometricUnit().get(e) + " ");
				out.write("\n");

				out.flush();
				numCreated++;
			}
			
			out.close();
		} catch (IOException x) {
			System.out.println("ERROR: RandomSearchMethod.doOpt(): " + x.getMessage());
			System.exit(-1);
		}
		
		// return the best cell
		Cell best = relaxedCell;
		for(Cell c : cells.keySet())
			if (cells.get(c) > cells.get(best))
				best = c;
			
		return best;
	}
	
	public static void usage() {
	      System.out.println("Usage: randsearch [OPTIONS] \n" +
	//		"  --d		Print debugging info\n" +
			"  --h		Print this help info and exit\n" +
	//		"  --I <file>	Input file in POSCAR format\n" +
			"  --outputTempFiles	Write intermediate files\n" +
			"  --compositionSpace <numElements> <Sym>* (<Num>*)*          System composition \n" +
			"  --objfcn gulpepa <potentialFile> <timeLimit(sec)> <cautious>  Objective function \n" +
			"  --outDir <directory> Write verbose output to directory \n" +
			"  --minLatticeLength <Num>			 In Angstroms \n" +
			"  --maxLatticeLength <Num>			 In Angstroms \n" +
			"  --minLatticeAngle <Num>			 In Degrees \n" +
			"  --maxLatticeAngle <Num>			 In Degrees \n" +
			"  --minNumAtoms <Num>		\n" +
			"  --maxNuAtoms <Num>	 \n" +
			"  --numStructures <Num>	 \n"
		);
	}

	
	public static void main(String args[]) {
		ArgumentParser aparser = new ArgumentParser(args);
		if (aparser.hasOption("h") || args.length == 0) { /* Have to have an input POSCAR */
			usage();
			System.exit(0);
		}
		
		OptiSystem sys = new OptiSystem(args);
		
		// make the output directory (throw error if it already exists)
		if (sys.getOutDir() != null) {
			if (! new File(sys.getOutDir()).mkdir()) {
				System.out.println("ERROR: Failed to make output directory " + sys.getOutDir());
				System.exit(-1);
			}
		}
		
		RandomSearchMethod method = new RandomSearchMethod(sys, sys.getOutDir());
		Cell optimalCell = method.doOpt();
		
		System.out.println(optimalCell);
	}
}
