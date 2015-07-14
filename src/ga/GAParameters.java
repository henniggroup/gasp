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

package ga;

import java.util.*;
import java.io.*;

import optimization.PDObjFcn;
import pdvisual.*;
import chemistry.CompositionSpace;
import chemistry.Element;
import utility.ArgumentParser;
import utility.Pair;
import utility.Utility;
import utility.Triplet;
import vasp.VaspIn;
import vasp.VaspOut;
import crystallography.Cell;
import crystallography.Isotropy;


// GAParameters is a singleton class which is responsible for parsing most of the 
// commandline and input file options.  It passes the ones which are particular
// to a certain algorithm on to that algorithm and stores the rest.
//
// The main algorithms and the Population are also owned by GAParameters since GAParameters
// has all the info about what particular algorithms should be used.  The result of this 
// is that the algorithms and Population are kind of pseudo-singletons in that they should
// only need to be accessed through GAParameters's getters.
//
// Argument-parsing is hierarchical.  For example, a command line option might be:
//   --objectiveFunction epa gulp gulppotl_file true  .
// The "else if" chain below recognizes the first token, '--objectiveFunction'.  It
// reads the next token, 'epa,' to decide that it should make an EnergyPerAtom 
// objective function.  The rest of the arguments, 'gulp gulppotl_file true,' are
// passed to the EnergyPerAtom constructor.  The EnergyPerAtom constructor reads the
// first token, 'gulp,' to see that it should use a GULP energy function, which
// happens to be called GulpEnergy.  The rest of the arguments, 'gulppotl_file true,'
// are passed on to GulpEnergy, which parses and stores the rest of the arguments.
// In this way, it is not necessary for GAParameters to be able to deal with all of the
// parameters of the whole program, but only the "top-level" ones.  
//
// For the sake of this file, parameters fall into one of two categories.  They are
// either always needed by the GA program (e.g. number of members in the population), 
// or they are only needed when using a particular algorithm (e.g. the gulp potential
// filename).  The first type are generally stored in a private variable in the
// GAParameters class.  The rest are stored in the class containing the relevant
// algorithm.
// 
// NB: In any of these cases, make sure to add an entry in the help statement in
// GAUtils.usage() as well as an entry in the toString() method below.
// Things to do when adding a new argument which is always necessary:
//   - make a private variable in which to store the parameter and set it 
//	   to a default value
//   - choose the commandline argument corresponding to that tag
//   - add an entry in the chain of 'else if's in setArgs().  make sure
//	   the tag starts with "--."  look at the others for examples.
//   - add a getter at the bottom of the file
//   - the option can be accessed where necessary by calling
//	   GAParameters.getParams().getMyNewOption();
// Things to do when adding an algorithm-specific argument:
//	 - These should be automatically passed on to the algorithm and handled there.
// Things to do when adding a new algorithm (say, a Variation):
//   - add a new if clause to the section of the "else-if" chain below which
//	   handles the "--variation" tag.  the clause should be similar to the 
//	   ones already there in that it creates a Variation object by passing
//	   the remaining arguments to the Variation's constructor and then adds
//	   that new Variation to the list.
//   - look to what's already here (e.g. StructureMut) to see how a Variation
//	   might want to handle its arguments.

public class GAParameters implements Serializable {
	static final long serialVersionUID = 1;

	private static GAParameters params;
	String inputFile;
	List<Pair<String,String[]>> argmap;
	
	// one Random for use by the whole program
	Random rand;
	
	// These are all the different parameters the program can take from the command line
	// along with their default values:
	// - verbosity
	private int verbosity = 4;
	// if dryRun is true, we don't save any output to disk
	private Boolean dryRun = false;
	private String runTitle = "default";
	// - hard constraints (in Angstroms and degrees):
	private double minInteratomicDistance = -1; 
	private List<Triplet<Element,Element,Double>> perSpeciesMIDs;
	private double maxLatticeLength = -1;
	private double minLatticeLength = -1;
	private double maxLatticeAngle = -1;
	private double minLatticeAngle = -1;
	private double maxCellHeight = Double.POSITIVE_INFINITY;
	private int maxNumAtoms = 1;
	private int minNumAtoms = 1;
	private int minNumSpecies = 1;
	private boolean doNonnegativityConstraint = true;
	private double dValue = 0;
	private List<Set<String>> notNearestNeighbors = null;
	private double maxNearestNeighborLength = 0;
	private boolean useNiggliReducedCell = false;
	private boolean use2DNiggliReducedCell = false;
	private boolean useSubstrate = false;
	private Cell substrate = null;
	private boolean writeHartkeFile = false;
	
	private boolean colorOutput = true;
	
	// default output dir is set in setArgs after we know the runTitle
	private String outDirName;
	private Boolean keepTempFiles = true;
	private boolean saveStateEachIter = false;
	
	// singleton
	private GARecord record = null;
	
	// These are the algorithms the program will use.
	private Selection sel = null;
	private Development dev = null;
	private Promotion pro = null;
	private List<String> objFcnArgs;
	private int minPopSize = 0;
	private int numCalcsInParallel = 1;
	// vars holds the variations that run on each generation to make the next generation.
	// each Variation is run the number of times stored in the corresponding location in
	// numToMake.
	private Vector<Variation> vars = null;
	private Vector<Double> initialVarProbs = null;
	private Vector<Double> endgameVarProbs = null;
	private Boolean inEndgame = false;
	private int endGameNumGens = 50;
	
	// population related options
	private int popSize = 40;
	private String redundancyGuardType = "none";
	private List<String> redundancyGuardArgs = null;
	private boolean useSurrogate = false;
	private List<String> surrogateArgs = null;
//	private static boolean unitsOnly = false;
	
	// constituents holds key-value pairs of the form (Atomic symbol, quantity) where the
	// quantity is the relative stoichiometry
	//private HashMap<String, Integer> constituents;
	//private Boolean fixStoich = false;
	private CompositionSpace compSpace;
	PDBuilder pdbuilder;
	// initialOrgCreators is a mapping between an OrgCreator and the (Integer) number of
	// organisms it should be used to create.
	private List<Pair<StructureOrgCreator,Integer>> initialOrgCreators = null;
	
	// convergence criteria
	Vector<ConvergenceCriterion> ccs = null;
	
	Boolean optimizeDensity = false;
	int numDenAdapt = 4;
	double weightDenAdapt = 0.5;
	
	private Cell[] seedStructures;
	
	// singleton
	private GAParameters() {
	    // initialize our prng
		rand = new Random();
		// initialize variations and convergence criteria vectors
		vars = new Vector<Variation>();
		initialVarProbs = new Vector<Double>();
		endgameVarProbs = new Vector<Double>();
		ccs = new Vector<ConvergenceCriterion>();
		// initialize constituents map
		//constituents = new HashMap<String, Integer>();
		// initialize initialOrgCreators map
		initialOrgCreators = new ArrayList<Pair<StructureOrgCreator,Integer>>();
		
		perSpeciesMIDs = new ArrayList<Triplet<Element,Element,Double>>();
	}
	
	// singleton
	public static GAParameters getParams() {
		if (params == null)
			params = new GAParameters();

		return params;
	}
	public static void setParams(GAParameters p) {
		params = p;
	}
	
	// print out a usage statement
	public static void usage(String errorStr, Boolean die) {
		System.out.println();
		GAOut.out().stdout("GASP: built " + utility.BuildDate.buildDate, GAOut.INFO);
		if (!errorStr.equals("")) {
			GAOut.out().stdout(errorStr, GAOut.CRITICAL);
			GAOut.out().stdout("", GAOut.CRITICAL);
		}
		System.out.println("Usage: ga <options> (--f <input-file> OR --r <resume-file>)");
		System.out.println("Arguments are case-insensitive.  All arguments can be passed in the input-file.");
		System.out.println("Optional flags [Current value]:");
		System.out.println("   --help : displays this message with parameters' default values");
		System.out.println("   --verbosity n : verbosity of output [" + GAParameters.getParams().getVerbosity() + "]");
		System.out.println("Genetic Algorithm");
		System.out.println("   --t runTitle : name the run");
		System.out.println("   --outDir : specify the output directory");
		System.out.println("   --dryRun : don't save any output to disk");
		System.out.println("   --keepTempFiles <true|false>");
		System.out.println("   --saveStateEachIter <true|false>");
		System.out.println("   --popSize <n> : use a non-initial population size of n");
		System.out.println("   --promotion <n> : promote n best structures (or the whole CH) to next gen");
//		System.out.println("   --constituents <fix stoichiometry?> <stoichiometry (e.g. Mn 1 O 2)>");
		System.out.println("   --compositionSpace <numElements> <Sym>*  :     System composition for phase diagram search");
		System.out.println("                   or <numElements> <Sym>* <amount of element1> <amount of element2> etc.");
		System.out.println("   --optimizeDensity <weight of adaptation> <num orgs to avg. over>");
		System.out.println("   --useRedundancyGuard <wholePopulation|perGeneration|both> <atomic misfit> <lattice misfit> <angle misfit> <use PBCs?>");
		System.out.println("   --useSurrogateModel <gulp header file> <potentials_specification file> <refit frequency>");
		System.out.println("   --endgameNumGens <n>");
		System.out.println("   --useNiggliReducedCell <true|false>");
		System.out.println("   --use2DNiggliReducedCell <true|false>");
		System.out.println("   --useSubstrate <true|false> <path to substrate POSCAR file");
		System.out.println("   --writeHartkeFile <boolean>");
		System.out.println("   --colorOutput <boolean>");
		System.out.println("Initial Population");
		System.out.println("   --initialPopulation <num> random givenVol <volumeperatom>");
		System.out.println("   --initialPopulation <num> random randomVol");
//		System.out.println("   --initialPopulation <num> resume <directory> <recalculate energies?>");
		System.out.println("   --initialPopulation <num> poscars <directory>");
		System.out.println("   --initialPopulation <num> manual");
		System.out.println("   --initialPopulation <num> units <numMols> <numAtoms_1>...<numAtoms_n> (<symbol_i> <x_i> <y_i> <z_i>)+ <numUnits_1>...<numUnits_n> <targetDensity> <densityTol> <unitsOnly?>");
		System.out.println("   --initialPopulation <num> supercell a b c maxll minll maxla minla maxh maxna minna <randomsocreator args>");
		System.out.println("Objective Functions");
		System.out.println("   --objectiveFunction cluster <padding length> <other obj fcn args from below...>");
		System.out.println("   --objectiveFunction surface <padding length> <other obj fcn args from below...>");
		System.out.println("   --objectiveFunction substrate <padding length> <other obj fcn args from below...>");
		System.out.println("   --objectiveFunction <epa/pd> gulp <gulp header file> <gulp potential file> <cautious?> <species needing a shell>");
		System.out.println("   --objectiveFunction <epa/pd> vasp <cautious?> <kpoints> <incar> <element potcar>+ ");
		System.out.println("   --objectiveFunction <epa/pd> ohmms <header> <footer> <cautious?>");
		System.out.println("   --objectiveFunction <epa/pd> lammps <potlFile> <units> <relax box?>");
		System.out.println("   --objectiveFunction <epa/pd> castep <cautious?> <kpointSpacing> <pressure> <paramFile> <element potcar>+ ");
		System.out.println("   --objectiveFunction <epa/pd> avogadro <avog header file>");
		System.out.println("   --objectiveFunction <epa/pd> dlpoly <loc> <potl>");
		System.out.println("   --objectiveFunction <epa/pd> mopac <execpath>");
		System.out.println("   --objectiveFunction <epa/pd> dftpp <dftpp_inputs> <cautious?> <element ppFile.fhi>*");
		System.out.println("   --objectiveFunction <epa/pd> generic");
		System.out.println("   --parallelize <numCalcsInParallel> <minPopSize>");
		System.out.println("Variation Algorithms");
		System.out.println("   --variation <percentage> <percentage> slicer <thicknessMean> <thicknessSigma> <majorShiftFrac> <minorShiftFrac> <maxAmplitude> <maxFreq> <growParents?> <doublingProb>");
		System.out.println("   --variation <percentage> <percentage> structureMut <rate> <sigmaAtoms> <sigmaLattice>");
		System.out.println("   --variation <percentage> <percentage> permutation <meanSwaps> <sigmaSwaps> <pairsToSwap (e.g. Mn-O)>");
		System.out.println("   --variation <percentage> <percentage> numStoichsMut <meanNumAtoms> <sigmaNumAtoms>");
//		System.out.println("   --variation <percentage> <percentage> supercell <re-relax children?>");
		System.out.println("Selection Algorithms");
		System.out.println("   --selection probDist <numParents> <selectionPower>");
		System.out.println("Convergence Criteria");
		System.out.println("   --convergenceCriterion maxFunctionEvals <n>");
		System.out.println("   --convergenceCriterion maxNumGens <n>");
		System.out.println("   --convergenceCriterion maxNumGensWOImpr <n> <dValue>");
		System.out.println("   --convergenceCriterion valueAchieved <maximum acceptable energy>");
		System.out.println("   --convergenceCriterion foundStructure <CIF filename> <rGuard misfits>");
		System.out.println("Hard Constraints");
		System.out.println("   --minInteratomicDistance d : minimum interatomic distance (Angstroms)");
		System.out.println("   --perSpeciesMID <symbol1> <symbol2> <distance>");
		System.out.println("   --maxLatticeLength d : maximum lattice vector length (Angstroms)");
		System.out.println("   --minLatticeLength d : minimum lattice vector length (Angstroms)");
		System.out.println("   --maxLatticeAngle d : maximum lattice angle (Degrees)");
		System.out.println("   --minLatticeAngle d : minimum lattice angle (Degrees)");
		System.out.println("   --maxCellHeight d : maximum height of cell in z-direction");
		System.out.println("   --maxNumAtoms n");
		System.out.println("   --minNumAtoms n");
		System.out.println("   --minNumSpecies n");
		System.out.println("   --doNonnegativityConstraint <boolean>");
		System.out.println("   --dValue x : discard organisms within a value of x of each other");
		
		if (die)
			System.exit(1);
	}
	
	// We combine the mappings created by parsing the input file and the
	// command-line arguments so that the command-line arguments take 
	// precedence over those in the input file.  thus we can
	// use the same code to set the private variables.
	public void setArgs(String[] args) {		
		
		// parse the input arguments
		ArgumentParser aParser;
		String inputFileName = null;
		for (int i = 0; i < args.length - 1; i++) {
			if (args[i].equals("--f"))
				inputFileName = args[i+1];
		}
		
		if(inputFileName == null) {
			if (verbosity >= 3)
				System.out.println("WARNING: No input file given.");
			aParser = new ArgumentParser(args);
		} else {
			aParser = new ArgumentParser(args, inputFileName);
		}
		
		// parse the combined set of arguments
		for (Pair<String,List<String>> p : aParser.getOptions()) {
			String flag = p.getFirst();
			List<String> arguments = p.getSecond();
			
			if (flag.equals("--help"))
				usage("", true);
			else if (flag.equalsIgnoreCase("verbosity"))
				verbosity = Integer.parseInt(arguments.get(0));
			else if (flag.equalsIgnoreCase("runTitle"))
				runTitle = arguments.get(0);
			else if (flag.equalsIgnoreCase("outDirName")) {
				outDirName = "";
				for (String s : arguments)
					outDirName += s;
			} else if (flag.equalsIgnoreCase("minInteratomicDistance"))
				minInteratomicDistance = Double.parseDouble(arguments.get(0));
			else if (flag.equalsIgnoreCase("perSpeciesMID")) {
				if (arguments.size() < 3)
					usage("Error: not enough arguments to perSpeciesMID", true);
				Element a = Element.getElemFromSymbol(arguments.get(0));
				Element b = Element.getElemFromSymbol(arguments.get(1));
				Double dist = Double.parseDouble(arguments.get(2));
				perSpeciesMIDs.add(new Triplet<Element,Element,Double>(a,b,dist));
			} else if (flag.equalsIgnoreCase("maxLatticeLength"))
				maxLatticeLength = Double.parseDouble(arguments.get(0));
			else if (flag.equalsIgnoreCase("minLatticeLength"))
				minLatticeLength = Double.parseDouble(arguments.get(0));
			else if (flag.equalsIgnoreCase("maxLatticeAngle"))
				maxLatticeAngle = Double.parseDouble(arguments.get(0));
			else if (flag.equalsIgnoreCase("minLatticeAngle"))
				minLatticeAngle = Double.parseDouble(arguments.get(0));
			else if (flag.equalsIgnoreCase("maxCellHeight"))
				maxCellHeight = Double.parseDouble(arguments.get(0));
			else if (flag.equalsIgnoreCase("maxNumAtoms"))
				maxNumAtoms = Integer.parseInt(arguments.get(0));
			else if (flag.equalsIgnoreCase("minNumAtoms"))
				minNumAtoms = Integer.parseInt(arguments.get(0));
			else if (flag.equalsIgnoreCase("minNumSpecies"))
				minNumSpecies = Integer.parseInt(arguments.get(0));
			else if (flag.equalsIgnoreCase("doNonnegativityConstraint"))
				doNonnegativityConstraint = Boolean.parseBoolean(arguments.get(0));
			else if (flag.equalsIgnoreCase("dValue"))
				dValue = Double.parseDouble(arguments.get(0));
			else if (flag.equalsIgnoreCase("endgameNumGens"))
				endGameNumGens = Integer.parseInt(arguments.get(0));
			else if (flag.equalsIgnoreCase("useNiggliReducedCell"))
				useNiggliReducedCell = Boolean.parseBoolean(arguments.get(0));
			else if (flag.equalsIgnoreCase("use2DNiggliReducedCell"))
				use2DNiggliReducedCell = Boolean.parseBoolean(arguments.get(0));
			else if (flag.equalsIgnoreCase("use2DNiggliReducedCell")) {
				useSubstrate = Boolean.parseBoolean(arguments.get(0));
				if (useSubstrate)
					substrate = VaspOut.getPOSCAR(arguments.get(1));
			} else if (flag.equalsIgnoreCase("popSize"))
				popSize = Integer.parseInt(arguments.get(0));
			else if (flag.equalsIgnoreCase("keepTempFiles"))
				keepTempFiles = Boolean.parseBoolean(arguments.get(0));
			else if (flag.equalsIgnoreCase("saveStateEachIter"))
				saveStateEachIter = Boolean.parseBoolean(arguments.get(0));
			else if (flag.equalsIgnoreCase("useRedundancyGuard")) {
				redundancyGuardType = arguments.get(0);
				redundancyGuardArgs = Utility.subList(arguments, 1);
			} else if (flag.equalsIgnoreCase("useSurrogateModel")) {
				useSurrogate = true;
				surrogateArgs =  arguments;
			} else if (flag.equalsIgnoreCase("optimizeDensity")) {
				optimizeDensity = true;
				if (optimizeDensity) {
					if (arguments.size() < 2)
						usage("Error: not enough arguments to optimizeDensity", true);
					weightDenAdapt = Double.parseDouble(arguments.get(0));
					numDenAdapt = Integer.parseInt(arguments.get(1));
				}
			}
			else if (flag.equalsIgnoreCase("notNearestNeighbors")) {
				maxNearestNeighborLength = Double.parseDouble(arguments.get(0));
				notNearestNeighbors = GAUtils.parsePairs(Utility.subList(arguments, 1));
			}
			else if (flag.equalsIgnoreCase("dryRun"))
				dryRun = Boolean.parseBoolean(arguments.get(0));
			else if (flag.equalsIgnoreCase("objectiveFunction")) {
				if (arguments.size() < 2)
					usage("Not enough parameters given to --objectiveFunction", true);
				objFcnArgs = arguments;
			}
			else if (flag.equalsIgnoreCase("parallelize")) {
				numCalcsInParallel = Integer.parseInt(arguments.get(0));
				minPopSize = Integer.parseInt(arguments.get(1));
			}
			else if (flag.equalsIgnoreCase("selection")) {
				if (arguments.get(0).equalsIgnoreCase("probDist"))
					sel = new ProbDistSelection( Utility.subList(arguments, 1));
				else if (arguments.get(0).equalsIgnoreCase("off")) 
					sel = null;
				else 
					usage("Unknown selection function " + arguments.get(0), true);
			}
			// kind of a hack. (it's because java collections don't support multiple identical keys)
			else if (flag.toLowerCase().startsWith("variation")) {
				if (arguments.size() < 3)
					usage("Malformed input to variation", true);
				String variation = arguments.get(2);
				initialVarProbs.add(Double.parseDouble(arguments.get(0)));
				endgameVarProbs.add(Double.parseDouble(arguments.get(1)));
				if (variation.equalsIgnoreCase("slicer"))
					vars.add(new Slicer(Utility.subList(arguments, 3)));
				else if (variation.equalsIgnoreCase("structureMut")) 
					vars.add(new StructureMut(Utility.subList(arguments, 3)));
				else if (variation.equalsIgnoreCase("permutation")) 
					vars.add(new Permutation(Utility.subList(arguments, 3)));
				else if (variation.equalsIgnoreCase("numStoichsMut")) 
					vars.add(new NumStoichsMut(Utility.subList(arguments, 3)));
		//		else if (variation.equalsIgnoreCase("supercell")) 
		//			vars.add(new SupercellVariation(Utility.subList(arguments, 3)));
				else 
					usage("Unknown variation function " + variation, true);
			}
			// similar to the variation handling above
			else if (flag.toLowerCase().startsWith("convergencecriterion")) {
				if (arguments.get(0).equalsIgnoreCase("maxFunctionEvals"))
					ccs.add(new NumFunctionEvalsCC( Utility.subList(arguments, 1) ));
				else if (arguments.get(0).equalsIgnoreCase("maxNumGens"))
					ccs.add(new NumGensCC( Utility.subList(arguments, 1)));
				else if (arguments.get(0).equalsIgnoreCase("maxNumGensWOImpr")) 
					ccs.add(new NumGensWOImprCC( Utility.subList(arguments, 1)));
				else if (arguments.get(0).equalsIgnoreCase("valueAchieved")) 
					ccs.add(new ValueAchievedCC( Utility.subList(arguments, 1)));
				else if (arguments.get(0).equalsIgnoreCase("foundStructure")) 
					ccs.add(new FoundStructureCC( Utility.subList(arguments, 1)));
				else 
					usage("Unknown convergence criterion " + arguments.get(0), true);
			}
			else if (flag.equalsIgnoreCase("promotion")) {
				pro = new Promotion(arguments);
			}
			else if (flag.toLowerCase().startsWith("initialpopulation")) {
				if (arguments.size() < 2)
					usage("Malformed input in initialpopulation", true);
				Integer numOrgs = new Integer(Integer.parseInt(arguments.get(0)));
				String creatorType = arguments.get(1);

				if (creatorType.equalsIgnoreCase("random")) {
					initialOrgCreators.add(new Pair<StructureOrgCreator,Integer>(new RandomSOCreator( Utility.subList(arguments, 2)), numOrgs));
//				} else if (creatorType.equalsIgnoreCase("resume")) {
//					initialOrgCreators.put(new ResumeSOCreator(GAUtils.subArray(values, 2)), numOrgs);
				} else if (creatorType.equalsIgnoreCase("manual")) {
					initialOrgCreators.add(new Pair<StructureOrgCreator,Integer>(new ManualSOCreator( Utility.subList(arguments, 2)), numOrgs));
				} else if (creatorType.equalsIgnoreCase("poscars")) {
					initialOrgCreators.add(new Pair<StructureOrgCreator,Integer>(new GivenSOCreator( Utility.subList(arguments, 2)), numOrgs));
				} else if (creatorType.equalsIgnoreCase("units")) {
					initialOrgCreators.add(new Pair<StructureOrgCreator,Integer>(new UnitsSOCreator( Utility.subList(arguments, 2)), numOrgs));
				} else if (creatorType.equalsIgnoreCase("supercell")) {
					initialOrgCreators.add(new Pair<StructureOrgCreator,Integer>(new SupercellSOCreator( Utility.subList(arguments, 2)), numOrgs));
				} else {
					usage("Unrecognized population type " + creatorType, true);
				}					
			} else if (flag.equalsIgnoreCase("compositionSpace")) {
				List<String> csArgs = new LinkedList<String>();
				for (String s : p.getSecond())
					csArgs.add(s);
				compSpace = new CompositionSpace(csArgs, false);
			} else if (flag.equalsIgnoreCase("writeHartkeFile")) {
				writeHartkeFile = Boolean.parseBoolean(arguments.get(0));
			} else if (flag.equalsIgnoreCase("colorOutput")) {
				colorOutput = Boolean.parseBoolean(arguments.get(0));
			}
			// we deal with the input file separately
			else if (!flag.equalsIgnoreCase("f") && verbosity >= 1)
				System.out.println("Ignoring unrecognized flag " + flag);
		}
		
		// make the development object
		dev = new StructureDev();

		// default output directory:
		if (outDirName == null)
			outDirName = new String("garun_" + runTitle);
		
		
		checkInputs();
	}
	
	public void checkInputs() {
		// some sanity checks, not systematic
		if (compSpace == null) 
			usage("No --compositionSpace option passed.",true);
		if (minLatticeLength > maxLatticeLength)
			usage("Error: minLatticeLength > maxLatticeLength.", true);
		if (minLatticeAngle > maxLatticeAngle)
			usage("Error: minLatticeAngle > maxLatticeAngle.", true);
		if (minNumAtoms > maxNumAtoms)
			usage("ERROR: minNumAtoms > maxNumAtoms.", true);
		if (minNumAtoms < 1)
			usage("ERROR: minNumAtoms should be at least 1.", true);
		if (!(redundancyGuardType.equalsIgnoreCase("both")
				|| redundancyGuardType.equalsIgnoreCase("none")
				|| redundancyGuardType.equalsIgnoreCase("wholePopulation")
				|| redundancyGuardType.equalsIgnoreCase("perGeneration")))
					usage("ERROR: Unrecognized RedundancyGuard type " + redundancyGuardType, true);
		if (objFcnArgs == null || objFcnArgs.size() < 2)
			usage("ERROR: Need an objective function.", true);
		String objFcnType = objFcnArgs.get(0);
		if (objFcnType.equalsIgnoreCase("pd") && this.getCompSpace().getNumDimensions() < 2)
			usage("ERROR: Can't use pd objFun w/ < 2 dimensions.", true);
		if (vars == null || vars.size() == 0)
			System.out.println("WARNING: Given no variation operators.");
	}
	
	public void setSeedGeneration(Cell[] initialPop) {
		seedStructures = initialPop;
	}

	
	public String toString() {
		StringBuilder result = new StringBuilder();
		String newline = GAUtils.newline();
		
		result.append("runTitle: " + runTitle + newline);
		result.append("inputFile: " + inputFile + newline);
		if (dryRun)
			result.append("dryRun" + newline);
		else 
			result.append("outDirName: " + outDirName + newline);
		result.append("verbosity: " + verbosity + newline);
		
		result.append(newline);
		result.append("minInteratomicDistance: " + minInteratomicDistance + newline);
		result.append("maxLatticeLength: " + maxLatticeLength + newline);
		result.append("minLatticeLength: " + minLatticeLength + newline);
		result.append("maxLatticeAngle: " + maxLatticeAngle + newline);
		result.append("minLatticeAngle: " + minLatticeAngle + newline);
		result.append("dValue: " + dValue + newline);

		result.append(newline);
		// actually should really not try to get an objective function until after we've already
		// initialized GAParameters
	//	if (getObjectiveFunctionInstance(null) != null)
	//		result.append(getObjectiveFunctionInstance(null).toString() + newline);
	//	else
	//		result.append("objectiveFunction: off" + newline);
		
		if (getSelection() != null)
			result.append(getSelection().toString() + newline);
		else
			result.append("selection: off" + newline);
		
		result.append(getPromotion() + newline);
		
		for (int i = 0; i < vars.size(); i++) {
			Variation v = vars.get(i);
			result.append(v + " (probability: " + params.getVarProb(i) + ") "+ newline);
		}
		
		if (getDevelopment() != null)
			result.append(getDevelopment().toString() + newline);
		else
			result.append("development: off" + newline);
		
/*		result.append("Constituents: ");
		Set<String> cons = constituents.keySet();
		Iterator<String> i = cons.iterator();
		while (i.hasNext()) {
			String symbol = i.next();
			result.append(symbol + " " + constituents.get(symbol) + " ");
		}
		result.append(newline);
		*/
		result.append("Composition space: \n");
		result.append(compSpace.toString() + "\n");
		
		result.append("StructureOrgCreators: " + newline);
		for (Pair<StructureOrgCreator,Integer> soc : initialOrgCreators)
			result.append(soc.getFirst() + newline);
		
		result.append("Convergence criteria: " + newline);
		for (ConvergenceCriterion cc : ccs)
			result.append(cc + newline);
		
		return result.toString();
	}
	
	public Generation makeEmptyGeneration() {
		return new Structures();
	}
	
	public void seedIDs(int s) {
		currentID = s;
	}
	
	public void setInEndgame(Boolean b) {
		inEndgame = b;
	}
	
	//singleton
	public Object clone() throws CloneNotSupportedException {
	    throw new CloneNotSupportedException(); 
	}
	
	public ObjectiveFunction getObjectiveFunctionInstance(Organism o) {
		return ObjFcnFactory.getObjectiveFunctionInstance(o, objFcnArgs);
	}
	
	public boolean doingPDRun() {
		for (String s : objFcnArgs)
			if (s.equalsIgnoreCase("pd"))
				return true;
		
		return false;
	}
	
	/********************* All the Getters go below here ********************/
	// All of these are either set during the commandline/input file parsing
	// or else they assume their default values set above.
	
	private int currentID = 0;
	public int getNewOrgID() {
		currentID++;
		return currentID;
	}

	public Selection getSelection() {
		return sel;
	}
	
	public Vector<Variation> getVariations() {
		return vars;
	}
	
	public Promotion getPromotion() {
		return pro;
	}
	
	public double getVarProb(int varNum) {
		if (inEndgame)
			return endgameVarProbs.get(varNum);
		else
			return initialVarProbs.get(varNum);
	}
	
	public Development getDevelopment() {
		return dev;
	}

	public double getMaxLatticeAngleDegrees() {
		return maxLatticeAngle;
	}

	public double getMaxLatticeLength() {
		return maxLatticeLength;
	}
	
	public double getMaxCellHeight() {
		return maxCellHeight;
	}

	public double getMinInteratomicDistance() {
		return minInteratomicDistance;
	}
	
	public List<Triplet<Element,Element,Double>> getPerSpeciesMIDs() {
		return perSpeciesMIDs;
	}

	public double getMinLatticeAngleDegrees() {
		return minLatticeAngle;
	}

	public double getMinLatticeLength() {
		return minLatticeLength;
	}
	
	public int getMaxNumAtoms() {
		return maxNumAtoms;
	}
	
	public int getMinNumAtoms() {
		return minNumAtoms;
	}
	
	public int getMinNumSpecies() {
		return minNumSpecies;
	}
	
	public boolean getDoNonnegativityConstraint() {
		return doNonnegativityConstraint;
	}
	
	public String getRunTitle() {
		return runTitle;
	}

	public double getDValue() {
		return dValue;
	}
	
	public List<Set<String>> getNotNearestNeighbors() {
		return notNearestNeighbors;
	}
	
	public double getMaxNearestNeighborLength() {
		return maxNearestNeighborLength;
	}

	public int getVerbosity() {
		return verbosity;
	}
	
	public Random getRandom() {
		return rand;
	}
	
	public Boolean getOptimizeDensity() {
		return optimizeDensity;
	}
	
	public int getNumDenAdapt() {
		return numDenAdapt;
	}
	
	public double getWeightDenAdapt() {
		return weightDenAdapt;
	}
	
	public String getOutDirName() {
		return outDirName;
	}
	
	public Boolean getDryRun() {
		return dryRun;
	}
	
	public Cell[] getSeedStructures() {
		return seedStructures;
	}
	
/*	public HashMap<String, Integer> getConstituents() {
		return constituents;
	}
	
	public Boolean getFixStoich() {
		return fixStoich;
	}
	*/
	
	public CompositionSpace getCompSpace() {
		return compSpace;
	}
	
	public Vector<ConvergenceCriterion> getCCs() {
		return ccs;
	}
	
	public int getPopSize() {
		return popSize;
	}
	
	public String getRedundancyGuardType() {
		return redundancyGuardType;
	}
	
	public List<String> getRedundancyGuardArgs() {
		return new ArrayList<String>(redundancyGuardArgs);
	}
	
	public boolean usingSurrogateModel() {
		return useSurrogate;
	}
	
	public List<String> getSurrogateArgs() {
		return new ArrayList<String>(surrogateArgs);
	}
	
	public String getTempDirName() {
		if (outDirName == null)
			return "/tmp";
		else
			return outDirName + "/temp";
	}
	
	public int getEndGameNumGens() {
		return endGameNumGens;
	}
	
	public boolean usingNiggliReducedCell() {
		return useNiggliReducedCell;
	}
	
	public boolean using2DNiggliReducedCell() {
		return use2DNiggliReducedCell;
	}
	
	public boolean usingSubstrate() {
		return useSubstrate;
	}
	
	public Cell getSubstrate() {
		return substrate;
	}
	
	public Boolean getKeepTempFiles() {
		return keepTempFiles;
	}
	
	public boolean getSaveStateEachIter() {
		return saveStateEachIter;
	}
	
	public int getMinPopSize() {
		return minPopSize;
	}
	
	public PDBuilder getPDBuilder() {
		if (pdbuilder == null) {
			Map<Element, Double> cps = new HashMap<Element,Double>();
			for (Element e : compSpace.getElements())
				cps.put(e, 0.0);
			pdbuilder = new PDBuilder(new LinkedList<IComputedEntry>(), compSpace.getElements(), cps );
		}
		
		return pdbuilder;
	}
	
	public int getNumCalcsInParallel() {
		return numCalcsInParallel;
	}
	
	public List<Pair<StructureOrgCreator,Integer>> getInitialOrgCreators() {
		return initialOrgCreators;
	}
	
	// singleton
	public GARecord getRecord() {
		if (record == null)
			record = new GARecord();
		return record;
	}
	
	public void setMinInteratomicDistance(double d) {
		minInteratomicDistance = d;
	}
	public void setMaxLatticeLength(double d) {
		maxLatticeLength = d;
	}
	public void setMinLatticeLength(double d) {
		minLatticeLength = d;
	}
	public void setMaxLatticeAngle(double d) {
		maxLatticeAngle = d;
	}
	public void setMinLatticeAngle(double d) {
		minLatticeAngle = d;
	}
	public void setMaxNumAtoms(int i) {
		maxNumAtoms = i;
	}
	public void setMinNumAtoms(int i) {
		minNumAtoms = i;
	}
	public void setCompSpace(CompositionSpace c) {
		compSpace = c;
	}
	public void setVerbosity(int v) {
		verbosity = v;
	}
		
	public class GARecord implements Serializable {
		
		static final long serialVersionUID = 1;
		
		private String newline = GAUtils.newline();
		
		// We put all our output in the directory outDirName.  This output consists of
		// all the CIF files of our organisms as well as a file which describes which
		// are in each generation and records their energies, outFileName.
		private String outDirName;
		private String outFileName = "index";
		private String paramFileName = "parameters";
		private String tempDirName;
		private File outFile;
		private File tempDir;
		private File outDir;
		private File paramFile;
		
		private int currentGenNum = 0;
		private Generation currentGen = null;
		
		//
		double bestDenEstimate;
		
		// singleton
		private GARecord() {		
			GAParameters params = GAParameters.getParams();
			
			String runTitle = params.getRunTitle();
			outDirName = params.getOutDirName();
			tempDirName = params.getTempDirName();

			// make the output directory
			outDir = new File(outDirName);
			if (!params.getDryRun() && !outDir.mkdir()) 
				GAParameters.usage("ERROR: Can't create directory " + outDirName, true);
			
			// make the output files
			outFile = new File(outDir, outFileName);
			
			// make the temp directory
			tempDir = new File(tempDirName);
			tempDir.mkdir();
			
			// save the run parameters to a file
			paramFile = new File(outDir, paramFileName);
			GAUtils.writeStringToFile(GAParameters.getParams().toString(), paramFile, false);
			
			if (verbosity >= 1) {
				System.out.print("Run " + runTitle + ": ");
				if (params.getDryRun())
					System.out.println("dry run");
				else
					System.out.println("outputting to directory " + outDirName);
			}
			
		}
		
		public String toString() {
			return "GARecord. outDir: " + outDirName;
		}
		
		// singleton
		public Object clone() throws CloneNotSupportedException {
		    throw new CloneNotSupportedException(); 
		}

		private String getGenSummary(Generation g) {
			StringBuilder result = new StringBuilder();
			
			result.append("Population summary (generation " + currentGenNum + "):" + newline);
			
			Iterator<Organism> i = g.iterator();
			if (!i.hasNext())
				return new String("Empty population?");		
			
			// report the lowest energy:
			Organism bestOrg = g.getNthBestOrganism(1);
			result.append("Best value: " + bestOrg.getValue() + " (organism " + bestOrg.getID() + ")");
			result.append(GAUtils.newline());
			result.append("Best density estimate: " + bestDenEstimate + " per atom");
			result.append(GAUtils.newline());
			result.append("Number of energy calculations thus far: " + ObjectiveFunction.getNumCalculations());
			
			return result.toString();
		}
		
		private String makePOSCARPath(StructureOrg s) {
			File f = new File(outDir, Integer.toString(s.getID()) + ".POSCAR");
			return f.getPath();
		}
		
		private String makeFindSymPath(StructureOrg s) {
			File f = new File(outDir, Integer.toString(s.getID()) + ".fso");
			return f.getPath();
		}
		
		// this writes status info to the screen, saves data to files, and updates
		// the generation counter.  it should be called by the main algorithm at the
		// end of each generation.
		public void finishGen(Generation g) {
			GAParameters params = GAParameters.getParams();
			
			// write the generation header (generation x N)
			GAUtils.writeStringToFile("generation " + Integer.toString(currentGenNum) + " " + g.getNumOrganisms() + newline, outFile, true);
			
			// save structures and findsym outputs
			for (Organism o : g) {
				// assume here that our organisms are StructureOrgs
				StructureOrg s = (StructureOrg)(o);
				// save the structure
				VaspIn.writePoscar(s.getCell(), makePOSCARPath(s), false);
				// save the findsym output
				File outFindSym = new File(makeFindSymPath(s));
				GAUtils.writeStringToFile(s.getCell().getFSOOutput(), outFindSym, false);
			}
			
			// write energy-sorted index (i.e. lowest to highest)
			for (Organism o : g.getOrganismsSorted()) {
				StructureOrg s = (StructureOrg) o;
				StringBuilder info = new StringBuilder();
				info.append(Integer.toString(s.getID()) + " ");
				info.append(Double.toString(s.getValue()) + " ");
				info.append(makePOSCARPath(s) + newline);
				GAUtils.writeStringToFile(info.toString(), outFile, true);
			}
			
			
			// update the best density info
			int nAdapt = params.getNumDenAdapt();
			double weight = params.getWeightDenAdapt();
			// find the average density of the nAdapt best organisms
			double avgDen = 0;
			for (int j = 1; j <= nAdapt; j++) {
				Cell currentStruct = ((StructureOrg)g.getNthBestOrganism(j)).getCell();
				avgDen += currentStruct.getVolume() / currentStruct.getNumSites();
			}
			avgDen /= nAdapt;
			// set the new best volume.  if this is the first generation, we don't use the weight
			if (getRecord().getGenNum() == 0)
				bestDenEstimate = avgDen;
			else {
				bestDenEstimate *= (1-weight);
				bestDenEstimate += weight * avgDen;
			}
			
			updateInEndgame(g);
			
			// update phase diagram stuff
			if (params.doingPDRun()) {
				String pdbuilder_out_fname = outDirName + "/gen" + currentGenNum + ".pdb.tgz";
				for (Organism o : g.organisms)
					if (! params.getPDBuilder().containsEntry((StructureOrg)o))
						params.getPDBuilder().addEntry((StructureOrg)o);
				Utility.writeSerializable((params.getPDBuilder().getPDData()), pdbuilder_out_fname);
			//	new TernPD3D(params.getPDBuilder().getPDData());
			}
			
			if (verbosity >= 1)
				System.out.println(getGenSummary(g));
			
			currentGenNum++;
			currentGen = g;
			
			if (params.getSaveStateEachIter()) {
				GAParameters.getParams().getInitialOrgCreators().clear();
				String save_out_fname = outDirName + "/gen" + currentGenNum + ".save.tgz";
				Utility.writeSerializable(params, save_out_fname);
				if (verbosity >= 3)
					System.out.println("Wrote save file to "+ save_out_fname);
			}
		}
		
		double bestEnergy = 0;
		int numGensWOImpr = 0;
		private void updateInEndgame(Generation g) {
			double newBest = g.getExtremeValues()[0];
			if (newBest < bestEnergy)
				numGensWOImpr = 0;
			else
				numGensWOImpr++;
			
			if (numGensWOImpr > GAParameters.getParams().getEndGameNumGens())
				GAParameters.getParams().setInEndgame(true);
			else
				GAParameters.getParams().setInEndgame(false);
			
			bestEnergy = newBest;
		}
		
		public int getGenNum() {
			return currentGenNum;
		}
		
		public File getOutFile() {
			return outFile;
		}
		
		public void setCurrentGenNum(int n) {
			currentGenNum = n;
		}
		
		public Generation getCurrentGen() {
			return currentGen;
		}
		
		public double getBestDensityEstimate() {
			return bestDenEstimate;
		}
	
		public void cleanup() {	
			if(!GAParameters.getParams().getKeepTempFiles()) {
				// delete temporary files
				if (tempDir != null && tempDir.exists())
					try {
						GAUtils.deleteDirectory(tempDir);
					} catch (IOException x) {
						GAOut.out().stdout("Error deleting temp directory: " + x.getMessage(), GAOut.WARNING);
					}
			}
		}

	}

	public boolean getWriteHartkeFile() {
		return writeHartkeFile;
	}

	public boolean getColorOutput() {
		return colorOutput;
	}
	
	public String getHartkeOutFile() {
		return outDirName + "/hartke.txt";
	}
}
