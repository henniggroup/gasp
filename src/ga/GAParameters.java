/* Genetic algorithm for crystal structure prediction.  Will Tipton.  Ceder Lab at MIT. Summer 2007 */

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
	private double maxLatticeLength = -1;
	private double minLatticeLength = -1;
	private double maxLatticeAngle = -1;
	private double minLatticeAngle = -1;
	private int maxNumAtoms = 1;
	private int minNumAtoms = 1;
	private boolean doNonnegativityConstraint = true;
	private double dValue = 0;
	private List<Set<String>> notNearestNeighbors = null;
	private double maxNearestNeighborLength = 0;
	private boolean useNiggliReducedCell = false;
	
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
	private String[] objFcnArgs;
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
	private String[] redundancyGuardArgs = null;
	private boolean useSurrogate = false;
	private String[] surrogateArgs = null;
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
		if (!errorStr.equals("")) {
			System.out.println(errorStr);
			System.out.println("");
		}
		System.out.println("Usage: GeneticAlgorithm <options> (--f <input-file> OR --r <resume-file>)");
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
//		System.out.println("   --constituents <fix stoichiometry?> <stoichiometry (e.g. Mn 1 O 2)>");
		System.out.println("   --compositionSpace <numElements> <Sym>* (<Num>*)*          System composition");
		System.out.println("   --optimizeDensity <weight of adaptation> <num orgs to avg. over>");
		System.out.println("   --useRedundancyGuard <wholePopulation|perGeneration|both> <atomic misfit> <lattice misfit> <angle misfit>");
		System.out.println("   --useSurrogateModel <gulp header file> <potentials_specification file> <refit frequency>");
		System.out.println("   --endgameNumGens <n>");
		System.out.println("   --useNiggliReducedCell <true|false>");
		System.out.println("Initial Population");
		System.out.println("   --initialPopulation <num> random givenVol <volumeperatom>");
		System.out.println("   --initialPopulation <num> random randomVol");
//		System.out.println("   --initialPopulation <num> resume <directory> <recalculate energies?>");
		System.out.println("   --initialPopulation <num> fromCifs <directory>");
		System.out.println("   --initialPopulation <num> manual");
		System.out.println("   --initialPopulation <num> units <symbol 1> <x 1> <y 1> <z 1> ... <true|false> <numUnits> <targetDensity> <densityTol>");
		System.out.println("Objective Functions");
		System.out.println("   --objectiveFunction <epa/pd> gulp <gulp header file> <gulp potential file> <cautious?> <species needing a shell>");
		System.out.println("   --objectiveFunction <epa/pd> vasp <cautious?> <kpoints> <incar> <element potcar>+ ");
		System.out.println("   --objectiveFunction <epa/pd> ohmms <header> <footer> <cautious?>");
		System.out.println("   --objectiveFunction <epa/pd> lammps <potlFile>");
		System.out.println("   --objectiveFunction <epa/pd> castep <cautious?> <kpointSpacing> <pressure> <paramFile> <element potcar>+ ");
		System.out.println("   --objectiveFunction <epa/pd> avogadro <avog header file>");
		System.out.println("   --parallelize <numCalcsInParallel> <minPopSize>");
		System.out.println("Variation Algorithms");
		System.out.println("   --variation1 <percentage> <percentage> slicer <thicknessMean> <thicknessSigma> <majorShiftFrac> <minorShiftFrac> <maxAmplitude> <maxFreq> <growParents?>");
		System.out.println("   --variation2 <percentage> <percentage> structureMut <rate> <sigmaAtoms> <sigmaLattice>");
		System.out.println("   --variation3 <percentage> <percentage> permutation <meanSwaps> <sigmaSwaps> <pairsToSwap (e.g. Mn-O)>");
		System.out.println("   --variation4 <percentage> <percentage> numStoichsMut <meanNumAtoms> <sigmaNumAtoms>");
		System.out.println("   --variation5 <percentage> <percentage> supercell");
		System.out.println("Selection Algorithms");
		System.out.println("   --selection probDist <numParents> <selectionPower>");
		System.out.println("Convergence Criteria");
		System.out.println("   --convergenceCriterion1 maxFunctionEvals <n>");
		System.out.println("   --convergenceCriterion2 maxNumGens <n>");
		System.out.println("   --convergenceCriterion3 maxNumGensWOImpr <n> <dValue>");
		System.out.println("   --convergenceCriterion4 valueAchieved <maximum acceptable energy>");
		System.out.println("   --convergenceCriterion5 foundStructure <CIF filename>");
		System.out.println("Hard Constraints");
		System.out.println("   --minInteratomicDistance d : minimum interatomic distance (Angstroms)");
		System.out.println("   --maxLatticeLength d : maximum lattice vector length (Angstroms)");
		System.out.println("   --minLatticeLength d : minimum lattice vector length (Angstroms)");
		System.out.println("   --maxLatticeAngle d : maximum lattice angle (Degrees)");
		System.out.println("   --minLatticeAngle d : minimum lattice angle (Degrees)");
		System.out.println("   --maxNumAtoms n");
		System.out.println("   --minNumAtoms n");
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
		// parse the commandline arguments
	//	ArgumentParser aParser = new ArgumentParser(args);
		argmap = (new ArgumentParser(args)).getCLParserData();
	
		// parse the input file arguments
		if(!isSet("--f")) {
			if (verbosity >= 2)
				System.out.println("Really no input file?");
		} else {
			inputFile = getValues("--f")[0];
			List<Pair<String,String[]>> inputFileArgs = parseInputFile();
			
			// combine the two maps such that the commandline args take precedence
			inputFileArgs.addAll(argmap);
			argmap = inputFileArgs;
		}
		
		if (argmap == null)
			usage("", true);
		
		// TODO: switch everything to use aParser
		if (!isSet("--compositionSpace")) 
			usage("No --compositionSpace option passed.",true);
		
		// parse the combined set of arguments
		for (Pair<String,String[]> p: argmap) {
			String flag = p.getFirst();
			if (flag.equals("--help"))
				usage("", true);
			else if (flag.equalsIgnoreCase("--verbosity"))
				verbosity = Integer.parseInt(getValues(flag)[0]);
			else if (flag.equalsIgnoreCase("--runTitle"))
				runTitle = getValues(flag)[0];
			else if (flag.equalsIgnoreCase("--outDirName")) {
				outDirName = "";
				String[] values = getValues(flag);
				for (String s : values)
					outDirName += s;
			}
			else if (flag.equalsIgnoreCase("--minInteratomicDistance"))
				minInteratomicDistance = Double.parseDouble(getValues(flag)[0]);
			else if (flag.equalsIgnoreCase("--maxLatticeLength"))
				maxLatticeLength = Double.parseDouble(getValues(flag)[0]);
			else if (flag.equalsIgnoreCase("--minLatticeLength"))
				minLatticeLength = Double.parseDouble(getValues(flag)[0]);
			else if (flag.equalsIgnoreCase("--maxLatticeAngle"))
				maxLatticeAngle = Double.parseDouble(getValues(flag)[0]);
			else if (flag.equalsIgnoreCase("--minLatticeAngle"))
				minLatticeAngle = Double.parseDouble(getValues(flag)[0]);
			else if (flag.equalsIgnoreCase("--maxNumAtoms"))
				maxNumAtoms = Integer.parseInt(getValues(flag)[0]);
			else if (flag.equalsIgnoreCase("--minNumAtoms"))
				minNumAtoms = Integer.parseInt(getValues(flag)[0]);
			else if (flag.equalsIgnoreCase("--doNonnegativityConstraint"))
				doNonnegativityConstraint = Boolean.parseBoolean(getValues(flag)[0]);
			else if (flag.equalsIgnoreCase("--dValue"))
				dValue = Double.parseDouble(getValues(flag)[0]);
			else if (flag.equalsIgnoreCase("--endgameNumGens"))
				endGameNumGens = Integer.parseInt(getValues(flag)[0]);
			else if (flag.equalsIgnoreCase("--useNiggliReducedCell"))
				useNiggliReducedCell = Boolean.parseBoolean(getValues(flag)[0]);
			else if (flag.equalsIgnoreCase("--popSize"))
				popSize = Integer.parseInt(getValues(flag)[0]);
			else if (flag.equalsIgnoreCase("--keepTempFiles"))
				keepTempFiles = Boolean.parseBoolean(getValues(flag)[0]);
			else if (flag.equalsIgnoreCase("--saveStateEachIter"))
				saveStateEachIter = Boolean.parseBoolean(getValues(flag)[0]);
			else if (flag.equalsIgnoreCase("--useRedundancyGuard")) {
				redundancyGuardType = getValues(flag)[0];
				redundancyGuardArgs = GAUtils.subArray(getValues(flag), 1);
			} else if (flag.equalsIgnoreCase("--useSurrogateModel")) {
				useSurrogate = true;
				surrogateArgs = getValues(flag);
			} else if (flag.equalsIgnoreCase("--optimizeDensity")) {
				String[] values = getValues(flag);
				optimizeDensity = true;
				if (optimizeDensity) {
					if (values.length < 2)
						usage("Error: not enough arguments to optimizeDensity", true);
					weightDenAdapt = Double.parseDouble(values[0]);
					numDenAdapt = Integer.parseInt(values[1]);
				}
			}
			else if (flag.equalsIgnoreCase("--notNearestNeighbors")) {
				String[] values = getValues(flag);
				maxNearestNeighborLength = Double.parseDouble(values[0]);
				notNearestNeighbors = GAUtils.parsePairs(GAUtils.subArray(values, 1));
			}
			else if (flag.equalsIgnoreCase("--dryRun"))
				dryRun = Boolean.parseBoolean(getValues(flag)[0]);
			else if (flag.equalsIgnoreCase("--objectiveFunction")) {
				String[] values = getValues(flag);
				if (values.length < 2)
					usage("Not enough parameters given to --objectiveFunction", true);
				objFcnArgs = values;
			}
			else if (flag.equalsIgnoreCase("--parallelize")) {
				numCalcsInParallel = Integer.parseInt(getValues(flag)[0]);
				minPopSize = Integer.parseInt(getValues(flag)[1]);
			}
			else if (flag.equalsIgnoreCase("--selection")) {
				String[] values = getValues(flag);
				if (values[0].equalsIgnoreCase("probDist"))
					sel = new ProbDistSelection(GAUtils.subArray(values, 1));
				else if (values[0].equalsIgnoreCase("off")) 
					sel = null;
				else 
					usage("Unknown selection function " + values[0], true);
			}
			// kind of a hack. (it's because java collections don't support multiple identical keys)
			else if (flag.toLowerCase().startsWith("--variation")) {
				String[] values = getValues(flag, 3);
				String variation = values[2];
				initialVarProbs.add(Double.parseDouble(values[0]));
				endgameVarProbs.add(Double.parseDouble(values[1]));
				if (variation.equalsIgnoreCase("slicer"))
					vars.add(new Slicer(GAUtils.subArray(values, 3)));
				else if (variation.equalsIgnoreCase("structureMut")) 
					vars.add(new StructureMut(GAUtils.subArray(values, 3)));
				else if (variation.equalsIgnoreCase("permutation")) 
					vars.add(new Permutation(GAUtils.subArray(values, 3)));
				else if (variation.equalsIgnoreCase("numStoichsMut")) 
					vars.add(new NumStoichsMut(GAUtils.subArray(values, 3)));
				else if (variation.equalsIgnoreCase("supercell")) 
					vars.add(new SupercellVariation(GAUtils.subArray(values, 3)));
				else 
					usage("Unknown variation function " + variation, true);
			}
			// similar to the variation handling above
			else if (flag.toLowerCase().startsWith("--convergencecriterion")) {
				String[] values = getValues(flag);
				if (values[0].equalsIgnoreCase("maxFunctionEvals"))
					ccs.add(new NumFunctionEvalsCC(GAUtils.subArray(values, 1)));
				else if (values[0].equalsIgnoreCase("maxNumGens"))
					ccs.add(new NumGensCC(GAUtils.subArray(values, 1)));
				else if (values[0].equalsIgnoreCase("maxNumGensWOImpr")) 
					ccs.add(new NumGensWOImprCC(GAUtils.subArray(values, 1)));
				else if (values[0].equalsIgnoreCase("valueAchieved")) 
					ccs.add(new ValueAchievedCC(GAUtils.subArray(values, 1)));
				else if (values[0].equalsIgnoreCase("foundStructure")) 
					ccs.add(new FoundStructureCC(GAUtils.subArray(values, 1)));
				else 
					usage("Unknown convergence criterion " + values[0], true);
			}
			else if (flag.equalsIgnoreCase("--promotion")) {
				String[] values = getValues(flag);
				pro = new Promotion(values);
			}
			else if (flag.toLowerCase().startsWith("--initialpopulation")) {
				String[] values = getValues(flag, 2);
				Integer numOrgs = new Integer(Integer.parseInt(values[0]));
				String creatorType = values[1];

				if (creatorType.equalsIgnoreCase("random")) {
					initialOrgCreators.add(new Pair<StructureOrgCreator,Integer>(new RandomSOCreator(GAUtils.subArray(values, 2)), numOrgs));
//				} else if (creatorType.equalsIgnoreCase("resume")) {
//					initialOrgCreators.put(new ResumeSOCreator(GAUtils.subArray(values, 2)), numOrgs);
				} else if (creatorType.equalsIgnoreCase("manual")) {
					initialOrgCreators.add(new Pair<StructureOrgCreator,Integer>(new ManualSOCreator(GAUtils.subArray(values, 2)), numOrgs));
				} else if (creatorType.equalsIgnoreCase("fromCifs")) {
					initialOrgCreators.add(new Pair<StructureOrgCreator,Integer>(new FromCifsSOCreator(GAUtils.subArray(values, 2)), numOrgs));
				} else if (creatorType.equalsIgnoreCase("units")) {
					initialOrgCreators.add(new Pair<StructureOrgCreator,Integer>(new UnitsSOCreator(GAUtils.subArray(values, 2)), numOrgs));
				} else {
					usage("Unrecognized population type " + creatorType, true);
				}					
			}
			else if (flag.equalsIgnoreCase("--constituents")) {
				System.out.println("WARNING: --constituents flag obsolete. Use --compositionSpace.");
		//		String[] values = getValues(flag);
		//		if (values.length % 2 != 1)
		//			usage("Error: Constituents requires an odd number of arguments", true);
		//		fixStoich = Boolean.parseBoolean(values[0]);
		//		for (int j = 1; j < values.length; j = j + 2) {
		//			constituents.put(values[j], Integer.parseInt(values[j+1]));
		//		}
			} else if (flag.equalsIgnoreCase("--compositionSpace")) {
				List<String> csArgs = new LinkedList<String>();
				for (String s : p.getSecond())
					csArgs.add(s);
				compSpace = new CompositionSpace(csArgs, false);
			}
			// we deal with the input file separately
			else if (!flag.equalsIgnoreCase("--f") && verbosity >= 1)
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
		if (objFcnArgs == null || objFcnArgs.length < 2)
			usage("ERROR: Need an objective function.", true);
		String objFcnType = objFcnArgs[0];
		if (objFcnType.equalsIgnoreCase("pd") && this.getCompSpace().getNumDimensions() < 2)
			usage("ERROR: Can't use pd objFun w/ < 2 dimensions.", true);
	}
	
	public void setSeedGeneration(Cell[] initialPop) {
		seedStructures = initialPop;
	}
	
	// returns a mapping of key-value pairs parsed from the input file
	private List<Pair<String,String[]>>  parseInputFile() {
		int verbosity = GAParameters.getParams().getVerbosity();
		List<Pair<String,String[]>> result = new ArrayList<Pair<String,String[]>>();
		String line = null;
		
		// make sure we have a file to parse
		if (inputFile == null || inputFile.equals("")) 
			usage("Give an input file.", true);
		
		// open the input file
		try {
			BufferedReader in = new BufferedReader(new FileReader(inputFile));
			while ((line = in.readLine()) != null) {
    			StringTokenizer t = new StringTokenizer(line);
    			// skip empty lines:
    			if (!t.hasMoreTokens())
    				continue;
    			String key = t.nextToken();
    			// allow for comment lines starting with #
    			if (key.charAt(0) == '#')
    				continue;
    			// read in all the values associated with the key
    			ArrayList<String> valuesList = new ArrayList<String>();
    			while (t.hasMoreTokens())
    				valuesList.add(t.nextToken());
    			// add the key-value pair to the mapping (add the -- to the front of the key)
    			key = (new String("--")).concat(key);
    			String[] values = new String[valuesList.size()];
    			values = valuesList.toArray(values);
    			result.add(new Pair<String,String[]>(key,values));
    			//result.put(key, values);
    			if (verbosity >= 5) {
    				System.out.print("Parsed from input file: " + key + " ");
    				for (int i = 0; i < values.length; i++)
    					System.out.print(values[i] + " ");
    				System.out.println("");
    			}
			}
		} catch(IOException x) {
			usage("GAParameters.parseInputFile: " + x.getMessage(), true);
		}
		
		return result;	
	}
	
	// returns the first value associated with the given flag.
	// if no such value is available, we print the usage statement and exit.
	private String[] getValues(String flag) {		
		String[] result = null;
		for (Pair<String,String[]> p : argmap) {
			if (p.getFirst().compareToIgnoreCase(flag) == 0) {
				result = p.getSecond();
				break;
			}
		}
		
		if(!isSet(flag) || result == null || result.length == 0)
			usage("Improper usage of " + flag, true);
		
	    return result;
	}
	// overloaded getValues:
	// make sure we have at least n values following the given flag
	private String[] getValues (String flag, int n) {
		String[] answer = getValues(flag);
		if (answer.length < n)
			usage("Improper usage of " + flag + ". Requires " + n + " arguments.", true);
		return answer;
	}
	
	// boolean function to determine whether or not a particular flag
	// was passed on the command line
	private Boolean isSet(String flag)
	{
		if (argmap == null)
			return false;
		
		for (Pair<String,String[]> p : argmap) {
			if (p.getFirst().compareToIgnoreCase(flag) == 0)
				return true;
		}
		return false;
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
	
	private PDBuilder makePDBuilder() {
		Map<Element, Double> cps = new HashMap<Element,Double>();
		for (Element e : compSpace.getElements())
			cps.put(e, 0.0);
			return new PDBuilder(new LinkedList<IComputedEntry>(), compSpace.getElements(), cps );
	}
	
	public ObjectiveFunction getObjectiveFunctionInstance(Organism o) {
		if (objFcnArgs == null || objFcnArgs.length < 1)
			usage("Not enough parameters passed to --ObjectiveFunction", true);
		
		ObjectiveFunction obj = null;
				
		String objFcnType = objFcnArgs[0];
		if (objFcnType.equalsIgnoreCase("epa")) 
			obj = new EnergyPerAtom(GAUtils.subArray(objFcnArgs, 1), o);
		else if (objFcnType.equalsIgnoreCase("pd")) {
			if (pdbuilder == null)
				pdbuilder = makePDBuilder();
			if (getRecord().getGenNum() < 1) {
				obj = new EnergyPerAtom(GAUtils.subArray(objFcnArgs, 1), o);
			} else {
				obj = new PDObjFcn(GAUtils.subArray(objFcnArgs, 1), o, pdbuilder);
			}
		} else
			usage("Unknown objective function " + objFcnType, true);	
		
		return obj;
	}
	
	public boolean doingPDRun() {
		String objFcnType = objFcnArgs[0];
		return objFcnType.equalsIgnoreCase("pd");
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

	public double getMaxLatticeAngle() {
		return maxLatticeAngle;
	}

	public double getMaxLatticeLength() {
		return maxLatticeLength;
	}

	public double getMinInteratomicDistance() {
		return minInteratomicDistance;
	}

	public double getMinLatticeAngle() {
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
	
	public String[] getRedundancyGuardArgs() {
		return redundancyGuardArgs;
	}
	
	public boolean usingSurrogateModel() {
		return useSurrogate;
	}
	
	public String[] getSurrogateArgs() {
		return surrogateArgs;
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
	
	public boolean getUseNiggliReducedCell() {
		return useNiggliReducedCell;
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
		private String energySorted = "index_sorted";
		private File outFile;
		private File tempDir;
		private File outDir;
		private File paramFile;
		private File energyFile;
		
		private int currentGenNum = 0;
		private Generation currentGen = null;
		
		//
		double bestDenEstimate;
		
		// singleton
		private GARecord() {		
			GAParameters params = GAParameters.getParams();
			int verbosity = params.getVerbosity();
			
			String runTitle = params.getRunTitle();
			outDirName = params.getOutDirName();
			tempDirName = params.getTempDirName();

			// make the output directory
			outDir = new File(outDirName);
			if (!params.getDryRun() && !outDir.mkdir()) 
				GAParameters.usage("Can't create directory " + outDirName, true);
			
			// make the output files
			outFile = new File(outDir, outFileName);
			energyFile = new File(outDir, energySorted);
			
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
		
		private String makeCIFPath(StructureOrg s) {
			File f = new File(outDir, Integer.toString(s.getID()) + ".cif");
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
			int verbosity = params.getVerbosity();
			
			// write the generation header (generation x N)
			GAUtils.writeStringToFile("generation " + Integer.toString(currentGenNum) + " " + g.getNumOrganisms() + newline, outFile, true);
			
			// write info on all the members of the population and save their CIFs
			Iterator<Organism> i = g.iterator();
			while (i.hasNext()) {
				// assume here that our organisms are StructureOrgs
				StructureOrg s = (StructureOrg)(i.next());
				// build and write the index info
				StringBuilder info = new StringBuilder();
				info.append(Integer.toString(s.getID()) + " ");
				info.append(Double.toString(s.getValue()) + " ");
				info.append(makeCIFPath(s) + newline);
				GAUtils.writeStringToFile(info.toString(), outFile, true);
				// save the cif
				s.getCell().writeCIF(makeCIFPath(s));
				// save the findsym output
				File outFindSym = new File(makeFindSymPath(s));
				GAUtils.writeStringToFile(Isotropy.getFindsymOutput(s.getCell()), outFindSym, false);
			}
			
			// write energy-sorted list (i.e. lowest to highest)
			GAUtils.writeStringToFile("generation " + Integer.toString(currentGenNum) + " " + g.getNumOrganisms() + newline, energyFile, true);
			List<Integer> energyList = new LinkedList<Integer>();
			for (int k=0; k<g.getNumOrganisms(); k++) {
				double lowest = Double.POSITIVE_INFINITY;
				int index = 0;
				for (Organism r: g) {
					double thisEnergy = r.getValue();
					if (thisEnergy < lowest && !energyList.contains(r.getID())) {
						index = r.getID();
						lowest = thisEnergy;						
					}
				}
				energyList.add(k,index);
			}			
			for (Integer p: energyList) {
				Organism m = g.getOrgByID(p);
				StructureOrg s = (StructureOrg)(m);
				StringBuilder info = new StringBuilder();
				info.append(Integer.toString(s.getID()) + " ");
				info.append(Double.toString(s.getValue()) + " ");
				info.append(makeCIFPath(s) + newline);
				GAUtils.writeStringToFile(info.toString(), energyFile, true);
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
						if (GAParameters.getParams().getVerbosity() >= 1)
							System.out.println("Error deleting temp directory: " + x.getMessage());
					}
			}
		}

	}
}
