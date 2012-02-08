package vasp;

import crystallography.Cell;
import ga.GAOut;
import ga.GAParameters;

import java.io.*;
import java.util.*;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import javax.xml.parsers.ParserConfigurationException;

import org.xml.sax.SAXException;

import utility.*;
import crystallography.Site;
import chemistry.Element;

/* Reads files in a vasp directory into appropriate data structures
 * 
 * Will Tipton
 * 
 * TODO:
 *  - make use of other vasp input datas?
 *  - In the future, this should really just do the parsing and then
 *    create a VaspData object which holds all the parameters, etc. that
 *    go into a vasp run.  Then VaspOut can take a VaspData file and 
 *    write vasp input files to disk.
 */

//NB: Convention: when runs fail, energy should evaluate to Double.POSITIVE_INFINITY

public class VaspOut {
	
	static final String vaspSuccessString = "reached required accuracy";

	public static Cell getPOSCAR(String poscarInFile) {
		String description = null;	/* Will hold first line of the POSCAR */
		
		/**************** Parse me a POSCAR *******************/
		List<Vect> latticeVectors = new LinkedList<Vect>();
		List<Site> basis = new LinkedList<Site>();
		
		StringTokenizer tokens = null;
		try {
			BufferedReader poscarReader = new BufferedReader(new FileReader(poscarInFile));
			/* Read the comment line */
			description = poscarReader.readLine(); 
			/* Get the scaling factor */
			double scalingFactor = Double.parseDouble(poscarReader.readLine()); 
			/* Get lattice vectors */
			for (int i = 0; i < Constants.numDimensions; i++) {
				tokens = new StringTokenizer(poscarReader.readLine());
				Double x = Double.parseDouble(tokens.nextToken()) * scalingFactor;
				Double y = Double.parseDouble(tokens.nextToken()) * scalingFactor;
				Double z = Double.parseDouble(tokens.nextToken()) * scalingFactor;
				latticeVectors.add(new Vect(x,y,z));
			}
			/* Get components */
			tokens = new StringTokenizer(poscarReader.readLine());
			List<Element> elements = new LinkedList<Element>();
			while (tokens.hasMoreTokens()) 
				elements.add(Element.getElemFromSymbol(tokens.nextToken()));
			
			/* Get numbers of components */
			tokens = new StringTokenizer(poscarReader.readLine());
			List<Integer> numsOfComponents = new LinkedList<Integer>();
			
			while (tokens.hasMoreTokens()) 
				numsOfComponents.add(Integer.parseInt(tokens.nextToken()));
			
			/* Skip Selective Dynamics and get Direct or Cartesian */
			String line = poscarReader.readLine();
			if (line.startsWith("S") || line.startsWith("s"))
				line = poscarReader.readLine();
			Boolean usesCartesianCoords = line.startsWith("C") || line.startsWith("c") 
											|| line.startsWith("k") || line.startsWith("K");
			/* Get site positions */
			for (int elemNum = 0; elemNum < numsOfComponents.size(); elemNum++) {
				for (int i = 0; i < numsOfComponents.get(elemNum); i++) {
					tokens = new StringTokenizer(poscarReader.readLine());
					Double x = Double.parseDouble(tokens.nextToken());
					Double y = Double.parseDouble(tokens.nextToken());
					Double z = Double.parseDouble(tokens.nextToken());
					Vect siteLoc = null;
					if (usesCartesianCoords)
						siteLoc = new Vect(x,y,z);
					else
						siteLoc = new Vect(x,y,z,latticeVectors);
					basis.add(new Site(elements.get(elemNum), siteLoc));
				}
			}
		} catch (Exception x) {
			GAOut.out().stdout("Warning: VaspOut.getCell() failed: " + x.getMessage(), GAOut.NOTICE);

			return null;
		} 
		
		/* Make the VaspData object which for now just has a Cell */
		return new Cell(latticeVectors, basis, null);

	}
	
	public static double getFinalEnergy(String outcarFile, boolean cautious) {
		String outcar = Utility.readStringFromFile(outcarFile);
		
		double energy = Double.POSITIVE_INFINITY;
		if (cautious && !vaspRunConverged(outcar))
			return energy;
		
		// get energy: (kinda a hack.)
		try {
			String energyRegexp = "TOTEN *= *[-\\.0-9]* eV";
			Pattern energyPattern = Pattern.compile(energyRegexp);
			Matcher energyMatcher = energyPattern.matcher(outcar);
			StringTokenizer energyLineTok = null;
			while (energyMatcher.find())
				energyLineTok = new StringTokenizer(energyMatcher.group(energyMatcher.groupCount()));
			energyLineTok.nextToken();energyLineTok.nextToken();
			energy = Double.parseDouble(energyLineTok.nextToken());
		} catch (Exception x) {
			GAOut.out().stdout("Warning: VaspOut.getFinalEnergy() failed: " + x.getMessage(), GAOut.NOTICE);
		}
		
		return energy;
	}
	
	public static boolean vaspRunConverged(String outcar) {
		return outcar.contains(vaspSuccessString);
	}

}
