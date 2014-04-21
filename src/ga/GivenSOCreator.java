/* Genetic algorithm for crystal structure prediction.  Will Tipton.  Ceder Lab at MIT. Summer 2007 */

package ga;

import java.io.*;
import java.util.List;

import vasp.VaspOut;

import crystallography.Cell;

// FromCifsSOCreator implements StructureOrgCreator.  It is given a directory, and it
// creates StructureOrgs from all POSCAR files in the directory.
// These files' filenames must end in 'POSCAR'.

public class GivenSOCreator implements StructureOrgCreator {
	
	String dirName;
	
	File[] poscarFiles;

	public GivenSOCreator(List<String> args) {
		if (args.size() < 1)
			GAParameters.usage("Not enough arguments passed to GivenSOCreator.", true);
		dirName = new String(args.get(0));		
	}
	
	public String toString() {
		return "FromCifsSOCreator: directory = " + dirName;
	}
	
	public StructureOrg makeOrganism(Generation g) {
		GAParameters params = GAParameters.getParams();
		
		// get a list of the CIF files in the directory
		if (poscarFiles == null) {
			File dir = new File(dirName);
			poscarFiles = dir.listFiles(new PoscarFileFilter());
			if (poscarFiles == null || poscarFiles.length == 0) {
				GAOut.out().stdout("WARNING: No POSCAR files found in " + dirName, GAOut.WARNING);
				GAOut.out().stdout("POSCAR files must end in 'POSCAR'.", GAOut.WARNING);

			}
		}
		
		// make and return an organism
		for (int i = 0; i < poscarFiles.length; i++)
			if (poscarFiles[i] != null) {
				StructureOrg result = new StructureOrg(null);
				GAOut.out().stdout("Making new StructureOrg from " + poscarFiles[i].getPath(), GAOut.NOTICE, result.getID());
				result.setCell(VaspOut.getPOSCAR(poscarFiles[i].getPath()));
				poscarFiles[i] = null;
				return result;
			}
		
		// if we get this far we've used all the poscars in the directory
		return null;
	}
	
	private class PoscarFileFilter implements FilenameFilter {
		public boolean accept(File f, String name) {
			return name.endsWith("POSCAR");
		}
	}
}
