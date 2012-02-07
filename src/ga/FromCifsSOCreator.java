/* Genetic algorithm for crystal structure prediction.  Will Tipton.  Ceder Lab at MIT. Summer 2007 */

package ga;

import java.io.*;
import java.util.List;

import crystallography.Cell;

// FromCifsSOCreator implements StructureOrgCreator.  It is given a directory, and it
// creates StructureOrgs from all Cif files in the directory.

public class FromCifsSOCreator implements StructureOrgCreator {
	
	String dirName;
	
	File[] cifFiles;

	public FromCifsSOCreator(List<String> args) {
		if (args.size() < 1)
			GAParameters.usage("Not enough arguments passed to FromCifsSOCreator.", true);
		dirName = new String(args.get(0));		
	}
	
	public String toString() {
		return "FromCifsSOCreator: directory = " + dirName;
	}
	
	public StructureOrg makeOrganism(Generation g) {
		GAParameters params = GAParameters.getParams();
		int verbosity = params.getVerbosity();
		
		// get a list of the CIF files in the directory
		if (cifFiles == null) {
			File dir = new File(dirName);
			cifFiles = dir.listFiles(new CifFileFilter());
		}
		
		if ((cifFiles == null || cifFiles.length == 0) && verbosity >= 2)
			System.out.println("Warning: no CIF files found in directory " + dirName);
		
		// make and return an organism
		for (int i = 0; i < cifFiles.length; i++)
			if (cifFiles[i] != null) {
				if (verbosity >= 3)
					System.out.println("Making new StructureOrg from " + cifFiles[i].getPath());
				StructureOrg result = new StructureOrg(Cell.parseCif(cifFiles[i]));
				cifFiles[i] = null;
				return result;
			}
		
		// if we get this far we've used all the Cif's in the directory
		return null;
	}
	
	private class CifFileFilter implements FilenameFilter {
		public boolean accept(File f, String name) {
			return name.endsWith(".cif");
		}
	}
}
