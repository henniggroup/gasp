/* Genetic algorithm for crystal structure prediction.  Will Tipton.  Ceder Lab at MIT. Summer 2007 */

package ga;

import java.io.File;
import java.util.LinkedList;
import java.util.List;
import java.util.Random;

import Jama.Matrix;

import utility.Constants;
import utility.RandomNumbers;
import utility.Vect;

import crystallography.Cell;
import crystallography.Site;

// StructureMut is a Variation operation which acts on StructureOrg organisms.
// It modifies the lattice of the StructureOrg's Structure by multiplication
// by a strain matrix (I+E_{ij}) where the E_{ij} are zero-mean, normally-distributed
// random variables with standard deviation sigmaLattice except that they are forced 
// to be between $-1$ and $1$.  Each of the atoms is
// perturbed with probability mutRate by addition to each of its real-space coords
// of a zero-mean, normally-distributed random variable with standard dev. sigmaAtoms.

public final class StructureMut implements Variation {
	private double mutRate;
	private double sigmaAtoms;
	private double sigmaLattice;
	
	public StructureMut(String[] args) {
		if (args == null || args.length < 3)
			GAParameters.usage("Not enough parameters given to StructureMut", true);
		
		mutRate = Double.parseDouble(args[0]);
		sigmaAtoms = Double.parseDouble(args[1]);
		sigmaLattice = Double.parseDouble(args[2]);
	}
	
	public String toString() {
		StringBuilder result = new StringBuilder();
		
		result.append("StructureMut variation");
		result.append(" (mutRate = " + mutRate + ", sigmaAtoms = " + sigmaAtoms
				+ ", sigmaLattice = " + sigmaLattice + ")");
		
		return result.toString();
	}
	
	public Organism doVariation(Generation parents, Generation offspring, Selection sel) {
		GAParameters params = GAParameters.getParams();
		Random rand = params.getRandom();
		int verbosity = params.getVerbosity();
		
		// pick a parent randomly 
		StructureOrg p = (StructureOrg)(sel.doSelection(parents, 1)[0]);
		Cell pStruct = p.getCell();
		
		// copy the parent's data into new structures
		List<Site> pSites = pStruct.getSites();
		List<Site> newSites = new LinkedList<Site>();

		// perturbing the lattice parameters
		double[][] strain = new double[3][3];
		// initialize the strain matrix to the identity
		// then add some zero mean Gaussian random variables
		for (int i = 0; i < strain.length; i++)
			for (int j = 0; j < strain.length; j++) {
				strain[i][j] = rand.nextGaussian() * sigmaLattice;
				// allow strain matrix perturbations to be only between -1 and 1
				if (strain[i][j] > 1)
					strain[i][j] = 1;
				if (strain[i][j] < -1)
					strain[i][j] = -1;
				// add the identity
				if (i == j)
					strain[i][j] += 1;		
			}

		// modify the basis to make newBasis
		double[][] newVectsD = new double[Constants.numDimensions][Constants.numDimensions];
		for (int i = 0; i < Constants.numDimensions; i++)
			for (int j = 0; j < Constants.numDimensions; j++)
				newVectsD[i][j] = pStruct.getLatticeVectors().get(i).getCartesianComponents().get(j);
		//newVects = MatrixMath.SquareMatrixMult(strain, newVects);
		newVectsD = (new Matrix(strain)).times(new Matrix(newVectsD)).getArray();
		
		List<Vect> newVects = new LinkedList<Vect>();
		for (int i = 0; i < Constants.numDimensions; i++)
				newVects.add(new Vect(newVectsD[i]));
		
		// perturb atomic positions
		for (int i = 0; i < pSites.size(); i++) {
			// perturb each site's location with probability mutRate
			List<Double> fracCoords = pSites.get(i).getCoords().getComponentsWRTBasis(pStruct.getLatticeVectors());
			if (rand.nextDouble() < mutRate) {
				// perturb by a Gaussian of stddev mutRadius along each axis
				for (int k = 0; k < Constants.numDimensions; k++) 
					fracCoords.set(k, GAUtils.renormZeroOne(fracCoords.get(k) + rand.nextGaussian() * sigmaAtoms));
			} 
			newSites.add(new Site(pSites.get(i).getElement(), new Vect(fracCoords, newVects)));
		}
		
				
		// make the new offspring
		StructureOrg result = new StructureOrg(new Cell(newVects, newSites));
		
		if (verbosity >= 5) {
			System.out.println("StructureMut created new StructureOrg:");
			System.out.println(result);
		}			
		
		return result;
	}	
	
	// just for testing
	public static void main(String[] args) {
		StructureOrg s1 = new StructureOrg(Cell.parseCell("/home/wtipton/cifs/143.cif", "cif"));
		s1.setFitness(-1);
		
		String[] smArgs = {"0.5", "0.5", "0.5"};
		StructureMut s = new StructureMut(smArgs);
		String[] selArgs = {"1", "0"};
		Selection sel = new ProbDistSelection(selArgs);

		Generation parents = new Structures();
		parents.addOrganism(s1);
		
		System.out.println(s1);
		StructureOrg o = (StructureOrg)s.doVariation(parents, null, sel);
	//	GAUtils.writeStringToFile(o.getCIF(), new File("offspring.cif"), false);
		System.out.println(o);
	}

}
