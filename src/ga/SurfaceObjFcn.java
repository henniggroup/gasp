package ga;

import java.util.ArrayList;

import java.util.List;

import java.util.Arrays;   // Ben added this

import chemistry.Element;
import crystallography.Cell;
import crystallography.Site;
import pdvisual.PDAnalyzer;
import utility.Triplet;
import utility.Utility;
import utility.Vect;
import vasp.VaspOut;

public class SurfaceObjFcn extends ObjectiveFunction {
	
	List<String> objFcnArgs;
	StructureOrg org;
	double padding;
	
	public SurfaceObjFcn(List<String> subArray, Organism o) {
		padding = Double.parseDouble(subArray.get(0));
		objFcnArgs = Utility.subList(subArray, 1);
		org = (StructureOrg)o;
	}
	
	/**
	 * Rotates the cell into the principle directions and then adds padding to the magnitude 
	 * of the c lattice vector (which is now aligned vertically). Takes care of atomic sites 
	 * so that the relative positions of the atoms are the same in the padded cell as the original.
	 */
	private void padOrg() {
		Cell oldCell = org.getCell();
		
		oldCell.rotatedIntoPrincDirs(); 
		
		List<Vect> basis = new ArrayList<Vect>();
		basis.add(oldCell.getLatticeVectors().get(0));
		basis.add(oldCell.getLatticeVectors().get(1));
		double c = oldCell.getCellLengths()[2];   
		basis.add(new Vect(0.0, 0.0, padding + c));
		
		List<Site> newSites = new ArrayList<Site>();
		for (Site s : oldCell.getSites())
			newSites.add(new Site(s.getElement(), s.getCoords().plus(new Vect(0.0, 0.0, (padding + c)/2))));	
		
		Cell newCell = new Cell(basis, newSites, oldCell.getLabel());
		org.setCell(newCell, false);
	}
	
	/**
	 * Rotates the cell into the principle directions and then removes vertical padding 
	 * from the cell without changing the relative positions of the atoms. The magnitude 
	 * of the vertical cell vector of the unpadded cell is equal to the vertical distance 
	 * between highest and lowest atoms in the cell plus the minimum interatomic distance 
	 * (mid). The atoms are arranged in the unpadded cell such that highest atom is mid/2 
	 * from the top of cell, and the lowest atom is mid/2 from the bottom of the cell.   
	 */
	private void unpadOrg() {
		Cell oldCell = org.getCell();
		if (oldCell == null)
			return;
		
		oldCell.rotatedIntoPrincDirs();
		double mid; // minimum interatomic distance
		
		// get the largest minimum interatomic distance
		if (GAParameters.getParams().getMinInteratomicDistance() == -1) {
			mid =  maxMID();	
		} else {
			mid = GAParameters.getParams().getMinInteratomicDistance();
		}

		// get a bounding box for the atoms
		double minz = Double.MAX_VALUE;
		double maxz = Double.MIN_VALUE;
		for (Site s : oldCell.getSites()) {
			List<Double> cartComps = s.getCoords().getCartesianComponents();
			minz = Math.min(minz, cartComps.get(2));
			maxz = Math.max(maxz, cartComps.get(2));
		}
		double zlen = maxz - minz;
		
		// make new list of sites where we subtract (0,0,minz-mid/2) off all the old ones
		List<Site> newSites = new ArrayList<Site>();
		Vect minv = new Vect(0.0, 0.0, minz - mid/2);
		for (Site s : oldCell.getSites())
			newSites.add(new Site(s.getElement(), s.getCoords().subtract(minv)));
		
		// make a new box which is xlen x ylen x zlen plus the min interatomic distance
		List<Vect> newBasis = new ArrayList<Vect>();
		newBasis.add(oldCell.getLatticeVectors().get(0));
		newBasis.add(oldCell.getLatticeVectors().get(1));
		newBasis.add(new Vect(0.0, 0.0, zlen + mid));
		
		Cell newCell = new Cell(newBasis, newSites, oldCell.getLabel());
		org.setCell(newCell, false);
		
	}
	
	
	/**
	 * Returns the largest per species minimum interatomic distance
	 * Precondition: the perSpeciesMID option has been used
	 */
	public double maxMID() {
		// Get an array of all the minimum interatomic distances
		List<Triplet<Element, Element, Double>> tripletList = GAParameters.getParams().getPerSpeciesMIDs();
		double[] MIDs = new double[tripletList.size()];
		for (int i = 0; i < tripletList.size(); i++) {
			MIDs[i] = tripletList.get(i).getThird();
		}
		
		// Sort through the array and find the max value
		double mid = 0;
		for (int j = 0; j < MIDs.length; j++) {
			if (MIDs[j] > mid) {
				mid = MIDs[j];
			}
		}
		return mid;
	}
	
	
	public void run() {
		
		padOrg();
		
		ObjectiveFunction energyFcn = ObjFcnFactory.getObjectiveFunctionInstance(org, objFcnArgs);

		// relax the cell - have to wait for it to finish before using results
		Thread t = energyFcn.evaluate();
		try {
			if (t != null)
				t.join();
		} catch (InterruptedException x) {
			GAOut.out().stdout("InterruptedException in energy calc thread in SurfaceObjFcn: " + x.getMessage(), GAOut.WARNING, org.getID());
		} catch (Exception x) {
			GAOut.out().stdout("ERROR: exception in SurfaceObjFcn:run", GAOut.CRITICAL);
			x.printStackTrace();
		}
		
		// updating structure w/ relaxed version is this done in EPA already
		//org.setValue((new PDAnalyzer(pdbuilder.getPDData())).getEnergyPerAtomAboveHull(org));
		
		unpadOrg();
	}

	// dont update numCalculations when implementing this.. the underlying objfcn will do it
	public Thread evaluate() {		
		// short circuit here if we've done the calculation already
		if (org.knowsValue())
			return null;
		
		// start the calculation and return the Thread
		Thread t = new Thread(this);
		t.start();
		return t;
	}
	
	// an ObjectiveFunction should also overload toString();
	public String toString() {
		return "SurfaceObjFcn"; // TODO: more? 
	}
	
	
	// for testing
	
/*	public static void main(String args[]) {
		
		
		String arg[] = {"20", "hi"};
		
		StructureOrg c = new StructureOrg(VaspOut.getPOSCAR("/Users/benjaminrevard/GA/padding_testing/POSCAR"));
		
		// Convert array of strings to list of strings for the constructor
		List<String> larg = Arrays.asList(arg);
		
		SurfaceObjFcn cof = new SurfaceObjFcn(larg, c);
		
		cof.padOrg();
		
		c.getCell().writeCIF("/Users/benjaminrevard/GA/padding_testing/POSCAR.padded.cif");
		
		cof.unpadOrg();
		
		c.getCell().writeCIF("/Users/benjaminrevard/GA/padding_testing/POSCAR.unpadded.cif");
		
		System.out.println(GAParameters.getParams().getMinInteratomicDistance());

	}*/
}
