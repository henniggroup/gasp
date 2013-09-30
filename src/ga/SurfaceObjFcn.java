package ga;

import java.util.ArrayList;
import java.util.List;

import crystallography.Cell;
import crystallography.Site;
import pdvisual.PDAnalyzer;
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
	
	private void padOrg() {
		Cell oldCell = org.getCell();
		
		List<Vect> basis = new ArrayList<Vect>();
		basis.add(oldCell.getLatticeVectors().get(0));
		basis.add(oldCell.getLatticeVectors().get(1));
		basis.add(new Vect(0.0, 0.0, padding));
		
		List<Site> newSites = new ArrayList<Site>();
		for (Site s : oldCell.getSites())
			newSites.add(new Site(s.getElement(), s.getCoords().plus(new Vect(0.0, 0.0, padding/2))));	
		
		Cell newCell = new Cell(basis, newSites, oldCell.getLabel());
		org.setCell(newCell, false);
	}
	
	private void unpadOrg() {
		Cell oldCell = org.getCell();
		if (oldCell == null)
			return;
		
		oldCell = oldCell.getCellRotatedIntoPrincDirs();
		
		double mid = GAParameters.getParams().getMinInteratomicDistance();

		// get a bounding box for the atoms
		double minz = Double.MAX_VALUE;
		double maxz = Double.MIN_VALUE;
		for (Site s : oldCell.getSites()) {
			List<Double> cartComps = s.getCoords().getCartesianComponents();
			minz = Math.min(minz, cartComps.get(2));
			maxz = Math.max(maxz, cartComps.get(2));
		}
		double zlen = maxz - minz;
		
		// make new list of sites where we subtract (minz) off all the old ones
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
	/*
	public static void main(String args[]) {
		
		
		String arg[] = {"20", "hi"};
		
		StructureOrg c = new StructureOrg(VaspOut.getPOSCAR("/home/wtipton/POSCAR"));
		
		SurfaceObjFcn cof = new SurfaceObjFcn(arg, c);
		
		cof.padOrg();
		
		c.getCell().writeCIF("/home/wtipton/POSCAR.padded.cif");
		
		cof.unpadOrg();
		
		c.getCell().writeCIF("/home/wtipton/POSCAR.unpadded.cif");

	}*/
}
