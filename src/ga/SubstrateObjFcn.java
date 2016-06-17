package ga;

import java.util.ArrayList;
import java.util.List;

import utility.Vect;
import vasp.VaspIn;
import vasp.VaspOut;
import crystallography.Cell;
import crystallography.Site;

public class SubstrateObjFcn extends SurfaceObjFcn {
	
	double substrateDistance;
	Cell substrate;

	public SubstrateObjFcn(List<String> subArray, Organism o) {
		super(subArray, o);
		substrateDistance = Double.parseDouble(objFcnArgs.get(0));
		substrate = VaspOut.getPOSCAR("/home/ay256/Desktop/2.vasp");
		//substrate = GAParameters.getParams().getSubstrate();
		substrate.rotatedIntoPrincDirs();
		substrate = removeWhiteSpace(substrate);
	}

	@Override
	public void padOrg() {
		Cell oldCell = org.getCell();
		
		oldCell.rotatedIntoPrincDirs();
		oldCell = removeWhiteSpace(oldCell);
		
		List<Vect> newBasis = new ArrayList<Vect>();
		newBasis.add(oldCell.getLatticeVectors().get(0));
		newBasis.add(oldCell.getLatticeVectors().get(1));
		double c = oldCell.getCellLengths()[2] + substrate.getCellLengths()[2];   
		newBasis.add(new Vect(0.0, 0.0, c));
		
		List<Site> newSites = new ArrayList<Site>();
		for (Site s : substrate.getSites())
			newSites.add(s);
		for (Site s : oldCell.getSites())
			newSites.add(new Site (s.getElement(), s.getCoords().subtract(
					new Vect (0.0, 0.0, substrate.getCellLengths()[2] + substrateDistance * 2))));
		
		Cell newCell = new Cell(newBasis, newSites, oldCell.getLabel());
		org.setCell(newCell, false);
		//super.padOrg();
	}

	public Cell removeWhiteSpace(Cell oldCell) {
		
		double [] bounds = oldCell.getAtomBoxZ();
		double zlen = bounds[1] - bounds[0];
		
		// make a new box
		List<Vect> newBasis = new ArrayList<Vect>();
		newBasis.add(oldCell.getLatticeVectors().get(0));
		newBasis.add(oldCell.getLatticeVectors().get(1));
		double c = bounds[1] - bounds[0] + substrateDistance;  
		newBasis.add(new Vect(0.0, 0.0, c));
		
		// get the new sites
		List<Site> newSites = new ArrayList<Site>();
		for (Site s : oldCell.getSites())
			newSites.add(new Site (s.getElement(), s.getCoords().subtract(
					new Vect (0.0, 0.0, bounds[0] - substrateDistance / 2))));
		return new Cell (newBasis, newSites, oldCell.getLabel());
	}

	@Override
	public void unpadOrg() {
		Cell oldCell = org.getCell();
		double [] bounds = oldCell.getAtomBoxZ();
		List<Site> newSites = new ArrayList<Site>();
		for (Site s : oldCell.getSites()) {
			if (s.getCoords().getCartesianComponents().get(2) > bounds[0] + substrate.getCellLengths()[2] + substrateDistance / 2)
				newSites.add(s);
		}
		Cell newCell = new Cell(oldCell.getLatticeVectors(), newSites, oldCell.getLabel());
		org.setCell(newCell, false);
	}

	@Override
	public String toString() {
		return "SubstrateObjFcn";
	}
	
	public static void main(String[] args) {
		Cell oldCell = VaspOut.getPOSCAR("/home/ay256/Desktop/23.vasp");
		Cell substrate = VaspOut.getPOSCAR("/home/ay256/Desktop/2.vasp");
		List<Vect> newBasis = new ArrayList<Vect>();
		newBasis.add(substrate.getLatticeVectors().get(0));
		newBasis.add(substrate.getLatticeVectors().get(1));
		newBasis.add(oldCell.getLatticeVectors().get(2));
		List<Site> newSites = new ArrayList<Site>();
		for (Site site : oldCell.getSites()) {
			Vect v = site.getCoords().changeBasis(newBasis);
			newSites.add(new Site(site.getElement(), v));
		}
		Cell newCell = new Cell(newBasis, newSites, oldCell.getLabel());
		
		List<String> arguments = new ArrayList<String>();
		arguments.add("10.0");
		arguments.add("2.0");
		Organism o = new StructureOrg(newCell);
		SubstrateObjFcn fcn = new SubstrateObjFcn(arguments, o);
		Cell whiteCell = fcn.removeWhiteSpace(oldCell);
		VaspIn.writePoscar(whiteCell, "/home/ay256/Desktop/whitespace.vasp", true);
		//fcn.padOrg();
		//VaspIn.writePoscar(fcn.org.getCell(), "/home/ay256/Desktop/padOrg.vasp", true);
		
	}

}
