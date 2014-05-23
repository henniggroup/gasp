package ga;

import java.util.ArrayList;
import java.util.List;

import utility.Vect;
import crystallography.Cell;
import crystallography.Site;

public class SubstrateObjFcn extends SurfaceObjFcn {
	
	double substrateDistance;
	Cell substrate;

	public SubstrateObjFcn(List<String> subArray, Organism o) {
		super(subArray, o);
		substrateDistance = Double.parseDouble(objFcnArgs.get(0));
		substrate = GAParameters.getParams().getSubstrate();
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
		double c = oldCell.getCellLengths()[2] + substrate.getCellLengths()[2] + substrateDistance;   
		newBasis.add(new Vect(0.0, 0.0, c));
		
		List<Site> newSites = new ArrayList<Site>();
		for (Site s : substrate.getSites())
			newSites.add(s);
		for (Site s : oldCell.getSites())
			newSites.add(new Site (s.getElement(), s.getCoords().plus(
					new Vect (0.0, 0.0, substrate.getCellLengths()[2] + substrateDistance))));
		
		Cell newCell = new Cell(newBasis, newSites, oldCell.getLabel());
		org.setCell(newCell, false);
		super.padOrg();
	}

	private Cell removeWhiteSpace(Cell oldCell) {
		double [] bounds = super.getAtomBox(oldCell);
		List<Vect> newBasis = new ArrayList<Vect>();
		newBasis.add(oldCell.getLatticeVectors().get(0));
		newBasis.add(oldCell.getLatticeVectors().get(1));
		double c = bounds[1] - bounds[0];  
		newBasis.add(new Vect(0.0, 0.0, c));
		
		List<Site> newSites = new ArrayList<Site>();
		for (Site s : oldCell.getSites())
			newSites.add(new Site (s.getElement(), s.getCoords().plus(
					new Vect (0.0, 0.0, bounds[0]))));
		return new Cell (newBasis, newSites, oldCell.getLabel());
	}

	@Override
	public void unpadOrg() {
		Cell oldCell = org.getCell();
		double [] bounds = super.getAtomBox(oldCell);
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

}
