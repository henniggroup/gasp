package pdvisual;

import chemistry.Composition;
import crystallography.Cell;

public class ManualComputedEntry implements IComputedEntry {
	
	double totalEnergy;
	Cell cell;
	
	public ManualComputedEntry(Cell c, double e) {
		cell = c;
		totalEnergy = e;
	}

	public double getEnergyPerAtom() {
		return totalEnergy / cell.getNumSites();
	}
	
	public Composition getComposition() {
		return cell.getComposition();
	}
	
	public Cell getCell() {
		return cell;
	}
	
	public double getTotalEnergy() {
		return totalEnergy;
	}
	
	public String getLabel() {
		return cell.getLabel();
	}
	
}
