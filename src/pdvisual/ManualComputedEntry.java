package pdvisual;

import ga.GAUtils;

import java.io.Serializable;

import chemistry.Composition;
import crystallography.Cell;

public class ManualComputedEntry implements IComputedEntry,Serializable {
	static final long serialVersionUID = 1;
	
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
	
	public String toString() {
		StringBuilder result = new StringBuilder();
		String newline = GAUtils.newline();
		
		result.append("ManualComputedEntry: " + getLabel() + ":" + newline);
		result.append("Total Energy: " + totalEnergy +  newline);
		result.append(cell.toString());
		
		return result.toString();
	}
	
}
