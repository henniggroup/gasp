package pdvisual;

import crystallography.*;
import chemistry.*;

public interface IComputedEntry {

	public double getEnergyPerAtom();
	
	public Composition getComposition();
	
	public Cell getCell();
	
	public double getTotalEnergy();
	
	public String getLabel();
	
}
