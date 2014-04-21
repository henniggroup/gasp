package pdvisual;

import chemistry.Element;
import java.util.*;

public interface IPDAnalyzer {

	
	public double getEnergyPerAtomAboveHull(int i);
//	public double getRelativeChemicalPotential(Element element, List<Integer> facet);
    public double getChemicalPotential(Element element, List<Integer> facet);
	
}
