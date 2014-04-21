package pdvisual;

import chemistry.*;
import crystallography.*;
import java.util.List;

public class PseudoPDAnalyzer implements IPDAnalyzer {

	PseudoPDData ppdd;
	
	public PseudoPDAnalyzer(PseudoPDData _ppdd) {
		ppdd = _ppdd;
	}
	
	public double getEnergyPerAtomAboveHull(int i) {
		return getEnergyPerAtomAboveHull(ppdd.getAllPhases().get(i));
	}
	
	public double getEnergyPerAtomAboveHull(MixedPhase m) {
	//	return m.getEnergy() - getEnergyPerAtomOnHull(m.getComposition());
		// TODO
		return 0;
	}
	
	public double getRelativeChemicalPotential(Element element, List<Integer> f) {
		// TODO
		return 0;
	}

    public double getChemicalPotential(Element element, List<Integer> f) {
		// TODO
		return 0;
	}
}
