/* Copyright 2011-2014 Will Tipton, Richard Hennig, Ben Revard, Stewart Wenner

This file is part of the Genetic Algorithm for Structure and Phase Prediction (GASP).

    GASP is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    GASP is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with GASP.  If not, see <http://www.gnu.org/licenses/>.
    
    
    */

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
