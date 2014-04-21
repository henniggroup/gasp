package castep;

import crystallography.*;
import java.text.NumberFormat;
import java.util.*;
import utility.*;

import chemistry.Element;
// writes internal data structures into CASTEP input files

public class CastepIn {
		
	private Cell cell;
	private List<Element> elements; /* Holds the elements involved in this
	 								 * CASTEP calculation in the order they come
	 								 * in the POTCAR file. */
	
	double kpointSpacing;
	double pressure;
	String paramFile;
	Map<Element,String> ppFileMap;

	public CastepIn(Cell _cell, double _kpointSpacing, double _pressure, String _paramFile, Map<Element,String> _ppFileMap) {
		cell = _cell;
		elements = cell.getComposition().getElements();
		kpointSpacing = _kpointSpacing;
		pressure = _pressure;
		ppFileMap = _ppFileMap;
		paramFile = _paramFile;
	}
	
	public Cell getCell() {
		return cell;
	}
	
	public List<Element> getElements() {
		List<Element> result = new LinkedList<Element>();
		result.addAll(elements);
		return result;
	}
	
	public String makeCellFileString() {
		/* Note that this will overwrite outPoscar. If this is not desired behavior,
		 * the calling routine should check if the output file already exists */
		
		StringBuilder result = new StringBuilder();
		
		result.append("%BLOCK LATTICE_CART\n");
	
		NumberFormat nf = NumberFormat.getInstance();
		nf.setMaximumFractionDigits(8);
		nf.setMinimumFractionDigits(8);
		
		List<Vect> latticeVectors = cell.getLatticeVectors();
		for (Vect vec : latticeVectors) {
			List<Double> comps = vec.getCartesianComponents();
			for (int i = 0; i < Constants.numDimensions; i++) {
				String coordStr = nf.format(comps.get(i));
				result.append(coordStr + " ");
			}
			result.append("\n");
		}
		
		result.append("%ENDBLOCK LATTICE_CART\n\n");
		
		result.append("%BLOCK POSITIONS_FRAC\n");
		for (Site s : cell.getSites()) {
			result.append(s.getElement().getSymbol() + " ");
			for (int i = 0; i < Constants.numDimensions; i++)
				result.append(nf.format(s.getCoords().getComponentsWRTBasis(cell.getLatticeVectors()).get(i)) + " ");
			result.append("\n");
		}
		
		result.append("%ENDBLOCK POSITIONS_FRAC\n\n");
		
		result.append("KPOINTS_MP_SPACING " + kpointSpacing + "\n\n");
		
		result.append("SNAP_TO_SYMMETRY\n\n");
		
		result.append("%BLOCK SPECIES_POT\n");
		for (Element e : cell.getComposition().getElements())
			result.append(e.getSymbol() + " " + ppFileMap.get(e) + "\n");
		result.append("%ENDBLOCK SPECIES_POT\n\n");
		
		result.append("%BLOCK EXTERNAL_PRESSURE\n");
		result.append("GPa\n");
		result.append(pressure + " 0.0 0.0\n");
		result.append(pressure + " 0.0 \n");
		result.append(pressure + " \n");
		result.append("%ENDBLOCK EXTERNAL_PRESSURE\n");



		return result.toString();
	}

	public void makeParam(String directory) {
		String paramStr = Utility.readStringFromFile(paramFile);
		Utility.writeStringToFile(paramStr, directory + "/"+CastepEnergy.castepPrefix+".param");
		
	}

	public void makeCell(String directory) {
		Utility.writeStringToFile(makeCellFileString(), directory + "/"+CastepEnergy.castepPrefix+".cell");
	}

}
