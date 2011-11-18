package vasp;

import crystallography.*;
import java.io.*;
import java.text.NumberFormat;
import java.util.*;
import utility.*;

import chemistry.Element;
// writes internal data structures into VASP input files

public class VaspIn {
		
	private Cell cell;
	private List<Element> elements; /* Holds the elements involved in this
	 								 * VASP calculation in the order they come
	 								 * in the POTCAR file. */
	String description;	/* Holds the comment in the first line of the POSCAR file. */
	
	String kpointsFile;
	String incarFile;
	Map<Element,String> potcarFileMap;
	
	public VaspIn(Cell _cell,  String kpointsFile, String incarFile, Map<Element,String> potcarFileMap) {
		this(_cell, "Created by VaspData", kpointsFile, incarFile, potcarFileMap);
	}
	
	public VaspIn(Cell _cell, String _description, String _kpointsFile, String _incarFile, Map<Element,String> _potcarFileMap) {
		cell = _cell;
		elements = cell.getComposition().getElements();
		description = _description;
		kpointsFile = _kpointsFile;
		potcarFileMap = _potcarFileMap;
		incarFile = _incarFile;
	}
	
	public Cell getCell() {
		return cell;
	}
	
	public List<Element> getElements() {
		List<Element> result = new LinkedList<Element>();
		result.addAll(elements);
		return result;
	}
	
	public String getDescription()	{
		return description;
	}
		
	public void writePoscar(String outPoscar, Boolean useCartesianCoords) {
		/* Note that this will overwrite outPoscar. If this is not desired behavior,
		 * the calling routine should check if the output file already exists */
		
		/* Write the output */
		try {
			NumberFormat nf = NumberFormat.getInstance();
			nf.setMaximumFractionDigits(8);
			nf.setMinimumFractionDigits(8);
			
			BufferedWriter writer = new BufferedWriter(new FileWriter(outPoscar));
			writer.write(description + "\n");
			writer.write("1.0 \n");
			List<Vect> latticeVectors = cell.getLatticeVectors();
			for (Vect vec : latticeVectors) {
				List<Double> comps = vec.getCartesianComponents();
				for (int i = 0; i < Constants.numDimensions; i++) {
					String coordStr = nf.format(comps.get(i));
					writer.write(coordStr + " ");
				}
				writer.write("\n");
			}
			/* Write elements */
			for (Element e : elements)
				writer.write(e.getSymbol() + " ");
			writer.write("\n");
			/* Write number of each type of site */
			for (Element e : elements)
				writer.write(cell.getNumSitesWithElement(e) + " ");
			writer.write("\n");
			if (useCartesianCoords)
				writer.write("Cartesian\n");
			else
				writer.write("Direct\n");
			/* Make sure w're printing these out in the right order */
			List<Site> basis = cell.getSites();	
			for (Element e : elements)
				for (Site s : basis) {
					if (s.getElement() == e) {
						List<Double> coords = null;
						if (useCartesianCoords)
							coords = s.getCoords().getCartesianComponents();
						else
							coords = s.getCoords().getComponentsWRTBasis(cell.getLatticeVectors());
						for (int i = 0; i < Constants.numDimensions; i++) {
							String coordStr = nf.format(coords.get(i));
							writer.write(coordStr + " ");
						}
						writer.write("\n");
					}
				}
			
			/* Write it all to disk */
			writer.flush();
		} catch (IOException x) {
			System.out.println("IOException in VaspOut.writePoscar(): " + x.getMessage());
		}
	}

	public void makeINCAR(String directory) {
		String incarStr = Utility.readStringFromFile(incarFile);
		Utility.writeStringToFile(incarStr, directory + "/INCAR");
		
	}

	public void makePOSCAR(String directory) {
		writePoscar(directory + "/POSCAR", false);
	}

	public void makeKPOINTS(String directory) {
		String kpointsStr = Utility.readStringFromFile(kpointsFile);
		Utility.writeStringToFile(kpointsStr, directory + "/KPOINTS");
	}

	public void makePOTCAR(String directory) {
		StringBuilder potcar = new StringBuilder();
		
		for (Element e : elements) {
//	System.out.println("potcar for " + e.getSymbol() + " : " + potcarFileMap.get(e));
			String potcarStr = Utility.readStringFromFile(potcarFileMap.get(e));
			potcar.append(potcarStr);
		}
		
		Utility.writeStringToFile(potcar.toString(), directory + "/POTCAR");
	}

}
