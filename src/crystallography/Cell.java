package crystallography;

import ga.GAOut;
import ga.GAParameters;
import ga.GAUtils;
import ga.Generation;
import ga.RandomSOCreator;
import ga.StructureOrg;
import ga.Structures;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.Serializable;
import java.util.*;

/*import org.openbabel.OBAtom;
import org.openbabel.OBBase;
import org.openbabel.OBConversion;
import org.openbabel.OBMol;
import org.openbabel.OBUnitCell;
import org.openbabel.openbabel_java;
import org.openbabel.openbabel_javaConstants;
import org.openbabel.vector3;
import org.openbabel.vectorVector3;*/

import Jama.Matrix;

import utility.*;
import vasp.VaspOut;
import chemistry.*;

//import org.openbabel.*;
/* Represents a crystallographic cell: lattice parameters and a basis.
 * IS IMMUTABLE!
 * 
 * Will Tipton
 * 
 * TODO:
 * 	- accessors.  make sure to pass deep copies of our internal objects
 */

public class Cell implements Serializable {
	
	static final long serialVersionUID = 1l;
	
	private List<Vect> latticeVectors;
	private List<Site> basis;
	private String label;
	
/*	static {
		System.loadLibrary("openbabel_java");
	}
*/	
	public Cell (Double _a, Double _b, Double _c, Double _alpha, Double _beta, Double _gamma, List<Site> _basis, String _label){
		this(getVectorsfromLParamsDegrees(_a, _b, _c, _alpha, _beta, _gamma), _basis, _label);
	}
	
	public Cell (List<Vect> _vectors, List<Site> _basis) {
		this(_vectors,_basis,null);
	}
	
	public Cell (List<Vect> _vectors, List<Site> _basis, String _label) {
		/* Initialize and do a copy of the vectors list.
		 * Shallow copy is ok since Vects are immutable. */
		/* Make sure we're given 3 vectors as we're necessarily in 3 dimensions	*/
		if (_vectors.size() != Constants.numDimensions)
			throw new IllegalArgumentException("Cell constructor given other than " + Constants.numDimensions + " vectors.");
		latticeVectors = new LinkedList<Vect>();
		latticeVectors.addAll(_vectors);
		
		/* do a deep copy of the basis list */
		// TODO: should probably be doing a clone of Sites here once Site implements Cloneable
		basis = new LinkedList<Site>();
		if (_basis != null)
			for (Site s : _basis)
				basis.add(s);
		
		// copy label
		label = _label;
	}
	
	public List<Vect> getLatticeVectors() {
		List<Vect> result = new LinkedList<Vect>();
		result.addAll(latticeVectors);
		return result;
	}
	
	/* Returns a list of length 6 containing
	 * a,b,c,alpha,beta,gamma
	 * in that order, where gamma is the angle between a and b, etc,
	 * and the angles are in degrees.
	 */
	public List<Double> getLatticeParametersDegrees() {
		List<Double> result = new LinkedList<Double>();
		for (Vect v : latticeVectors)
			result.add(v.length());
		result.add(latticeVectors.get(1).angleToInDegrees(latticeVectors.get(2)));
		result.add(latticeVectors.get(0).angleToInDegrees(latticeVectors.get(2)));
		result.add(latticeVectors.get(0).angleToInDegrees(latticeVectors.get(1)));
		return result;
	}
	
	public double[] getCellLengths() {
		double[] result = new double[Constants.numDimensions];
		for (int i = 0; i < Constants.numDimensions; i++)
			result[i] = latticeVectors.get(i).length();
		return result;
	}
	
	public double[] getCellAnglesDegrees() {
		double[] result = new double[Constants.numDimensions];
		result[0] = (latticeVectors.get(1).angleToInDegrees(latticeVectors.get(2)));
		result[1] = (latticeVectors.get(0).angleToInDegrees(latticeVectors.get(2)));
		result[2] = (latticeVectors.get(0).angleToInDegrees(latticeVectors.get(1)));
		return result;
	}
	
	/* Returns a supercell of the current cell with new lattice vectors
	 * linear combinations of those of the current cell specified by coefficents.
	 */
	public static Cell getSupercell(Cell c, List<List<Integer>> coefficients) {	
		// NB: important that all atoms are in cell for this to work.
		Cell cell = c.getCellWithAllAtomsInCell();
		
		/* Make sure we're making the right number of new vectors */
		if (coefficients.size() != Constants.numDimensions) 
			throw new IllegalArgumentException("Cell.getSupercell() given coefficients with size=" + coefficients.size());
		
//		System.out.println(coefficients.get(0).get(0) + " " + coefficients.get(0).get(1) + " " + coefficients.get(0).get(2) + " ");
//		System.out.println(coefficients.get(1).get(0) + " " + coefficients.get(1).get(1) + " " + coefficients.get(1).get(2) + " ");
//		System.out.println(coefficients.get(2).get(0) + " " + coefficients.get(2).get(1) + " " + coefficients.get(2).get(2) + " ");

		/* Make the new vectors by taking linear combinations of the old ones */
		List<Vect> newVectors = new LinkedList<Vect>();
		for (int i = 0; i < Constants.numDimensions; i++) {
			List<Integer> coefs = coefficients.get(i);
			/* Make sure we have the right number of coefficients */
			if (coefs.size() != Constants.numDimensions) 
				throw new IllegalArgumentException("Cell.getSupercell() given coefs with size=" + coefs.size());
			/* Make the new lattice vector and add it to the newVectors list */
			// TODO: let's be honest, we could do this better.  ffs, use some jama.
			Double x = coefs.get(0)*(cell.latticeVectors.get(0)).getCartesianComponents().get(0)
						+ coefs.get(1)*(cell.latticeVectors.get(1)).getCartesianComponents().get(0)
						+ coefs.get(2)*(cell.latticeVectors.get(2)).getCartesianComponents().get(0);
			Double y = coefs.get(0)*(cell.latticeVectors.get(0)).getCartesianComponents().get(1)
						+ coefs.get(1)*(cell.latticeVectors.get(1)).getCartesianComponents().get(1)
						+ coefs.get(2)*(cell.latticeVectors.get(2)).getCartesianComponents().get(1);
			Double z = coefs.get(0)*(cell.latticeVectors.get(0)).getCartesianComponents().get(2)
						+ coefs.get(1)*(cell.latticeVectors.get(1)).getCartesianComponents().get(2)
						+ coefs.get(2)*(cell.latticeVectors.get(2)).getCartesianComponents().get(2);
			newVectors.add(new Vect(x, y, z));
		}
		/* Strategy for finding all the sites to make the basis of the supercell
		 * TODO:document me :)
		 */
		int minx = 0, maxx = 0, miny = 0, maxy = 0, minz = 0, maxz = 0;
	//	int minx = -1, maxx = 1, miny = -1, maxy = 1, minz = -1, maxz = 1;
		for (int i = 0; i <= 1; i++)
			for (int j = 0; j <= 1; j++)
				for (int k = 0; k <= 1; k++) {
					int combox = i * coefficients.get(0).get(0) + j * coefficients.get(1).get(0) + k * coefficients.get(2).get(0);
					minx = Math.min(minx, combox);
					maxx = Math.max(maxx, combox);
					int comboy = i * coefficients.get(0).get(1) + j * coefficients.get(1).get(1) + k * coefficients.get(2).get(1);
					miny = Math.min(miny, comboy);
					maxy = Math.max(maxy, comboy);
					int comboz = i * coefficients.get(0).get(2) + j * coefficients.get(1).get(2) + k * coefficients.get(2).get(2);
					minz = Math.min(minz, comboz);
					maxz = Math.max(maxz, comboz);
				}
		
	//	System.out.println("minx " + minx + " maxx " + maxx + " miny " + miny + " maxy " + maxy + " minz " + minz + " maxz " + maxz);
		List<Site> newBasisCandidates = new LinkedList<Site>();
		for (Site s : cell.basis) {
			for(int i = minx; i <= maxx; i++) {
				for(int j = miny; j <= maxy; j++) {
					for (int k = minz; k <= maxz; k++) {
						Vect candidateSiteLoc = s.getCoords().plus(new Vect(new Double(i),new Double(j),new Double(k),cell.latticeVectors));
						List<Double> newFracCoords = candidateSiteLoc.getComponentsWRTBasis(newVectors);
						// if (candidate is for-real in the new cell)
				//		if (newFracCoords.get(0) > -Constants.epsilon && newFracCoords.get(0) <= 1 - Constants.epsilon
				//				&& newFracCoords.get(1) > -Constants.epsilon && newFracCoords.get(1) <= 1 - Constants.epsilon
				//				&& newFracCoords.get(2) > -Constants.epsilon && newFracCoords.get(2) <= 1 - Constants.epsilon)
				//		if (newFracCoords.get(0) >= 0  && newFracCoords.get(0) < 1
				//				&& newFracCoords.get(1) >= 0 && newFracCoords.get(1) < 1 
				//				&& newFracCoords.get(2) >= 0 && newFracCoords.get(2) < 1 )
							newBasisCandidates.add((new Site(s.getElement(),candidateSiteLoc.getVectShiftedIntoBasis(newVectors))));
					}
				}
			}
		}
		
		for (int i = 0; i < newBasisCandidates.size(); i++) {
			Cell cWithRedundantSites = new Cell(newVectors, newBasisCandidates, cell.getLabel());
			Site s = cWithRedundantSites.getSite(i);
			for (Site r : cWithRedundantSites.getAtomsInSphereSorted(s.getCoords(), Constants.epsilon)) {
				if (r != s) {					
					newBasisCandidates.remove(r);
				}
			}			
		}
		
		return new Cell(newVectors, newBasisCandidates, cell.getLabel());
	}
	
	public List<Site> getSites() {
		return basis;
	}
	
	public int	getBasisSize() {
		return basis.size();
	}
	
	/* TODO: revise this if Sites ever support partial occupancies or non-element decorations */
	public int getNumSitesWithElement(Element e) {
		int result = 0;
		for (Site s : basis) {
			if (s.getElement().equals(e))
				result++;
		}
		return result;
	}
	
	public List<Site> getSitesWithElement(Element e) {
		List<Site> result = new ArrayList<Site>();
		for (Site s : basis) {
			if (s.getElement().equals(e))
				result.add(s);
		}
		return result;
	}
	
	public double[] getLatticeVectorsArray() {
		double result[] = new double[9];
		result[0] = latticeVectors.get(0).getCartesianComponents().get(0);
		result[1] = latticeVectors.get(0).getCartesianComponents().get(1);
		result[2] = latticeVectors.get(0).getCartesianComponents().get(2);
		result[3] = latticeVectors.get(1).getCartesianComponents().get(0);
		result[4] = latticeVectors.get(1).getCartesianComponents().get(1);
		result[5] = latticeVectors.get(1).getCartesianComponents().get(2);
		result[6] = latticeVectors.get(2).getCartesianComponents().get(0);
		result[7] = latticeVectors.get(2).getCartesianComponents().get(1);
		result[8] = latticeVectors.get(2).getCartesianComponents().get(2);
		return result;
	}
	
	public static List<Vect> getVectorsfromLParamsDegrees(double a, double b, double c, double alpha, double beta, double gamma) {
		// make new axes:
		List<Vect> newCellVectors = new LinkedList<Vect>();
		// make \vec{a}:
		newCellVectors.add(new Vect(a,0.0,0.0));
		// make \vec{b}:
		newCellVectors.add(new Vect(b*Math.cos(gamma*Math.PI/180),b*Math.sin(gamma*Math.PI/180),0.0));
		// make \vec{c}:
		double c1 = c * Math.cos(beta*Math.PI/180);
		double c2 = (c*Math.cos(alpha*Math.PI/180) - Math.cos(gamma*Math.PI/180)*c*Math.cos(beta*Math.PI/180))/(Math.sin(gamma*Math.PI/180));
		double c3 = Math.sqrt(c*c - c1*c1 - c2*c2);
		newCellVectors.add(new Vect(c1,c2,c3));
		
		return newCellVectors;
	}
	
	public boolean isMalformed() {
		Vect v = new Vect(0.0, 0.0, 0.0);
		
		// check if we have a degenerate lattice
		try {
			v.getComponentsWRTBasis(this.getLatticeVectors());
		} catch( IllegalArgumentException x ) {
			return true;
		}
		
		return false;
	}
	
	public Cell getCellRotatedIntoPrincDirs() {
		// so we keep all the sites with the same fractional coordinates
		// and just make new axes with the same lengths/angles but in
		// better directions
		
		List<Double> lParams = getLatticeParametersDegrees();
		double a = lParams.get(0);
		double b = lParams.get(1);
		double c = lParams.get(2);
		double alpha = lParams.get(3);
		double beta = lParams.get(4);
		double gamma = lParams.get(5);
		List<Vect> newCellVectors = getVectorsfromLParamsDegrees(a, b, c, alpha, beta, gamma);
		
		// make new sites with the same fractional coords in the new cell
		List<Site> newBasis = new LinkedList<Site>();
		for (Site s : basis) {
			List<Double> fracCoords = s.getCoords().getComponentsWRTBasis(latticeVectors);
			newBasis.add(new Site(s.getElement(), new Vect(fracCoords.get(0), fracCoords.get(1), fracCoords.get(2), newCellVectors)));
		}
		
		return new Cell (newCellVectors, newBasis, getLabel());
	}
	
	public double getVolume() {
		double a[] = getLatticeVectorsArray();
		return a[0]*a[4]*a[8] + a[1]*a[5]*a[6] + a[2]*a[3]*a[7]
		       - a[2]*a[4]*a[6] - a[1]*a[3]*a[8] - a[0]*a[5]*a[7];
	}
	
	
	public void writeCIF(String outFile) {
		// Note that this will overwrite outFile. If this is not desired behavior,
		// the calling routine should check if the output file already exists 
		
		//TODO: option to write cartesian as opposed to fractional coords
		
		// Write the output 
		Utility.writeStringToFile(getCIF(), outFile);
	}
	
	public String getCIF() {		
		//TODO: option to write cartesian as opposed to fractional coords?
		
		StringBuilder result = new StringBuilder();
		
		// Write header stuff 
		result.append("data_cif\n\n");
		result.append("_audit_creation_method 'Autocreated by VaspOut'\n\n"); //TODO: something better from VaspData
		result.append("_symmetry_space_group_name_H-M 'P 1'\n");
		result.append("_symmetry_Int_Tables_number 1\n");
		result.append("_symmetry_cell_setting triclinic\n\n");
		// Write lattice params 
		List<Double> latticeParams = getLatticeParametersDegrees();
		result.append("_cell_length_a " + latticeParams.get(0) + "\n");
		result.append("_cell_length_b " + latticeParams.get(1) + "\n");
		result.append("_cell_length_c " + latticeParams.get(2) + "\n");
		result.append("_cell_angle_alpha " + latticeParams.get(3) + "\n");
		result.append("_cell_angle_beta " + latticeParams.get(4) + "\n");
		result.append("_cell_angle_gamma " + latticeParams.get(5) + "\n\n");
		// Write sites 
		result.append("loop_\n");
		result.append("_atom_site_label\n");
		result.append("_atom_site_fract_x\n");
		result.append("_atom_site_fract_y\n");
		result.append("_atom_site_fract_z\n");
		result.append("_atom_site_occupancy\n");
		List<Site> basis = getSites();
		for (Element e : this.getComposition().getElements()) {
			int counter = 1;
			for (Site s : this.getSitesWithElement(e)) {
				result.append(s.getElement().getSymbol() + counter + " ");	
				List<Double> coords = s.getCoords().getComponentsWRTBasis(getLatticeVectors());
				for (int i = 0; i < Constants.numDimensions; i++)
					result.append(coords.get(i) + " ");
				result.append("1.0000"); // TODO: fixme: real occupancy 
				result.append("\n");
				counter ++;
			}
		}
		
		return result.toString();
	}
	
	
	public String toStringJustVectors() {
		StringBuilder output = new StringBuilder();
		double a[] = getLatticeVectorsArray();
		output.append(new PrintfFormat("  Lattice vectors: %10.6f").sprintf(a[0]));
		output.append(new PrintfFormat(" %10.6f").sprintf(a[1]));
		output.append(new PrintfFormat(" %10.6f\n                   ").sprintf(a[2]));
		output.append(new PrintfFormat("%10.6f ").sprintf(a[3]));
		output.append(new PrintfFormat("%10.6f ").sprintf(a[4]));
		output.append(new PrintfFormat("%10.6f\n                   ").sprintf(a[5]));
		output.append(new PrintfFormat("%10.6f ").sprintf(a[6]));
		output.append(new PrintfFormat("%10.6f ").sprintf(a[7]));
		output.append(new PrintfFormat("%10.6f\n").sprintf(a[8]));
		output.append("  Volume = " + getVolume());
		return output.toString();
	}
	
	public String getLabel() {
		return label;
	}
	
	public Cell getCellWithAllAtomsInCell() {
		
		List<Site> newBasis = new LinkedList<Site>();
		
		for (Site s : basis)
			newBasis.add(new Site(s.getElement(), s.getCoords().getVectShiftedIntoBasis(latticeVectors)));
		
		
		return new Cell(latticeVectors, newBasis, label);
	}
	
	public Composition getComposition() {
		Map <Element, Integer> components = new HashMap <Element, Integer>();
		
		for (Site s : basis)
			if (components.containsKey(s.getElement()))
				components.put(s.getElement(), components.get(s.getElement()) + 1);
			else
				components.put(s.getElement(), 1);
		
		return new Composition(components, false);
	}
	
	public String toString() {
		StringBuilder output = new StringBuilder();
		
		output.append(toStringJustVectors() + "\n");
		
		output.append("  " + getBasisSize()  + " sites: \n");
		for (Site s : basis) {
			output.append("  " + s.getElement().getSymbol());
			for (double d : s.getCoords().getComponentsWRTBasis(getLatticeVectors()))
				output.append(" " + d);
			output.append("\n");
		}
			//output.append("  " + s.toString() + "\n");
		
		return output.toString();
	}

	public int getNumSites() {
		return basis.size();
	}
	
    private static double SignNoZero(double in){
        return in >= 0 ? 1 : -1; 
    }
    
    List<Vect> getRecipLVects() {
    	Vect a1 = latticeVectors.get(0);
    	Vect a2 = latticeVectors.get(1);
    	Vect a3 = latticeVectors.get(2);
    	
    	List<Vect> result = new ArrayList<Vect>(3);
    	double oneoverdenom = 1/(a1.dot(a2.cross(a3)));
    	result.add(a2.cross(a3).scalarMult(2 * Math.PI * oneoverdenom));
    	result.add(a3.cross(a1).scalarMult(2 * Math.PI * oneoverdenom));
    	result.add(a1.cross(a2).scalarMult(2 * Math.PI * oneoverdenom));
    	return result;
    }

    /**
     * method returning the Niggli reduced cell representation for a given structure.
     * this code is adapted from convasp 5.6.
     * Warning: this code doesn't check if the cell is actually primitive
     */
    private Cell niggliReducedCell;
    public Cell getNigliReducedCell() {
    	if (niggliReducedCell != null)
    		return niggliReducedCell;
    
//    	if (true)
//throw new RuntimeException("nrc");
    	    	
        double TOL = 1e-8;

        // Initialize matrices for tranformations (3x3).
        double[][] m1 = new double[3][3];
        m1[0][1] = 1.0; m1[1][0] = 1.0; m1[2][2] = -1.0;
        double[][] m2 = new double[3][3];
        m2[0][0] = -1.0; m2[1][2] = 1.0; m2[2][1] = 1.0;
        double[][] m3 = new double[3][3];
        double[][] m4 = new double[3][3];
        double[][] m5 = Matrix.identity(3,3).getArrayCopy();
        double[][] m6 = Matrix.identity(3,3).getArrayCopy();
        double[][] m7 = Matrix.identity(3,3).getArrayCopy();
        double[][] m8 = Matrix.identity(3,3).getArrayCopy();
        m8[0][2] = 1.0;
        m8[1][2] = 1.0;

        // Initialize a, b, c, ksi, eta, zeta
        double[] indat = new double[6];
        indat[0] = this.getCellLengths()[0];
        indat[1] = this.getCellLengths()[1];
        indat[2] = this.getCellLengths()[2];
        indat[3] = latticeVectors.get(1).angleToInDegrees(latticeVectors.get(2)) * Math.PI / 180;
        indat[4] = latticeVectors.get(0).angleToInDegrees(latticeVectors.get(2)) * Math.PI / 180;
        indat[5] = latticeVectors.get(0).angleToInDegrees(latticeVectors.get(1)) * Math.PI / 180;
        double a = indat[0] * indat[0];
        double b = indat[1] * indat[1];
        double c = indat[2] * indat[2];
        double ksi = 2.0 * indat[1] * indat[2] * Math.cos(indat[3]);
        double eta = 2.0 * indat[0] * indat[2] * Math.cos(indat[4]);
        double zeta = 2.0 * indat[0] * indat[1] * Math.cos(indat[5]);
        // Dummy variables
        double temp, temp1, temp2, temp3;
        // Initialize tranformation matrix
        double[][] P = Matrix.identity(3,3).getArrayCopy();

        // start loop
        int maxNumIters = 100000; // just some big number
        int numIters = 0;
        while (true) {
        	// break if we're in an infinite loop which can happen apparently.
        	// but rlyrly shouldnt happen much.
        	if (numIters > maxNumIters) {
        		GAOut.out().stdout("WARNING: breaking out of infinite loop in Cell.getNigliReducedCell", GAOut.WARNING);
        		break;
        	}
        	numIters++;
            // Step 1
            if (((a - b) > TOL)
                    || ((Math.abs(a - b) < TOL) && (Math.abs(ksi) > Math
                            .abs(eta)))) {
                temp = a;
                a = b;
                b = temp;
                temp = -ksi;
                ksi = -eta;
                eta = temp;
                P = (new Matrix(P)).times(new Matrix(m1)).getArray();
       //         P = MatrixMath.SquareMatrixMult(P, m1);

            }

            // Step 2
            if (((b - c) > TOL)
                    || ((Math.abs(b - c) < TOL) && (Math.abs(eta) > Math
                            .abs(zeta)))) {
                temp = b;
                b = c;
                c = temp;
                temp = -eta;
                eta = -zeta;
                zeta = temp;
                P = (new Matrix(P)).times(new Matrix(m2)).getArray();
           //     P = MatrixMath.SquareMatrixMult(P, m2);

                continue;
            }

            // Step 3
            if ((ksi * eta * zeta) > 0.0) {
                m3[0][0] = SignNoZero(ksi);
                m3[1][1] = SignNoZero(eta);
                m3[2][2] = SignNoZero(zeta);
                P = (new Matrix(P)).times(new Matrix(m3)).getArray();
          //      P = MatrixMath.SquareMatrixMult(P, m3);
                ksi = Math.abs(ksi);
                eta = Math.abs(eta);
                zeta = Math.abs(zeta);
            }

            // Step 4

            if (ksi * eta * zeta <= 0) {
                m4[0][0] = -(SignNoZero(ksi));
                m4[1][1] = -(SignNoZero(eta));
                m4[2][2] = -(SignNoZero(zeta));
                if (Math.abs(ksi) < TOL)
                    m4[0][0] = m4[1][1] * m4[2][2];
                if (Math.abs(eta) < TOL)
                    m4[1][1] = m4[0][0] * m4[2][2];
                if (Math.abs(zeta) < TOL)
                    m4[2][2] = m4[0][0] * m4[1][1];
                P = (new Matrix(P)).times(new Matrix(m4)).getArray();
        //        P = MatrixMath.SquareMatrixMult(P, m4);
                ksi = -Math.abs(ksi);
                eta = -Math.abs(eta);
                zeta = -Math.abs(zeta);

            }

            // Step 5
            if (((Math.abs(ksi) - b) >= TOL)
                    || ((Math.abs(ksi - b) < TOL) && (2.0 * eta < zeta))
                    || ((Math.abs(ksi + b) < TOL) && (zeta < 0))) {
                m5[1][2] = -(SignNoZero(ksi));
                P = (new Matrix(P)).times(new Matrix(m5)).getArray();
          //      P = MatrixMath.SquareMatrixMult(P, m5);
                temp1 = b + c - ksi * SignNoZero(ksi);
                temp2 = eta - zeta * SignNoZero(ksi);
                temp3 = ksi - 2.0 * b * SignNoZero(ksi);
                c = temp1;
                eta = temp2;
                ksi = temp3;

                continue;
            }

            // Step 6
            if (((Math.abs(eta) - a) > TOL)
                    || ((Math.abs(eta - a) < TOL) && (2.0 * ksi < zeta))
                    || ((Math.abs(eta + a) < TOL) && (zeta < 0))) {
                m6[0][2] = -SignNoZero(eta);
                P = (new Matrix(P)).times(new Matrix(m6)).getArray();
          //      P = MatrixMath.SquareMatrixMult(P, m6);

                temp1 = a + c - eta * SignNoZero(eta);
                temp2 = ksi - zeta * SignNoZero(eta);
                temp3 = eta - 2.0 * a * SignNoZero(eta);
                c = temp1;
                ksi = temp2;
                eta = temp3;

            }

            // Step 7
            if (((Math.abs(zeta) - a) > TOL)
                    || ((Math.abs(zeta - a) < TOL) && (2.0 * ksi < eta))
                    || ((Math.abs(zeta + a) < TOL) && (eta < 0))) {
                m7[0][1] = -SignNoZero(zeta);
                P = (new Matrix(P)).times(new Matrix(m7)).getArray();
         //       P = MatrixMath.SquareMatrixMult(P, m7);
                temp1 = a + b - zeta * SignNoZero(zeta);
                temp2 = ksi - eta * SignNoZero(zeta);
                temp3 = zeta - 2.0 * a * SignNoZero(zeta);
                b = temp1;
                ksi = temp2;
                zeta = temp3;
                continue;
            }

            // Step 8
            if ((ksi + eta + zeta + a + b < 0)
                    || ((ksi + eta + zeta + a + b < 0) && (2.0 * (a + eta)
                            + zeta > 0))) {
            	P = (new Matrix(P)).times(new Matrix(m8)).getArray();
    //            P = MatrixMath.SquareMatrixMult(P, m8);
                temp1 = a + b + c + ksi + eta + zeta;
                temp2 = 2.0 * b + ksi + zeta;
                temp3 = 2.0 * a + eta + zeta;
                c = temp1;
                ksi = temp2;
                eta = temp3;

                continue;
            }
            
            break;
        }
        // 
        double[] outdat = new double[6];
        outdat[0] = Math.sqrt(a);
        outdat[1] = Math.sqrt(b);
        outdat[2] = Math.sqrt(c);
        outdat[3] = Math.acos(ksi / 2.0 / outdat[1] / outdat[2]);
        outdat[4] = Math.acos(eta / 2.0 / outdat[0] / outdat[2]);
        outdat[5] = Math.acos(zeta / 2.0 / outdat[0] / outdat[1]);
        
        // Get Niggli cell
        List<List<Integer>> coefficients = new LinkedList<List<Integer>>();
        for (int i = 0; i < Constants.numDimensions; i++) {
        	List<Integer> col = new LinkedList<Integer>();
        	for (int j = 0; j < Constants.numDimensions; j++) {
        		int entry = Math.round(Math.round((P[j][i])));
        		if (entry != P[j][i]) // P is double[][] but damn well better hold only integer values
        			System.out.println("ERROR: entry != P[i][j] in getNiggliReducedCell");
        		col.add(entry);
        	}
        	coefficients.add(col);
        }
        Cell result = getSupercell(this, coefficients);
        
        // Checks
        // Make sure that a,b,c,alpha,beta,gamma are the same from
        // direct calculation and from using P to get niggli_lat.
        double[] poutdat = new double[6];
        poutdat[0] = result.getCellLengths()[0];
        poutdat[1] = result.getCellLengths()[1];
        poutdat[2] = result.getCellLengths()[2];
        poutdat[3] = result.getCellAnglesDegrees()[0] * Math.PI / 180;
        poutdat[4] = result.getCellAnglesDegrees()[1] * Math.PI / 180;
        poutdat[5] = result.getCellAnglesDegrees()[2] * Math.PI / 180;
            		
        int flag = 0;
        for (int i = 0; i < 6; i++) {
            if (Math.abs(poutdat[i] - outdat[i]) > 2 * TOL) {
                flag = 1;
            }
        }
        if (flag == 1) {
                      System.out.println("ERROR: CellReduceFuncs/GetNiggliCell");
                      System.out.println("ERROR: Lattice parameters/angles as calculated");
                      System.out.println("ERROR: with Niggli algorithm and P do not match." );
                      System.out.println("ERROR: a,b,c,alpha,beta,gamma from direct algorithm: ");
                      for (int j = 0; j < 6; j++)
                          System.out.println(outdat[j]);
                      System.out.println("ERROR: a,b,c,alpha,beta,gamma from P matrix: ");
                      for (int j = 0; j < 6; j++)
                          System.out.println(poutdat[j]);
        }

  //      result.setDescription(this.getDescription());
  /*      for (int i = 0; i < this.m_sites.length; i++) {
            Vect p = this.m_sites[i].m_loc;
            Vect scaledPoint = p.inBasis(nigglibasis);
            scaledPoint = result.translateIntoCell(scaledPoint);
            result.addSiteChecked(new Site(scaledPoint, this.m_sites[i].m_sp,
                    this.m_sites[i].m_occ));
        } */
        
        if (this.getBasisSize() != result.getBasisSize()) {
        	GAOut.out().stdout("ERROR: Niggli cell reduction gained or lost atoms", GAOut.CRITICAL);
        	GAOut.out().stdout(this.toString(), GAOut.CRITICAL);
        	GAOut.out().stdout(result.toString(), GAOut.CRITICAL);
        	(new Exception()).printStackTrace();
        }
        
        niggliReducedCell = result;
        return result; 
    }


	public Site getSite(int i) {
		return basis.get(i);
	}
	
	public Cell scaleTo(double givenVol) {
		double currentVol = getVolume();
		double scaleFactor = Math.pow((givenVol / currentVol) , (1.0/Constants.numDimensions));
		
		List<Vect> newVects = new LinkedList<Vect>();
		for (Vect v : latticeVectors)
			newVects.add(v.scalarMult(scaleFactor));
		
		List<Site> newSites = new LinkedList<Site>();
		for (Site s : basis)
			newSites.add(new Site(s.getElement(), s.getCoords().scalarMult(scaleFactor)));
		
		return new Cell(newVects, newSites);
	}
	
	private boolean pointIsInParallelopiped(Vect corner, List<Vect>latticeVectors, int num) {
		for (Double d : corner.getComponentsWRTBasis(latticeVectors))
			if (d > num || d < -num)
				return false;
		return true;
	}
	
	private boolean parallelopipedDefEncompassesSphere(List<Vect> latticeVectors, int num, Vect center, double dist) {
		// get corners of cube that encompasses sphere
		List<Vect> corners = new LinkedList<Vect>();
		for (int i = -1; i <= 1; i = i + 2)
			for (int j = -1; j <= 1; j = j + 2)
				for (int k = -1; k <= 1; k = k + 2)
					corners.add(
							center.plus( (Vect.xHat().scalarMult((double)i)
										 .plus(Vect.yHat().scalarMult((double)j))
										 .plus(Vect.zHat().scalarMult((double)k))).scalarMult(dist) )
								);
		
		// check that all of the corners are in the parallelogram
		for (Vect corner : corners)
			if (! pointIsInParallelopiped(corner, latticeVectors, num))
				return false;
		return true;
		
	}
	
	private String fsoOutput;
	public String getFSOOutput() {
		if (fsoOutput == null)
			fsoOutput = Isotropy.getFindsymOut(this);
		return fsoOutput;
	}


	private Cell wyckoffCell;
	public Cell getWyckoffCell() {
		if (wyckoffCell == null)
			wyckoffCell = Isotropy.getWyckoffCell(this);
		return wyckoffCell;
	}
	
	// average position of all the sites.  doesnt necessarily make sense when using PBCs!
	public Vect getCentroid() {
		Vect result = new Vect(0.0, 0.0, 0.0);
		
		for (Site s : this.getSites()) {
			result = result.plus(s.getCoords());
		}
		
		result = result.scalarMult(1.0 / this.getNumSites());
		
		return result;
	}
	
	// use PBCs by default
	public boolean matchesCell(Cell other, double atomicMisfit, double lengthMisfit, double angleMisfit) {
		return matchesCellPBCs(other, atomicMisfit, lengthMisfit, angleMisfit);
	}
	
	
	// we want to : 
	// - translate the centroid to the origin
	// - apply rotations so that site s1 is on the (+) axis and site s2 is in the positive quadrange of the x-y plane.
	private Cell getCellRotatedWithSitesInPrincDirections(Vect v1, Vect v2) {
		Vect minusCentroid = getCentroid().scalarMult(-1.0);
		
		Vect p1 = v1.plus(minusCentroid);

		// get rotation matrix to rotate p1 onto the x axis
		// this is a rotation of angle cos^-1 ( (p1 . x-axis) / (2 * |p1| * |xaxis|  ) 
		//     about the vector (p1 x xaxis)   (<- x = cross product)
		double p1_rot_angle = p1.angleToInDegrees(Vect.xHat());
		Vect p1_rot_vect = p1.cross(Vect.xHat());
		Matrix p1_rot_mat = Utility.getRotationMatrixDegrees(p1_rot_angle, p1_rot_vect);
		
		// make p2 by v1 translated by centroid and by the first rotation
		// get rotation matrix to then rotate p2 onto the x-y plane
		// this is a rotation about the x-axis of angle found by:
		// - project p2 into the y-z plane, and find the angle between the y-axis and that projection
		Vect p2 = v2.plus(minusCentroid).leftMultByMatrix(p1_rot_mat);
		Vect p2proj = new Vect(0.0, p2.getCartesianComponents().get(1), p2.getCartesianComponents().get(2));
		double p2_rot_angle_deg = p2proj.angleToInDegrees(Vect.yHat());
		
		Matrix p2_rot_mat;
		if (p2.getCartesianComponents().get(2) < 0) // make sure we're rotating into the POSITIVE quadrant of the x-y plane
			p2_rot_mat = Utility.getRotationMatrixDegrees(p2_rot_angle_deg, Vect.xHat());
		else
			p2_rot_mat = Utility.getRotationMatrixDegrees(360.0 - p2_rot_angle_deg, Vect.xHat());
		
		Matrix total_rot = p2_rot_mat.times(p1_rot_mat);
		
		List<Vect> newVects = new ArrayList<Vect>();
		for (Vect b : this.getLatticeVectors())
			newVects.add(b.leftMultByMatrix(total_rot));
		
		List<Site> newSites = new ArrayList<Site>();
		for (Site s : this.getSites())
			newSites.add(new Site(s.getElement(), s.getCoords().plus(minusCentroid).leftMultByMatrix(total_rot)));
		
		return new Cell(newVects, newSites, this.getLabel());
	}
	
	/* Comparing cluster geometries:
	 *  - let c1, c2 be the centers of both clusters
	 *  - select p11 and p12 two arbitrary sites from cluster 1
	 *  - search over all cluster 2 sites to find those, p21, st
	 *    d(c1, p11) = d(c2, p21) (up to a tolerance) (compare elements also)
	 *  - for each of these, search over all cluster 2 sites p22 st
	 *    d(c1, p21) = d(c2, p22) _and_ d(p12,p11) = d(p22,p21) (compare elements also)
	 *  - this fixes an orientation of the two clusters
	 *  - so, make a rotated + shifted copy of cluster 2 so that c1,c2  p21,211  and p22,p12 all match up
	 *  - if all the other points line up too, we have a match
	 */
	public boolean matchesCellNoPBCs(Cell other, double atomicMisfit) {
		
		if (this.getNumSites() != other.getNumSites())
			return false;
		if (! this.getComposition().equals(other.getComposition()))
			return false;
		if (this.getNumSites() < 2)
			return true;
		
		Vect c1 = this.getCentroid();
		Vect c2 = other.getCentroid();
		
		Site p11 = this.getSite(0);
		Site p12;
		// try to find a p12 such that p12, p11, and c1 are not collinear
		int p12indx = 1;
		do {
			p12 = this.getSite(p12indx);
			p12indx++;
		} while (p12indx < this.getNumSites() && Vect.pointsAreCollinear(c1,p12.getCoords(),p11.getCoords(),atomicMisfit));
		
		
		Cell thisStd = this.getCellRotatedWithSitesInPrincDirections(p11.getCoords(), p12.getCoords());
		
		for (Site p21 : other.getSites()) {
			if (Math.abs( c1.getCartDistanceTo(p11.getCoords()) - c2.getCartDistanceTo(p21.getCoords()) ) > atomicMisfit)
				continue;
			if (p11.getElement() != p21.getElement())
				continue;
			
			for (Site p22 : other.getSites()) {
				if (p21 == p22)
					continue;
				if (Math.abs( c1.getCartDistanceTo(p12.getCoords()) - c2.getCartDistanceTo(p22.getCoords()) ) > atomicMisfit)
					continue;
				if (Math.abs( p11.getCoords().getCartDistanceTo(p12.getCoords()) - p21.getCoords().getCartDistanceTo(p22.getCoords()) ) > atomicMisfit)
					continue;
				if (p12.getElement() != p22.getElement())
					continue;
				
				// at this point, the correspondence c1<->c2, p11<->p21, p12<->p22 are 3 points that match between
				// the two clusters.  we can now imagine rotate/shifting one so that these 3 points lie on top of
				// each other, and then check that all the other points line on top of each other too to establish
				// the match.
				// the idea:
				// - generate transforms for both clusters such that:
				//   - centroid is at the origin
				//   - the first points (i.e. p11 and p21) are rotated onto the x axis
				//   - the second points (i.e. p12 and p22) are rotated onto the positive quadrante of the x-y plane
				Cell otherCandStd = other.getCellRotatedWithSitesInPrincDirections(p21.getCoords(), p22.getCoords());
				boolean successfulMatchOfWholeStructs = true;
				for (Site s : otherCandStd.getSites()) {
					boolean foundMatchForSite = false;
					for (Site t : thisStd.getAtomsInSphereSortedIgnoringPBCs(s.getCoords(), atomicMisfit)) {
						if (s.getElement() == t.getElement()) {
							foundMatchForSite = true;
							break;
						}
					}
					if (!foundMatchForSite) {
						successfulMatchOfWholeStructs = false;
						break;
					}
				}
				if (successfulMatchOfWholeStructs)
					return true;
			}
		}
		
		return false;
	}
	
	// angle misfit given in degrees
	public boolean matchesCellPBCs(Cell other, double atomicMisfit, double lengthMisfit, double angleMisfit) {
		// get Niggli reduced cell of this and other
		// and shift this cell so that one atom is at the origin
		/*
		Cell n1 = this.getNigliReducedCell().getWyckoffCell();
		Cell n2 = null;
		
		if (other.getBasisSize() <= 0)
			n2 = other.getNigliReducedCell().getWyckoffCell();
		else
			n2 = other.getNigliReducedCell().getWyckoffCell().getCellWithSiteIShiftedToOrigin(0);
			*/
		Cell n1 = this.getWyckoffCell();
		Cell n2 = null;
		
		if (other.getBasisSize() <= 0)
			n2 = other.getWyckoffCell();
		else
			n2 = other.getWyckoffCell().getCellWithSiteIShiftedToOrigin(0);
		
		if (n1 == null || n2 == null)
			return false;
		
		// if they have a different number of wyckoff positions, they don't match
		if (n1.getBasisSize() != n2.getBasisSize())
			return false;
		int numWSites = n1.getBasisSize();
		
		// also check that the particular elements match
		for (Element e : n1.getComposition().getElements())
			if (n2.getNumSitesWithElement(e) != n1.getNumSitesWithElement(e))
				return false;
		
		//for (each of the choices of others lattice vectors being considered a, b, and c)
		for (Cell n1rotated : n1.getCellWithAlternateAxesLabeling()) { 
			List<Double> parms1 = n1rotated.getLatticeParametersDegrees();
			List<Double> parms2 = n2.getLatticeParametersDegrees();
			//	if (lattices dont match up to tolerance), continue
			if (Math.abs(parms1.get(0) - parms2.get(0)) > lengthMisfit
				|| Math.abs(parms1.get(1) - parms2.get(1)) > lengthMisfit
				|| Math.abs(parms1.get(2) - parms2.get(2)) > lengthMisfit
				|| Math.abs(parms1.get(3) - parms2.get(3)) > angleMisfit
				|| Math.abs(parms1.get(4) - parms2.get(4)) > angleMisfit
				|| Math.abs(parms1.get(5) - parms2.get(5)) > angleMisfit)
				continue;
			
			// for (each choice of atom in other cell being at the origin) 
			for(int i = 0; i < numWSites; i++) {
				Cell n1shifted = n1rotated.getCellWithSiteIShiftedToOrigin(i);
				boolean allSitesLineUp = true;
				for (int j = 0; j < numWSites; j++) {
					//if (site j in n1s does not line up with a site in n2) 
					if (!sitesArrayContainsElement(n1shifted.getAtomsInSphereSorted(n2.getSite(j).getCoords(), atomicMisfit), n2.getSite(j).getElement())) {
						allSitesLineUp = false;
						break;
					}
				}
				if (allSitesLineUp)
					return true;
			}			
		}

		return false;
	}
	
	private boolean sitesArrayContainsElement(List<Site> sites, Element e) {
		for (Site s : sites)
			if (s.getElement() == e)
				return true;
		return false;
	}
	
	private List<Cell> getCellWithAlternateAxesLabeling() {
		List<Cell> result = new LinkedList<Cell>();
		
		// add unaltered cell
		result.add(this);
		
		// can't think of a good way to do this in a loop...
		List<Vect> vects;
		List<Site> sites = new LinkedList<Site>();
		double a = this.getLatticeParametersDegrees().get(0);
		double b = this.getLatticeParametersDegrees().get(1);
		double c = this.getLatticeParametersDegrees().get(2);
		double alpha = this.getLatticeParametersDegrees().get(3);
		double beta = this.getLatticeParametersDegrees().get(4);
		double gamma = this.getLatticeParametersDegrees().get(5);
		
		// make a new cell by switching a,c
		vects = getVectorsfromLParamsDegrees(c, b, a, gamma, beta, alpha);
		sites.clear();
		for (Site s : this.getSites()) {
			List<Double> coords = s.getCoords().getComponentsWRTBasis(this.getLatticeVectors());
			double temp = coords.get(0);
			coords.set(0, coords.get(2));
			coords.set(2, temp);
			sites.add(new Site(s.getElement(), new Vect(coords, vects)));
		}
		result.add(new Cell(vects, sites));
		
		// make a new cell by switching a,b
		vects = getVectorsfromLParamsDegrees(b, a, c, beta, alpha, gamma);
		sites.clear();
		for (Site s : this.getSites()) {
			List<Double> coords = s.getCoords().getComponentsWRTBasis(this.getLatticeVectors());
			double temp = coords.get(0);
			coords.set(0, coords.get(1));
			coords.set(1, temp);
			sites.add(new Site(s.getElement(), new Vect(coords, vects)));
		}
		result.add(new Cell(vects, sites));
		
		// make a new cell by switching b,c
		vects = getVectorsfromLParamsDegrees(a, c, b, alpha, gamma, beta);
		sites.clear();
		for (Site s : this.getSites()) {
			List<Double> coords = s.getCoords().getComponentsWRTBasis(this.getLatticeVectors());
			double temp = coords.get(1);
			coords.set(1, coords.get(2));
			coords.set(2, temp);
			sites.add(new Site(s.getElement(), new Vect(coords, vects)));
		}
		result.add(new Cell(vects, sites));
		
		return result;
	}
	
	public Cell getCellWithSiteIShiftedToOrigin(int i) {
		List<Site> newSites = new LinkedList<Site>();
		
		List<Vect> basisVects = this.getLatticeVectors();
		
		Site siteI = this.getSite(i);
		double iFracX = siteI.getCoords().getComponentsWRTBasis(basisVects).get(0);
		double iFracY = siteI.getCoords().getComponentsWRTBasis(basisVects).get(1);
		double iFracZ = siteI.getCoords().getComponentsWRTBasis(basisVects).get(2);
		
		for (Site s : this.getSites()) {
			
			double newFracX = s.getCoords().getComponentsWRTBasis(basisVects).get(0) - iFracX;
			double newFracY = s.getCoords().getComponentsWRTBasis(basisVects).get(1) - iFracY;
			double newFracZ = s.getCoords().getComponentsWRTBasis(basisVects).get(2) - iFracZ;
			
			if (newFracX < 0) newFracX += 1;
			if (newFracY < 0) newFracY += 1;
			if (newFracZ < 0) newFracZ += 1;
			
			newSites.add(new Site(s.getElement(), new Vect(newFracX, newFracY, newFracZ, basisVects)));
		}
		
		return new Cell(this.getLatticeVectors(), newSites, this.getLabel());
	}
	
	public List<Site> getAtomsInSphereSortedIgnoringPBCs(final Vect center, double dist) {
		List<Site> result = new ArrayList<Site>();
		
		for (Site s : basis)
			if (center.getCartDistanceTo(s.getCoords()) < dist)
				result.add(s);
		
		return result;
	}

	// return nearest neighborslist sorted by distance from center
	public List<Site> getAtomsInSphereSorted(final Vect center, double dist) {
		return getAtomsInSphereSorted(center, dist, false);
	}

		
	public List<Site> getAtomsInSphereSorted(final Vect center, double dist, boolean newSites) {
		//Map<Site, Double> resultsMap = new HashMap<Site,Double>();
		List<Site> result = new ArrayList<Site>();
		
		// get multiple of cell large enough to encompass sphere: kind of a hack
		// 
	//	int num = 1;
	//	while(!parallelopipedDefEncompassesSphere(latticeVectors, num, center, dist))
	//		num++;
		
		List<Vect> recipL = this.getRecipLVects();
		
		List<Double> cFracCoords = center.getComponentsWRTBasis(this.getLatticeVectors());
		int maxx = (int)Math.ceil(cFracCoords.get(0) + dist * recipL.get(0).length() / (2 * Math.PI));
		int minx = (int)Math.floor(cFracCoords.get(0) - dist * recipL.get(0).length() / (2 * Math.PI));
		int maxy = (int)Math.ceil(cFracCoords.get(1) + dist * recipL.get(1).length() / (2 * Math.PI));
		int miny = (int)Math.floor(cFracCoords.get(1) - dist * recipL.get(1).length() / (2 * Math.PI));
		int maxz = (int)Math.ceil(cFracCoords.get(2) + dist * recipL.get(2).length() / (2 * Math.PI));
		int minz = (int)Math.floor(cFracCoords.get(2) - dist * recipL.get(2).length() / (2 * Math.PI));
		
		for (int i = minx; i <= maxx; i++) {
			for (int j = miny; j <= maxy; j++) {
				for (int k = minz; k <= maxz; k++) {
					for (Site s : basis) {
						Vect trialCoords = s.getCoords().plus(latticeVectors.get(0).scalarMult((double)i))
										    .plus(latticeVectors.get(1).scalarMult((double)j))
										    .plus(latticeVectors.get(2).scalarMult((double)k));
						if (center.getCartDistanceTo(trialCoords) <= dist) {
							if (newSites)
								result.add(new Site(s.getElement(), trialCoords));
							else
								result.add(s);
							continue;
						}
						//	result.add(new Site(s.getElement(), trialCoords));
					}
				}
			}
		}
		
		Comparator<Site> c = new Comparator<Site>() {
			public int compare(Site s1, Site s2) {
				if (s1.getCoords().getCartDistanceTo(center) > s2.getCoords().getCartDistanceTo(center))
					return 1;
				else if (s1.getCoords().getCartDistanceTo(center) == s2.getCoords().getCartDistanceTo(center))
					return 0;
				else 
					return -1;
			}
		};
		
		/*
		for (int i = 0; i < result.size(); i++) {
			System.out.println(result.get(i) + " : " + result.get(i).getCoords().getCartDistanceTo(center));
			for (int j = 0; j < result.size(); j++) {
				System.out.println(result.get(j) + " : " + result.get(j).getCoords().getCartDistanceTo(center) 
						+ " : " + c.compare(result.get(i), result.get(j)));
			}
		}*/
		
		Collections.sort(result, c);
		
		return result;
	}

	public Site[] getAtomsInSphereSortedSlow(Vect center, double dist) {
		Map<Site, Double> resultsMap = new HashMap<Site,Double>();
		
		// get multiple of cell large enough to encompass sphere: kind of a hack
		// 
		int num = 1;
		while(!parallelopipedDefEncompassesSphere(latticeVectors, num, center, dist))
			num++;
		
		List<Site> candidateSites = new LinkedList<Site>();
		for (int i = -num; i <= num; i++) {
			for (int j = -num; j <= num; j++) {
				for (int k = -num; k <= num; k++) {
					for (Site s : basis) {
						candidateSites.add(new Site(s.getElement(), 
								s.getCoords().plus(latticeVectors.get(0).scalarMult((double)i))
											 .plus(latticeVectors.get(1).scalarMult((double)j))
											 .plus(latticeVectors.get(2).scalarMult((double)k)) ));
					}
				}
			}
		}
		

		for (Site candidateSite : candidateSites)
			if (center.getCartDistanceTo(candidateSite.getCoords()) < dist)
				resultsMap.put(candidateSite, center.getCartDistanceTo(candidateSite.getCoords()));
		
		// sort them by distance
		int numResults = resultsMap.keySet().size();
		Site result[] = new Site[numResults];
		
		// fill out the results array, removing results from resultsMap as we go
		for (int i = 0; i < numResults; i++) {
			Iterator<Site> it = resultsMap.keySet().iterator();
			// find the nearest remaining site
			double shortestDist = Double.POSITIVE_INFINITY;
			Site nearestSite = null;
			do {
				Site s = it.next();
				double sDist = resultsMap.get(s);
				if (sDist < shortestDist) {
					shortestDist = sDist;
					nearestSite = s;
				}
			} while (it.hasNext());
			result[i] = nearestSite;
			resultsMap.remove(nearestSite);
		}
		
		return result;
	}
	
	/*
	public String getCIF() {
		StringBuilder result = new StringBuilder();
		String newline = GAUtils.newline();
		
		// get the necessary data
		double[] lengths = getCellLengths();
		double[] angles = GAUtils.anglesToDegrees(getCellAngles());
		List<Site> sites = getSites();
		
		// build the cif file
		result.append("data_cif" + newline + newline);
		result.append("_audit_creation_method 'generated by StructureOrg.java'" + newline + newline);
		result.append("_symmetry_space_group_name_H-M 'P 1'" + newline);
		result.append("_symmetry_Int_Tables_number 1" + newline);
		result.append("_symmetry_cell_setting triclinic" + newline + newline);
		result.append("_cell_length_a " + lengths[0] + newline);
		result.append("_cell_length_b " + lengths[1] + newline);
		result.append("_cell_length_c " + lengths[2] + newline);
		result.append("_cell_angle_alpha " + angles[0] + newline);
		result.append("_cell_angle_beta " + angles[1] + newline);
		result.append("_cell_angle_gamma " + angles[2] + newline + newline);
		result.append("loop_" + newline);
		result.append("_atom_site_label" + newline);
		result.append("_atom_site_fract_x" + newline);
		result.append("_atom_site_fract_y" + newline);
		result.append("_atom_site_fract_z" + newline);
		result.append("_atom_site_occupancy" + newline);
		for (int i = 0; i < sites.size(); i++) {
			List<Double> loc = sites.get(i).getCoords().getComponentsWRTBasis(getLatticeVectors());
			result.append(sites.get(i).getElement().getSymbol() + " ");
			result.append(loc.get(0) + " " + loc.get(1) + " " + loc.get(2) + " 1.0000" + newline);
		}
		
		return result.toString();
	} */
	
/*
	public static Cell parseCell(String fName, String format) {
		OBMol mol = new OBMol();
		OBConversion conv = new OBConversion();
		conv.SetInFormat(format);
		if (!conv.ReadFile(mol, fName)) {
			System.out.println("Cell.parseCell failed to read " + fName + " with format " + format + ".");
			return null;
		}
		
		OBUnitCell c = openbabel_java.toUnitCell((mol.GetData(openbabel_javaConstants.UnitCell)));
		
		List<Site> sites = new ArrayList<Site>();
		for (int i = 1; i <= mol.NumAtoms(); i++) {
			Element e = Element.getElemFromZ((int) mol.GetAtom(i).GetAtomicNum());
			sites.add(new Site(e, new Vect( mol.GetAtom(i).GetX(), 
											mol.GetAtom(i).GetY(),
											mol.GetAtom(i).GetZ())));
		}
		
		// convert format of vectors
		List<Vect> vs = new ArrayList<Vect>();
		vectorVector3 vects = c.GetCellVectors();
		for (int i = 0; i < Constants.numDimensions; i++) 
			vs.add(new Vect(vects.get(i).GetX(), vects.get(i).GetY(), vects.get(i).GetZ()));
		
		String label = "Parsed by parseCell.";
		return new Cell(vs, sites, label);
	}
	
	public String getFormattedCell(String format) {
		OBMol mol = new OBMol();
		OBConversion conv = new OBConversion();
		conv.SetOutFormat(format);
		
		// make the unit cell
		OBUnitCell c = new OBUnitCell();
		vector3 v1 = new vector3(this.latticeVectors.get(0).getCartesianComponents().get(0),
								 this.latticeVectors.get(0).getCartesianComponents().get(1),
								 this.latticeVectors.get(0).getCartesianComponents().get(2));
		vector3 v2 = new vector3(this.latticeVectors.get(1).getCartesianComponents().get(0),
								 this.latticeVectors.get(1).getCartesianComponents().get(1),
								 this.latticeVectors.get(1).getCartesianComponents().get(2));
		vector3 v3 = new vector3(this.latticeVectors.get(2).getCartesianComponents().get(0),
								 this.latticeVectors.get(2).getCartesianComponents().get(1),
								 this.latticeVectors.get(2).getCartesianComponents().get(2));
		c.SetData(v1, v2, v3);
		/*
		openbabel_javaConstants.SetData;
		openbabel_java.SetData(mol, c);
		mol.SetData(c); 
		
		return "TODO";
		// add a
	}
	
*/

	// parses a CIF file (of the format written by both this file and by GULP)
	// and returns a structure
	
	public static Cell parseCif(File cifFile) {
		String line = null;
		double la = 0, lb = 0, lc = 0, a = 0, b = 0, g = 0;

		List<Vect> latVects = null;
		List<Site> sites = new LinkedList<Site>();

		// We're going to assume that every line we see after the
		// _atom_site_occupancy
		// tag describes a site like, e.g. C 0.46502 0.01282 0.11016 1.0000
		Boolean readingInAtoms = false;

		try {
			BufferedReader r = new BufferedReader(new FileReader(cifFile));
			while ((line = r.readLine()) != null) {
				StringTokenizer t = new StringTokenizer(line);
				// skip empty lines:
				if (!t.hasMoreTokens())
					continue;
				if (readingInAtoms) {
					String siteStr = t.nextToken();
					String symbol = null;
					if (siteStr.length() < 2)
						symbol = siteStr;
					else
					    symbol = siteStr.substring(1, 2).matches("[A-Za-z]") ? siteStr.substring(0, 2) : siteStr.substring(0, 1); 		
					
					double x = Double.parseDouble(t.nextToken());
					double y = Double.parseDouble(t.nextToken());
					double z = Double.parseDouble(t.nextToken());
					// String occupancy = t.nextToken();
					sites.add(new Site(Element.getElemFromSymbol(symbol), new Vect(x,y,z, latVects)));
				} else {
					String firstWord = t.nextToken();
					if (firstWord.equals("_cell_length_a"))
						la = Double.parseDouble(t.nextToken());
					else if (firstWord.equals("_cell_length_b"))
						lb = Double.parseDouble(t.nextToken());
					else if (firstWord.equals("_cell_length_c"))
						lc = Double.parseDouble(t.nextToken());
					else if (firstWord.equals("_cell_angle_alpha"))
						a = Double.parseDouble(t.nextToken());
					else if (firstWord.equals("_cell_angle_beta"))
						b = Double.parseDouble(t.nextToken());
					else if (firstWord.equals("_cell_angle_gamma"))
						g = Double.parseDouble(t.nextToken());
					else if (firstWord.equals("_atom_site_occupancy")) {
						readingInAtoms = true;
						latVects = (new Cell(la,lb,lc,a,b,g,null,null)).getLatticeVectors();
					}
				}
			}
		} catch (FileNotFoundException x) {
			System.out.println("parseGulpCif: CIF file not found: " + cifFile.getAbsolutePath());
			return null;
		} catch (IOException y) {
			System.out.println("parseGulpCif: problem reading CIF file(?)" + cifFile.getAbsolutePath());
		}
		
		// maybe GULP died before spitting out a CIF
		if (sites.size() == 0)
			return null;

		Cell answer = new Cell(latVects, sites);
		
		return answer;
	} 
	
	public boolean satisfiesPerSpeciesMIDs(List<Triplet<Element,Element,Double>> mids) {
		for (Triplet<Element,Element,Double> mid : mids) {
			Element a = mid.getFirst();
			Element b = mid.getSecond();
			Double minid = mid.getThird();
			
			for (int i = 0; i < getNumSites(); i++) {
				List<Site> sitesInSphere = getAtomsInSphereSorted(getSite(i).getCoords(), minid);
				boolean sphereHasA = false;
				boolean sphereHasB = false;
				for (Site s : sitesInSphere) {
					if (!sphereHasA && s.getElement().equals(a)) {
						sphereHasA = true;
						continue;
					}
					if (s.getElement().equals(b))
						sphereHasB = true;
				}
				if (sphereHasA && sphereHasB)
					return false;
			}
		}
		return true;
	}

	public boolean satisfiesMinInteratomicDistance(double minid) {
		for (int i = 0; i < getNumSites(); i++) {
			// it's no good if there are any other atoms in the minimum radius sphere
			List<Site> sitesInSphere = getAtomsInSphereSorted(getSite(i).getCoords(), minid);
			if (sitesInSphere.size() > 1) {
				return false;
			}
		}
		return true;
	}
	
	// parses a CIF file (of the format written by openbabel)
	// and returns a structure
	public static Cell parseAvogCif(File cifFile) {
		String line = null;
		double la = 0, lb = 0, lc = 0, a = 0, b = 0, g = 0;

		List<Vect> latVects = null;
		List<Site> sites = new LinkedList<Site>();

		Boolean readingInAtoms = false;

		try {
			BufferedReader r = new BufferedReader(new FileReader(cifFile));
			while ((line = r.readLine()) != null) {
				StringTokenizer t = new StringTokenizer(line);
				// skip empty lines:
				if (!t.hasMoreTokens())
					continue;
				if (readingInAtoms) {
					String symbol = t.nextToken();
					t.nextToken();
					double x = Double.parseDouble(t.nextToken());
					double y = Double.parseDouble(t.nextToken());
					double z = Double.parseDouble(t.nextToken());
					// String occupancy = t.nextToken();
					sites.add(new Site(Element.getElemFromSymbol(symbol), new Vect(x,y,z)));
				} else {
					String firstWord = t.nextToken();
					if (firstWord.equals("_cell_length_a"))
						la = Double.parseDouble(t.nextToken());
					else if (firstWord.equals("_cell_length_b"))
						lb = Double.parseDouble(t.nextToken());
					else if (firstWord.equals("_cell_length_c"))
						lc = Double.parseDouble(t.nextToken());
					else if (firstWord.equals("_cell_angle_alpha"))
						a = Double.parseDouble(t.nextToken());
					else if (firstWord.equals("_cell_angle_beta"))
						b = Double.parseDouble(t.nextToken());
					else if (firstWord.equals("_cell_angle_gamma"))
						g = Double.parseDouble(t.nextToken());
					else if (firstWord.equals("_atom_site_Cartn_z")) {
						readingInAtoms = true;
						latVects = (new Cell(la,lb,lc,a,b,g,null,null)).getLatticeVectors();
					}
				}
			}
		} catch (FileNotFoundException x) {
			System.out.println("parseAvogCif: CIF file not found: " + cifFile.getAbsolutePath());
			return null;
		} catch (IOException y) {
			System.out.println("parseAvogCif: problem reading CIF file(?)" + cifFile.getAbsolutePath());
		}
		
	//	List<Vect> stBasis = new LinkedList<Vect>();
	//	Vect x = new Vect(1.0,0.0,0.0); Vect y = new Vect(0.0,1.0,0.0); Vect z = 
		
		// maybe openbabel died before spitting out a CIF
		if (sites.size() == 0)
			return null;

		Cell answer = new Cell(latVects, sites);
		
		return answer;
		
	} 
	

	public double getHeight() {
		List<Vect>lVects = getLatticeVectors();
		Vect a = lVects.get(0);
		Vect b = lVects.get(1);
		Vect c = lVects.get(2);
		Vect normalToXandY = a.cross(b);
		return normalToXandY.scalarMult(c.dot(normalToXandY)/normalToXandY.dot(normalToXandY)).length();
	}
	
/*
	// Convert from one type of file to another
	// Can be used to standardize CIF files (i.e. CIF to CIF conversion)
	public static void convertCell(String input, String formatIn, String output, String formatOut) {
		// Initialize
		System.loadLibrary("openbabel_java");
		OBConversion conv = new OBConversion();
		OBMol mol = new OBMol();
		
		// Set formats
		conv.SetInFormat(formatIn);
		conv.SetOutFormat(formatOut);
		
		conv.ReadFile(mol, input);
		conv.WriteFile(mol, output);
	}
*/
	
	// just for testing
	public static void main(String args[]) {
/*		Cell a = Cell.parseCif(new File("/home/wtipton/projects/ga_for_crystals/test_runs/lj_cluster_test/garun_lj_cluster/4363.cif"));
	//	Cell c = VaspOut.getPOSCAR("/home/wtipton/cifs/POSCAR_HCP");
		Cell b = Cell.parseCif(new File("/home/wtipton/projects/ga_for_crystals/test_runs/lj_cluster_test/garun_lj_cluster/4387.cif"));

		Cell arot = a.getCellRotatedWithSitesInPrincDirections(a.getSite(1).getCoords(), a.getSite(0).getCoords());
		System.out.println(a.getSite(1).getCoords());
		System.out.println(arot.getSite(1).getCoords());

		System.out.println(a.matchesCellNoPBCs(b, .1));
*/
		//Cell c2 = StructureOrg.parseCif(new File("/home/wtipton/cifs/2.cif"));
		/*Cell a = Cell.parseCif(new File("/home/wtipton/andysexamples/205.cif"));
		Cell b = Cell.parseCif(new File("/home/wtipton/andysexamples/29.cif"));
		Cell c = Cell.parseCif(new File("/home/wtipton/andysexamples/3851.cif"));
	//	System.out.println(b);
//		System.out.println(a);
		System.out.println(a.matchesCell(b, 0.5, 0.5, 0.05));
		System.out.println(b.matchesCell(a, 0.5, 0.5, 0.05));
		System.out.println(a.matchesCell(c, 0.5, 0.5, 0.05));
		System.out.println(c.matchesCell(a, 0.5, 0.5, 0.05));
		System.out.println(c.matchesCell(b, 0.5, 0.5, 0.05));
		System.out.println(b.matchesCell(c, 0.5, 0.5, 0.05));
		
		System.out.println(a.matchesCell(a, 0.1, 0.05, 0.05));
		System.out.println(b.matchesCell(b, 0.1, 0.05, 0.05));
		System.out.println(c.matchesCell(c, 0.1, 0.05, 0.05)); */
		
		
	//	Cell a = VaspOut.getPOSCAR("/Users/benjaminrevard/GA/1075.unrelaxed.POSCAR");
	//	System.out.println(a.satisfiesMinInteratomicDistance(1.075));
	//	System.out.println(a.getHeight());
	/*	Cell a = VaspOut.getPOSCAR("/Users/benjaminrevard/GA/materials/POSCAR");
		System.out.println(a.toString());
		System.out.println(a.getHeight()); */
		List<Triplet<Element,Element,Double>> perSpeciesMIDs = new ArrayList<Triplet<Element,Element,Double>>();
		perSpeciesMIDs.add(new Triplet<Element,Element,Double>(Element.getElemFromSymbol("Al"),
				Element.getElemFromSymbol("Al"),0.02));
		Cell c = VaspOut.getPOSCAR("/home/wtipton/POSCAR");
		System.out.println(c.satisfiesPerSpeciesMIDs(perSpeciesMIDs));

	//	a.getNigliReducedCell();
	/*	
		for (int i = 0; i < a.getNumSites(); i++) {
			// it's no good if there are any other atoms in the minimum radius sphere
			System.out.println(a.getSite(i).toString());
			List<Site> sitesList = a.getAtomsInSphereSorted(a.getSite(i).getCoords(), .1);
			if (sitesList.size() > 1) {
				System.out.println("Organism failed minimum interatomic distance constraint: " + sitesList.size());
			}
		} */

		
//		a.writeCIF("test");
//		System.out.println(Cell.parseCell("test", "cif").toString());
	//	c.writeCIF("/home/wtipton/cifs/1.cif");
	//	c.writeCIF("/home/wtipton/cifs/2.cif");
		
		/*
		Generation g = new Structures();
		GAParameters.getParams().setMaxNumAtoms(10);
		GAParameters.getParams().setMinNumAtoms(1);
		GAParameters.getParams().setMaxLatticeAngle(140);
		GAParameters.getParams().setMinLatticeAngle(40);
		GAParameters.getParams().setMaxLatticeLength(15);
		GAParameters.getParams().setMinLatticeLength(2);
		List<Composition> comps = new LinkedList<Composition>();
		Element els[] = {Element.getElemFromSymbol("Pt")};
		Double cs[] = {1.0};
		comps.add(new Composition(els,cs));
		GAParameters.getParams().setCompSpace(new CompositionSpace(comps));
		String arg[] = {"randomVol"};
		RandomSOCreator soc = new RandomSOCreator(arg);
		
	//	System.out.println(c.getNigliReducedCell());
	//	System.out.println(c.getCellWithAllAtomsInCell());
		for (int i = 0; i < 10000; i++) {
			StructureOrg s = soc.makeOrganism(g);
			//System.out.println(c.getNigliReducedCell());
			s.getCell().getNigliReducedCell();
			s.getCell().getCellWithSiteIShiftedToOrigin(0).getNigliReducedCell();
		}
		*/
		
/*		for (Site s : c.getAtomsInSphereSorted(new Vect(2.0,0.0,0.0), 3))
			System.out.println(s);
		System.out.println("------");
		for (Site s : c.getAtomsInSphereSortedSlow(new Vect(2.0,0.0,0.0), 3))
			System.out.println(s);
		*/
	//	System.out.println(c1.matchesCell(c2, 0.1, 0.2, 5*Math.PI/180));
		
		/*
		List<Site> sites = new LinkedList<Site>();
		sites.add(new Site(Element.getElemFromSymbol("C"), new Vect(0.0, 0.0, 0.0)));
		
		Cell c = new Cell(5.0,5.0,5.0,90*Math.PI/180,20*Math.PI/180,90*Math.PI/180, sites, "testcell");
		List<Vect> vects = new LinkedList<Vect>();
		vects.add(new Vect(new double[]{11.287795,   0.000000,   0.000000}));
    	vects.add(new Vect(new double[]{11.409892 , 10.459132 ,  0.000000}));
    	vects.add(new Vect(new double[]{-1.502377  ,-2.150436 , 13.552370}));
		Cell bob = new Cell(vects, sites);
		
    	System.out.println(bob);
    	System.out.println("---------");    	
    	System.out.println(bob.getNigliReducedCell());
    	*/
	
	/*	System.out.println(c);
		Vect center = new Vect(0.0,0.0,5.0);
		Site[] nns = c.getAtomsInSphereSorted(center, 2);
		
		c.writeCIF("/home/wtipton/bob.cif");
		
		for (Site s : nns) {
			System.out.println("-------------------");
			System.out.println(s.getCoords().toString());
			System.out.println(s.getCoords().toString(c.getLatticeVectors()));
			System.out.println(s.getCoords().getCartDistanceTo(center));
		}
		*/
		
	}

	
}