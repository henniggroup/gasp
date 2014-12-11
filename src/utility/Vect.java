/*
 * Copyright 2011-2014 Will Tipton, Richard Hennig, Ben Revard, Stewart Wenner

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
package utility;

import java.io.Serializable;
import java.util.*;

import vasp.VaspIn;
import vasp.VaspOut;
import crystallography.Cell;
import crystallography.Site;
import Jama.*;
/*  Represents a vector.
 *  Immutable!
 * 
 *  Represent everything in cartesian coords under the covers.
 * 
 *  TODO:
 *  - maybe refactor this stuff to use Jama more extensively under the covers
 */

public class Vect implements Serializable {
	
	static final long serialVersionUID = 1l;
	
	/* fractional coordinates */
	private List<Double> frac;
	/* cartesian coordinates */
	private List<Double> cart;
	/* basis given in cartesian coordinates */
	private List<Vect> basis;
	
	public Vect(List<Double> cartComps) {
		this(cartComps, null);
	}
	
	public Vect (List<Double> cartComps, List<Vect> _basis) {
		
		if (_basis != null) {
			basis = _basis;
			frac = cartComps;
			setCartCoordToMatchFracCoord();
		} else
			cart = new ArrayList<Double>(cartComps);
	}
	public Vect(Double _u, Double _v, Double _w) {
		this(_u,_v,_w,null);	
	}
	
	public Vect(Double _u, Double _v, Double _w, List<Vect> _basis) {
		if (_basis != null) {
			basis = _basis;
			frac = new ArrayList<Double>();
			frac.add(_u); frac.add(_v); frac.add(_w);
			// A null basis indicates that u,v,and w are cartesian coords,
			// i.e. the basis is the standard unit vectors 
			// Make sure our basis has three vectors 
			if (_basis.size() != 3)
				throw new IllegalArgumentException("ThreeVector given basis with other than 3 vectors.");
			setCartCoordToMatchFracCoord();
		} else {
			cart = new ArrayList<Double>();
			cart.add(_u); cart.add(_v); cart.add(_w);
		}
	}
	
	public Vect(double[] ds) {
		this(ds[0], ds[1], ds[2]);
	}
	
	public Vect leftMultByMatrix(Matrix m) {
		if (m.getColumnDimension() != this.getDimension())
			throw new RuntimeException("dims dont match in Vect.leftMultByMatrix");
		
		return new Vect(m.times(new Matrix(this.getPackedCartComps())).getColumnPackedCopy());
	}
	
	public double[][] getPackedCartComps() {
		List<Double> ccomps = getCartesianComponents();
		double result[][] = new double[ccomps.size()][1];
		for (int i = 0; i < ccomps.size(); i++)
			result[i][0] = ccomps.get(i);
		return result;
	}

	public List<Double> getCartesianComponents() {
		return new ArrayList<Double>(cart);
	}
	
	private void setCartCoordToMatchFracCoord() {
		if (basis.isEmpty() || frac.size() != basis.get(0).getDimension())
			throw new RuntimeException("Vect.getCartCompsFromNoncartComps");
		
		int dim = basis.get(0).getDimension();
		List<Double> result = new ArrayList<Double>();
		for (int i = 0; i < dim; i++){
			double entry = 0;
			for (int j = 0; j < dim; j++)
				entry += frac.get(j) * basis.get(j).getCartesianComponents().get(i);
			result.add(entry);
		}
		cart = result;
	}
	
	public Vect plusWRT(Vect shift, List<Vect> basis) {
		Vect result = this.plus(shift);
		
		List<Double> fractComps = result.getComponentsWRTBasis(basis);
		for (int i = 0; i < Constants.numDimensions; i++) {
			double numResid = Math.floor(fractComps.get(i));
			result = result.plus(basis.get(i).scalarMult(-numResid));
		}

		return result;
	}
	
	
	// returns this x other
	public Vect cross(Vect other) {
		List<Double> comps1 = this.getCartesianComponents();
		double a1 = comps1.get(0);
		double a2 = comps1.get(1);
		double a3 = comps1.get(2);
		List<Double> comps2 = other.getCartesianComponents();
		double b1 = comps2.get(0);
		double b2 = comps2.get(1);
		double b3 = comps2.get(2);
		
		return new Vect(a2*b3 - a3*b2, a3*b1 - a1*b3, a1*b2 - a2*b1);
	}
	
	public Vect plus(Vect v) {
		List<Double> vComps = v.getCartesianComponents();
		List<Double> ourComps = getCartesianComponents();
		List<Double> newComps = new ArrayList<Double>();
		
		if (vComps.size() != ourComps.size())
			throw new RuntimeException("Vector: trying to add vectors of different lengths");
		
		for (int i = 0; i < vComps.size(); i++)
			newComps.add(vComps.get(i) + ourComps.get(i));

		return new Vect(newComps);
	}
	
	public Vect subtract(Vect v) {
		List<Double> vComps = v.getCartesianComponents();
		List<Double> ourComps = getCartesianComponents();
		List<Double> newComps = new ArrayList<Double>();
		
		if (vComps.size() != ourComps.size())
			throw new RuntimeException("Vector: trying to subtract vectors of different lengths");
		
		for (int i = 0; i < vComps.size(); i++)
			newComps.add(vComps.get(i) - ourComps.get(i));
		
		return new Vect(newComps);
	}
	
	public Vect scalarMult(Double s) {
		List<Double> newV = new ArrayList<Double>();
		for (int i = 0; i < getDimension(); i++)
			newV.add(s * cart.get(i));

		return new Vect(newV);
	}
	
	private static Vect xhat = null;
	public static Vect xHat() {
		if (xhat == null)
			xhat = new Vect(1.0,0.0,0.0);
		return xhat;
	}
	
	private static Vect yhat = null;
	public static Vect yHat() {
		if (yhat == null)
			yhat = new Vect(0.0,1.0,0.0);
		return yhat;
	}
	
	private static Vect zhat = null;
	public static Vect zHat() {
		if (zhat == null)
			zhat = new Vect(0.0,0.0,1.0);
		return zhat;
	}
	
	public Vect changeBasis(List<Vect> basis) {
		return new Vect(this.frac, basis);
	}

	/* Find the coordinates of our ThreeVector with respect to an alternate basis.
	 * We have:
	 *           X = B * A
	 * where X is our ThreeVector's cartesian coords, B is the new basis, and 
	 * A is what we want.  So, we just find
	 *           A = B^-1 * X .
	 * TODO:
	 *  - see if running Jama's least-squares routine actually gives better results
	 *    than solving the system "manually".
	 */
	public List<Double> getComponentsWRTBasis(List<Vect> basis) {
		int dim = getDimension();
		
		if (basis.size() != dim)
			throw new RuntimeException("Vector: trying get components wrt basis of wrong size.");

		List<Double> cartComps = getCartesianComponents();
		Matrix X = new Matrix(dim,1);
		for (int i = 0; i < dim; i++)
			X.set(i, 0, cartComps.get(i));

		Matrix B = new Matrix(dim,dim);
		for (int i = 0; i < dim; i++)
			for (int j = 0; j < dim; j++)
				B.set(i, j, basis.get(j).getCartesianComponents().get(i));
		
		if (Math.abs(B.det()) < Constants.epsilon)
			throw new IllegalArgumentException("ThreeVector.getComponentsWRTBasis given degenerate basis.");
		
		Matrix A = B.inverse().times(X);
		
		List<Double> result = new ArrayList<Double>();
		for (int i = 0; i < dim; i++)
			result.add(A.get(i, 0));
		
		return result;
	}
	
	public double getCartDistanceTo(Vect other) {
		List<Double> otherCCs = other.getCartesianComponents();
		List<Double> thisCCs = getCartesianComponents();
		
		if (otherCCs.size() != thisCCs.size())
			throw new RuntimeException("Vector: trying to get distance between vectors of different lengths.");
		
		double sum = 0.0;
		for (int i = 0; i < getDimension(); i++)
			sum += (otherCCs.get(i)-thisCCs.get(i))*(otherCCs.get(i)-thisCCs.get(i));
		
		return Math.sqrt(sum);
	}
	
	public String toString() {
		StringBuilder result = new StringBuilder();
		List<Double> cartCoords = getCartesianComponents();
		result.append(getDimension());
		result.append("-Vector: ( ");
		for (int i = 0; i < getDimension(); i++) {
			result.append(cartCoords.get(i));
			result.append(" ");
		}
		result.append(")");
		
		return result.toString();
	}
	
	public String toString(List<Vect> lattice) {
		StringBuilder result = new StringBuilder();
		List<Double> coords = this.getComponentsWRTBasis(lattice);
		result.append(getDimension());
		result.append("-Vector: ( ");
		for (int i = 0; i < getDimension(); i++) {
			result.append(coords.get(i));
			result.append(" ");
		}
		result.append(")");
		
		return result.toString();
	}
	
	public Double angleToInDegrees(Vect another) {
		return Math.acos((this.dot(another))/(this.length()*another.length()))
				* 180.0 / Math.PI;
	}
	
	public Vect getVectShiftedIntoBasis(List<Vect> b) {
		List<Double> fracComps = this.getComponentsWRTBasis(b);
		
		List<Double> newFracComps = new ArrayList<Double>();
		
		for (Double c : fracComps)
			newFracComps.add(Utility.getPositiveFracPart(c));
		
		return new Vect(newFracComps, b);
	}
	
	public Double dot(Vect v) {
		List<Double> vComps = v.getCartesianComponents();
		List<Double> ourComps = getCartesianComponents();
		
		if (vComps.size() != ourComps.size())
			throw new RuntimeException("Vector: trying to dot vectors of different lengths");
		
		Double result = 0.0;
		for (int i = 0; i < getDimension(); i++)
			result += vComps.get(i)*ourComps.get(i);
		return result;
	}
	
	public int getDimension() {
		return cart.size();
	}
	
	public Double length() {
		List<Double> cartComps = getCartesianComponents();
		Double resultSquared = 0.0;
		for (Double c : cartComps)
			resultSquared += c*c;
		return Math.sqrt(resultSquared);
	}
	
	// Magnitude of projection of w onto given vector
	public Double magnitudeOfProjection(Vect w) {
		double dotprod = Math.abs(this.dot(w));
		double l = this.length();
		
		return dotprod/l;
	}
	
	public Vect getNormalized() {
		List<Double> cartComps = this.getCartesianComponents();
		Double l = this.length();
		for (int i=0; i<Constants.numDimensions; i++)
			cartComps.set(i, cartComps.get(i)/l);
		
		return new Vect(cartComps.get(0), cartComps.get(1), cartComps.get(2));
	}
	
	public Vect getSumNormalized(Double n) {
		List<Double> cartComps = getCartesianComponents();
		double sum = 0;
		for (Double d : cartComps)
			sum += d;
		return this.scalarMult(n / sum);
	} 
	
	public static Vect getNullVect(int numDims) {
		List<Double> comps = new ArrayList<Double>();
		for (int i = 0; i < numDims; i++)
			comps.add(0.0);
		return new Vect(comps);
	}
	
	public static boolean pointsAreCollinear(Vect p1, Vect p2, Vect p3, double epsilon) {
		Vect v1 = p2.subtract(p1);
		Vect v2 = p3.subtract(p1);
		return v1.cross(v2).length() < epsilon;
	}
	
	/*
	 * For the sake of testing the ThreeVector class
	 */
	public static void main(String[] args) {
//		Vect bob = new Vect(1.0, 1.0, 1.0);
//		
//		Vect e1 = new Vect(3.0, 0.0, 0.0);
//		Vect e2 = new Vect(0.0, 3.0, 0.0);
//		Vect e3 = new Vect(0.0, 0.0, 3.0);
//		List<Vect> basis = new ArrayList<Vect>();
//		basis.add(e1);
//		basis.add(e2);
//		basis.add(e3);
//		
//		Vect shift = new Vect(-2.0,-1.0,-2.0);
//		
//		bob = bob.plusWRT(shift, basis);
//		
//		System.out.println(bob);
//		
//		List<Double> components = bob.getComponentsWRTBasis(basis);
//		System.out.println(components.get(0) + " " + components.get(1) + " " + components.get(2));
//		
//		System.out.println(e1.cross(e2));
//		System.out.println(Vect.pointsAreCollinear(e1, e2, e3, 0.0001));
		
		
		//Cell substrate = VaspOut.getPOSCAR("/home/ay256/Desktop/23.POSCAR");
//		List<Vect> basis = new ArrayList<Vect>();
//		basis.add(new Vect (1.0, 0.0, 0.0));
//		basis.add(new Vect (0.0, 1.0, 0.0));
//		basis.add(new Vect (0.0, 0.0, 1.0));
//		System.out.println("Making Vect v");
//		Vect v = new Vect(0.5, 0.5, 0.5, basis);
//		System.out.println(v.toString());
//		List<Vect> newBasis = new ArrayList<Vect>();
//		newBasis.add(new Vect (2.0, 0.0, 0.0));
//		newBasis.add(new Vect (0.0, 2.0, 0.0));
//		newBasis.add(new Vect (0.0, 0.0, 2.0));
//		v.changeBasis(newBasis);
//		assert v.basis == newBasis;
//		for (Double d : v.frac)
//			System.out.println("fractional coord" + d);
//		System.out.println("After changing basis: " + v.toString());
		
		Cell oldCell = VaspOut.getPOSCAR("/home/ay256/Desktop/2.vasp");
		Cell substrate = VaspOut.getPOSCAR("/home/ay256/Desktop/23.POSCAR");
		List<Vect> newBasis = new ArrayList<Vect>();
		newBasis.add(substrate.getLatticeVectors().get(0));
		newBasis.add(substrate.getLatticeVectors().get(1));
		newBasis.add(oldCell.getLatticeVectors().get(2));
		List<Site> newSites = new ArrayList<Site>();
		for (Site site : oldCell.getSites()) {
			Vect v = site.getCoords().changeBasis(newBasis);
			newSites.add(new Site(site.getElement(), v));
		}
		Cell newCell = new Cell(newBasis, newSites, oldCell.getLabel());
		VaspIn.writePoscar(newCell, "/home/ay256/Desktop/newCell.vasp", true);
		VaspIn.writePoscar(oldCell, "/home/ay256/Desktop/oldCell.vasp", true);
		
	}




}
