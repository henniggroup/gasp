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

package pdvisual;

import java.util.*;

import utility.Utility;
import crystallography.*;
import chemistry.*;

import Jama.*;


public class PseudoPDBuilder {
	
	PDData pdd;
	
	static double EPSILON = 1E-4;
	
	public PseudoPDBuilder(PDData _pdd) {
		pdd = _pdd;
	}
	
	// this is where we build the pseudopddata.
	// takes a list of points in the coordinate-space that specify the hyperplane
	// of the crosssection corresponding to our pseudopd
	public PseudoPDData getPseudoPDData(List<Composition> subspace) {
		
		// Convert the problem to an algebraic one: want to find intersection of
		// pdd's phase space (#elems dimensions) and the hyperplane in that
		// space defined by the reaction coordinates.
		
		
		// solve for intersections along the reaction coords
		List<Composition> comps = getIntersections(subspace);
		// make sure the list of interesting compositions includes the "endpoints"
		comps.addAll(subspace);
		
		// get the mixedphases associated w/ all the vertices in the problem
		// from the original PDData
		PDAnalyzer pda = new PDAnalyzer(pdd);
			
		// and the list of MixedPhases from the compositions of intersection
		List<MixedPhase> phases = new LinkedList<MixedPhase>();
		for (Composition c : comps) {
			List<Integer> f = pda.getFacet(c);
			// ick
			if (f == null) {
				System.out.println("WARNING: could not find facet for composition" + c);
				phases.add(new MixedPhase(c,pdd.getBasisEntries()));
			}
			else
				phases.add(new MixedPhase(c,pdd.getEntriesOnFacet(pdd.getIndxFacets().indexOf(f))));
		
		}
		
		// make the list of basis indices
		List<Integer> basisIndxs = new LinkedList<Integer>();
		for (Composition c : subspace)
			basisIndxs.add(comps.indexOf(c));
		
		/*
		 * This is a quick add to fix the bug where some phases on the edge of the
		 * phase diagrams weren't represented in the pseudo
		 * TODO: come with something better that would fix the problem in the intersection
		 * contruction above...
		 * Geoffroy
		 */
		
		//for any stable phase on the convex hull we check if it's in the subspace
		//if yes, we check if there is already a mixed phase containing this pure phase
		// 			if yes, we don't add it of course
		//			otherwise, we do add it...
		
		List<IComputedEntry> list=pdd.getStableEntries();
		for(IComputedEntry e:list)
		{
			Matrix a;
			try {
				 a = doDecomposition(pdd.getElements(), subspace, e.getComposition());
			} catch (RuntimeException x) {
				// if there was no good solution, move on to the next one:
				//  - in this case, there was no solution to the linear system
				continue;
			}
			//  - in this case, it was outside the subdomain of interest
			boolean outside = false;
			for(Double d : a.getRowPackedCopy())
				if (d < - EPSILON)
					outside = true;
			if (outside)
				continue;
		
			// we have an entry in the subspace: check if already present in phases
			boolean already=false;
			for(MixedPhase m:phases)
			{
				if(m.getPhases().containsKey(e))
				{
					if(m.getPhases().get(e)>1-EPSILON)
						already=true;
				}
			}
			if(!already)
			{
				// we have a really new entry add it to the pseudopd
				Map<IComputedEntry, Double> eMap = new HashMap<IComputedEntry, Double>();
				eMap.put(e, 1.0);
				phases.add(new MixedPhase(eMap));
			}
		}
		
		
		
		// now add MixedPhases corresponding to known structures (IComputedEntries) by checking
		// each of the unstable IComputedEntries in the original phase diagram to see
		// if it lies in the subspace
		List<Integer> indxUnstableEntries = pdd.getIndxUnstableEntries();
		for (Integer indx : indxUnstableEntries) {
			IComputedEntry e = pdd.getEntry(indx);
		
			Matrix a;
			try {
				 a = doDecomposition(pdd.getElements(), subspace, e.getComposition());
			} catch (RuntimeException x) {
				// if there was no good solution, move on to the next one:
				//  - in this case, there was no solution to the linear system
				continue;
			}
			//  - in this case, it was outside the subdomain of interest
			boolean outside = false;
			for(Double d : a.getRowPackedCopy())
				if (d < - EPSILON)
					outside = true;
			if (outside)
				continue;
		
			// we have an entry in the subspace: add it to the pseudopd
			Map<IComputedEntry, Double> eMap = new HashMap<IComputedEntry, Double>();
			eMap.put(e, 1.0);
			phases.add(new MixedPhase(eMap));
		}
		
		
		PseudoPDData newPseudoPDData = new PseudoPDData(phases, basisIndxs, pdd);
		// find the convex hull
		if (comps.size() > 2) {
			PDConvexable convexableData = new PDConvexable(newPseudoPDData);
			Convexhull chull = new Convexhull(convexableData, true);
			List<List<Integer>> facets = chull.getFaLists(); 
			
			/**
			 * one last check -- remove facets where the formation energy of
			 * one of the points on the facet is positive (i.e., it is part of the
			 * convex hull, but in the Eform > 0 space)
			 */
			Iterator<List<Integer>> iFacet = facets.iterator();
			List<Double> eforms = newPseudoPDData.getFormEnergiesPerAtom(); 
			while(iFacet.hasNext()){
			    List<Integer> facet = iFacet.next();
			    boolean removeFacet = false;
			    for(Integer ivert : facet){
			        if(eforms.get(ivert) > PDBuilder.EFORM_TOL){
			            removeFacet = true;
			            break; 
			        }
			    }
			    if(removeFacet){
	                iFacet.remove();		        
			    }
			}
			
			newPseudoPDData.setFacets(facets);
		} else {
			newPseudoPDData.setFacets(new LinkedList<List<Integer>>());
		}

		// put it all together
		return newPseudoPDData;
	}
	
	public static Matrix doDecomposition(List<Element> units, List<Composition> basis, Composition c) {
		// to see if e is in the subspace, we're going to try to solve a linear system.
		// make the composition matrix, Cn
		double Cn[][] = new double[units.size()][basis.size()];
		for (int i = 0; i < units.size(); i++)
			for (int j = 0 ; j < basis.size(); j++)
				Cn[i][j] = basis.get(j).getFractionalCompo(units.get(i));
		Matrix MCn = new Matrix(Cn);
		// make a matrix C of the input composition (e.getComposition())
		double C[][] = new double[units.size()][1];
		for (int i = 0; i < units.size(); i++)
			C[i][0] = c.getFractionalCompo(units.get(i));
		Matrix MC = new Matrix(C);
		// try to solve it
		Matrix a = MCn.solve(MC);
		// check that it was a good solution
		double resid = MCn.times(a).minus(MC).normInf(); //picks off maximum row
		if(resid > EPSILON){
		    throw new RuntimeException("No good decomposition possible");
		}
		return a;
	}
	
	private List<Composition> getIntersections(List<Composition> subspace) {	
		List<Composition> result = new LinkedList<Composition>();
		
		// for each pair of ADJACENT facets
		List<List<Integer>> adjacentFacets = pdd.getAdjacencyList();
		for(List<Integer> pair : adjacentFacets) {
			int indxFacet1 = pair.get(0);
			int indxFacet2 = pair.get(1);
			// solve the problem to get a point 
			Composition min = getIntersection(indxFacet1, indxFacet2, subspace, false);
			if (min != null) {
				// check for uniqueness
				boolean skipit = false;
				for (Composition c : result)
					if (c.equals(min)) 
						skipit = true;
				if (skipit)
					continue;
				// check that we're in the right subdomain
				//  - for the constraints
				if (!linearComboIsPositive(pdd.getElements(), subspace, min))
					continue;
				//  - for the first facet
				List<Composition> facet1subspace = new LinkedList<Composition>();
				for (Integer i : pdd.getIndxFacets().get(indxFacet1))
					facet1subspace.add(pdd.getEntry(i).getComposition());
				if (!linearComboIsPositive(pdd.getElements(), facet1subspace, min))
					continue;
				//  - for the second facet
				List<Composition> facet2subspace = new LinkedList<Composition>();
				for (Integer i : pdd.getIndxFacets().get(indxFacet2))
					facet2subspace.add(pdd.getEntry(i).getComposition());
				if (!linearComboIsPositive(pdd.getElements(), facet2subspace, min))
					continue;	
				// ok, we're good. add it
				result.add(min);
			}	
		}
		
		return result;
	}
	
	// skip compositions that are outside the desired range.  more precisely,
	// we currently have a list of compositions along an entire slice of the
	// original phase diagram, whereas we only want those between the
	// two (or more) selected points.  if it's outside the given subdomain,
	// we won't be able to make a linear combination of the bases with nonnegative
	// coefficients.
	// and we wanna do this check for the subspace defined by the cut
	// and by both facets. 
	private boolean linearComboIsPositive(List<Element> units, List<Composition> basis, Composition c) {
		Matrix a;
		try {
			 a = doDecomposition(units, basis, c);
		} catch (RuntimeException x) {
			// we should never get here, all of these should be in the subspace
			throw x;
		}
		//  - in this case, it was outside the subdomain of interest
		boolean allPositiveCoefs = true;
		for(Double d : a.getRowPackedCopy())
			if (d < - EPSILON)
				allPositiveCoefs = false;
		
		return allPositiveCoefs;
	}
	
	// set up and solve an LP problem to find the point along the intersection
	// of two facets with the highest or lowest energy.  return the corresponding
	// Composition
	private Composition getIntersection(int indxFacet1, int indxFacet2, List<Composition> subspace, boolean maximize) {
		// our linear programming problem has several components:
		//	structural variables(columns) - one for each constituent ELEMENT in our system
		//	auxillary variables(rows) - each equal 1
		//	optimization function: energy
		//	other constraints: compositions must sum to 1 and each lie between 0 and 1
		
		  int dim = pdd.getDimension();
		  // structural vars (constituent fractions and energy)
		  int numCols = dim + 1;
		  // get elements in an array
		  Element elements[] = new Element[dim];
		  elements = pdd.getElements().toArray(elements);
	      
	      // constraints come from two facets and the cut
	      // facet 1
	      List<Double[]> facet1Points = new LinkedList<Double[]>();
	      List<Integer> facet1 = pdd.getIndxFacets().get(indxFacet1);
	      for (Integer f : facet1) {
	    	  Composition c = pdd.getEntry(f).getComposition();
	    	  Double newPoint[] = new Double[numCols];
	    	  for (int i = 0; i < dim; i++)
	    		  newPoint[i] = c.getFractionalCompo(elements[i]);
	    	  newPoint[numCols - 1] = pdd.getEntry(f).getTotalEnergy();
	    	  facet1Points.add(newPoint);
	      }
	      double[][] facet1Constraints = getRowsFromPoints(facet1Points);     
	      // facet 2
	      List<Double[]> facet2Points = new LinkedList<Double[]>();
	      List<Integer> facet2 = pdd.getIndxFacets().get(indxFacet2);
	      for (Integer f : facet2) {
	    	  Composition c = pdd.getEntry(f).getComposition();
	    	  Double newPoint[] = new Double[numCols];
	    	  for (int i = 0; i < dim; i++)
	    		  newPoint[i] = c.getFractionalCompo(elements[i]);
	    	  newPoint[numCols - 1] = pdd.getEntry(f).getTotalEnergy();
	    	  facet2Points.add(newPoint);
	      }
	      double[][] facet2Constraints = getRowsFromPoints(facet2Points); 
	      // desired subspace 
	      // TODO: no reason to do this computation over and over. move it 
	      List<Double[]> subspacePoints = new LinkedList<Double[]>();
	      for (Composition c : subspace) {
	    	  Double newPoint[] = new Double[numCols];
	    	  for (int i = 0; i < dim; i++)
	    		  newPoint[i] = c.getFractionalCompo(elements[i]);
	    	  newPoint[numCols - 1] = 0.0; // energy coord doesnt matter
	    	  subspacePoints.add(newPoint);
	      }
	      double[][] subspaceConstraints = getRowsFromPoints(subspacePoints); 
	      // get rid of all dependence on the energy coordinate
	      for (int i = 0; i < subspaceConstraints.length; i++)
	    	  subspaceConstraints[i][numCols - 1] = 0;
	      
	      // auxillary vars (2 facets and 1 rxncoord and 1 physical constraints)
	      int numRows = facet1Constraints.length
	      				+ facet2Constraints.length
	      				+ subspaceConstraints.length
	      				+ 1;
	      Matrix RHS = new Matrix(numRows,1);
	      
	      // fill in the matrix
	      int counter = 1;
	      int ia[] = new int[numCols*numRows + 1]; //
	      int ja[] = new int[numCols*numRows + 1];
	      double ar[] = new double[numCols*numRows + 1];
	      // fill in from facet 1:
	      int rowsOffset = 0;
	      for (int i = 0; i < facet1Constraints.length; i++) {
	    	  double currentRow[] = facet1Constraints[i];
	    	  int rowNum = i + 1 + rowsOffset;
	    	  for (int j = 0; j < currentRow.length; j++) {
	    		  ia[counter] = rowNum;
	    		  ja[counter] = j + 1;
	    		  ar[counter] = currentRow[j];
	    		  if (ar[counter] != 0) { // don't count zero-valued entries
	    			  counter++;
	    		  }
    	    	  // first row should equal 1; rest equal 0
    	    	  if (i == 0)
    	    		  RHS.set(rowNum - 1, 0, 1);
    	    	  else
    	    		  RHS.set(rowNum - 1, 0, 0);
	    	  }

	      }
	      // fill in from facet 2:
	      rowsOffset += facet1Constraints.length;
	      for (int i = 0; i < facet2Constraints.length; i++) {
	    	  double currentRow[] = facet2Constraints[i];
	    	  int rowNum = i + 1 + rowsOffset;
	    	  for (int j = 0; j < currentRow.length; j++) {
	    		  ia[counter] = rowNum;
	    		  ja[counter] = j + 1;
	    		  ar[counter] = currentRow[j];
	    		  if (ar[counter] != 0) { // don't count zero-valued entries
	    			  counter++;
	    		  }
	    	  }
	    	  // first row should equal 1; rest equal 0
	    	  if (i == 0)
	    		  RHS.set(rowNum - 1, 0, 1);
	    	  else
	    		  RHS.set(rowNum - 1, 0, 0);
	      }
	      // fill in from subspace:
	      rowsOffset += facet2Constraints.length;
	      for (int i = 0; i < subspaceConstraints.length; i++) {
	    	  double currentRow[] = subspaceConstraints[i];
	    	  int rowNum = i + 1 + rowsOffset;
	    	  for (int j = 0; j < currentRow.length; j++) {
	    		  ia[counter] = rowNum;
	    		  ja[counter] = j + 1;
	    		  ar[counter] = currentRow[j];
	    		  if (ar[counter] != 0) { // don't count zero-valued entries
	    			  counter++;
	    		  }
	    	  }
	    	  // first row should equal 1; rest equal 0
	    	  if (i == 0)
	    		  RHS.set(rowNum - 1, 0, 1);
	    	  else
	    		  RHS.set(rowNum - 1, 0, 0);
	      }

	      // physical constraint
	      rowsOffset += subspaceConstraints.length;
	      int rowNum = 1 + rowsOffset;
	      for (int j = 1; j <= dim; j++) {
	    	  ia[counter] = rowNum;
	    	  ja[counter] = j;
	    	  ar[counter] = 1;
	    	  counter++;
	      }
    	  RHS.set(rowNum - 1, 0, 1);
	      
	      // debugging: let's take a look at the matrix:
	      double bob[][] = new double[numRows][numCols];
	      for (int i = 1; i <= counter - 1; i++)
	    	  bob[ia[i]-1][ja[i]-1] = ar[i];
	      Matrix C = new Matrix(bob);
	      
	      Matrix CTC = C.transpose().times(C);
	      if (Math.abs(CTC.det()) < Math.pow(EPSILON,CTC.getColumnDimension()))
	    	  return null;

	      Matrix CInverse;
	//      try {
	    	  CInverse = CTC.inverse().times(C.transpose());
	  //    } catch (RuntimeException x) {
	    	  // no solution!
	    //	  return null;
	     // }
	      
	      Matrix MResult = CInverse.times(RHS);
	      
	      Double amounts[] = new Double[pdd.getDimension()];
	      for (int i = 0; i < pdd.getDimension(); i++) {
	    	  amounts[i] = MResult.get(i, 0);
	    	  // remove unphysical solutions
	    	  if (amounts[i] <  -EPSILON || amounts[i] > 1 + EPSILON) 
	    		  return null;
	      }
	      
		    
		  return new Composition(elements, amounts);

	}
	
	// returns a matrix, each row of which hold coefficients of a constraint.
	// RHS of the first constraint is 1; rest are 0;
	private static double[][] getRowsFromPoints(List<Double[]> origPoints) {
		
		boolean successful = false;
		double[][] result = null;
		List<List<Double[]>> pointsPowerSet = null;;
		
		List<Double[]> points = new LinkedList<Double[]>(origPoints);
		
		while (successful == false) {
			try {
				// make F*:
				Double[] f1 = points.get(0);
				// d = numrows, m = numcols = num constraints
				int d = f1.length;
				int m = points.size();
				Matrix FStar = new Matrix(d, m);
				// fill in column 1
				for (int i = 0; i < d; i++)
					FStar.set(i, 0, f1[i]);
				// fill in the rest
				for (int i = 0; i < d; i++)
					for (int j = 1; j < m; j++)
						FStar.set(i, j, points.get(j)[i]-f1[i]);
		
				// get the SVD.U
				Matrix U = FStar.svd().getU();
				// normalize the columns so that they each sum to 1
				for (int j = 0; j < m; j++) {
					double sum = 0;
					for (int i = 0; i < d; i++)
						sum += U.get(i, j);
					if (sum != 1.0)
						for (int i = 0; i < d; i++) {
							double value = U.get(i, j);
							U.set(i, j, value/sum);
						}
				}
				
				// make a new matrix, FStar extended to be invertible:
				Matrix FStarI = new Matrix(d,d);
				FStarI.setMatrix(0, d-1, 0, m-1, FStar);
				FStarI.setMatrix(0, d-1, m, d-1, U.getMatrix(0, d-1, m, d-1));
				
				// take an inverse
				double[][] R = FStarI.inverse().getArray();
				
				// result is the inverse without rows 2 through (num input points)
				result = new double[d-m+1][];
				result[0] = R[0];
				for (int i = 1; i < result.length; i++)
					result[i] = R[m-1+i];
				
				return result;
			} catch(RuntimeException x) {
				// so if we get here, the input vectors were linearly dependent.
				// we need to remove one or more of them until it's not		
				if (pointsPowerSet == null) {
					pointsPowerSet = Utility.getPowerSet(origPoints);
				}
				
				if (pointsPowerSet.isEmpty())
					throw x; // we should never get here :)
				
				points = pointsPowerSet.get(0);
				pointsPowerSet.remove(0);
			}
		}
		return result;
	}
	
/*	private double[] getEnergyCoefs(List<Integer> origFacet, Element[] basis) throws RuntimeException {
		
		boolean successful = false;
		double[] result = new double[basis.length];
		
		List<Integer> facet = new LinkedList<Integer>(origFacet);
		
		// might have to repeat this calculation with a subset of the given basis if we
		// find that it is linearly dependent.
		List<List<Integer>> facetPowerSet = null;
		while (successful == false) {
			try {
				// get \vec{E}
			    double basisEnergies[][] = new double[facet.size()][1];
			    for (int i = 0; i < facet.size(); i++)
			    	basisEnergies[i][0] = pdd.getEntry(facet.get(i)).getEnergy();
			    Matrix ME = new Matrix(basisEnergies);
				
			    // get F
				double F[][] = new double[basis.length][facet.size()];
				// fill in F
				for (int i = 0; i < basis.length; i++)
					for (int j = 0 ; j < facet.size(); j++)
						F[i][j] = pdd.getEntry(facet.get(j)).getComposition().getNormAmount(basis[i]);
				// make it a matrix
				Matrix MF = new Matrix(F);
				
				// do it
				Matrix bob = MF.transpose().times(MF).inverse().times(MF.transpose());
				Matrix MResult = bob.transpose().times(ME);
		
				// return the first row
				for (int i = 0 ; i < basis.length; i++)
						result[i] = MResult.get(i, 0);
				successful = true;
			} catch (RuntimeException x) {
				// so if we get here, the input vectors were linearly dependent.
				// we need to remove one or more of them until it's not		
				if (facetPowerSet == null) {
					facetPowerSet = new LinkedList<List<Integer>>();
					Utility.getPowerSet(origFacet, facetPowerSet);
				}
				
				if (facetPowerSet.isEmpty())
					throw x; // we should never get here :)
				
				facet = facetPowerSet.get(0);
				facetPowerSet.remove(0);
				
			}
		}
		return result;
	}*/

}
