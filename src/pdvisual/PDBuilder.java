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

import Jama.Matrix;

import chemistry.*;
import crystallography.*;

import ga.GAOut;
import ga.GAParameters;
import ga.Organism;

import java.io.Serializable;
import java.util.*;


/**
 * PDBuilder.java
 * 
 * Takes composition/source info and creates a PDData object.
 * Plan on this being a Builder patterned class to deal with 
 * different types of data sources.
 * 
 * External interface consists of two calls:
 *  - construct a PDBuilder object - give it various info about the system
 *  - call getPDData() - to get the finished PDData object
 * Most actual computation will be done during the getPDData call.
 * 
 * @author wtipton
 *
 */
public class PDBuilder implements Serializable {
	final static long serialVersionUID = 1;

    public static double DET_TOL = 0.00001;
    /**
     * if a vertex of the convex hull lies above this number in formation energy
     * then it shouldn't be a facet ! 
     */
    public static double EFORM_TOL = 0.00001;
    protected List<IComputedEntry> entries;
    protected List<Element> elements;
    protected Map<Element, Double> chempots;

    public PDBuilder(List<IComputedEntry> _entries, List<Element> _elements, Map<Element, Double> _chempots) {
        entries = new LinkedList<IComputedEntry>(_entries);
        elements = new LinkedList<Element>(_elements);
        chempots = new HashMap<Element, Double>(_chempots);
    }
    
    public void addEntry(IComputedEntry c) {
    	entries.add(c);
    }

    public PDData getPDData() {

        PDData newPDData = new PDData(entries, elements, chempots);
        
        /*
         * qconvex doesn't like when we have a number of equals to the dimension of the convex hull space
         * Ex: when we have a binary (dim 2) with only 2 entries basically a trivially immiscible system
         * qconvex would throw something like: qhull   error: not enough points (2) to construct initial simplex (need 3)
         * 
         * In this case we just return directly the PDData and we don't compute the convexhull
         * 
         * Geoffroy, 21 october 2008
         */

        if (entries.size() == elements.size()) {
            List<List<Integer>> facets = new LinkedList<List<Integer>>();
            List<Integer> list = new LinkedList<Integer>();
            for (int i = 0; i < entries.size(); i++) {
                list.add(i);
            }
            facets.add(list);
            newPDData.setFacets(facets);
            return newPDData;
        }

        // find the convex hull
        PDConvexable convexableData = new PDConvexable(newPDData); 
        Convexhull chull = new Convexhull(convexableData, false);
        List<List<Integer>> facets = chull.getFaLists();
        
        newPDData.setCHullVolume(chull.getVolume());
        
        // remove vertical facets:
        // note that the PDk option in qconvex is not good for this as it is susceptible
        // to precision error
        /**
         * Fri Jul 18 01:59:29 EDT 2008 Chris
         * One, of probably several ways to remove vertical facets is to determine when the
         * matrix c[i][j] == amount of i^th element in composition of j^th vertex and
         * dim(c) = numElement x numElement has a determinant of zero (i.e., 
         * it lies in a subspace of the composition space
         * Example: facets in Li-Ti-O include some which lie along the Ti-O axis
         *          and are parallel to the energy axis 
         */
        Iterator<List<Integer>> iFacet = facets.iterator();
        while (iFacet.hasNext()) {
            List<Integer> facet = iFacet.next();
            boolean needToRemove = false;
            //build array 
            double[][] C = new double[elements.size()][elements.size()];

            for (int i = 0; i < facet.size(); i++) {
                Composition ci = entries.get(facet.get(i)).getComposition();
                for (int j = 0; j < elements.size(); j++) {
                    C[j][i] = ci.getFractionalCompo(elements.get(j));
                }
            }
            Matrix m = new Matrix(C);
            if (Math.abs(m.det()) < DET_TOL) {
                needToRemove = true;
            }
            if (needToRemove) {
    			GAOut.out().stdout("removing facet (all vertices on vertical face of hull):", GAOut.DEBUG);
                for (int i = 0; i < facet.size(); i++) {
                	GAOut.out().stdout(newPDData.getEntry(facet.get(i)).toString(), GAOut.DEBUG);
                }
                iFacet.remove();
            }
        }
        /**
         * remove facets where the formation energy of
         * one of the points on the facet is positive (i.e., it is part of the
         * convex hull, but in the Eform > 0 space)
         * 
         * and lastly, in the case where points with negative formation energy exist,
         * and no points with positive formation energy exist, we have the "flat"
         * facet whose endpoints are the elemental phases, and we want to remove this
         */
        iFacet = facets.iterator();
        List<Double> eforms = newPDData.getFormEnergiesPerAtom();
        
        boolean haveNegEformEntries = false;
        for (Double e : eforms) {
        	if (e < 0) {
        		haveNegEformEntries = true;
        		break;
        	}
        }
        
        while (iFacet.hasNext()) {
            List<Integer> facet = iFacet.next();
            boolean removeFacet = false;
            for (Integer ivert : facet) {
                if (eforms.get(ivert) > EFORM_TOL) {
                    removeFacet = true;
                    break;
                }
            }
            if (!removeFacet && haveNegEformEntries) {
            	removeFacet = true;
                for (Integer ivert : facet) {
                	// remove facet if eforms of all endpoints are essentially 0
                	removeFacet &= Math.abs(eforms.get(ivert)) < EFORM_TOL;
                }
            }
            if (removeFacet) {
                iFacet.remove();
            }
        }

        newPDData.setFacets(facets);

        // compute the adjacent facets
        // build the adjacency list
        List<List<Integer>> adjacencyList = new LinkedList<List<Integer>>();
        // get two facets and see if they are adjacent.  they are adjacent if they have
        // (dim - 1) vertices in common.
        for (int i = 0; i < facets.size(); i++) {
            for (int j = 0; j < facets.size(); j++) {
                if (i == j) {
                    continue;
                }
                List<Integer> fi = facets.get(i);
                List<Integer> fj = facets.get(j);
                int numVerticesInCommon = 0;
                for (int m = 0; m < fi.size(); m++) {
                    for (int n = 0; n < fj.size(); n++) {
                        if (fi.get(m).equals(fj.get(n))) {
                            numVerticesInCommon++;
                        }
                    }
                }
                if (numVerticesInCommon == elements.size() - 1) {
                    List<Integer> newPair = new LinkedList<Integer>();
                    newPair.add(i);
                    newPair.add(j);
                    adjacencyList.add(newPair);
                }
            }
        }
        newPDData.setAdjacencyList(adjacencyList);

        return newPDData;
    }

	public boolean containsEntry(IComputedEntry o) {
		for (IComputedEntry e : entries)
			if (o == e) // just compare references?
				return true;
		return false;
	}
}
