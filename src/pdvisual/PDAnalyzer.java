package pdvisual;

import Jama.Matrix;

import chemistry.*;
import crystallography.*;

import java.util.*;

/**
 * PDAnalyzer.java
 * 
 * Takes a PDData object which we can then query for information
 * about them.
 * 
 * @author wtipton
 *
 */
public class PDAnalyzer implements IPDAnalyzer {

    /**
     * Fri Jul 18 02:27:07 EDT 2008 Chris
     * this is the composition tolerance used to determine if a composition lies
     * on a vertex in the phaes diagram 
     */
    public static double COMP_TOL = 0.000001;
    
    public static double DET_TOL = 0.0001;
    /**
     * Fri Jul 18 02:27:57 EDT 2008 Chris
     * this is the composition tolerance used to determine when a given composition
     * does not lie in a particular facet @see compIsInFacet()
     * 
     */
    public static double NEG_COMP_TOL = -0.00000001;
    PDData pdd;

    public PDAnalyzer(PDData _pdd) {
        pdd = _pdd;
    }

    // Returns the distance above the convex hull of a given entry
    // which may or may not be in the PDData's allEntries list.  In this
    // case, if this function returns a negative number, it indicates
    // that the phase diagram was incomplete.
    public double getEnergyPerAtomAboveHull(IComputedEntry e) {
        return e.getEnergyPerAtom() - getEnergyPerAtomOnHull(e.getComposition());
    }
    // returns the distance above the convex hull of the entry specified
    // by the given index into the pdd.allEntries list

    public double getEnergyPerAtomAboveHull(int i) {
        return getEnergyPerAtomAboveHull(pdd.getEntry(i));
    }

    // Returns a PDData object representing the phase diagram without
    // the specified entry.
    public PDData getPDwoEntry(IComputedEntry e) throws IllegalArgumentException {
        // get the index of the entry to remove
        int indexOfE = pdd.getAllEntries().indexOf(e);
        if (indexOfE == -1) {
            throw new IllegalArgumentException("Given entry not in the phase diagram");
        }

        // copy the list and remove e
        List<IComputedEntry> entries = new LinkedList<IComputedEntry>(pdd.getAllEntries());
        ;
        entries.remove(indexOfE);

        // make a new PDData
        PDBuilder pdb = new PDBuilder(entries, pdd.getElements(), pdd.getChemPots());
        return pdb.getPDData();
    }
    // returns a PDData object representing the phase diagram without
    // the entry specified by the given index into the pdd.allEntries list

    public PDData getPDwoEntry(int i) {
        return getPDwoEntry(pdd.getEntry(i));
    }

    /**
     * gives the phase diagram without the entry and without the phases not
     * on the convexhull
     * @param e Entry
     * @return
     *@author geoffroy
     *@date Sep 13, 2008
     */
    public PDData getPDwoEntryandUnstable(IComputedEntry e) {
        // get the index of the entry to remove
        int indexOfE = pdd.getStableEntries().indexOf(e);
        if (indexOfE == -1) {
            throw new IllegalArgumentException("Given entry not in the phase diagram and/or stable");
        }

        // copy the list and remove e
        List<IComputedEntry> entries = new LinkedList<IComputedEntry>(pdd.getStableEntries());
        ;
        entries.remove(indexOfE);

        // make a new PDData
        PDBuilder pdb = new PDBuilder(entries, pdd.getElements(), pdd.getChemPots());
        return pdb.getPDData();
    }

    /**
     * Shyue Ping Ong
     * Returns range of stable elemental chemical potentials for a phase.
     * @param entry Entry of interest
     * @param el Element of interest
     * @return double[2] array with [0] element having the min. chem. pot and [1] element having the max.
     */
    /*
    public double[] getStableChemPotRange(ComputedEntry entry, Element el) {
        double[] range = new double[2];
        range[0] = +10000;
        range[1] = -10000;
        List<List<Integer>> facets = pdd.getIndxFacets();
        int entryIndex = pdd.getAllEntries().indexOf(entry);

        for (int i = 0, n = facets.size(); i < n; i++) {
            List<Integer> facet = facets.get(i);
            if (facet.contains(entryIndex)) {
                double chempot = getChemicalPotential(el, facet);
                range[0] = (range[0] < chempot) ? range[0] : chempot;
                range[1] = (range[1] > chempot) ? range[1] : chempot;
            }
        }


        // Shyue Ping Ong:
        //   If upper stability chempot is the same as the ref chem pot, this phase
        //  is assumed to be stable up to +ve infinity (not strictly true, since
        //  additional phases appearing at chem pot above the ref chem pot can
        //  supersede current phases, but good enough for our purposes now)
         
        if (Math.abs(range[1] - pdd.getRefStates().get(el).getEnergyPerAtom()) < 1e-4) {
            range[1] = Float.POSITIVE_INFINITY;
        }

         // If phase does not contain element, the stability is bounded by
         // -infinity.  This should always be true.  The phase without the
         // element is basically the most stable "reduced" form.  It should not be
         // possible to get anything more stable than that.
         
        if (entry.getComposition().getOrigAmount(el) == 0) {
            range[0] = Float.NEGATIVE_INFINITY;
        }
        return range;
    } */

    // returns the MixedPhase on the phase diagram with given composition
    public MixedPhase getPhase(Composition c) {
        List<Integer> f = getFacet(c);
        List<IComputedEntry> entries = new LinkedList<IComputedEntry>();
        for (Integer i : f) {
            entries.add(pdd.getEntry(i));
        }
        return new MixedPhase(c, entries);
    }

    // Returns an index into the pdd.facets list indicating the facet in which
    // the given composition lies.  In the degenerate case (e.g., the composition
    // is on a vertex), the result is undefined.
    // TODO: more smartly!!!1 (point location problem)
    public List<Integer> getFacet(Composition c) {

        for (List<Integer> facet : pdd.getIndxFacets()) {
            if (compoIsInFacet(facet, c)) {
                return facet;
            }
        }

        System.out.println("Warning: no facet found for composition " + c);

        return null;
    }

    /**
     * Thu Jul 17 16:14:48 EDT 2008 Chris
     * Will, I'm not sure this method is working quite right. I've found in the
     * Li-Ti-O system that for (Li,Ti,O)=(0,1/3,2/3) this method gives the wrong facet
     *
     * Fri Jul 18 02:29:37 EDT 2008 Chris
     * Bug was in PDBuilder which did not sanitize the set of facets by removing
     * vertical facets. PDBuilder has been updated as well as this method
     * (see additional comments below)
     *
     * @param indxFacet
     * @param c
     * @return
     */
    private boolean compoIsInFacet(List<Integer> facet, Composition c) {

        /**
         * check to see if this composition corresponds to a vertex
         * Chris: I think we want the tolerance on the composition comparison to
         *        be quite small.
         */
        for (Integer f : facet) {
            if (pdd.getEntry(f).getComposition().equals(c)) {
                return true;
            }
        }
        //
        int numElems = pdd.getElements().size();

        /**
         * Thu Jul 17 16:18:04 EDT 2008
         * if compoIsInFacet == true then we have
         * c = \sum_i a_i c_i where c_i is the composition of the vertex
         * and a_i >= 0 \forall i
         *
         * take care to embed every composition into the possibly higher dimensional
         * space of the phase diagram
         *
         */
        /**
         * the following should be true for all facets of the phase diagram
         */
        if (facet.size() != numElems) {
            throw new RuntimeException("Internal error !! List of IComputedEntries \n" +
                    " in facet has length smaller than dimension of phase diagram !");
        }

        double[][] comp = new double[numElems][1];
        /**
         * vC[i][j] == amount of i^th element in composition of j^th vertex
         */
        double[][] vC = new double[numElems][numElems];
        for (int i = 0; i < numElems; i++) {
            Composition ci = this.pdd.getEntry(facet.get(i)).getComposition();
            comp[i][0] = c.getFractionalCompo(this.pdd.getElements().get(i)); //i^th row of comp vector
            for (int j = 0; j < numElems; j++) {
                vC[j][i] = ci.getFractionalCompo(this.pdd.getElements().get(j));
            }
        }

        Matrix jamavC = new Matrix(vC);
        Matrix jamacomp = new Matrix(comp);
        /**
         * Thu Jul 17 18:07:20 EDT 2008
         * ok, i think the bug has been found -- facet list includes stuff that is on
         * a vertical 'face' of a convex hull (i.e., some faces lie parallel to the energy axis !!)
         *
         * I will fix by filtering out these facets when the PDData is built in the PDBuilder
         */
        if (Math.abs(jamavC.det()) < PDBuilder.DET_TOL) {
            throw new RuntimeException("Internal Error !! Facet in PDData lies " +
                    "parallel to energy axis. ");
        }
        Matrix sol = jamavC.solve(jamacomp);
        /**
         * check to ensure that a_i >= 0
         */
        for (int i = 0; i < sol.getRowDimension(); i++) {
            if (sol.get(i, 0) < NEG_COMP_TOL) {
                return false;
            }
        }
        if (true) {
            return true;
        }
        /**
         * @TODO the code below isn't reachable -- clean after testing !!
         */
        // this stores n+1 determinant values
        double[] det = new double[c.getElements().size() + 1];
        double[][] M0 = new double[numElems][numElems];

        for (int i = 0; i < M0.length; i++) {
            for (int j = 0; j < M0[i].length - 1; j++) {
                M0[i][j] = pdd.getEntry(facet.get(i)).getComposition().getFractionalCompo(pdd.getElements().get(j));
            }
            M0[i][numElems - 1] = 1;
        }
        Matrix MJama = new Matrix(M0);
        det[0] = MJama.det();

        for (int i = 1; i < det.length; i++) {
            double[][] Mi = new double[numElems][numElems];
            // create Mi, it's a copy of M0 but with the ith line replaced
            // by the point compo ccordinates
            for (int j = 0; j < M0.length; j++) {
                for (int k = 0; k < M0[j].length; k++) {
                    if (j == i - 1) {
                        if (k != numElems - 1) {
                            Mi[j][k] = c.getFractionalCompo(pdd.getElements().get(k));
                        } else {
                            Mi[j][k] = 1;
                        }
                    } else {
                        Mi[j][k] = M0[j][k];
                    }
                }
            }
            Matrix MJama2 = new Matrix(Mi);
            det[i] = MJama2.det();
        }

        // check if the dets are all the same sign...
        int start = 0;
        if (det[0] == 0) {
            start++;
        }
        if (det[start] < 0) {
            for (int i = start + 1; i < det.length; i++) {
                if (det[i] > 0) {
                    return false;
                }
            }
        } else {
            for (int i = start + 1; i < det.length; i++) {
                if (det[i] < 0) {
                    return false;
                }
            }
        }

        return true;
    }

    /**
     * given a composition c, determine the lowest energy facet corresponding
     * to c and determine the chemical potentials of each element in facet.
     * @todo IMPLEMENT ! this just returns null
     * @param c
     * @return
     */
    public Map<Element, Double> getChemPots(Composition c) {
        Map<Element, Double> result = new HashMap<Element, Double>();

        List<Integer> f = getFacet(c);
        for (Element el : pdd.getElements()) {
            result.put(el, getChemicalPotential(el, f));
        }

        return result;
    }

    /**
     * This method gives the chemical potential of a given element for
     * a given facet
     * @param el
     * @param facet 
     * @return
     *
     */
    
    public double getChemicalPotential(Element el, List<Integer> facet) {
        //first of all get the equation of the facet
        List<IComputedEntry> vertices = new LinkedList<IComputedEntry>();
        for (Integer i : facet) {
            vertices.add(pdd.getEntry(i));
        }

        // create matrix M with compo in each line
        // if (x1,y1,z1), (x2,y2,z2) are the coordinates
        //we have M(1,:)=(x1,y1,1) etc...
        double[][] M = new double[facet.size()][facet.size()];
        for (int i = 0; i < facet.size(); i++) {
            for (int j = 0; j < facet.size() - 1; j++) {
                M[i][j] = vertices.get(i).getComposition().getFractionalCompo(pdd.getElements().get(j));
            }
            M[i][facet.size() - 1] = 1;
        }

        //let's invert M


        Matrix MJama = new Matrix(M);
        if (Math.abs(MJama.det()) < DET_TOL) {
            return Double.NEGATIVE_INFINITY;
        }
        Matrix MinvJama = MJama.inverse();

        double[][] Minv = MinvJama.getArray();
        //Minv*t(-E1,_E2,_E3,-E4...) is the solution for the
        //eqcoeff in a column vector
        double[] eqcoeff = new double[facet.size()];
        //make a -E1 -E2 -E3 -E4.. vector
        double[] Evect = new double[facet.size()];
        for (int i = 0; i < facet.size(); i++) {
            Evect[i] = -pdd.getEnergiesPerAtom().get(facet.get(i));
        }

        Matrix MinvM = new Matrix(Minv);
        double[][] Evect2 = new double[1][Evect.length];
        for (int i = 0; i < Evect.length; i++)
        	Evect2[0][i] = Evect[i];
        Matrix EvectM = new Matrix(Evect2);
        eqcoeff = MinvM.times(EvectM).getColumnPackedCopy();

        //now, we have the equation for this facet
        //we look then at the distance in E from ce to the hyperplane

        //compute the energy of the point at composition on the target=1.0 axis


       // for(int i=0;i<this.elem.length;i++)
       // {
       // if(this.elem[i].equals(target))
       // return eqcoeff[i]+eqcoeff[eqcoeff.length-1];
       // }

        double[] x = new double[facet.size()];
        int elIndex = pdd.getElements().indexOf(el);
        x[elIndex] = 1;

        x[facet.size() - 1] = 1.0;

        double result = 0.0;

        for (int j = 0; j < facet.size(); j++) {
            result += x[j] * eqcoeff[j];
        }

        return (-result);
    } 

    /*
    public double getRelativeChemicalPotential(Element el, List<Integer> facet) {
        double refEnergy = this.pdd.getRefStates().get(el).getEnergyPerAtom();
        
        //  Check if MUREFS has any reference state for the element
        //  if yes use that if no, use the element reference state energy in the pd
         

        return getChemicalPotential(el, facet) - (refEnergy + PDDefaults.getChemPotShift(el, refEnergy));
    } */

    /*
     * Shyue Ping Ong Nov 1 2008
     * Returns list of critical chemical potentials for Element in phase diagram,
    duplicates removed.
     * */
    /*
    public List<Double> getCriticalChemicalPotentials(Element el) {
        double TOLERANCE = 1e-4;
        List<Double> chemPotList = new ArrayList<Double>();
        for (List<Integer> facet : pdd.getIndxFacets()) {
            double chempotFacet = getChemicalPotential(el, facet);
            boolean found = false;
            for (double check : chemPotList) {
                if (Math.abs(check - chempotFacet) < TOLERANCE) {
                    found = true;
                    break;
                }

            }
            if (!found) {
                chemPotList.add(chempotFacet);
            }
        }
        return chemPotList;
    } */

    /*
    public Map<ComputedEntry, double[]> getRelativeChemicalPotentialRange(Element target) {

        Map<ComputedEntry, double[]> result = new HashMap<ComputedEntry, double[]>();

        HashMap<List<Integer>, Double> facetmap = new HashMap<List<Integer>, Double>();
        for (List<Integer> fac : pdd.getIndxFacets()) {
            double mu = getRelativeChemicalPotential(target, fac);
            facetmap.put(fac, mu);
        }


        //then look for each VaspEntry in the chull. what are the different chemical pot
        Iterator<Integer> iter = pdd.getIndxStableEntries().iterator();
        while (iter.hasNext()) {
            ComputedEntry ve = pdd.getEntry(iter.next());

            Iterator<List<Integer>> iterfac = pdd.getIndxFacets().iterator();
            double[] values = new double[0];
            while (iterfac.hasNext()) {
                List<Integer> fac = iterfac.next();
                if (fac.contains(pdd.getAllEntries().indexOf(ve))) {
                    values = Utility.appendArray(values, facetmap.get(fac));
                }
            }
            Arrays.sort(values);
            double[] newvalue = new double[2];
            newvalue[0] = values[0];
            newvalue[1] = values[values.length - 1];
            //result.put(ve, values[values.length-1]-values[0]);
            result.put(ve, newvalue);
        }
        return result;

    } */

    /**
     * returns the energy/atom on the hull
     *
     * @param c
     * @return
     */
    public double getEnergyPerAtomOnHull(Composition c) {

        // for numerical reasons it's best to do this in two steps:
        //  1) search to see if there is a stable vertex at the given composition

        List<IComputedEntry> entries = pdd.getAllEntries();
        for (int i = 0; i < entries.size(); i++) {
            if (pdd.isStable(i) && entries.get(i).getComposition() == c) {
                return entries.get(i).getEnergyPerAtom();
            }
        }

        //	2) otherwise it's probably at the interior of a facet so things work out ok
        //		if we look for the facet it's in and do the linear combination thing

        // get the facet in which the composition lies

        List<Integer> facet = getFacet(c);

        // make a list of Compositions from which the given Composition will be made
        List<IComputedEntry> basis = new LinkedList<IComputedEntry>();
        for (Integer i : facet) {
            basis.add(pdd.getEntry(i));
        }

        // get coefficient of each constituent
        MixedPhase mphase = new MixedPhase(c, basis);
        return mphase.getEnergyPerAtom();
    }


    // NB: doing traces like this kinda defeats the purpose of having done
    // the pseudophasediagrams right in the first place.  just make a pseudopddata
    // and work with it's mixedphases.
    Map<Composition, Double> getEnergyTrace(Composition a, Composition b) {
        //TODO
		/*
        int numSteps = 10;
        for (int i = 0; i < numSteps; i++) {
        double cFracs[] = new double[3];
        for (int j = 0; j < 3; j++)
        cFracs[j] = comp1[j]*(i/numSteps) + comp2[j]*((numSteps-i)/numSteps);
        Composition co = new Composition(elems, cFracs);
        //MixedPhase m = new MixedPhase(co,pdd.getEntriesOnFacet(pd));
        System.out.println(co);
        System.out.println(pda.getEnergyOnHull(co));
        }*/
        return null;
    }

    /**
     * returns the volume equivalent to a linear combination of stable phases
     * on the hull for a given composition.
     *
     * NOTE: returns an extensive quantity. In other words, if you ask for the
     *       volume of "A1 B2" and this method returns X, then you ask for the
     *       volume of "A2 B4" this method will return 2X.
     *
     * @param c
     * @return
     */
    public double getVolumeOnHull(Composition c) {
        MixedPhase m = this.getPhase(c);
        return m.getVolume();
    }

    /**
     * DOES NOT WORK !!!! returns null
     * @param a
     * @param b
     * @return
     */
    public Map<Composition, Double> getVolumeTrace(Composition a, Composition b) {
        //TODO
        return null;
    }

    /**
     * Note that the definition for the inverse distance to the hull
     * is now:
     * 1)remove all unstable phases
     * 2) remove the given phase
     * 3) rebuild the convex hull
     * 4) look how far the removed phase is from the new hull
     * @param i
     * @return
     *@author geoffroy
     *@date Sep 14, 2008
     */
    public double getInvDistToHullPerAtom(int i) {
        return getInvDistToHullPerAtom(pdd.getAllEntries().get(i));
    }

    // is positive if given entry is on hull.
    public double getInvDistToHullPerAtom(IComputedEntry e) {
        PDAnalyzer pda = new PDAnalyzer(this.getPDwoEntryandUnstable(e));
        return pda.getEnergyPerAtomOnHull(e.getComposition()) - e.getEnergyPerAtom();
    }

    Map<Composition, MixedPhase> getCompositionTrace(Composition a, Composition b) {
        //TODO
        return null;
    }

    /**
     * just spits out a report with all the info you wanted about the phases
     * in a phase Diagram: formation energy, distance to the hull, chemical pot etc...
     * Geoffroy
     *
     * @return a String of report
     */
    /*
    public String getReport(Element chempotTarget, int nbofelems) {
        String results = "";
        Map<ComputedEntry, double[]> murange = this.getRelativeChemicalPotentialRange(chempotTarget);
        int count = 0;
        for (ComputedEntry e : this.pdd.getAllEntries()) {
            if (e.getComposition().getElements().length == nbofelems) {
                String formulanice = Utility.removeChar(Utility.getNormalizedFormula(e.getComposition().getformulasum()), ' ');
                if (this.pdd.getIndxStableEntries().contains(this.pdd.getAllEntries().indexOf(e))) {
                    results += e.getStructureID() + " " + formulanice + " " + e.getEnergyPerAtom() + " " + this.getEnergyPerAtomAboveHull(e) + " " + murange.get(e)[0] + " " + murange.get(e)[1] + " " + this.getInvDistToHullPerAtom(count) + "\n";
                } else {
                    results += e.getStructureID() + " " + formulanice + " " + e.getEnergyPerAtom() + " " + this.getEnergyPerAtomAboveHull(e) + " null null " + " " + "null" + "\n";
                }
            }
            count++;
        }
        return results;
    } */
}
