package pdvisual;

import java.util.*;
import java.io.*;
import chemistry.*;

/**
 * PDData.java
 * 
 * Holds data representing a computed phase diagram, the stable as
 * well as the unstable compounds:
 *  - List<ComputedEntry> allICEs - all ICE's we've considered for this PD
 *  - List<Element> elements - constituent elements
 *  - Map <Element, Double> elements - map from elements to applied chemical potentials
 *  - Set<Set<Integer>> indxFacets - indices into everything of Entries making up the facets
 * 
 * @author wtipton
 * 
 * @author Chris Fischer <ccfisch@gmail.com>
 * CHANGES: Fri Aug  8 01:45:48 EDT 2008
 * 1) the method getEnergies() should return the *total* energy of each entry
 * 2) added method getEnergiesPerAtom() that can be used to construct phase 
 *    diagram
 * 3) handle treatment of input chemical potentials. For the purpose of a phase
 *    diagram the chemical potential of element 'i' is the energy at which 'i' 
 *    is supplied to the system -- so it is the reference state of the end
 *    members. To incorporate this effect we now analyze the supplied element
 *    energies in @see setRefStates(). In this method each input chemical
 *    potential is compared to the lowest energy/atom element state and these
 *    element states are shifted accordingly so that the lowest energy/atom
 *    entry has an energy corresponding to the input chemical potential (through
 *    the use of GrandPotEntry objects).  
 *    
 * CHANGES: Sat Aug  9 14:11:55 EDT 2008
 * 1) add input parameter check to ensure that at least one entry is supplied
 *    for each reference state (element)
 *
 */
public class PDData implements Serializable, IPDPlottable {
    // version identifier
    static final long serialVersionUID = 1;    // Keep an ordered list of all entries we've considered in building this PDData.
    // For the sake of storage space, if we're going to be storing lots and lots of
    // PDDatas somewhere, we should probably not repeatedly store the entries, but
    // some sort of reference to them.
    private List<IComputedEntry> allICEs;
    private final Map<Element, IComputedEntry> refStates;    // Keep a list of the elements that bound this phase diagram.
    private List<Element> elements;
    private List<Integer> basisIndxs;
    private List<Integer> indxUnstableEntries;
    /*
     *  Keep a mapping to tell us the chemical potential of each of the elements.
     *  
     *  Tue Aug  5 02:14:15 EDT 2008: Chris
     *  Note we don't need chempots for every element, just those which you want
     *  to constrain to be a particular value. For example in the Ca-O system, I
     *  may want the O energy to be -4.251 eV, this is higher than the energy of
     *  O in GGA, but the Ca energy should be the same
     */
    private Map<Element, Double> chempots;    // Keep a list of lists to represent the facets in the phase diagram.  
    // In particular, each List contained in indxFacets is a list of integers
    // which are the indices into allICEs of the entries on that facet.
    private List<List<Integer>> indxFacets;
    // Keep a record of which facets are adjacent.
    private List<List<Integer>> adjacencyList;

    public PDData(List<IComputedEntry> all, List<Element> els, Map<Element, Double> cps, List<List<Integer>> f) {
        checkInputParameters(all, els);
        chempots = cps;
        elements = els;
        indxFacets = f;
        allICEs = all;
        refStates = setRefStates();
    }

    public PDData(List<IComputedEntry> all, List<Element> els, Map<Element, Double> cps) {
        this(all, els, cps, null);
    }

    private void checkInputParameters(List<IComputedEntry> ents, List<Element> el) {
        for (Element e : el) {
            boolean refPresent = false;
            for (IComputedEntry ent : ents) {
                List<Element> entEls = ent.getComposition().getElements();
                if (entEls.contains(e)) {
                    refPresent = true;
                }
            }
            if (!refPresent) {
                throw new IllegalArgumentException("ERROR: in the list of " +
                        "entries supplied, no reference is available for the element " +
                        e.getSymbol() + " !");
            }
        }
    }

    /**
     * Last Update: Tue Aug  5 02:57:08 EDT 2008 by Chris Fischer 
     * 
     * this method picks off the elemental entry with the lowest energy/atom or
     * sets the reference state to that given in the chempots input map
     * for all elements in the phase diagram
     * 
     * the tricky part is we need to shift the energies of all the element states
     * (up or down) so that the lowest "energy" structure has the same value as the
     * input chemical potential. I'm not sure this is the cleanest solution, but let's
     * try it for now ... 
     * 
     * @return
     */
    
    
    private Map<Element, IComputedEntry> setRefStates() {
        Map<Element, IComputedEntry> lowE = new HashMap<Element, IComputedEntry>();
        // create the list of indices of the basis elemental entries here, also
        basisIndxs = new LinkedList<Integer>();

        for (Element e : elements) {
            double lowestE = Double.MAX_VALUE;
            for (IComputedEntry ent : allICEs) {
                List<Element> el = ent.getComposition().getElements();
                if (el.size() != 1 || !el.get(0).equals(e)) {
                    continue;
                }
                double ePerAtom = ent.getEnergyPerAtom();
                if (ePerAtom < lowestE) {
                    lowE.put(e, ent);
                    lowestE = ePerAtom;
                }
            }
            basisIndxs.add(allICEs.indexOf(lowE.get(e)));
        }

        
        //  now we know the lowest energy structure for each element. apply a shift
        //  to each pure element structure so that the reference states match the given
        //  chemical potentials 
         /*
        for (int idx = 0; idx < allICEs.size(); idx++) {
            ComputedEntry ent = allICEs.get(idx);
            List<Element> els = ent.getComposition().getElements();
            if (els.size() != 1 || !chempots.containsKey(els.get(0))) {
                continue;
            }
            // mu_i = lowEperAtom - shift_i
            ComputedEntry lowEent = lowE.get(els.get(0));
            
             // if the lowEent has already been updated to a GrandPotEntry then
             // we need to get the original energy to find the shift
             
            double epaLow = (lowEent instanceof GrandPotEntry ? ((GrandPotEntry) lowEent).getOrigEntry().getEnergyPerAtom() : lowEent.getEnergyPerAtom());

            double shift = epaLow - chempots.get(els.get(0));
            Map<Element, Double> appliedMu = new HashMap<Element, Double>();
            appliedMu.put(els.get(0), shift);
            ComputedEntry shifted = new GrandPotEntry(ent, appliedMu);

            allICEs.remove(idx);
            allICEs.add(idx, shifted);

            
             // replace reference state if this entry is the lowest energy entry
             
            if (ent.equals(lowEent)) {
                lowE.put(els.get(0), shifted);
            }
        }
        */
        return lowE; 
    } 

    /**
     * returns a copy of the ref states for this PDData object
     * @return
     */
    /*
    public Map<Element, ComputedEntry> getRefStates() {
        return new HashMap<Element, ComputedEntry>(refStates);
    } */
    // return a copy of the list
    public List<IComputedEntry> getAllEntries() {
        return new LinkedList<IComputedEntry>(allICEs);
    }

    public IComputedEntry getEntry(int i) {
        return allICEs.get(i);
    }

    public List<IComputedEntry> getEntriesOnFacet(int indxFacet) {
        List<IComputedEntry> result = new LinkedList<IComputedEntry>();

        for (Integer i : indxFacets.get(indxFacet)) {
            result.add(allICEs.get(i));
        }
        return result;
    }

    public List<List<Integer>> getIndxFacets() {
        return indxFacets;
    }

    public Map<Element, Double> getChemPots() {
        return chempots;
    }

    public void setFacets(List<List<Integer>> f) {
        indxFacets = f;
    }

    public Set<Integer> getIndxStableEntries() {
        Set<Integer> indxStable = new HashSet<Integer>();

        for (List<Integer> facet : indxFacets) {
            indxStable.addAll(facet);
        }
        return indxStable;
    }

    /**
     * return the entries stable
     * @return
     *@author geoffroy
     *@date Sep 13, 2008
     */
    public List<IComputedEntry> getStableEntries() {
        List<IComputedEntry> list = new LinkedList<IComputedEntry>();
        Set<Integer> set = this.getIndxStableEntries();
        for (Integer i : set) {
            list.add(this.getEntry(i));
        }
        return list;
    }

        /**
     * return the stable entries by name
     * @return
     *@author geoffroy
     *@date Sep 13, 2008
     */
    public List<String> getStableEntriesName() {
        List<IComputedEntry> list = getStableEntries();
        List<String> namelist = new ArrayList<String>(list.size());
        for (IComputedEntry entry : list) {
            namelist.add(entry.getLabel());
        }
        return namelist;
    }

    public int getDimension() {
        // If this were to change (to generalize for pseudo-pd's etc..) then also
        // consider changing the definition of the 'dim' variable in MyConvexable.
        return elements.size();
    }
    // returns a list of the coordinates of all entries
    public List<Double[]> getCompositionFractions() {
        List<Double[]> coords = new LinkedList<Double[]>();
        Iterator<IComputedEntry> iter = getAllEntries().iterator();

        while (iter.hasNext()) {
            coords.add(getCompositionFraction(iter.next()));
        }
        return coords;
    }

    public List<Integer> getIndxUnstableEntries() {
        // start out with a list of all entries, and remove ones which lie on a facet
        if (indxUnstableEntries == null) {
            indxUnstableEntries = new LinkedList<Integer>();

            for (int i = 0; i < allICEs.size(); i++) {
                indxUnstableEntries.add(i);
            }
            for (List<Integer> f : indxFacets) {
                for (Integer i : f) {
                    indxUnstableEntries.remove(i);
                }
            }
        }

        return indxUnstableEntries;
    }

    public Double[] getCompositionFraction(IComputedEntry e) {
        Double[] newVertex = new Double[getDimension()];

        for (int i = 0; i < getDimension(); i++) {
            newVertex[i] = e.getComposition().getFractionalCompo(getElements().get(i));
        }
        return newVertex;
    }

    public Element[] getElementsArray() {
        Element[] result = new Element[elements.size()];
        result = elements.toArray(result);
        return result;
    }

    public List<Element> getElements() {
        return elements;
    }

    /**
     * @return the energy of each entry (note these are TOTAL energies, not normalized !)
     *
     * @see com.cmcweb.analysis.IPDPlottable#getEnergies()
     */
    public List<Double> getEnergies() {
        List<Double> result = new LinkedList<Double>();

        for (IComputedEntry e : allICEs) {
            result.add(e.getTotalEnergy());
        }
        return result;
    }

    /**
     * @return the energy of each entry normalized per atom
     * @return
     */
    public List<Double> getEnergiesPerAtom() {
        List<Double> result = new LinkedList<Double>();

        for (IComputedEntry e : allICEs) {
            /** Shyue Ping Ong : I modified the line below to make use of the getEnergyPerAtom()
             * which is part of the ComputedEntry interface.  The original formulation
             * result.add(e.getEnergy() / e.getStructure().getNumSites());
             * causes catastrophic errors when the entry is a modified one in which the relevant
             * sites do not conform to the sites in the structure!!  e.g. in a grand canonical construct.
             * */
            result.add(e.getEnergyPerAtom());
        }
        return result;
    }

    /**
     * returns the energy of the i^th entry 
     * @param i
     * @return
     */
    public double getEnergy(int i) {
        return allICEs.get(i).getTotalEnergy();
    }

    /**
     * returns the formation energies (per atom) of all entries
     * 
     * @return
     *
     * @see com.cmcweb.analysis.IPDPlottable#getFormEnergiesPerAtom()
     */
    
    public List<Double> getFormEnergiesPerAtom() {
        List<Double> eforms = new ArrayList<Double>(allICEs.size());
        for (IComputedEntry ice : allICEs) {
            eforms.add(getEformPerAtom(ice));
        }
        return eforms;
    }

    private double getEformPerAtom(IComputedEntry i) {
        double energy = i.getEnergyPerAtom();
        for (Element e : elements) {
          //  energy -= chempots.get(e) * i.getComposition().getFractionalCompo(e);
        	if (!refStates.containsKey(e)) {
        		System.out.println("WARNING: No reference state for " + e.getSymbol() + " in PDData.getEformPerAtom().  Using 0.0.");
        		continue;
        	}
            energy -= refStates.get(e).getEnergyPerAtom() * i.getComposition().getFractionalCompo(e);
        }
        return energy;
    }
    

    public void setAdjacencyList(List<List<Integer>> adjList) {
        adjacencyList = adjList;
    }

    public List<List<Integer>> getAdjacencyList() {
        return adjacencyList;
    }

    public List<String> getVertexLabels() {
        List<String> result = new LinkedList<String>();

        for (IComputedEntry e : allICEs) {
            result.add(e.getLabel());
        //result.add(e.getComposition().toStringNice());
        }
        return result;
    }

    public List<String> getAxisLabels() {
        List<String> result = new LinkedList<String>();

        Iterator<Element> i = elements.iterator();
        while (i.hasNext()) {
            Element e = i.next();
            result.add(e.getSymbol());
        }

        return result;
    }

    public String toString() {
        StringBuilder sb = new StringBuilder();
        sb.append("Stable Entries: \n");

        for (int i : getIndxStableEntries()) {
            sb.append(allICEs.get(i).toString());
        }
        return sb.toString();
    }

    public boolean isPurePhase(int i) {
        return true;
    }

    public List<Integer> getPhasesWSameCompAs(int indxEntry) {
        List<Integer> result = new LinkedList<Integer>();

        for (int i = 0; i < allICEs.size(); i++) {
            if (allICEs.get(i).getComposition().equals(allICEs.get(indxEntry).getComposition())) {
                result.add(i);
            //sameCompoList.add(allICEs.get(i));
            //for now we'll just sort the results by energy per atom
            }
        }
        List<Integer> sorted = new LinkedList<Integer>();

        Integer chosenOne = null;
        while (result.size() != 0) {
            double energymin = Double.MAX_VALUE;
            for (Integer i : result) {
                if (energymin > allICEs.get(i).getEnergyPerAtom()) {
                    energymin = allICEs.get(i).getEnergyPerAtom();
                    chosenOne = i;
                //    System.out.println(chosenOne);
                }
            }
            result.remove(chosenOne);
            sorted.add(chosenOne);
        }

        return sorted;
    }

    public boolean isStable(int i) {
        return !getIndxUnstableEntries().contains(i);
    }

    public Composition getComposition(int i) {
        return getEntry(i).getComposition();
    }

    public int getNumEntries() {
        return allICEs.size();
    }

    public List<Composition> getCompositions() {
        List<Composition> result = new LinkedList<Composition>();

        for (IComputedEntry e : allICEs) {
            result.add(e.getComposition());
        }
        return result;
    }

    public List<Integer> getBasisEntryIndxs() {
        return basisIndxs;
    }

    public List<IComputedEntry> getBasisEntries() {
        List<IComputedEntry> result = new LinkedList<IComputedEntry>();

        for (Integer i : basisIndxs) {
            result.add(allICEs.get(i));
        }
        return result;
    }
}
