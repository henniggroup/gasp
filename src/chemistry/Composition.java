package chemistry;

import java.io.Serializable;
import java.util.*;

public class Composition implements Serializable {
	
	static final long serialVersionUID = 1l;

	Map <Element, Double> components;
	
	public static double COMP_TOL = 0.001;
	
	public Composition(Map <Element, Double> _components) {
		components = _components;
	}
	
	public Composition(Map <Element, Integer> _components, boolean f) {
		components = new HashMap<Element,Double>();
		for (Element e : _components.keySet())
				components.put(e, (double)(_components.get(e)));
	}
	
	public Composition(Element[] elementsArray, Integer[] values) {
		if (elementsArray.length != values.length)
			System.out.println("ERROR: elementsArray.length != values.length in Composition()");
		
		components = new HashMap<Element,Double>();
		for (int i = 0; i < elementsArray.length; i++)
			components.put(elementsArray[i], values[i].doubleValue());
	}
	
	public Composition(Element[] elementsArray, Double[] values) {
		if (elementsArray.length != values.length)
			System.out.println("ERROR: elementsArray.length != values.length in Composition()");
		
		components = new HashMap<Element,Double>();
		for (int i = 0; i < elementsArray.length; i++)
			components.put(elementsArray[i], values[i]);
	}

	public double getNumComponents() {
		double sum = 0;
		for (Double i : components.values())
			sum += i;
		return sum;
	}
	
	public int getNumElements() {
		return components.keySet().size();
	}
	
	public List<Element> getElements() {
		List<Element> result = new LinkedList<Element>();
		for (Element e : components.keySet())
			if (components.get(e) > COMP_TOL)
				result.add(e);
		return result;
	}
	
	public double getFractionalCompo(Element e) {
		if (! components.containsKey(e))
			return 0;
		return components.get(e) / getNumComponents();
	}
	
	public double getOrigAmount(Element e) {
		if (! components.containsKey(e))
			return 0.0;
		
		return components.get(e);
	}
	
	public String toString() {
		StringBuilder ans = new StringBuilder();
		
		for (Element e : components.keySet()) 
			ans.append(e.getSymbol() + " " + components.get(e) + " ");
		
		return ans.toString();
	}
	
	public boolean isEmpty() {
		return getNumComponents() == 0;
	}
	
	public Map<Element,Double> getStoichiometricUnit() {
		// TODO: return lowest common denominator stoichiometric unit ??
		
		// return a copy of our internal structure
		Map<Element,Double> result = new HashMap<Element,Double>();
		for (Element e : components.keySet())
			result.put(e, components.get(e));
		return result;
	}
	
	public Map<Element,Double> getFractionalComposition() {
		Map<Element,Double> result = new HashMap<Element,Double>();
		for (Element e : components.keySet())
			result.put(e, getFractionalCompo(e));
		return result;
	}
	
    public boolean equals(Composition o) {
        return equals(o, COMP_TOL); 
    }    
    
    public boolean equals (Composition o, double tol) {
        Composition me = scrubSmallAmounts(this, tol, false);
        Composition so = scrubSmallAmounts(o, tol, false); 
        
        // make sure they have the same elements
    	for (Element e : me.getElements())
    		if (! so.getElements().contains(e))
    			return false;
    	for (Element e : so.getElements())
    		if (! me.getElements().contains(e))
    			return false;

    	for (Element e : me.getElements())
    		if ( Math.abs(so.getFractionalCompo(e) - me.getFractionalCompo(e)) > tol )
    			return false;

        
        return true;

    }
    
    /**
     * This method will remove components of a composition object that have an
     * amount less than the given @param tol. When @param use_orig_amounts == true
     * then the scrubbing is done on the original amount array 
     * 
     * Example: 
     * 
     * let tol = 0.01 && use_orig_amounts == false
     * scrubbing A0.999B0.001 will yield A 
     * 
     * scrubbing AB0.001 with use_orig_amounts==true will yield A
     * 
     * @param x
     * @param tol
     * @param use_orig_amounts
     * @return
     */
    public static Composition scrubSmallAmounts(Composition x, double tol, boolean use_orig_amounts){
        List<Element> scrubbedE = new ArrayList<Element>(x.getNumElements());
        List<Double>  scrubbedD = new ArrayList<Double>(x.getNumElements()); 
        for(int i=0; i<x.getNumElements(); i++){
        	Element e = x.getElements().get(i);
            double amt = (use_orig_amounts ? x.getOrigAmount(e) : x.getFractionalCompo(e)); 
            if(amt > tol){
                scrubbedE.add(x.getElements().get(i));
                scrubbedD.add(amt); 
            }
        }
        Element[] els = new Element[scrubbedE.size()];
        scrubbedE.toArray(els); 
        Double[] amts = new Double[scrubbedD.size()];
        scrubbedD.toArray(amts); 
        return new Composition(els, amts); 
    }


}
