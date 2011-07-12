package chemistry;

import ga.GAParameters;

import java.io.Serializable;
import java.util.*;

import Jama.Matrix;

import pdvisual.PseudoPDBuilder;

import utility.*;

public class CompositionSpace implements Serializable {
	static final long serialVersionUID = 1;
	
	private List<Element> elements;
	private List<Vect> basis;
	private List<Composition> endpoints;
	
	public CompositionSpace(List<String> compTokens, boolean hi) {
		if (compTokens.size() == 0)
			throw new IllegalArgumentException("Malformed compositionSpace passed.");
		Iterator<String> compStrs = compTokens.iterator();
		//  now, read the number of elements
		int numElems = Integer.parseInt(compStrs.next());
		if ((compTokens.size() - 1) % numElems != 0)
			throw new IllegalArgumentException("Malformed compositionSpace passed.");
		// read off the elements
		List<Element> elements = new LinkedList<Element>();
		for (int i = 0; i < numElems; i++)
			elements.add(Element.getElemFromSymbol(compStrs.next()));
		// read off the compositions
		List<Composition> endpoints = new LinkedList<Composition>();
		while (compStrs.hasNext()) {
			Map<Element,Integer> compMap = new HashMap<Element,Integer>();
			for (int i = 0; i < numElems; i++)
				compMap.put(elements.get(i), Integer.parseInt(compStrs.next()));
			endpoints.add(new Composition(compMap, false));
		}	
		setElementsAndBasis(endpoints);
	}
	
	public CompositionSpace(List<Composition> _endpoints) {
		setElementsAndBasis(_endpoints);
	}
	
	private void setElementsAndBasis(List<Composition> _endpoints) {
		// make the list of elements
		elements = new LinkedList<Element>();
		for (Composition c : _endpoints) 
			for (Element e : c.getElements())
				if (!elements.contains(e))
					elements.add(e);
		
		// make the basis
		basis = new LinkedList<Vect>();
		for (Composition c : _endpoints) {
			List<Double> compFracs = new LinkedList<Double>();
			for (Element e : elements) 
				compFracs.add(c.getFractionalCompo(e));
			basis.add((new Vect(compFracs)).getSumNormalized(1.0));
		}
		
		endpoints = new LinkedList<Composition>(_endpoints);
	}
	
	/*
	// return a random composition from this space
	public Composition getRandomCompInSpace() {
		// get a linear combination of the basis vectors with random, positive coefficients (normalized)
		Vect linearComp = Vect.getNullVect(elements.size());
		
		for (Vect b : basis)
			linearComp = linearComp.plus(b.scalarMult(Math.random()));
		
		Map<Element, Integer> m = new HashMap<Element,Integer>();
		for (int i = 0; i < linearComp.getDimension(); i++)
			m.put(elements.get(i), linearComp.getCartesianComponents().get(i));
		
		return new Composition(m); 
	} */
	
	// return a random integer-coefficient composition from this space
	// TODO: not sure if this "max" thing is a great idea
/*	public Composition getRandomIntegerCompInSpace(int max, int minNumAtoms, int maxNumAtoms) {
		// get a linear combination of the basis vectors with random, positive coefficients 
		
		Map<Element, Integer> m = new HashMap<Element,Integer>();
		for (Composition c : endpoints) {
			int numOfThisComp = RandomNumbers.getUniformIntBetween(0, max);
			for (Element e : elements) {
				int numOfElem = Math.round(Math.round((c.getStoichiometricUnit().get(e)))) * numOfThisComp;
				if (m.keySet().contains(e))
					m.put(e, m.get(e) + numOfElem);
				else
					m.put(e, numOfElem);
			}
		}
		
		// divide out the gcd
		int gcd = Utility.gcd(m.values());
		for (Element e : m.keySet())
			m.put(e, m.get(e) / gcd);

		return new Composition(m, false);
	}*/
	
	public Composition getRandomIntegerCompInSpace(int minNumAtoms, int maxNumAtoms) {
		// get a linear combination of the basis vectors with random, positive coefficients 
		
		Map<Element, Double> m = new HashMap<Element,Double>();
		double sumAtoms = 0.0;
		for (Composition c : endpoints) {
			double fracOfThisComp = GAParameters.getParams().getRandom().nextDouble();
			for (Element e : elements) {
				double numOfElem = c.getStoichiometricUnit().get(e) * fracOfThisComp;
				sumAtoms += numOfElem;
				if (m.keySet().contains(e))
					m.put(e, m.get(e) + numOfElem);
				else
					m.put(e, numOfElem);
			}
		}
		
		int targetNumAtoms = RandomNumbers.getUniformIntBetween(minNumAtoms, maxNumAtoms);
		
		Map<Element, Integer> c = new HashMap<Element,Integer>();

		for (Element e : m.keySet()) 
			c.put(e, Math.round(Math.round(m.get(e) * targetNumAtoms / sumAtoms)));
		
		System.out.println("composition: " + c);

		return new Composition(c, false);
	}
	
	public boolean contains(Composition c) {
		// make sure c doesnt have extra elements
		for (Element ce : c.getElements())
			if (! this.getElements().contains(ce))
				return false;		
		
		// make a vector to represent the composition
		List<Double> fractComps = new LinkedList<Double>();
		for (Element e : this.getElements())
			fractComps.add(c.getFractionalCompo(e));
		Vect cVect1 = new Vect(fractComps);
		
		//	 use jama's least squares solver to solve for coefs of cVect wrt basis
		//       then take that linear combo of the basis and see if it differs from cVect
		//       by more than some threshold.
		
		Matrix a = null;
		try {
			 a = PseudoPDBuilder.doDecomposition(elements, endpoints, c);
		} catch (RuntimeException x) {
			return false;
		}
		for(Double d : a.getRowPackedCopy())
			if (d < - 0.00001)
				return false;
		
		return true;
		/*
		double[] coefs = PseudoPDBuilder.doDecomposition(elements, endpoints, c).getRowPackedCopy();
		Vect cVect2 = Vect.getNullVect(Constants.numDimensions);
		for (int i = 0; i < coefs.length; i++)
			cVect2 = cVect2.plus(coefs[i] * )
		
		Vect cVect2 = new Vect(coefs.getRowPackedCopy());
		
		return (cVect1.plus(cVect2.scalarMult(-1.0))).length() < 0.0001;
		*/
	}
	
	public List<Element> getElements() {
		return elements;
	}
	
	public String toString() {
		return "Elements " + elements + " with endpoints " + endpoints;
	}

}
