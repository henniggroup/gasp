package ga;

import java.util.*;

import chemistry.Composition;
import chemistry.CompositionSpace;
import chemistry.Element;

import crystallography.Cell;

public class SupercellSOCreator implements StructureOrgCreator {
	
	private RandomSOCreator creator;
	private int a, b, c; // supercell dimensions
	
	public SupercellSOCreator(List<String> args) {
		
		// parse the options passed in
		if (args == null || args.size() < 11)
			GAParameters.usage("Not enough parameters given to SupercellSOCreator", true);
		
		// parse supercell dimensions
		a = Integer.parseInt(args.get(0));
		b = Integer.parseInt(args.get(1));
		c = Integer.parseInt(args.get(2));
		
		// parse underlying cell constraints
		double maxll = Double.parseDouble(args.get(3));
		double minll = Double.parseDouble(args.get(4));
		double maxla = Double.parseDouble(args.get(5));
		double minla = Double.parseDouble(args.get(6));
		double maxh = Double.parseDouble(args.get(7));
		int maxna = Integer.parseInt(args.get(8));
		int minna = Integer.parseInt(args.get(9));

		creator = new RandomSOCreator(args.subList(10, args.size()), maxll, minll, maxla, minla, maxh, maxna, minna);

	}

	@Override
	public StructureOrg makeOrganism(Generation g) {
		
		Cell smallcell = creator.makeRandomOrg().getCell();
		
		System.out.println(smallcell);
		
		List<List<Integer>> coefficients = new ArrayList<List<Integer>>();
		List<Integer> l1 = new ArrayList<Integer>();
		l1.add(a); l1.add(0); l1.add(0);
		List<Integer> l2 = new ArrayList<Integer>();
		l2.add(0); l2.add(b); l2.add(0);
		List<Integer> l3 = new ArrayList<Integer>();
		l3.add(0); l3.add(0); l3.add(c);
		
		coefficients.add(l1);
		coefficients.add(l2);
		coefficients.add(l3);
		
		return new StructureOrg(Cell.getSupercell(smallcell, coefficients));
	}
	
	// for testing
	public static void main(String args[]) {
		GAParameters p = GAParameters.getParams();
		p.setMinInteratomicDistance(1.0);
		p.setMinNumAtoms(2);
		p.setMaxNumAtoms(40);
		Composition comp = new Composition(Element.getElemFromZ(6));
		List<Composition> comps = new ArrayList<Composition>();
		comps.add(comp);
		CompositionSpace compSpace = new CompositionSpace(comps);
		p.setCompSpace(compSpace);
		List<String> a = new ArrayList<String>();
		a.add("3");
		a.add("3");
		a.add("1");
		a.add("5");
		a.add("4");
		a.add("120");
		a.add("60");
		a.add("2");
		a.add("7");
		a.add("9");
		a.add("givenVol");
		a.add("20");
		SupercellSOCreator soc = new SupercellSOCreator(a);
		Generation g = new Structures();
		System.out.println(soc.makeOrganism(g));
	}

}
