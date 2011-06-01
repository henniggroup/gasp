package gulp;

import java.io.Serializable;
import java.util.ArrayList;
import java.util.List;
import chemistry.Element;

public class GulpPotential implements Serializable {
	private static final long serialVersionUID = 1l;

	private String type;
	private List<Element> species;
	private List<Double> parameters;
	private List<Double> parameters_min;
	private List<Double> parameters_max;
	private List<Integer> opt_orders;
	private double cutoff_min, cutoff_max;
	
	public GulpPotential(String _type, List<Element> _species, List<Double> _parameters, List<Double> _parameters_min, List<Double> _parameters_max, List<Integer> _optOrders, double _cutoff_min, double _cutoff_max) {
		setParams(_type, _species, _parameters, _parameters_min, _parameters_max, _optOrders, _cutoff_min, _cutoff_max);
	}
	
	public GulpPotential(GulpPotential oldPot, List<Double> newParams, List<Boolean> fitParams) {
		if (GulpSurrogate.getNumFittedParams(fitParams) != newParams.size())
			System.out.println("ERROR: oldPot.getNumFittedParams() != newParams.size() in GulpPotential.");
		
		newParams = new ArrayList<Double>(newParams);
		List<Double> params = new ArrayList<Double>();
		for (int i = 0; i < fitParams.size(); i++)
			if (fitParams.get(i))
				params.add(newParams.remove(0));
			else
				params.add(oldPot.parameters.get(i));
		setParams(oldPot.type, oldPot.species, params, oldPot.parameters_min, oldPot.parameters_max, oldPot.opt_orders, oldPot.cutoff_min, oldPot.cutoff_max);
	}
	
	private void setParams(String _type, List<Element> _species, List<Double> _parameters, List<Double> _parameters_min, List<Double> _parameters_max, List<Integer> _optOrders, double _cutoff_min, double _cutoff_max) {
		type = _type;
		species = _species;
		parameters = _parameters;
		parameters_min = _parameters_min;
		parameters_max = _parameters_max;
		opt_orders = _optOrders;
		cutoff_min = _cutoff_min;
		cutoff_max = _cutoff_max;
	}
	
	public boolean paramIsPinned(int i) {
		return parameters.get(i) <= parameters_min.get(i) || parameters.get(i) >= parameters_max.get(i);
	}
	
	public String toGulpPotlStr(List<Boolean> fitParams) {
		StringBuilder result = new StringBuilder();
		String nl = System.getProperty("line.separator");
		
		result.append(type + nl);
		for (Element e : species)
			result.append(e.getSymbol() + " core ");
		for (int i = 0; i < parameters.size(); i++) 
			result.append(parameters.get(i) + " ");
		
		result.append(cutoff_min + " " + cutoff_max + " ");
		if (fitParams != null)
			for (Boolean b : fitParams)
				result.append(b ? "1 " : "0 ");
		
		return result.toString();
	}
	
	public int getNumParameters() {
		return parameters.size();
	}
	
	public void enforceConstraints() {
		for (int i = 0; i < parameters.size(); i++) 
			parameters.set(i, Math.min(Math.max(parameters.get(i), parameters_min.get(i)), parameters_max.get(i)));
	}
	
	public boolean isPinned(int i) {
		return (parameters.get(i) >= parameters_max.get(i) || parameters.get(i) <= parameters_min.get(i));
	}
	
	public boolean violatesConstraints() {
		boolean result = false;
		for (int i = 0; i < parameters.size(); i++)
			if (parameters.get(i) > parameters_max.get(i) || parameters.get(i) < parameters_min.get(i))
				result = true;
		return result;
	}
	
	public List<Integer> getOptIndexes() {
		return new ArrayList<Integer>(opt_orders);
	}
}
