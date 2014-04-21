<<<<<<< HEAD
=======
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

>>>>>>> 0e3189c40547bbd59ea42c4f91890d7511fb7797
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
