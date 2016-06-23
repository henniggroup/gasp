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

package crystallography;

import java.io.Serializable;

import utility.Vect;
import chemistry.Element;

/*
 * TODO:
 *  - String description is temporary.  Need an actual way to specify an element located
 *    on a Site, not to mention partial occupancies or molecules.
 */

public class Site implements Serializable{
	
	static final long serialVersionUID = 1l;
	
	private Vect coords;
	private Element element;
	private boolean relax; // indicated whether this site should be allowed to relax in the energy calculation
	
	public Site(Element _element, Vect _v) {
		element = _element;
		coords = _v;
		relax = true; // default to true - we generally want to allow atoms to relax
	}
	
	public Site(Element _element, Vect _v, boolean _relax) {
		element = _element;
		coords = _v;
		relax = _relax;
	}
	
	public Vect getCoords() {
		return coords;
	}
	
	// toString using fractional coords
	public String toString() {
		StringBuilder output = new StringBuilder();
		
		output.append(element.getSymbol() + " ");
		for (int i = 0; i < 3; i++)
			output.append(coords.getCartesianComponents().get(i) + " ");
		
		return output.toString();
	}
	
	/* TODO: revise this if Sites ever support partial occupancies or non-element decorations */
	public Element getElement() {
		return element;
	}
	
	public boolean getRelax() {
		return relax;
	}
}
