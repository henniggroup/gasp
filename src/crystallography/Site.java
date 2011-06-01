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
	
	public Site(Element _element, Vect _v) {
		element = _element;
		coords = _v;
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
}
