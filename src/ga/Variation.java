<<<<<<< HEAD
/* Genetic algorithm for crystal structure prediction.  Will Tipton.  Ceder Lab at MIT. Summer 2007 */
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

package ga;

import java.io.Serializable;

// A Variation operation takes a pair of generations and uses them to create an organism to
// be added to the offspring generation.  Neither should be modified.  Application of Variations,
// along with Promotion, are how a new Generation is created from a parent Generation.

// When implementing, keep in mind that Structure, Organism, etc. mostly don't have copy constructors.

public interface Variation extends Serializable {
	
	public Organism doVariation(Generation parents, Generation offspring, Selection sel);
	
}
