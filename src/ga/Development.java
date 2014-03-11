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

package ga;

// The Development interface is implemented by methods which oversee
// the "growing up" of an organism by implementing doDevelop.  If the
// organism is obviously unfit (hard constraints), doDevelop may return
// false, and the organism should not be considered any further.  
// doDevelop may modify the Organism (e.g. structure relaxation), but
// should not modify the Generation.

public interface Development {
	
	public Boolean doDevelop(Generation g, Organism o);

}
