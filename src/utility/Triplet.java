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
package utility;

import java.io.Serializable;

public class Triplet<A,B,C> implements Serializable {
    static final long serialVersionUID = 1;

    protected A first;
    protected B second;
    protected C third;

    public Triplet(A a, B b, C c) {
        first = a;
        second = b;
        third = c;
    }

    public A getFirst() {
        return first;
    }

    public B getSecond() {
        return second;
    }

    public C getThird() {
        return third;

    }
}
