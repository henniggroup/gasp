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
