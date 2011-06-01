package utility;

import java.io.Serializable;

public class Pair<T, S> implements Serializable
{
  public static final long serialVersionUID = 1l;
  
  public Pair(T f, S s)
  { 
    first = f;
    second = s;   
  }

  public T getFirst()
  {
    return first;
  }

  public S getSecond() 
  {
    return second;
  }

  public String toString()
  { 
    return "(" + first.toString() + ", " + second.toString() + ")"; 
  }

  private T first;
  private S second;
}
