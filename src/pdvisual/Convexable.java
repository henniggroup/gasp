package pdvisual;
/**
 * This is the interface you have to implement for
 * objects from which you want to compute a convex hull.
 * 
 * @author geoffroy
 * Date: 22 October 2007
 */


public interface Convexable {
	/**
	 * getAllverticesforConvexHull should return an array of double
	 * each line is a vertice and each column the coord of the vertice.
	 * @return
	 */
	public double[][] getAllVerticesforConvexHull();
	/**
	 * returns the dimension in which you're building your convex hull
	 * 
	 */
	public int getdimension();
	

}
