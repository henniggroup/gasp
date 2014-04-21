package pdvisual;

import java.util.Iterator;
import java.util.List;

public class PDConvexable implements Convexable {
	
	int dim;
	double[][] allVertices;
	IPDPlottable pdd;
	
	public PDConvexable(IPDPlottable _pdd) {
		// make the allVertices list
		// we use the elements order in this!
		
		pdd = _pdd;
		
		List<Double[]> coords = pdd.getCompositionFractions();
		//List<Double> energies = pdd.getFormEnergiesPerAtom();
		List<Double> energies = pdd.getEnergiesPerAtom(); 

             
		/**
		 * Thu Jul 24 14:36:26 EDT 2008 Chris BUGFIX
		 * chem potential information was not used ! now updated in PDData
		 */
		allVertices = new double[coords.size()][pdd.getDimension()];
		Iterator<Double[]> iterv = coords.iterator();
		for(int i = 0; i < allVertices.length ; i++) {
			Double[] v = iterv.next();
			for (int j = 0; j < v.length; j++)
				allVertices[i][j] = v[j].doubleValue();
		}
		for(int i = 0; i < energies.size(); i++)
			allVertices[i][pdd.getDimension() - 1] = energies.get(i);
	}
	
	// getAllverticesforConvexHull should return an array of double
	// each row is a vertex and each column the coord of the vertex.
	public double[][] getAllVerticesforConvexHull() {
		return allVertices;
	}
	
	// returns the dimension in which you're building your convex hull
	public int getdimension() {
		return pdd.getDimension();
	}
}