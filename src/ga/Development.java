/* Genetic algorithm for crystal structure prediction.  Will Tipton.  Ceder Lab at MIT. Summer 2007 */

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
