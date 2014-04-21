/* Genetic algorithm for crystal structure prediction.  Will Tipton.  Ceder Lab at MIT. Summer 2007 */

package ga;

// StructureOrgCreator is an interface to objects which create StructureOrgs,
// generally for use in the initial population.

public interface StructureOrgCreator {
	
	public StructureOrg makeOrganism(Generation g);

}
