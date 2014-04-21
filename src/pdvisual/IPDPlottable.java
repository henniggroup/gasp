package pdvisual;

import java.util.List;

import chemistry.Composition;
import chemistry.Element;


public interface IPDPlottable {
	
	// returns the dimension of the plottable pd object
	public int getDimension();
	public int getNumEntries();
    public List<Element> getElements();
	public List<Integer> getIndxUnstableEntries();
	public List<Double> getEnergiesPerAtom(); 
	public List<Double> getFormEnergiesPerAtom();
	public List<List<Integer>> getIndxFacets();
	public List<String> getVertexLabels();
	public List<String> getAxisLabels();
	public boolean isPurePhase(int i);
	// returns the indices of the "basis" entries, e.g. the elements in a normal phase diagram
	public List<Integer> getBasisEntryIndxs();
	// get (REAL) entrie with same composition as entry i
	public List<Integer> getPhasesWSameCompAs(int i);
	public Composition getComposition(int i);
	public List<Composition> getCompositions();
	public List<Double[]> getCompositionFractions();
	

}
