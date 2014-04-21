package pdvisual;

public class PDAnalyzerFactory {

	public static IPDAnalyzer getIPDAnalyzer(IPDPlottable p) {
		if (p.getClass() == PDData.class) {
			return new PDAnalyzer((PDData)p);
		} else if (p.getClass() == PseudoPDData.class) {
			return new PseudoPDAnalyzer((PseudoPDData)p);
		} else {
			return null;
		}
	}
}
