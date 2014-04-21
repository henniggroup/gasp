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

package chemistry;

import java.io.Serializable;

/*
 * IMMUTABLE!
 * 
 */

public class Element implements Comparable, Serializable {
	
	static final long serialVersionUID = 1l;
	
	private int Z;
	
	private static Element[] elements;
	// from http://www.science.co.il/PTelements.asp
	// 		format:
	// 		Z mass name symbol density group ionization-energy
	static int ZEntryLoc = 0;
	static int massEntryLoc = 1;
	static int nameEntryLoc = 2;
	static int symbolEntryLoc = 3;
	static int densityEntryLoc = 4;
	static String[][] elementData = {{"1","1.0079","Hydrogen","H","0.00008988","1","13.5984"},
		{"2","4.0026","Helium","He","0.0001785","18","24.5874"},
		{"3","6.9410","Lithium","Li","0.530","1","5.3917"},
		{"4","9.0122","Beryllium","Be","1.850","2","9.3227"},
		{"5","10.8110","Boron","B","2.340","13","8.2980"},
		{"6","12.0107","Carbon","C","2.260","14","11.2603"},
		{"7","14.0067","Nitrogen","N","0.0012506","15","14.5341"},
		{"8","15.9994","Oxygen","O","0.001429","16","13.6181"},
		{"9","18.9984","Fluorine","F","0.001696","17","17.4228"},
		{"10","20.1797","Neon","Ne","0.0008999","18","21.5645"},
		{"11","22.9897","Sodium","Na","0.970","1","5.1391"},
		{"12","24.3050","Magnesium","Mg","1.740","2","7.6462"},
		{"13","26.9815","Aluminum","Al","2.700","13","5.9858"},
		{"14","28.0855","Silicon","Si","2.330","14","8.1517"},
		{"15","30.9738","Phosphorus","P","1.820","15","10.4867"},
		{"16","32.0650","Sulfur","S","2.070","16","10.3600"},
		{"17","35.4530","Chlorine","Cl","0.003214","17","12.9676"},
		{"18","39.9480","Argon","Ar","0.0017837","18","15.7596"},
		{"19","39.0983","Potassium","K","0.860","1","4.3407"},
		{"20","40.0780","Calcium","Ca","1.550","2","6.1132"},
		{"21","44.9559","Scandium","Sc","2.990","3","6.5615"},
		{"22","47.8670","Titanium","Ti","4.540","4","6.8281"},
		{"23","50.9415","Vanadium","V","6.110","5","6.7462"},
		{"24","51.9961","Chromium","Cr","7.190","6","6.7665"},
		{"25","54.9380","Manganese","Mn","7.430","7","7.4340"},
		{"26","55.8450","Iron","Fe","7.870","8","7.9024"},
		{"27","58.9332","Cobalt","Co","8.900","9","7.8810"},
		{"28","58.6934","Nickel","Ni","8.900","10","7.6398"},
		{"29","63.5460","Copper","Cu","8.960","11","7.7264"},
		{"30","65.3900","Zinc","Zn","7.130","12","9.3942"},
		{"31","69.7230","Gallium","Ga","5.910","13","5.9993"},
		{"32","72.6400","Germanium","Ge","5.320","14","7.8994"},
		{"33","74.9216","Arsenic","As","5.720","15","9.7886"},
		{"34","78.9600","Selenium","Se","4.790","16","9.7524"},
		{"35","79.9040","Bromine","Br","3.120","17","11.8138"},
		{"36","83.8000","Krypton","Kr","0.003733","18","13.9996"},
		{"37","85.4678","Rubidium","Rb","1.630","1","4.1771"},
		{"38","87.6200","Strontium","Sr","2.540","2","5.6949"},
		{"39","88.9059","Yttrium","Y","4.470","3","6.2173"},
		{"40","91.2240","Zirconium","Zr","6.510","4","6.6339"},
		{"41","92.9064","Niobium","Nb","8.570","5","6.7589"},
		{"42","95.9400","Molybdenum","Mo","10.220","6","7.0924"},
		{"43","98.0000","Technetium","Tc","11.500","7","7.2800"},
		{"44","101.0700","Ruthenium","Ru","12.370","8","7.3605"},
		{"45","102.9055","Rhodium","Rh","12.410","9","7.4589"},
		{"46","106.4200","Palladium","Pd","12.020","10","8.3369"},
		{"47","107.8682","Silver","Ag","10.500","11","7.5762"},
		{"48","112.4110","Cadmium","Cd","8.650","12","8.9938"},
		{"49","114.8180","Indium","In","7.310","13","5.7864"},
		{"50","118.7100","Tin","Sn","7.310","14","7.3439"},
		{"51","121.7600","Antimony","Sb","6.680","15","8.6084"},
		{"52","127.6000","Tellurium","Te","6.240","16","9.0096"},
		{"53","126.9045","Iodine","I","4.930","17","10.4513"},
		{"54","131.2930","Xenon","Xe","0.005887","18","12.1298"},
		{"55","132.9055","Cesium","Cs","1.870","1","3.8939"},
		{"56","137.3270","Barium","Ba","3.590","2","5.2117"},
		{"57","138.9055","Lanthanum","La","6.150","3","5.5769"},
		{"58","140.1160","Cerium","Ce","6.770","101","5.5387"},
		{"59","140.9077","Praseodymium","Pr","6.770","101","5.4730"},
		{"60","144.2400","Neodymium","Nd","7.010","101","5.5250"},
		{"61","145.0000","Promethium","Pm","7.300","101","5.5820"},
		{"62","150.3600","Samarium","Sm","7.520","101","5.6437"},
		{"63","151.9640","Europium","Eu","5.240","101","5.6704"},
		{"64","157.2500","Gadolinium","Gd","7.900","101","6.1501"},
		{"65","158.9253","Terbium","Tb","8.230","101","5.8638"},
		{"66","162.5000","Dysprosium","Dy","8.550","101","5.9389"},
		{"67","164.9303","Holmium","Ho","8.800","101","6.0215"},
		{"68","167.2590","Erbium","Er","9.070","101","6.1077"},
		{"69","168.9342","Thulium","Tm","9.320","101","6.1843"},
		{"70","173.0400","Ytterbium","Yb","6.900","101","6.2542"},
		{"71","174.9670","Lutetium","Lu","9.840","101","5.4259"},
		{"72","178.4900","Hafnium","Hf","13.310","4","6.8251"},
		{"73","180.9479","Tantalum","Ta","16.650","5","7.5496"},
		{"74","183.8400","Tungsten","W","19.350","6","7.8640"},
		{"75","186.2070","Rhenium","Re","21.040","7","7.8335"},
		{"76","190.2300","Osmium","Os","22.600","8","8.4382"},
		{"77","192.2170","Iridium","Ir","22.400","9","8.9670"},
		{"78","195.0780","Platinum","Pt","21.450","10","8.9587"},
		{"79","196.9665","Gold","Au","19.320","11","9.2255"},
		{"80","200.5900","Mercury","Hg","13.550","12","10.4375"},
		{"81","204.3833","Thallium","Tl","11.850","13","6.1082"},
		{"82","207.2000","Lead","Pb","11.350","14","7.4167"},
		{"83","208.9804","Bismuth","Bi","9.750","15","7.2856"},
		{"84","209.0000","Polonium","Po","9.300","16","8.4170"},
		{"85","210.0000","Astatine","At","0.000","17","9.3000"},
		{"86","222.0000","Radon","Rn","0.00973","18","10.7485"},
		{"87","223.0000","Francium","Fr","0.000","1","4.0727"},
		{"88","226.0000","Radium","Ra","5.500","2","5.2784"},
		{"89","227.0000","Actinium","Ac","10.070","3","5.1700"},
		{"90","232.0381","Thorium","Th","11.720","102","6.3067"},
		{"91","231.0359","Protactinium","Pa","15.400","102","5.8900"},
		{"92","238.0289","Uranium","U","18.950","102","6.1941"},
		{"93","237.0000","Neptunium","Np","20.200","102","6.2657"},
		{"94","244.0000","Plutonium","Pu","19.840","102","6.0262"},
		{"95","243.0000","Americium","Am","13.670","102","5.9738"},
		{"96","247.0000","Curium","Cm","13.500","102","5.9915"},
		{"97","247.0000","Berkelium","Bk","14.780","102","6.1979"},
		{"98","251.0000","Californium","Cf","15.100","102","6.2817"},
		{"99","252.0000","Einsteinium","Es","0.000","102","6.4200"},
		{"100","257.0000","Fermium","Fm","0.000","102","6.5000"},
		{"101","258.0000","Mendelevium","Md","0.000","102","6.5800"},
		{"102","259.0000","Nobelium","No","0.000","102","6.6500"},
		{"103","262.0000","Lawrencium","Lr","0.000","102","4.9000"},
		{"104","261.0000","Rutherfordium","Rf","0.000","4","0.0000"},
		{"105","262.0000","Dubnium","Db","0.000","5","0.0000"},
		{"106","266.0000","Seaborgium","Sg","0.000","6","0.0000"},
		{"107","264.0000","Bohrium","Bh","0.000","7","0.0000"},
		{"108","277.0000","Hassium","Hs","0.000","8","0.0000"},
		{"109","268.0000","Meitnerium","Mt","0.000","9","0.0000"}};
	
	private Element (int _Z) {
		Z = _Z;
	}
	
	private Element() {
		// no public constructors
	}
	
	private String[] getEntryFromZ(int _Z) {
		for (String[] entry : elementData) {
			if (Integer.parseInt(entry[ZEntryLoc]) == _Z)
				return entry;
		}
		throw new RuntimeException("Element.getEntryFromZ() tried to get nonexistent Z: " + _Z);
	}
	
	private static int getZFromSymbol(String symb) {
		for (String[] entry : elementData) {
			if (symb.compareToIgnoreCase(entry[symbolEntryLoc]) == 0)
				return Integer.parseInt(entry[ZEntryLoc]);
		}
		throw new RuntimeException("Element.getZFromSymbol() tried to get nonexistent symbol: " + symb);
	}
	
	private static int getZFromName(String name) {
		for (String[] entry : elementData) {
			if (name.compareToIgnoreCase(entry[nameEntryLoc]) == 0)
				return Integer.parseInt(entry[ZEntryLoc]);
		}
		throw new RuntimeException("Element.getZFromName() tried to get nonexistent name: " + name);
	}
	
	public String getSymbol() {
		return getEntryFromZ(Z)[symbolEntryLoc];
	}
	
	public double getAtomicMass() {
		return Double.parseDouble(getEntryFromZ(Z)[massEntryLoc]);
	}
	
	public double getDensity() {
		return Double.parseDouble(getEntryFromZ(Z)[densityEntryLoc]);
	}
	
	public String getName() {
		return getEntryFromZ(Z)[nameEntryLoc];
	}
	
	public int getZ() {
		return Z;
	}
	
	public String toString() {
		return getEntryFromZ(Z)[nameEntryLoc];
	}
	
	private static int getEntryNumFromZ(int _Z) {
		return _Z - 1;
	}
	
	private static void initializeElementsArray() {
		elements = new Element[elementData.length];
		for (int i = 0; i < elementData.length; i++)
			elements[i] = new Element(Integer.parseInt(elementData[i][ZEntryLoc]));
	}
	
	public static Element getElemFromSymbol(String sym) {
		if (elements == null)
			initializeElementsArray();
		
		return elements[getEntryNumFromZ(getZFromSymbol(sym))];
	}
	
	public static Element getElemFromZ(int Z) {
		if (elements == null)
			initializeElementsArray();
		
		return elements[getEntryNumFromZ(Z)];
	}
	
	public static Element getElemFromName(String name) {
		if (elements == null)
			initializeElementsArray();
		
		return elements[getEntryNumFromZ(getZFromName(name))];
	}
	
	@Override
	public int compareTo(Object o) {
		return Z - ((Element)o).Z;
	}
	
	public boolean equals(Object aThat) {
		if ( this == aThat ) 
			return true;
		if ( !(aThat instanceof Element) ) 
			return false;
		Element that = (Element)aThat;
		return that.Z == this.Z;
	}
	
	public int hashCode() {
		return Z;
	}

}
