/* Genetic algorithm for crystal structure prediction.  Will Tipton.  Ceder Lab at MIT. Summer 2007 */

package ga;

import java.io.*;
import java.util.*;
import java.lang.reflect.Array;

import chemistry.Element;

// GAUtils holds some general utility functions 

public class GAUtils {
	
	// returns a String containing a newline character
	public static String newline() {
		return System.getProperty("line.separator");
	}
	
	// returns a subarray of an array starting at start and of length length
	public static  <T> T[] subArray(T[] src, int start, int length) {
		// the warning is harmless:
		T[] dest = (T[])Array.newInstance(src.getClass().getComponentType(),length);
		
		System.arraycopy(src, start, dest, 0, length);
		return dest;
	}
	
	// returns a subarray of an array of strings starting at start and going until the end
	public static <T> T[] subArray(T[] src, int start) {
		int length = src.length - start;
		return subArray(src, start, length);
	}

	// returns the nth largest double in an Arraylist (the largest is for n=1).  
	// modifies the order of the list.  could be done better.
	public static double doubleSelect(List<Double> a, int n) {
		Collections.sort(a);
		Collections.reverse(a);
		return ((Double)(a.get(n-1))).doubleValue();
	}
	
	public static <T,S> void putIntoMap(Map<T,S> m, T[] keys, S[] values) {
		if (keys.length != values.length) 
			System.err.println("GAUtils.putIntoMap: Error: keys.length != values.length");
		
		for (int i = 0; i < keys.length; i++)
			m.put(keys[i], values[i]);
	}
	
	// modifies its input:
	public static Double[] normalizeArray(Double[] a) {
		double sum = 0;
		for (int i = 0; i < a.length; i++)
			sum += a[i];
		for (int i = 0; i < a.length; i++)
			a[i] /= sum;
		return a;
	}
	
	// writes a string to a file, optionally appending
	public static void writeStringToFile(String s, File f, Boolean append) {
		if (GAParameters.getParams().getDryRun())
			return;
		try {
			BufferedWriter out = new BufferedWriter(new FileWriter(f, append));
			out.append(s);
			out.close();
		} catch (IOException x) {
			System.out.println("GARecord.writeStringToFile(): " + x.getMessage());
		}
	}
	
	// reads a file into a string
	public static String readStringFromFile(File f) {
		StringBuilder result = new StringBuilder();
		String s = null;
		
		try {
			BufferedReader b = new BufferedReader(new FileReader(f));
            while ((s = b.readLine()) != null) 
            	result.append(s + newline());				
		} catch (FileNotFoundException x) {
			System.out.println("GAUtils.readStringFromFile: " + x.getMessage());
			System.exit(1);
		} catch (IOException x) {
			System.out.println("GAUtils.readStringFromFile: " + x.getMessage());
			System.exit(1);
		}

		return result.toString();
	}
	
	public static int strArrayContains(String[] a, String s) {
		for (int i = 0; i < a.length; i++)
			if (s.equalsIgnoreCase(a[i]))
				return i;
		return -1;
	}
	
	// take an array of angles in radians and returns those angles in degrees
	public static double[] anglesToDegrees(double[] anglesRad) {
		double[] anglesDeg = new double[anglesRad.length];
		
		for (int i = 0; i < anglesRad.length; i++)
			anglesDeg[i] = anglesRad[i] * 180.0 / Math.PI;
		
		return anglesDeg;
	}
	
	// heh only use this for small numbers
	public static double renormZeroOne(double d) {
		if (d > 1)
			return renormZeroOne(d - 1);
		if (d < 0)
			return renormZeroOne(d + 1);
		return d;
	}
	
	// checks whether or not three numbers satisfy the Euclidean triangle inequality
	public static Boolean satisfiesTriangleInequality(double a, double b, double c) {
		return (a + b >= c && b + c >= a && a + c >= b);
	}
	
	// takes in an array of strings of the form Mg-O and returns list of Sets of String,
	// where each set contains the two strings (e.g. Mg and O)
	public static List<Set<String>> parsePairs(List<String> args) {
		GAParameters params = GAParameters.getParams();
		// parse the pairs of atoms which can be exchanged
		ArrayList<Set<String>> result = new ArrayList<Set<String>>();
		for (int i = 0; i < args.size(); i++) {
			// get the pair
			String[] symbols = args.get(i).split("-");
			// make sure we've got a pair of species that are in the system
			if (symbols.length != 2)
				GAParameters.usage("Malformed species pair " + args.get(i), true);
			if (!params.getCompSpace().getElements().contains(Element.getElemFromSymbol(symbols[0]))
					|| !params.getCompSpace().getElements().contains(Element.getElemFromSymbol(symbols[1])))
				GAParameters.usage("Species not in the system: " + args.get(i), true);
			// TODO: how should we handle above contingency in PD obj fcn case (when it happens,
			// to, say, select an elemental phase)?
			// store the pair
			Set<String> newPair = new HashSet<String>();
			newPair.add(symbols[0]);
			newPair.add(symbols[1]);
			result.add(newPair);
		}

		return result;
	}
	
	public static Double sumDoubles(Double[] ds) {
		Double sum = 0.0;
		for (int i = 0; i < ds.length; i++)
			sum += ds[i];
		return sum;
	}
	
	// ignores case
	public static Boolean listHasString(Collection<String> a, String s) {
		Iterator<String> i = a.iterator();
		String t;
		while (i.hasNext()) {
			t = i.next();
			if (t.equalsIgnoreCase(s))
				return true;
		}
		return false;
	}
	
	// This function returns the sign of x, 1 when x=0;
	public static int SignNoZero(double x){
	  if (x == 0)
		  System.out.println("poo");
	  if (x >= 0) 
	    return 1;
	  return -1;
	}
	
	public static void printArray(Object[] a) {
		for (int i = 0; i < a.length; i++)
			System.out.println(a[i]);
	}
	
	// rounds r to the nearest d
	public static double roundToNearest(double r, double d) {
		if (d == 0)
			return r;		
		double temp = Math.round(Math.round(r/d));
		return temp*d;
	}
	
	// recursively deletes a directory
	public static void deleteDirectory(File dir) throws IOException
	{
	   File[] files = dir.listFiles();

	   for(File f : files) {
	      if(f.isDirectory())
	         deleteDirectory(f);

	      f.delete();
	   }
	}
	
	// for testing
	public static void main(String[] args) {
		List<Double> bob = new ArrayList<Double>();
		bob.add(1.0);
		bob.add(2.0);
		bob.add(3.0);
		bob.add(4.0);
		bob.add(5.0);
		for (Double i : bob)
			System.out.println(i);
		System.out.println("selected: " + doubleSelect(bob,1));
		for (Double i : bob)
			System.out.println(i);
	}
}
