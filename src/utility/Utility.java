package utility;

import java.io.*;
import java.lang.reflect.Array;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import java.util.zip.GZIPInputStream;
import java.util.zip.GZIPOutputStream;
import java.util.*;

public class Utility {
	
	public static void writeStringToFile(String s, String path) {
		FileWriter fwriter = null;
		BufferedWriter writer = null;
		try {
			fwriter = new FileWriter(path);
			writer = new BufferedWriter(fwriter);
			writer.write(s);
			writer.flush();
			fwriter.close();
			writer.close();
		} catch (IOException x) {
			System.out.println("IOException in writeStringToFile(): " + x.getMessage());
		} finally {
			if (fwriter != null)
				try { fwriter.close(); } catch (IOException e) {  }//ignore
		}
	}
	
	  /**
	  * Fetch the entire contents of a text file, and return it in a String.
	  * This style of implementation does not throw Exceptions to the caller.
	  *
	  * from: http://www.javapractices.com/topic/TopicAction.do?Id=42
	  * 
	  * @param aFile is a file which already exists and can be read.
	  */
	  public static String readStringFromFile(String path) {
	    //...checks on aFile are elided
	    StringBuilder contents = new StringBuilder();
	    
	    File aFile = new File(path);
	    
	    try {
	      //use buffering, reading one line at a time
	      //FileReader always assumes default encoding is OK!
	      BufferedReader input =  new BufferedReader(new FileReader(aFile));
	      try {
	        String line = null; //not declared within while loop
	        /*
	        * readLine is a bit quirky :
	        * it returns the content of a line MINUS the newline.
	        * it returns null only for the END of the stream.
	        * it returns an empty String if two newlines appear in a row.
	        */
	        while (( line = input.readLine()) != null){
	          contents.append(line);
	          contents.append(System.getProperty("line.separator"));
	        }
	      } finally {
	        input.close();
	      }
	    }
	    catch (IOException ex){
	      ex.printStackTrace();
	    }
	    
	    return contents.toString();
	  }
	
	public static boolean stringContains(String superStr, String subStr) {
		Pattern pattern = Pattern.compile(subStr);
		Matcher energyMatcher = pattern.matcher(superStr);
		return energyMatcher.find();
	}
	
	public static int maxInt(Collection<Integer> nums) {

		int answer = Integer.MIN_VALUE;
		
		for (Integer i : nums)
			if (answer < i)
				answer = i;
		
		return answer;
	}
	
    /**
     * see Tim Mueller's ArrayUtils....
     * @param orig
     * @param obj
     * @return
     */

    public static <T> T[] appendArray(T[] orig, T obj){
        T[] appa =  (T[])Array.newInstance(orig.getClass().getComponentType(), orig.length + 1);
        System.arraycopy(orig, 0, appa, 0, orig.length);
        appa[orig.length] = obj;
        return appa;
    }
	
	// find greatest common divisor of a list of numbers
	// could obv be done better (euclid) if speed is important.
	public static int gcd(Collection<Integer> nums) {
		int max = (int) maxInt(nums);
		
		int maxDivisor = 1;
		
		for (int i = 1; i < max; i++) {
			boolean iDividesAllNums = true;
			for (int j : nums)
				if (j % i != 0)
					iDividesAllNums = false;
			if (iDividesAllNums)
				maxDivisor = i;
		}
			
		return maxDivisor;
	}
    // takes a list of objects and returns a list of all subsets of that object.
    // guarantees the larger subsets come first
    public static <T> List<List<T>> getPowerSet(List<T> list) {
            List<List<T>> result = new LinkedList<List<T>>();
            
            if (list.size() == 0) {
                    result.add(new LinkedList<T>());
                    return result;
            }
            
            T firstElem = list.get(0);
            List<T> listWOFirstElem = new LinkedList<T>(list);
            listWOFirstElem.remove(0);
            
            List<List<T>> subProblemResult = getPowerSet(listWOFirstElem);
            List<List<T>> subProblemResult2 = new LinkedList<List<T>>();
            for (List<T> L : subProblemResult) {
                    List<T> firstElemConsL = new LinkedList<T>(L);
                    firstElemConsL.add(0, firstElem);
                    subProblemResult2.add(firstElemConsL);
            }
            
            // combine them in order
            Iterator<List<T>> i = subProblemResult.iterator();
            Iterator<List<T>> j = subProblemResult2.iterator();

            List<T> iCurr = null;
            List<T> jCurr = null;
            if (i.hasNext())
                    iCurr = i.next();
            if (j.hasNext())
                    jCurr = j.next();
            while (iCurr != null && jCurr != null) {
                            if (iCurr.size() > jCurr.size()) {
                                    result.add(iCurr);
                                    if (i.hasNext())
                                            iCurr = i.next();
                                    else 
                                            iCurr = null;
                            } else {
                                    result.add(jCurr);
                                    if (j.hasNext())
                                            jCurr = j.next();
                                    else 
                                            jCurr = null;
                            }
            }
            if (jCurr != null)
                    result.add(jCurr);
            if (iCurr != null)
                    result.add(iCurr);
            while (j.hasNext())
                    result.add(j.next());
            while (i.hasNext())
                    result.add(i.next());
            
            return result;
    }
    
    public static double[] todouble(Double[] Ds) {
    	double[] result = new double[Ds.length];
    	for (int i = 0; i < Ds.length; i++)
    		result[i] = Ds[i];
    	return result;
    }
    
    public static double getPositiveFracPart(double d) {
    	long l = (long) d;
    	double frac = d - l;
    	if (frac < 0)
    		return frac + 1;
    	else
    		return frac;
    }
    
	// see http://java.sun.com/developer/technicalArticles/Programming/compression/
	public static void writeSerializable(Serializable s, String fName) {
      // serialize the objects 
		FileOutputStream fos = null;
		try {
	      fos = new FileOutputStream(fName);
	      GZIPOutputStream gz = new GZIPOutputStream(fos);
	      ObjectOutputStream oos = new ObjectOutputStream(gz);
	      oos.writeObject(s);
	      oos.flush();
	      oos.close();
		} catch (FileNotFoundException x) {
		  System.out.println("FileNotFoundException in writeSerializable:" + x.getMessage());
		  x.printStackTrace();
		} catch (IOException x) {
		  System.out.println("IOException in writeSerializable:" + x.getMessage());
		  x.printStackTrace();
		} finally {
			if (fos != null)
				try { fos.close(); } catch (IOException e) {  }//ignore
		}
	}
	
	public static Serializable readSerializable(String fName) {
		Serializable result = null;
		FileInputStream fis = null;
		try {
	      fis = new FileInputStream(fName);
	      GZIPInputStream gs = new GZIPInputStream(fis);
	      ObjectInputStream ois = new ObjectInputStream(gs);
	      result = (Serializable) ois.readObject();
	      ois.close();
		} catch (FileNotFoundException x) {
			  System.out.println(x.getMessage());
		} catch (IOException x) {
		  System.out.println(x.getMessage());
		} catch (ClassNotFoundException x) {
		  System.out.println(x.getMessage());
		} finally {
			if (fis != null)
				try { fis.close(); } catch (Exception x) { } // ignore
		}
	
		return result;
	}

}
