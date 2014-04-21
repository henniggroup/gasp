package pdvisual;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.PrintStream;
import java.util.LinkedList;
import java.util.List;
import java.util.StringTokenizer;
import java.util.*;

/**
 * this class creates a convexhull.
 * It needs an object implementing the interface Convexable.
 * It uses qhull to compute the convex hull (http://www.qhull.org/).
 * Put the qconvex binary in the pathToQhull directory.
 * You also need to create a script in this directory: callqconvexfacets
 * The script should contain only one line with: qconvex Fx <$1
 * 
 * Thu May 29 12:49:08 EDT 2008 <ccfisch@mit.edu>
 *   Updated this class to use the class com.cmcweb.extcallers.QHullCaller since
 *   the if(WINDOWS) then X else if (MAC) then Y, etc stuff was starting to 
 *   propagate into several different parts of the code. 
 *
 * 
 * @author geoffroy
 * Date: 22 October 2007
 */



public class Convexhull {
	
	Convexable allVertices;
	int dim;
	/* facets of the convex hull
	the line are facets and columns are the vertices in the facet
	there is a nb of facet lines and a dim column */
	int[][] facets;
	double volume;

	public Convexhull(Convexable allVertices, boolean removeTopFaces)
	{
		this.allVertices = allVertices;
		this.dim = allVertices.getdimension();

		QHullCaller qhcaller;
		if (removeTopFaces)
			qhcaller = new QHullCaller(QHullCaller.QHullExe.QCVXFACETS, "Fn FS PD" + (dim - 1) + " ");  
		else
			qhcaller = new QHullCaller(QHullCaller.QHullExe.QCVXFACETS, "Fn FS");  
		

		//format the input for qconvex
		writeInputFile(qhcaller.getInputAsPrintStream()); 
		
		//run qconvex on this input file
		qhcaller.runProcess();

		//fetch the output and process 
		processOutput(qhcaller.getProcessOutput()); 
		
		//clean up open filehandles and process stuff
		qhcaller.cleanUp(); 
	}
	
	/**
	 * process output of qconvex
	 * @param br
	 */
	private void processOutput(BufferedReader br){
	    
        String output;
        try{ 
        	output = br.readLine();
       // 	System.out.println(output+"***@@");
        	int numFacets = Integer.parseInt(output);
            facets = new int[numFacets][dim];
            // parse facets
            for (int i = 0; i < numFacets; i++)
            {   
            	output = br.readLine();
                //ccfisch: is this used anywhere ? I'm commenting out now but maybe Geoffroy can respond if needed
                //System.out.println(output+"***@@");
            
                StringTokenizer token = new StringTokenizer(output);
                for (int j = 0; j < dim; j++)
                    facets[i][j] = Integer.parseInt(token.nextToken());
            }
            // parse adjacency info:
            // actually, we're going to compute this ourselves later since we have to
            // remove some facets first
            
            // skip adjacency list info
            numFacets = Integer.parseInt(br.readLine());
            for (int i = 0; i < numFacets; i++) 
            	output = br.readLine();
            
    /*    	numFacets = Integer.parseInt(br.readLine());
            adjacencyList = new LinkedList<List<Integer>>();
            // parse facets
            for (int i = 0; i < numFacets; i++)
            {   
            	output = br.readLine();
                //ccfisch: is this used anywhere ? I'm commenting out now but maybe Geoffroy can respond if needed
                //System.out.println(output+"***@@");
            
                StringTokenizer token = new StringTokenizer(output);
                int numNeighbors = Integer.parseInt(token.nextToken());
                for (int j = 0; j < numNeighbors; j++) {
                	int neighbor = Integer.parseInt(token.nextToken());
                	if (i < neighbor) {
                		List<Integer> newPair = new LinkedList<Integer>();
                		newPair.add(i);
                		newPair.add(neighbor);
                		adjacencyList.add(newPair);
                	}
                }
            }
            
      */
            
            // read volume info
        	output = br.readLine();
        	output = br.readLine();
            StringTokenizer token = new StringTokenizer(output);
            token.nextToken();token.nextToken();
            volume = Double.parseDouble(token.nextToken());

        }catch(IOException ioe){
            System.err.println("IOException in Convexhull.processOutput:" + ioe.getLocalizedMessage());            
        }
	    
	}
	
	/**
	 * format input for qconvex 
	 * @param ps
	 */
	private void writeInputFile(PrintStream ps){

	    double[][] coord = this.allVertices.getAllVerticesforConvexHull(); 
	    ps.println(dim);
        ps.println(coord.length);
        
        for(int i=0;i<coord.length;i++)
        {
            for(int j=0;j<coord[i].length;j++)
                ps.print(coord[i][j]+" ");
            ps.print("\n");
        }
        ps.close();
	    
       
	}
	
	public double getVolume() {
		return volume;
	}

	public int[][] getFacets()
	{
		return this.facets;
	}
	
	public List<List<Integer>> getFaLists() {
	    List<List<Integer>> ret = new LinkedList<List<Integer>>();
	    for(int[] facet : facets) {
	        List<Integer> fset = new LinkedList<Integer>();
	        for(int vtex : facet)
	        	fset.add(vtex); 
	        ret.add(fset);
	    }
	    return ret; 
	}
	
}
