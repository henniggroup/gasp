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

package crystallography;

import ga.GAOut;
import ga.GAParameters;
import ga.StructureOrg;

import java.io.*;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.LinkedList;
import java.util.List;

import utility.Constants;
import utility.Vect;
import vasp.VaspOut;

public class Isotropy {
	
	//  in case program bombs, try again w/ different accuracy
	public static double accuracies[] = {.1,0.04,0.01,0.005,0.001,0.0005}; 
	
	private static String getFindsymInput(Cell c, double accuracy) {
		StringBuilder result = new StringBuilder();
		result.append("title \n");
		result.append(accuracy + "  accuracy\n");
		result.append("1  lattice vectors given as vectors\n");
		for (int i = 0; i < Constants.numDimensions; i++) {
			for (int j = 0; j < Constants.numDimensions; j++)
				result.append(c.getLatticeVectors().get(i).getCartesianComponents().get(j) + " ");
			result.append("\n");
		}
		result.append("1  form of primitive lattice vectors\n");
		result.append("1 0 0\n");
		result.append("0 1 0\n");
		result.append("0 0 1\n");
		result.append(c.getNumSites() + "\n");
		// seems to be important that the integers specifying the types start at 1.
		// findsym crashes otherwise when trying to free some memory
		for (int i = 0; i < c.getNumSites(); i++)
			result.append(c.getComposition().getElements().indexOf(c.getSite(i).getElement())+1 + " ");
		result.append("\n");
		for (int i = 0; i < c.getNumSites(); i++) {
			for (double d : c.getSite(i).getCoords().getComponentsWRTBasis(c.getLatticeVectors()))
				result.append(d + " ");
			result.append("\n");
		}
		
		return result.toString();
	}
	
	private static String runFindsym(String input) {
		StringBuilder outBuilder = new StringBuilder();

		Process p = null;
		try {
	        String line;
	      //  String env[] = new String[1];
	    //    env[0] = "ISODATA=" + findsymDir;
	     //   Process p = Runtime.getRuntime().exec(findsymDir + "/" + findsymBin, env);
	        p = Runtime.getRuntime().exec("callfindsym");
	        
	        OutputStreamWriter inw = new OutputStreamWriter(p.getOutputStream());
	        BufferedWriter in = new BufferedWriter(inw);
	        try {
	        	in.write(input);
	        	in.flush();
	        } finally {
		        in.close();
		        inw.close();
	        }
	        
	        InputStreamReader outw = new InputStreamReader(p.getInputStream());
	        BufferedReader out = new BufferedReader(outw);
	        try {
		        while ((line = out.readLine()) != null) {
		        	outBuilder.append(line + "\n");
		        }
	        } finally {
		        out.close(); 
		        outw.close();
	        }
	        
	        // flush the error stream. apparently if all the streams arent flushed
	        // the java runtime can leak file descriptors and it'll run out and
	        // crash eventually
	        InputStreamReader eoutw = new InputStreamReader(p.getErrorStream());
	        BufferedReader eout = new BufferedReader(eoutw);
	        try {
		        while ((line = eout.readLine()) != null) 
		        	;
	        } finally {
		        eout.close(); 
		        eoutw.close();
	        }
	        
	      } catch (Exception err) {
	        err.printStackTrace();
	      } finally {
	    	  try {
	    		  p.getErrorStream().close();
	    	  } catch (IOException e) {
	    		  e.printStackTrace();
	    	  }
	    	  try {
	    		  p.getOutputStream().close();
	    	  } catch (IOException e) {
	    		  e.printStackTrace();
	    	  }
	    	  try {
	    		  p.getOutputStream().close();
	    		  p.getInputStream().close();
	    	  } catch (IOException e) {
	    		  e.printStackTrace();
	    	  }
	      }

	      return outBuilder.toString();
	}
	
	private static int getSpaceGroupIntFromFSOutput(String output) {
		
		String lines[] = output.split("\\r?\\n");

		for (String line : lines) {
			if (line.startsWith("_symmetry_Int_Tables_number")) {
				return Integer.parseInt(line.split(" ")[1]);
			}
		}
		
		System.out.println("ERROR: couldn't find spacegroup in findsym output: " + output);
		return -1;
	}
	
	// NB: call Cell.getFSOOutput instead to use cached version
	public static String getFindsymOut(Cell cell) {
		
		// the idea here is to run findsym w/ a variety of different tolerance parameters
		// and use the result that finds the highest symmetry
		String output = null;
		int bestSG = 0;
		String bestOutput = null;
		for (int i = 0; i < accuracies.length; i++){
			output = runFindsym(getFindsymInput(cell, accuracies[i]));
			if (output.contains("bombed")) 
				continue;
			
			int sg = getSpaceGroupIntFromFSOutput(output);
			if (sg >= bestSG) {
				bestSG = sg;
				bestOutput = output;
			}
		}
		
		
		return bestOutput;
	}
	
	private static Cell parseWyckoffCell(String output, Cell origCell) throws Exception {
		// get lattice parameters and site types
		boolean nextLineIsLatticeParams = false;
		String latticeParms = null;
		boolean nextLineIsTypes = false;

		List<Integer> types = new LinkedList<Integer>();
		for (String line : output.split("\n")) {
		//	System.out.print (line);
			if (nextLineIsLatticeParams)
				latticeParms = line;
			nextLineIsLatticeParams = line.startsWith("Values of a,b,c,alpha,beta,gamma");
			if (line.startsWith("Position of each atom"))
				nextLineIsTypes = false;
			if (nextLineIsTypes) {
				String typeStrs[] = line.split(" +");
                // remove the first element of the list if it's an empty string (sometimes happens with older versions of findsym)
                if (typeStrs[0].equalsIgnoreCase("")) {
                    String trimmedTypeStrs[] = Arrays.copyOfRange(typeStrs, 1, typeStrs.length); 
                    typeStrs = trimmedTypeStrs;
                }
				for (int i = 0; i < typeStrs.length; i++) {								
					String typeGroup[] = typeStrs[i].split("\\*"); 						
					if (typeGroup.length == 2) {								
						for (int j = 1; j <= Integer.parseInt(typeGroup[0]); j++) { 	
							types.add(Integer.parseInt(typeGroup[1]));					
						}																
					} else {															
						types.add(Integer.parseInt(typeGroup[0]));						
					}																				
				}
			}
			if (!nextLineIsTypes)
				nextLineIsTypes = line.startsWith("Type of each atom");	
		}
		
		/*
		if (latticeParms == null) {
			System.out.println("ERROR in Isotropy.getWyckoffCell: findsym output probably malformed");
			System.out.println(cell);
			System.out.println(getFindsymInput(cell));
			System.out.println(output);
			return null;
		} */
		
		// Make the lattice parameters
		String lParms[] = latticeParms.split("  *");
		double a = Double.parseDouble(lParms[1]);
		double b = Double.parseDouble(lParms[2]);
		double c = Double.parseDouble(lParms[3]);
		double alpha = Double.parseDouble(lParms[4]);
		double beta = Double.parseDouble(lParms[5]);
		double gamma = Double.parseDouble(lParms[6]);
		List<Vect> basis = Cell.getVectorsfromLParamsDegrees(a,b,c,alpha,beta,gamma);		
		
		// get wyckoff positions
		boolean positionsStarted = false;
		List<Site> sites = new LinkedList<Site>();
		for (String line : output.split("\n")) {
			if (!positionsStarted && line.startsWith("Wyckoff position"))
				positionsStarted = true;
			if (positionsStarted && line.startsWith("------"))
				break;
			if (positionsStarted) {
				if (!line.startsWith("Wyckoff position")) {
					String sitesStrs[] = line.split("  *");
					int typeIndx = Integer.parseInt(sitesStrs[1]) - 1;
					double xf = Double.parseDouble(sitesStrs[2]);
					double yf = Double.parseDouble(sitesStrs[3]);
					double zf = Double.parseDouble(sitesStrs[4]);
					int type = types.get(typeIndx);
					sites.add(new Site(origCell.getComposition().getElements().get(type - 1), new Vect(xf,yf,zf,basis)));
				}
			}
		}
		 
		return new Cell(basis, sites);
	}
	
	public static Cell getWyckoffCell(Cell cell) {
		Cell answer = null;

		String output = cell.getFSOOutput();
		try {
			answer = parseWyckoffCell(output,cell);
		} catch (Exception x) {
			GAOut.out().stdout(x.getLocalizedMessage(), GAOut.INFO);
			GAOut.out().stdout("Warning: Exception in parsing output from findsym for all accuracies tried for cell and input:", GAOut.INFO);
			GAOut.out().stdout(cell.toString(), GAOut.INFO);
			GAOut.out().stdout(getFindsymInput(cell, accuracies[0]), GAOut.INFO);
			GAOut.out().stdout((new Exception()).getStackTrace().toString(), GAOut.INFO);
			
			answer = cell;
		}
		
		return answer;
	}
	
	public static void main(String args[]) {
		//Cell c = VaspOut.getPOSCAR(args[0]);
	//	Cell c = VaspOut.getPOSCAR("/home/brevard/testing/S2Mo.POSCAR");
	//	Cell d = VaspOut.getPOSCAR("/home/bcr48/GA/2D/garun_C11/34257.POSCAR");
	//	Cell e = VaspOut.getPOSCAR("/home/bcr48/GA/2D/garun_C11/35690.POSCAR");
	//	Cell f = VaspOut.getPOSCAR("/home/bcr48/GA/2D/garun_C11/40495.POSCAR");
	//	Cell g = VaspOut.getPOSCAR("/home/bcr48/GA/2D/garun_C11/40914.POSCAR");
		
		//Cell c = StructureOrg.parseCif(new File("/home/wtipton/cifs/2.cif"));
	//	Cell a = Cell.parseCif(new File("/home/wtipton/10631.cif"));
	//	Cell b = Cell.parseCif(new File("/home/wtipton/11929.cif"));
		
	//	System.out.println(c.matchesCell(d, 0.1, 0.5, 0.3));
	//	System.out.println(c.matchesCell(e, 0.1, 0.5, 0.3));
	//	System.out.println(c.matchesCell(f, 0.1, 0.5, 0.3));
	//	System.out.println(d.matchesCell(e, 0.1, 0.5, 0.3));
	//	System.out.println(d.matchesCell(f, 0.1, 0.5, 0.3));
	//	System.out.println(e.matchesCell(f, 0.1, 0.5, 0.3));

		//String fsInput = utility.Utility.readStringFromFile("/home/wtipton/fstest");
		
	//	c = VaspOut.getPOSCAR("/home/wtipton/cifs/POSCAR");
		
	//	String fsoutput = getFindsymOut(c);
	//	System.out.println(fsoutput);	
	//	try {
	//		System.out.println(parseWyckoffCell(fsoutput, c));
	//	} catch (Exception x) {
	//		System.out.println(x);
	//	}
	//	System.out.println(d.getFSOOutput());
	//	System.out.println(e.getFSOOutput()); 
	//	System.out.println(f.getFSOOutput());
	//	System.out.println(g.getFSOOutput());
//		System.out.println(getFindsymOutput(b.getNigliReducedCell()));	
	//	b.getNigliReducedCell().writeCIF("/home/wtipton/11134.nig.cif");
	//	System.out.println(getWyckoffCell(c));
		
		//System.out.println(runFindsym(fsInput));
	}

}
