package crystallography;

import java.io.File;
import java.util.*;

import utility.ArgumentParser;
import utility.Vect;
import vasp.*;
import utility.*;

public class SupercellOptimizer {
	
	/* Returns a list of doubles which represent linear combinations
	 * of the input cell vectors of a supercell that, given input
	 * parameters, maximizes distance between any two parts of the cell
	 * with its image.
	 * Ported from Supercell utility written in C by Richard Hennig
	 */
	public static List<List<Integer>> getOptimalSupercell(Cell cell, Boolean deb, int N, int v, double rmax, boolean verbose) {
		  int[] sc = new int[9];            /* Supercell definition */
		  int[] L = new int[9];             /* Supercell with largest rws */
		  double[] x = new double[12];
		  double r, dmin, V, v0, dtot;
		  double[] c = new double[3]; 
		  double[] d = new double[3];
		  double[] a = new double[9];          /* Lattice vectors of primitive cell */

		  List<Vect> latticeVectors = cell.getLatticeVectors();
		  int counter = 0;
		  for (Vect vec : latticeVectors) {
			  for (Double comp : vec.getCartesianComponents()) {
				  a[counter] = comp;
				  counter++;
			  }
		  }

		  v0 = cell.getVolume();

		  int[] l = new int[3*(2*N+1)*(2*N+1)*(2*N+1)];

		  /* Generate list of lattice vectors shorter than rmax*/
		  if (deb) System.out.printf ("Generate list of lattice vectors shorter than %f\n", rmax);
		  rmax *= rmax;
		  int Nlist = 0;
		  for (sc[0]=-N; sc[0]<=N; sc[0]++) {
		    for (sc[1]=-N; sc[1]<=N; sc[1]++) {
		      for (sc[2]=-N; sc[2]<=N; sc[2]++) {
			x[0] = sc[0]*a[0] + sc[1]*a[3]+sc[2]*a[6];
			x[1] = sc[0]*a[1] + sc[1]*a[4]+sc[2]*a[7];
			x[2] = sc[0]*a[2] + sc[1]*a[5]+sc[2]*a[8];
			r = x[0]*x[0]+x[1]*x[1]+x[2]*x[2];
			if (r < rmax) {
			  /* Add lattice vector to list */
			  l[3*Nlist+0] = sc[0];
			  l[3*Nlist+1] = sc[1];
			  l[3*Nlist+2] = sc[2];
			  Nlist ++;
			}
		      }
		    }
		  }
		  if (deb) System.out.printf("Found %d lattice vectors shorter than %f.\n", Nlist, Math.sqrt(rmax));

		  /* Loop over all possible lattice vectors shorter than rmax and calculate Rws */
		  if (deb) System.out.printf("Loop over all possible lattice vectors shorter than rmax and calculate Rws\n");

		  dtot = 0.0;
		  for (int i=0; i<Nlist-2; i++) {
		    for (int j=i+1; j<Nlist-1; j++) {
		      for (int k=j+1; k<Nlist; k++) {
			/* Set lattice vectors */
			x[0] = l[3*i+0]*a[0] + l[3*i+1]*a[3]+l[3*i+2]*a[6];
			x[1] = l[3*i+0]*a[1] + l[3*i+1]*a[4]+l[3*i+2]*a[7];
			x[2] = l[3*i+0]*a[2] + l[3*i+1]*a[5]+l[3*i+2]*a[8];

			x[3] = l[3*j+0]*a[0] + l[3*j+1]*a[3]+l[3*j+2]*a[6];
			x[4] = l[3*j+0]*a[1] + l[3*j+1]*a[4]+l[3*j+2]*a[7];
			x[5] = l[3*j+0]*a[2] + l[3*j+1]*a[5]+l[3*j+2]*a[8];

			x[6] = l[3*k+0]*a[0] + l[3*k+1]*a[3]+l[3*k+2]*a[6];
			x[7] = l[3*k+0]*a[1] + l[3*k+1]*a[4]+l[3*k+2]*a[7];
			x[8] = l[3*k+0]*a[2] + l[3*k+1]*a[5]+l[3*k+2]*a[8];

		        x[9]  = Math.sqrt(x[0]*x[0] + x[1]*x[1] + x[2]*x[2]);
		        x[10] = Math.sqrt(x[3]*x[3] + x[4]*x[4] + x[5]*x[5]);
		        x[11] = Math.sqrt(x[6]*x[6] + x[7]*x[7] + x[8]*x[8]);

			/* Calculate Volume */
			V = x[0]*x[4]*x[8] + x[1]*x[5]*x[6] + x[2]*x[3]*x[7] 
			   -x[2]*x[4]*x[6] - x[1]*x[3]*x[8] - x[0]*x[5]*x[7];
			
			if (V/v0> 0.999*v && V/v0<1.001*v) {
			  /* Calculate Rws */
			
			  c[0] = 0.5*(x[0]+x[3]+x[6]);  /* Center of cell */
			  c[1] = 0.5*(x[1]+x[4]+x[7]);
			  c[2] = 0.5*(x[2]+x[5]+x[8]);
			  if (deb) System.out.printf("Cell center= (%12.6f %12.6f %12.6f)\n", c[0], c[1], c[2]);

			  d[0] = Math.abs(0.5*V/x[10]/x[11]);
			  d[1] = Math.abs(0.5*V/x[9]/x[11]);
			  d[2] = Math.abs(0.5*V/x[9]/x[10]);

			  dmin = d[0];
			  if (d[1]<dmin) dmin = d[1];
			  if (d[2]<dmin) dmin = d[2];

			  if (dmin>=dtot) {
			    dtot = dmin;
			    L[0]=l[3*i+0];
			    L[1]=l[3*i+1];
			    L[2]=l[3*i+2];
			    L[3]=l[3*j+0];
			    L[4]=l[3*j+1];
			    L[5]=l[3*j+2];
			    L[6]=l[3*k+0];
			    L[7]=l[3*k+1];
			    L[8]=l[3*k+2];
			    if (deb) System.out.printf("  N = (%d %d %d)\n      (%d %d %d)\n      (%d %d %d)\n  rmin = %f\n",
				   L[0], L[1], L[2], L[3], L[4], L[5], L[6], L[7], L[8], dtot);

			  }

			  if (deb) System.out.printf("For supercell (%d %d %d), (%d %d %d), (%d %d %d): rmin = %f\n",
			    l[3*i+0], l[3*i+1], l[3*i+2], l[3*j+0], l[3*j+1], l[3*j+2], l[3*k+0], l[3*k+1], l[3*k+2], dmin);
			}
		      }
		    }
		  }
		  if (verbose) {
			  System.out.printf("Found:\n");
			  System.out.printf("  rmin = %f\n  dmin = %f\n", dtot, 2*dtot);
		  }

		  List<List<Integer>> result = new LinkedList<List<Integer>>();
		  counter = 0;
		  for (int i = 0; i < Constants.numDimensions; i++) {
			  List<Integer> li = new LinkedList<Integer>();
			  for (int j = 0; j < Constants.numDimensions; j++) {
				  li.add(L[counter]);
				  counter++;
			  }
			  result.add(li);
		  }
		  return result;
	}
	
	/* Default argument values */
	private static Boolean debug = false;
	private static int n = 4;
	private static int v = 1;
	private static double rmax = 10.0;
	private static String inPoscar = null;
	private static String outFile = null;
	private static Boolean rotate = false;
	private static Boolean writeCart = false;
	private static Boolean niggliReduce = false;
	
	public static void usage() {
	      System.out.println("Usage: supercell [OPTIONS] --I <inFile>\n" +
			"  --d		Print debugging info\n" +
			"  --h		Print this help info and exit\n" +
			"  --I <file>	Input file in POSCAR format\n" +
			"  --O <file>	Write an output file\n" +
			"  --n <nmax>	Maximum value for supercell lattice entries [" + n + "]\n" +
			"  --r <rmax>	Maximum length for lattice vectors [" + rmax + "]\n" +
			"  --rot    	Rotate output cell into principle directions  \n" +
			"  --nig        Perform Niggli cell reduction on final cell \n" +
			"  --cart    	Write Cartesian (as opposed to fractional) coords in POSCAR \n" +
			"  --v <N>	Supercell volume in multiples of input cell [" + v + "]\n" +
			"  --c		Write output file in CIF format. (POSCAR otherwise.)\n" +
			"  --s <a11> <a12> <a13> <a21> <a22> <a23> <a31> <a32> <a33> \n" +
			"      Skip search for optimal vector combinations and use matrix: \n" +
			"        <a11> <a12> <a13> \n" +
			"        <a21> <a22> <a23> \n" +
			"        <a31> <a32> <a33> \n" +
			"      If no matrix given, assume the identity. \n"

		);
	}

	public static void main(String[] args) {
		
		/* Get input arguments */
		ArgumentParser aparser = new ArgumentParser(args);
		if (aparser.hasOption("h") || !aparser.hasOption("I")) { /* Have to have an input POSCAR */
			usage();
			System.exit(0);
		}
		inPoscar = aparser.getArgument("I");
		if (aparser.hasOption("d"))
			debug = true;
		if (aparser.hasOption("rot"))
			rotate = true;
		if (aparser.hasOption("nig"))
			niggliReduce = true;
		if (aparser.hasOption("cart"))
			writeCart = true;
		if (aparser.hasOption("O"))
			outFile = aparser.getArgument("O");
		if (aparser.hasOption("n"))
			n = Integer.parseInt(aparser.getArgument("n"));
		if (aparser.hasOption("v"))
			v = Integer.parseInt(aparser.getArgument("v"));
		if (aparser.hasOption("r"))
			rmax = Double.parseDouble(aparser.getArgument("r"));
		/* small sanity check */
		if (aparser.hasOption("v") && aparser.hasOption("s"))
			System.out.println("Probably shouldn't have both options -v and -s.");
		
		/* read input POSCAR into crystal */
		Cell inputCell = VaspOut.getPOSCAR(inPoscar);
		
		/* Print input cell parameters */
		System.out.printf("Primitive cell\n");
		System.out.println(inputCell.toStringJustVectors());
		
		/* Get supercell parameters */
		List<List<Integer>> cellCombos = null;
		Cell outputCell = null;
		if (!aparser.hasOption("s")) /* Get optimal supercell parameters */
			cellCombos = getOptimalSupercell(inputCell, debug, n , v, rmax, true);
		else if (aparser.hasArguments("s")) { /* Get supercell parameters from input */
			List<String> inputMatrixEntries = aparser.getArguments("s");
			if (inputMatrixEntries.size() != 9) {
				System.out.println("Input matrix wants 9 entries.");
				return;
			}
			cellCombos = new LinkedList<List<Integer>>();
			Iterator<String> entry = inputMatrixEntries.iterator();
			for (int i = 0; i < 3; i++) {
				List<Integer> vector = new LinkedList<Integer>();
				for (int j = 0; j < 3; j++) {
					vector.add(Integer.parseInt(entry.next()));
				}
				cellCombos.add(vector);
			}
		}
		
		if (cellCombos != null) { /* make the supercell */
			outputCell = Cell.getSupercell(inputCell, cellCombos);
		} else {
			outputCell = inputCell;
		}
		
		// maybe rotate the cell into principle directions
		if (rotate)
			outputCell = outputCell.getCellRotatedIntoPrincDirs();
		
		// maybe 
		if (niggliReduce)
			outputCell = outputCell.getNigliReducedCell();
		
		// Print output
		if (cellCombos != null) {
			System.out.printf("Supercell\n");
			System.out.printf("  N = (%d %d %d)\n      (%d %d %d)\n      (%d %d %d)\n",
				 cellCombos.get(0).get(0),cellCombos.get(0).get(1),cellCombos.get(0).get(2),
				 cellCombos.get(1).get(0),cellCombos.get(1).get(1),cellCombos.get(1).get(2),
				 cellCombos.get(2).get(0),cellCombos.get(2).get(1),cellCombos.get(2).get(2));
			System.out.println(outputCell.toStringJustVectors());
		}	
		
		// maybe write an output file
		if (outFile != null) {
		
			/* check to make sure the output file doesn't already exist */
			if ((new File(outFile)).exists()) {
				System.out.println("Error: Output file " + outFile + " exists.  Not writing output.");
				return;
			}
			
			/* write output poscar or cif */
			if (aparser.hasOption("c")) {
				outputCell.writeCIF(outFile);
				System.out.println("Wrote output CIF to: " + outFile + "\n");
			} else {
				VaspIn vaspin = new VaspIn(outputCell, "Autogenerated Supercell", null, null, null);
				vaspin.writePoscar(outFile, writeCart);
				System.out.println("Wrote output POSCAR to: " + outFile + "\n");
			}
			
			/* sanity check for good measure */
			if (aparser.hasOption("v") && v * inputCell.getBasisSize() != outputCell.getBasisSize()) {
				throw new RuntimeException("SupercellOptimizer: v * inputCell.getBasisSize() != outputCell.getBasisSize()");
			}
		}
	}
}
