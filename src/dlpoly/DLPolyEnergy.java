package dlpoly;

import ga.Energy;
import ga.GAParameters;
import ga.GAUtils;
import ga.StructureOrg;
import ga.UnitsSOCreator;
import gulp.GulpEnergy;

import java.io.*;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.List;
import java.util.Scanner;
import java.util.StringTokenizer;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import utility.Constants;
import utility.Utility;
import utility.Vect;

import crystallography.Cell;
import crystallography.Site;

import org.openbabel.*;

public class DLPolyEnergy implements Energy {

	public DLPolyEnergy(String[] args) {		
		
		if (args == null || args.length < 1)
			GAParameters.usage("Not enough parameters given to AvogadroEnergy", true);
		
		// read in the location of DLPoly executable
		String loc = args[0];
	}
	
	public double getEnergy(StructureOrg c) {
		
		runDLPoly(c);

		//TODO: this obviously will need to change
		return 0.0;
	}
	
	private static void writeConfig(StructureOrg c) {
		GAParameters params = GAParameters.getParams();
		
//		String outdir = params.getTempDirName() + "/CONFIG";
		String outdir = "/home/skw57/polyfiles/CONFIG";
		String title = "Structure " + c.getID() + "\n";
		// '0' -- coordinates only, '3' -- parallelepiped boundary conditions
		String record = "0         3" + "\n";
		
		
		List<Vect> latVects = c.getCell().getLatticeVectors();
		List<Site> sites = c.getCell().getSites();
		
		String lattice = "";
		for (Vect v: latVects) {
			List<Double> xyz = v.getCartesianComponents();
			for (int i=0; i<Constants.numDimensions; i++) {
				lattice = lattice + xyz.get(i) + "     ";
			}
			lattice = lattice + "\n";
		}
		
		String atoms = "";
		int n = 1;
		for (Site s: sites) {
			List<Double> coords = s.getCoords().getCartesianComponents();
			atoms = atoms + s.getElement().getSymbol() + "         " + n + "\n";
			for (int k=0; k<Constants.numDimensions; k++) {
				atoms = atoms + coords.get(k) + "      ";
			}
			atoms = atoms + "\n";
			n++;
		}
		String total = title + record + lattice + atoms;
		Utility.writeStringToFile(total, outdir);
	}
	
	//TODO: just need one CONTROL file for the whole run
	public static void writeControl() {
		GAParameters params = GAParameters.getParams();
		String newline = GAUtils.newline();
		try {
//			File f = new File(params.getTempDirName(),"CONTROL");
			File f = new File("/home/skw57/polyfiles/","CONTROL");

			// write the file
			BufferedWriter out = new BufferedWriter(new FileWriter(f));
			out.write(params.getRunTitle() + newline);		
			out.write("finish");
			out.close();
		} catch (IOException e) {
				System.out.println("DLPolyEnergy IOException: " + e.getMessage());
		}
	}
	
	//TODO: might be a problem running in parallel (again with the #units deal)
	//TODO: actually seems like a huge pain in the ass, combining UnitsSOCreator #s with potential (GULP-esque) information
	private static void writeField(StructureOrg c) {
		GAParameters params = GAParameters.getParams();
		int difUnits = UnitsSOCreator.getDifUnits();
		int[] numUnits = UnitsSOCreator.getTargets();
		int[] numAtoms = UnitsSOCreator.getNumAtoms();
		List<Site> sites = c.getCell().getSites();
		
//		String outdir = params.getTempDirName() + "/FIELD";
		String outdir = "/home/skw57/polyfiles/FIELD";
		String title = "Structure " + c.getID() + "\n";
		String units = "UNITS eV \n";
		
		String num = "MOLECULES " + difUnits + "\n";
		
		String molecules = "";
		int loc = 0;
		for (int k=0; k<difUnits; k++) {
			molecules = molecules + "mol" + (k+1) + "\n";
			molecules = molecules + "NUMMOLS " + numUnits[k] + "\n";
			molecules = molecules + "ATOMS " + numAtoms[k] + "\n";
			for (int l=0; l<numAtoms[k]; l++) {
				molecules = molecules + "    " + sites.get(loc+l).getElement().getSymbol() + "     " + sites.get(loc+l).getElement().getAtomicMass() + "\n";
			}
			molecules = molecules + "FINISH \n";
			loc = loc + numAtoms[k]*numUnits[k];
		}
		
		
		String total = title + units + num + molecules + "CLOSE";
		Utility.writeStringToFile(total, outdir);
	}
	
	public static void runDLPoly(StructureOrg c) {
		// Write input files
		writeConfig(c);
		writeControl();
		writeField(c);
		
		BufferedReader stdInput = null;
		BufferedReader stdError = null;
		try {
			// Execute DL_Poly
			//TODO: needs to be pointed towards the proper directory
			Process p = Runtime.getRuntime().exec("calldlpoly");

			stdInput = new BufferedReader(new InputStreamReader(
					p.getInputStream()));
			stdError = new BufferedReader(new InputStreamReader(
					p.getErrorStream()));
			
		} catch (IOException e) {
			System.out.println("IOException in DLPolyEnergy.runDLPoly: " + e.getMessage());
			System.exit(-1);
		} finally {
			if (stdInput != null) 
				try{ stdInput.close(); } catch (Exception x) { } //ignore
			if (stdError != null) 
				try{ stdError.close(); } catch (Exception x) { } //ignore
		}
	}


	//testing
/*	public static void main(String[] args) {
		String[] argsIn = {"/home/skw57/Downloads/dl_poly_4.02/execute/"};
		DLPolyEnergy bob = new DLPolyEnergy(argsIn);
		Cell c = Cell.parseCif(new File("/home/skw57/2.cif"));
		
		bob.getEnergy(new StructureOrg(c));
		
	//	String output = GulpEnergy.runGULP("mno2_poly.gin");
	//	System.out.println(output);
	//	System.out.println(GulpEnergy.parseFinalEnergy(output, bob.cautious));
	//	System.out.println(output.contains("failed"));
	}
*/	
}