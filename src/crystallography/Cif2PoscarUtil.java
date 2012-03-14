package crystallography;

import java.io.File;

import vasp.VaspIn;
import crystallography.Cell;

public class Cif2PoscarUtil {
	
	public static void main(String args[]) {
		if (args.length < 2) {
			System.out.println("Usage: cif2poscar file.cif out.POSCAR");
			System.exit(0);
		}
			
		(new VaspIn(Cell.parseCif(new File(args[0])), null, null, null)).writePoscar(args[1], false);
	}

}
