package castep;

import java.io.File;
import ga.StructureOrg;

import crystallography.Cell;

public class CastepTest {
	
	
	public static void main(String args[]) {
		
		Cell c = Cell.parseCell("/home/wtipton/castepTests/test.cif", "cif");
		
		String ceArgs[] = {"true", "0.1", "5.0", "/home/wtipton/castepTests/castepParam", "H", "/home/wtipton/castepTests/H_00PBE.usp", "O", "home/wtipton/castepTests/O_00PBE.usp"};
		
		CastepEnergy in = new CastepEnergy(ceArgs);
		
		in.getEnergy(new StructureOrg(c));
	}

}
