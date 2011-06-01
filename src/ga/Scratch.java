/* Genetic algorithm for crystal structure prediction.  Will Tipton.  Ceder Lab at MIT. Summer 2007 */

package ga;

import java.io.File;

import vasp.VaspOut;

import crystallography.Cell;
import crystallography.Site;

public class Scratch {

	/**
	 * @param args
	 */
	public static void main(String[] args) {
		
		Cell c = VaspOut.getPOSCAR("/home/wtipton/cifs/POSCAR_SC");
		String largs[] = {"/home/wtipton/projects/cacd6/potlfile"};
		LammpsEnergy leng = new LammpsEnergy(largs);
		
		System.out.println(leng.getEnergy(new StructureOrg(c)));
		
//		System.loadLibrary("example");
//	    System.out.println(example.fact(5));
//	     System.out.println(example.get_time());

	/*	System.out.println("Current working directory:");
		System.out.println(System.getProperty("user.dir"));
		System.out.println(new File(".").getAbsolutePath());

		System.out.println("Changing working directory...");
		System.setProperty("user.dir", "C:\\Java");

		System.out.println("Current working directory:");
		System.out.println(System.getProperty("user.dir"));
		System.out.println(new File(".").getAbsolutePath());
		*/
	/*	Cell structure = StructureOrg.parseCif(new File("281.relaxed.cif"));

		for (int i = 0; i < structure.getNumSites(); i++) {
			// it's no good if there are any other atoms in the minimum radius sphere
			Site[] sites = structure.getAtomsInSphereSorted(structure.getSite(i).getCoords(), 1.3);
			if (sites.length > 1) {
				System.out.println("Organism failed minimum interatomic distance constraint.");
			}
		}	
*/

		//java.sql.Connection c = com.cmcweb.db.postgresql.PostgresConnection.getPsqlConnection("jdbc:postgresql://atlov.mit.edu:5432/magichat", "anubhavj", "zK86R50r");
		//Cell str = StructureExtractor.getStructureFromICSDCode(c, "85276-ICSD");
		//StructureOrg s =  new StructureOrg(str);	
		
	//	System.out.println(s);
		
	//	DBQueryGUI db = new DBQueryGUI();
	//	db.setConnection("jdbc:postgresql://atlov.mit.edu:5432/magichat", "anubhavj", "zK86R50r");
		
	/*	Structure s1 = new CIFParser("/home/wtipton/icsd_4266.cif", System.out).getEntry(0).getStructure();
		Structure s2 = new CIFParser("/home/wtipton/icsd_412735.cif", System.out).getEntry(0).getStructure();
		
		FitterData fd = StructureFitter.getFitterDataInstance();
		fd.setAtomicMisfitTol(1.0);
		fd.setCellMisfitTol(0.001);
		new StructureFitter(s1, s2, fd, null);
		System.out.println(fd.getcellmisfit());
		System.out.println(fd.getatomicmisfit());
	*/	
		//GAUtils.writeStringToFile(s.getCIF(), new File("/home/wtipton/yba2cu4o8/85276.cif"), false);
		
		
		

	}

}
