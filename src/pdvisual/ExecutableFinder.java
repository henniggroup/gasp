/*
 * com.cmcweb.extcallers
 * Created on May 22, 2008
 * by ccfisch
 */
package pdvisual;

import java.io.File;

/**
 * 
 * simple class to find executables in a user's path
 * 
 * @author ccfisch, Chris Fischer <ccfisch@mit.edu>
 * Date: May 22, 2008
 * 
 */
public class ExecutableFinder {

    
    /**
     * returns the full path to an executable with given name if found, 
     * otherwise null
     * 
     * @param name
     * @return
     */
    public static String findExecutable(String name){

        String path = null;
        if(System.getProperty("os.name").toUpperCase().indexOf("WIND") != -1){
            path = System.getenv("path"); 
        }else {
            path = System.getenv("PATH");
        }
        String[] dirs = path.split(File.pathSeparator); 
        for(String pathDir : dirs){
            File test = new File(pathDir, name);
            if(test.isFile()){
                return test.getPath(); 
            }
        }
        return null; 
    }
    
    /**
     * @param args
     */
    public static void main(String[] args) {
        if(args == null || args.length == 0 || args[0].equals("-h")){
            System.out.println("ExecutableFinder [-h] [name]");
            System.out.println("   searches through PATH to find 'name' ");
            System.exit(0); 
        }
        String res = findExecutable(args[0]);
        if(res != null){
            System.out.println("found "+args[0]+" ... '"+res+"'");
        }else{
            System.out.println("could not find "+args[0]);
        }
    }

}
