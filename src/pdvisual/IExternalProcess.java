/*
 * com.cmcweb.extcallers
 * Created on May 22, 2008
 * by ccfisch
 */
package pdvisual;

import java.io.BufferedReader;
import java.io.FileOutputStream;
import java.io.PrintStream;

/**
 * This is supposed to be a generic interface for calling and retrieving the
 * output of an external process (e.g., qhull, gulp, ...)
 * 
 * the idea here is to try to localize all OS depending calling features rather
 * than scatter them throughout the rest of the java code. 
 * 
 * @author ccfisch, Chris Fischer <ccfisch@mit.edu>
 * Date: May 22, 2008
 * 
 */
public interface IExternalProcess {
    
    public void setInputFilename(String s);
    
    /**
     * this is an output stream (that you can write to) for the input
     * file of the process 
     * (use if your process takes a generic binary file as input)
     * 
     * @return
     */
    public FileOutputStream getInputAsStream(); 
    
    /**
     * this is a print stream corresponding to the input file of the process
     * (i.e., use if your process takes an ascii encoded file as input)
     * 
     * @return
     */
    public PrintStream getInputAsPrintStream(); 
    
    /**
     * execute process ! 
     */
    public void runProcess();
    
    /**
     * the output of the process (could be standard output, or perhaps a file
     * that is created by the process). use for post-processing. 
     * 
     * @return
     */
    public BufferedReader getProcessOutput();
 
}
