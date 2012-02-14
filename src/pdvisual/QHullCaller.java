/*
 * com.cmcweb.extcallers
 * Created on May 22, 2008
 * by ccfisch
 */
package pdvisual;


import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.PrintStream;


/**
 * this class will run qhull binaries with using a temporary file and
 * give you the output.
 * 
 * NOTES: 
 * on unix systems this assumes the following programs are available and can be
 * found in your path: bash, qconvex, qvoronoi 

 * on windows systems you need to have "qvoronoi.exe" and "qconvex.exe" in your 
 * path along with access to "cmd.exe"
 * 
 * Utilizes java facilities for creating temporary files and deleting them on exit
 * for example @see http://www.exampledepot.com/egs/java.io/CreateTempFile.html
 * 
 * @TODO anubhav and charles, when you have the chance please test.  
 * 
 * @author ccfisch, Chris Fischer <ccfisch@mit.edu>
 * Date: May 22, 2008
 * 
 */
public class QHullCaller implements IExternalProcess {

    public enum QHullExe {
        QVORONOI("qvoronoi.exe","qvoronoi", " Fv "),
        QVERTEX("qvoronoi.exe", "qvoronoi", " p "),
        QCVXFACETS("qconvex.exe", "qconvex", " i Qt ");
        
        QHullExe(String exew, String exen, String args){
            exe_win = exew;
            exe_nix = exen;
            qh_args = args; 
        }
        
        private final String exe_win; //exe for windows
        private final String exe_nix; //exe for linux
        private final String qh_args; //arguments passed to qhull programs 
    }
    
    private static String DEFAULT_INFILE = "inputqhull";

    private QHullExe m_qhexe; 
    private String m_infileName = DEFAULT_INFILE;
    private String m_exeToCall;
    private String m_extraArgs = ""; 
    
    private FileOutputStream m_inputfos;
    private PrintStream m_inputps;
    private File m_inputfile; 
    
    private BufferedReader m_output; 
    private Process m_p; 
    
    /**
     * you select the qhull program you want to run with the QHullExe enum
     * it is up to you to format the input properly with either @see getInputAsStream
     * or @see getInputAsPrintStream !!
     * 
     * @param q, otherargs
     */
    public QHullCaller(QHullExe q, String otherargs){
        this.m_qhexe = q; 
        if (otherargs != null)
        	this.m_extraArgs = otherargs; 
        setExecutable(); 
    }
    
    private void setExecutable(){
        
        if(System.getProperty("os.name").toUpperCase().indexOf("WIN") != -1){
            m_exeToCall = ExecutableFinder.findExecutable(this.m_qhexe.exe_win);
            if(m_exeToCall == null){
                throw new RuntimeException("Unable to find the executable '"+
                       this.m_qhexe.exe_win+"' in your path ! For this code \n"
                       +"to work it will need to be somewhere in your path shown here: "+
                       System.getenv("path"));
            }
        }else{
            m_exeToCall = ExecutableFinder.findExecutable(this.m_qhexe.exe_nix);
            if(m_exeToCall == null){
                throw new RuntimeException("Unable to find the executable '"+
                       this.m_qhexe.exe_nix+"' in your path ! For this code \n"
                       +"to work it will need to be somewhere in your path shown here: "+
                       System.getenv("PATH"));
            }            
        }
    }
    
    /**
     * NOTE: returns null if you've already called @see getInputAsPrintStream
     * @return
     *
     * @see com.cmcweb.extcallers.IExternalProcess#getInputOutputStream()
     */
    public FileOutputStream getInputAsStream() {
        if(m_inputps != null){ return null; }
        try {
            m_inputfile = File.createTempFile(m_infileName, "");
            m_inputfos = new FileOutputStream(m_inputfile.getPath());
        } catch (IOException e) {
            e.printStackTrace();
        }

        return m_inputfos; 
        
    }

    /**
     * 
     * NOTE: returns null if you've already called @see getInputAsPrintStream
     * @return
     *
     * @see com.cmcweb.extcallers.IExternalProcess#getInputPrintStream()
     */
    public PrintStream getInputAsPrintStream() {
        if(m_inputfos != null){ return null; }
        try {
            m_inputfile = File.createTempFile(m_infileName, "");
            m_inputps = new PrintStream(m_inputfile.getPath()); 
        } catch (IOException e) {
            e.printStackTrace();
        }
        
        return m_inputps; 
    }

    /**
     * @return
     *
     * @see com.cmcweb.extcallers.IExternalProcess#getProcessOutput()
     */
    public BufferedReader getProcessOutput() {
        return this.m_output; 
    }

    /**
     * 
     *
     * @see com.cmcweb.extcallers.IExternalProcess#runProcess()
     */
    public void runProcess() {
        if(m_inputps != null ){ m_inputps.close(); }
        if(m_inputfos != null ){ 
            try {
                m_inputfos.close();
            } catch (IOException e) {
                e.printStackTrace();
            } 
        }
        
        //now run process according to different os's catching the output
        if(System.getProperty("os.name").toUpperCase().indexOf("WIN") !=-1){
            this.m_output = runWin(); 
        }else{
            this.m_output = runNix(); 
        }
        m_inputfile.deleteOnExit();
    }
    
    /**
     * not the cleanest (requires bash), but should work for linux and macs
     * 
     * NOTE: we could remove the bash dependency by piping the input file to the
     *       standard input stream of the process.
     * 
     * @return
     */
    private BufferedReader runNix(){
        String bashPath = ExecutableFinder.findExecutable("bash");
        if(bashPath == null){ throw new RuntimeException("ERROR ! failed to find the bash interpreter in your path !"); }
        String[] cmd=new String[3];
        cmd[0]=bashPath;
        cmd[1]="-c";
        cmd[2]=m_exeToCall+this.m_qhexe.qh_args+" "+this.m_extraArgs+" < "+m_inputfile.getPath();
        try {
            m_p=Runtime.getRuntime().exec(cmd);
        } catch (IOException e1) {
            e1.printStackTrace();
            throw new RuntimeException("Error while executing command: ["+
                    cmd[0]+" "+cmd[1]+" "+cmd[2]+"]");            
        }
        BufferedReader r = new BufferedReader(new InputStreamReader(m_p.getInputStream()));
        try {
        	r.mark(100);
        	String firstline=r.readLine();
			if(firstline==null){
			    throw new RuntimeException(
				        "Error! QHull output missing after ["+
				        cmd[0] + " " + cmd[1] + " " +  cmd[2] + "]");
			}
			r.reset();
		} catch (IOException e) {
			e.printStackTrace();
		} 
		return r;
    }
        
    private BufferedReader runWin(){

        File tmp = null; 
        String[] cmd=new String[3];
        try {
            tmp = File.createTempFile("temp", "");
            cmd[0] = "cmd.exe" ;
            cmd[1] = "/C" ;
            cmd[2]= "\"\"" + this.m_exeToCall + "\"" + 
                   this.m_qhexe.qh_args +
                   " "+
                   this.m_extraArgs +
                   "< \""+ m_inputfile.getPath() + "\" > \""+ tmp.getPath() + "\"\"";
            m_p=Runtime.getRuntime().exec(cmd);
        } catch (IOException e1) {
            e1.printStackTrace();
            throw new RuntimeException("Error while executing command: ["+
                    cmd[0]+" "+cmd[1]+" "+cmd[2]+"]");
        } 
        
        try {
            m_p.waitFor();
        } catch (InterruptedException e) {
            e.printStackTrace();
        }
        
        BufferedReader r = null;        
        try {
            r = new BufferedReader(new InputStreamReader(new FileInputStream(tmp.getPath())));
            try {
            	r.mark(100);
            	String firstline=r.readLine();
    			if(firstline==null){
                    throw new RuntimeException(
                            "Error! QHull output missing after ["+
                            cmd[0] + " " + cmd[1] + " " +  cmd[2] + "]");    			    
    			}
    			r.reset();
    		} catch (IOException e) {
    			// TODO Auto-generated catch block
    			e.printStackTrace();
    		} 
        } catch (FileNotFoundException e) {
            e.printStackTrace();
        }
        return r; 
    }

    /**
     * @param s
     *
     * @see com.cmcweb.extcallers.IExternalProcess#setInputFilename(java.lang.String)
     */
    public void setInputFilename(String s) {
        this.m_infileName = s; 
    }

    public void cleanUp(){
        if(this.m_p != null) 
        	this.m_p.destroy();
        this.m_p = null;
        this.m_output = null; 
        if (m_inputfile.exists())
        	m_inputfile.delete();
    }
    
}
