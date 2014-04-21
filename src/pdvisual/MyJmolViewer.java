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

package pdvisual;
// see http://www.biojava.org/wiki/BioJava:CookBook:PDB:Jmol

import java.awt.Container;
import java.awt.Dimension;
import java.awt.Graphics;
import java.awt.Rectangle;
import java.awt.event.WindowAdapter;
import java.awt.event.WindowEvent;
import java.util.ArrayList;

import javax.swing.JFrame;
import javax.swing.JPanel;
import org.jmol.adapter.smarter.SmarterJmolAdapter;
import org.jmol.api.JmolAdapter;
import org.jmol.api.JmolSimpleViewer;

import com.sun.xml.internal.bind.v2.schemagen.xmlschema.List;

import crystallography.Cell;

public class MyJmolViewer {

	    JmolSimpleViewer viewer;
	    IComputedEntry structure; 
	 
	    JmolPanel jmolPanel;
	    JFrame frame ;
	 /*
	    public static void main(String[] args){
	        try {
	 
	            PDBFileReader pdbr = new PDBFileReader();          
	            pdbr.setPath("/Path/To/PDBFiles/");
	 
	            String pdbCode = "5pti";
	 
	            Structure struc = pdbr.getStructureById(pdbCode);
	 
	            SimpleJmolExample ex = new SimpleJmolExample();
	            ex.setStructure(struc);
	 
	 
	        } catch (Exception e){
	            e.printStackTrace();
	        }
	    }*/
	 
	 
	    public MyJmolViewer(IComputedEntry e) {
	    	
	        frame = new JFrame();
	        Container contentPane = frame.getContentPane();
	        jmolPanel = new JmolPanel();
	 
	        jmolPanel.setPreferredSize(new Dimension(200,200));
	        contentPane.add(jmolPanel);
	        
	        setStructure(e);
	 
	        frame.pack();
	        frame.setVisible(true); 
	 
	    }
	    public void setStructure(IComputedEntry s) {
	 
	    //    frame.setName(s.getPDBCode());
	 
	        // actually this is very simple
	        // just convert the structure to a PDB file
	 
	    //    String pdb = s.toPDB();
	 
	   //     structure = s;
	        JmolSimpleViewer viewer = jmolPanel.getViewer();
	 
	        // Jmol could also read the file directly from your file system
	        //viewer.openFile("/Path/To/PDB/1tim.pdb");
	 
	        // send the PDB file to Jmol.
	        // there are also other ways to interact with Jmol, but they require more
	        // code. See the link to SPICE above...
	        java.util.List<java.util.List<Integer>>coefs = new ArrayList<java.util.List<Integer> >();
	        for (int i = 0; i < 3; i++) {
	        	java.util.List<Integer> vect = new ArrayList<Integer>(3);
	        	vect.add(0);vect.add(0);vect.add(0);vect.set(i, 3);
	        	coefs.add(vect);
	        }
	        
	        viewer.openStringInline(Cell.getSupercell(s.getCell(), coefs).getCIF());
	        viewer.evalString("wireframe 0.15; spacefill 0.45;");
	     //   viewer.evalString("select *; spacefill off; wireframe off; backbone 0.4;  ");
	      //  viewer.evalString("color chain;  ");
	        this.viewer = viewer;
	 
	    }
	 
	    public void setTitle(String label){
	        frame.setTitle(label);
	    }
	 
	    public JmolSimpleViewer getViewer(){
	 
	        return jmolPanel.getViewer();
	    }
	 

	    static class JmolPanel extends JPanel {
	        /**
	         * 
	         */
	        private static final long serialVersionUID = -3661941083797644242L;
	        JmolSimpleViewer viewer;
	        JmolAdapter adapter;
	        JmolPanel() {
	            adapter = new SmarterJmolAdapter();
	            viewer = JmolSimpleViewer.allocateSimpleViewer(this, adapter);
	 
	        }
	 
	        public JmolSimpleViewer getViewer() {
	            return viewer;
	        }
	 
	        public void executeCmd(String rasmolScript){
	            viewer.evalString(rasmolScript);
	        }
	 
	 
	        final Dimension currentSize = new Dimension();
	        final Rectangle rectClip = new Rectangle();
	 
	        public void paint(Graphics g) {
	            getSize(currentSize);
	            g.getClipBounds(rectClip);
	            viewer.renderScreenImage(g, currentSize, rectClip);
	        }
	    }
	 
	}