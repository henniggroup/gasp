<<<<<<< HEAD
=======
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

>>>>>>> 0e3189c40547bbd59ea42c4f91890d7511fb7797
package pdvisual;

import java.awt.*;
import java.util.Iterator;
import java.util.LinkedList;
import jv.object.*;
import jv.project.PvDisplayIf;
import jv.viewer.PvViewer;
import chemistry.*;
import crystallography.*;
import javax.swing.*;
import ga.StructureOrg;

import java.util.List;

import utility.Constants;
import utility.Utility;
import vasp.VaspIn;

import java.awt.event.*;
import java.io.File;

public class TernPD3D extends jv.object.PsMainFrame implements ActionListener, ItemListener {
	/** 3D-viewer window for graphics output and which is embedded into the applet */
	protected	PvViewer			m_viewer;

	protected IPDPlottable pdd;
	
	TernPD3DProj m_project;
	
	//Where the GUI is created:
	MenuBar menuBar;
	Menu optionsMenu;
	MenuItem miMaxDistFromHull;
	Double maxDistFromHull = 1000.0;
	CheckboxMenuItem cbmiDisplayPosEnergyPoints;
	CheckboxMenuItem cbmiEnergyType;
	CheckboxMenuItem cbmiAutoCenter;
	Menu showLabelsMenu;
	CheckboxMenuItem cbmiShowLabelsNone;
	CheckboxMenuItem cbmiShowLabelsPseudo;
	CheckboxMenuItem cbmiShowLabelsReal;
	CheckboxMenuItem cbmiShowLabelsAll;
	CheckboxMenuItem cbmiShowAxes;
	CheckboxMenuItem cbmiBigMode;
	Menu colorByChempotsSubmenu;
	java.util.List<CheckboxMenuItem> colorByChempotChoices;
	Menu colorByInstabilitySubmenu;
	CheckboxMenuItem cbmiChempotNums;
	CheckboxMenuItem cbmiColorNone;
	CheckboxMenuItem cbmiColorBinary;
	CheckboxMenuItem cbmiColorSpectrum;
	CheckboxMenuItem cbmiColorType;
	
	Menu toolsMenu;
	MenuItem miZoomBox;
	CheckboxMenuItem cbmiMakePDD;
	MenuItem miGetStructureInfo;
	
	MenuItem miExportPDD;
	MenuItem miExportCHPDB;
	MenuItem miExportData;
	MenuItem miExportStructs;
	MenuItem miExportVoltageCurve;
	
	SelectPtsForPDD_IP selectPtsInfoPanel;

	public TernPD3D(IPDPlottable _pdd) {
		super("");
		pdd = _pdd;

		// Make a title
		java.util.List<String> axisLabels = pdd.getAxisLabels();
		StringBuilder title = new StringBuilder();
		title.append("PD: ");
		Iterator<String> i = axisLabels.iterator();
		while (i.hasNext()) {
			title.append(i.next());
			if (i.hasNext())
				title.append("-");
		}
		
		// Create top level window of application containing the applet
		this.setTitle(title.toString());
		this.pack();
		
		m_viewer = new PvViewer();
		
		// add a menu
		makeMenu();
		selectPtsInfoPanel = new SelectPtsForPDD_IP(this);
		
		init();
		
		this.setBounds(new Rectangle(420, 5, 660, 550));
		this.setVisible(true);
		
	}

	/**
	 * Configure and initialize the viewer, load system and user projects.
	 * One of the user projects must be selected here.
	 */
	public void init() {
		//if (m_frame.getComponentCount() > 0 && m_viewer != null)
		//	m_frame.remove((Component)(m_viewer.getDisplay()));
		//m_viewer = new PvViewer();
		
		// remove the old project (if extant)
	//	if (m_project != null)
	//		m_viewer.removeProject(m_project);
		// Create and load a new project
		m_project = new TernPD3DProj(pdd, this);
		m_viewer.addProject(m_project);
		m_viewer.selectProject(m_project);

		// Get 3d display from viewer and add it to applet
		PvDisplayIf disp = m_viewer.getDisplay();
		this.add((Component)disp);
		this.pack();
		
		// set background color
		m_viewer.getDisplay().setBackgroundColor(Color.WHITE);
		
		// add the control panel
		//m_frame.add(m_viewer.getPanel(PsViewerIf.PROJECT), BorderLayout.EAST);
		
		// Choose initial panel in control window (must press F1 inside the applet)
		m_viewer.showPanel(PsViewerIf.CAMERA);

		m_viewer.start();
	}

 //  public void addEntryVisualizer(EntryVisualizer visualizer)
  //  {
   //     m_project.addEntryVisualizer(visualizer);
    //}

    public void setNewPDD(IPDPlottable _pdd)
    {
        pdd = _pdd;
        init();
    }
    
	private void makeMenu() {
		// make the menubar
//		menuBar = new MenuBar();
		menuBar = m_viewer.newMenuBar(this);
		// make the options menu
		optionsMenu = new Menu("Options");
		// make max distance from hull menuitem
		miMaxDistFromHull = new MenuItem("Max dist from hull: " + maxDistFromHull, null);
		miMaxDistFromHull.addActionListener(this);
		optionsMenu.add(miMaxDistFromHull);
		//menu.addSeparator();
		// make display positive energy points menuitem
		cbmiDisplayPosEnergyPoints = new CheckboxMenuItem("Display E>0 points");
		cbmiDisplayPosEnergyPoints.addItemListener(this);
		cbmiDisplayPosEnergyPoints.setState(false);
		optionsMenu.add(cbmiDisplayPosEnergyPoints);
		// add checkboxmenuitem to toggle auto-fitting
		cbmiAutoCenter = new CheckboxMenuItem("Auto center/fit on change");
		cbmiAutoCenter.setState(false);
		cbmiAutoCenter.addItemListener(this);
		optionsMenu.add(cbmiAutoCenter);
		optionsMenu.addSeparator();
		// add item to select what to plot on the energy axis
		cbmiEnergyType = new CheckboxMenuItem("Plot formation energies");
		cbmiEnergyType.setState(true);
		cbmiEnergyType.addItemListener(this);
		optionsMenu.add(cbmiEnergyType);
		// add menu toggle showing vertex labels
		Menu showLabelsMenu = new Menu("Show vertex labels");
		cbmiShowLabelsNone = new CheckboxMenuItem("None");
		cbmiShowLabelsNone.setState(false);
		cbmiShowLabelsNone.addItemListener(this);
		showLabelsMenu.add(cbmiShowLabelsNone);
		cbmiShowLabelsPseudo = new CheckboxMenuItem("Pseudo only");
		cbmiShowLabelsPseudo.setState(false);
		cbmiShowLabelsPseudo.addItemListener(this);
		showLabelsMenu.add(cbmiShowLabelsPseudo);
		cbmiShowLabelsReal = new CheckboxMenuItem("Real only");
		cbmiShowLabelsReal.setState(false);
		cbmiShowLabelsReal.addItemListener(this);
		showLabelsMenu.add(cbmiShowLabelsReal);
		cbmiShowLabelsAll = new CheckboxMenuItem("All");
		cbmiShowLabelsAll.setState(true);
		cbmiShowLabelsAll.addItemListener(this);
		showLabelsMenu.add(cbmiShowLabelsAll);
		optionsMenu.add(showLabelsMenu);
		// add item to toggle showing axes
		cbmiShowAxes = new CheckboxMenuItem("Display axes");
			// default to displaying axes in 2D case
		cbmiShowAxes.setState(pdd.getDimension() == 2);
		cbmiShowAxes.addItemListener(this);
		optionsMenu.add(cbmiShowAxes);
		// make big/presentation mode option
		cbmiBigMode = new CheckboxMenuItem("Presentation mode");
		cbmiBigMode.setState(false);
		cbmiBigMode.addItemListener(this);
		optionsMenu.add(cbmiBigMode);
		// make visualize chempots menuitem
		optionsMenu.addSeparator();
		colorByChempotsSubmenu = new Menu("Visualize chemical potentials");
		cbmiChempotNums = new CheckboxMenuItem("Show chempot number");
		cbmiChempotNums.addItemListener(this);
		colorByChempotsSubmenu.add(cbmiChempotNums);
		colorByChempotsSubmenu.addSeparator();
		colorByChempotChoices = new LinkedList<CheckboxMenuItem>();
		CheckboxMenuItem mi = new CheckboxMenuItem("None");
		mi.setState(true);
		mi.addItemListener(this);
		colorByChempotChoices.add(mi);
		colorByChempotsSubmenu.add(mi);
		for (String e : pdd.getAxisLabels()) {
			mi = new CheckboxMenuItem(e);
			mi.addItemListener(this);
			colorByChempotChoices.add(mi);
			colorByChempotsSubmenu.add(mi);
		}
		optionsMenu.add(colorByChempotsSubmenu);
		// make color by instability submenu
		colorByInstabilitySubmenu = new Menu("Color vertices");
		cbmiColorNone = new CheckboxMenuItem("None");
		cbmiColorNone.setState(false);
		cbmiColorNone.addItemListener(this);
		colorByInstabilitySubmenu.add(cbmiColorNone);
		colorByInstabilitySubmenu.addSeparator();
		cbmiColorBinary = new CheckboxMenuItem("Stability: Binary");
		cbmiColorBinary.setState(true);
		cbmiColorBinary.addItemListener(this);
		colorByInstabilitySubmenu.add(cbmiColorBinary);
		cbmiColorSpectrum = new CheckboxMenuItem("Stability:Spectrum");
		cbmiColorSpectrum.setState(false);
		cbmiColorSpectrum.addItemListener(this);
		colorByInstabilitySubmenu.add(cbmiColorSpectrum);
		cbmiColorType = new CheckboxMenuItem("Entry type");
		cbmiColorType.setState(false);
		cbmiColorType.addItemListener(this);
		colorByInstabilitySubmenu.add(cbmiColorType);
		optionsMenu.add(colorByInstabilitySubmenu);
		
		// make the tools menu
		toolsMenu = new Menu("Tools");
		
		// make zoom box tool menu item
		miZoomBox = new MenuItem("Box zoom", null);
		miZoomBox.addActionListener(this);
		toolsMenu.add(miZoomBox);
		// make "Make pseudophasediagram" menu option
		cbmiMakePDD = new CheckboxMenuItem("View pseudopd tool");
		cbmiMakePDD.addItemListener(this);
		cbmiMakePDD.setState(false);
		toolsMenu.add(cbmiMakePDD);
		
		// make "Get structure info" menu option
		miGetStructureInfo = new MenuItem("Get structure info", null);
		miGetStructureInfo.addActionListener(this);
		toolsMenu.add(miGetStructureInfo);
		
		// reorder the menubar!
		Menu fileMenu = menuBar.getMenu(0);
		Menu inspectorMenu = menuBar.getMenu(1);
		Menu methodMenu = menuBar.getMenu(2);
		Menu windowMenu = menuBar.getMenu(3);
		Menu helpMenu = menuBar.getMenu(4);
		menuBar.remove(fileMenu);
		menuBar.remove(inspectorMenu);
		menuBar.remove(methodMenu);
		menuBar.remove(windowMenu);
		menuBar.remove(helpMenu);
		
		// add export options to file menu
		fileMenu.addSeparator();
		miExportPDD = new MenuItem("Export PDD", null);
		miExportPDD.addActionListener(this);
		fileMenu.add(miExportPDD);
		miExportData = new MenuItem("Export Data", null);
		miExportData.addActionListener(this);
		fileMenu.add(miExportData);
		miExportCHPDB = new MenuItem("Export Convex Hull", null);
		miExportCHPDB.addActionListener(this);
		fileMenu.add(miExportCHPDB);
		miExportStructs = new MenuItem("Export Structures", null);
		miExportStructs.addActionListener(this);
		fileMenu.add(miExportStructs);
		miExportVoltageCurve = new MenuItem("Export V Curve", null);
		miExportVoltageCurve.addActionListener(this);
		fileMenu.add(miExportVoltageCurve);
		
		menuBar.add(fileMenu);
		menuBar.add(optionsMenu);
		menuBar.add(toolsMenu);
		menuBar.add(inspectorMenu);
		menuBar.add(methodMenu);
		menuBar.add(windowMenu);
		menuBar.add(helpMenu);
		
		// add the menubar to the frame
		this.setMenuBar(menuBar);
	}
	
	public void resetGraph() {
	    m_project.init();
		m_project.m_vectors.update(m_project.m_vectors);
		if (cbmiAutoCenter.getState()) {
			m_project.getDisplay().center();
			m_project.getDisplay().fit();
		}
	}
	
	// for the menuitems
	public void actionPerformed(ActionEvent e) {
		Object source = e.getSource();
		
		if (source == miMaxDistFromHull) {
			// pop up a dialog to ask the max distance from hull to plot
            String s = (String)JOptionPane.showInputDialog(
                    this,
                    "Maximum energy above hull:",
                    "0");

			//If a string was returned, say so.
			if ((s != null) && (s.length() > 0)) {
			    try {
			    	maxDistFromHull = Double.parseDouble(s);
			    } catch (NumberFormatException x) {
			    	JOptionPane.showMessageDialog(this, "Ill formatted response: " + s);
			    	return;
			    }
			} 
			
			// update the menu
			miMaxDistFromHull.setLabel("Max dist from hull: " + maxDistFromHull);
			// update graph
			resetGraph();
		} else if (source == miZoomBox) {
			m_project.getDisplay().setMajorMode(PvDisplayIf.MODE_SCALE_RECT);
		} else if (source == miGetStructureInfo) {
			m_project.getDisplay().setMajorMode(PvDisplayIf.MODE_PICK);
		} else if (source == miExportPDD) {
			//Create a file chooser
			final JFileChooser fc = new JFileChooser();
	        int returnVal = fc.showOpenDialog(this);
	        if (returnVal == JFileChooser.APPROVE_OPTION) {
	            String fileName = fc.getSelectedFile().getAbsolutePath();
	            Utility.writeSerializable((PDData)pdd, fileName);
	        }
		} else if (source == miExportCHPDB) {
			//Create a file chooser
			final JFileChooser fc = new JFileChooser();
	        int returnVal = fc.showOpenDialog(this);
	        if (returnVal == JFileChooser.APPROVE_OPTION) {
	            String fileName = fc.getSelectedFile().getAbsolutePath();
	            //This is where a real application would open the file.
	            PDData p = (PDData)pdd;
	            Utility.writeSerializable((new PDBuilder(p.getStableEntries(), p.getElements(), p.getChemPots())).getPDData(), fileName);
	        }
		} else if (source == miExportData) {
			//Create a file chooser
			final JFileChooser fc = new JFileChooser();
	        int returnVal = fc.showOpenDialog(this);
	        if (returnVal == JFileChooser.APPROVE_OPTION) {
	            String fileName = fc.getSelectedFile().getAbsolutePath();
	            //This is where a real application would open the file.
	            StringBuilder result = new StringBuilder();
	            for (int i = 0; i < m_project.m_vectors.getNumVertices(); i++) {
	            	result.append(m_project.m_vectors.getVertex(i).getName() + " ");
	            	for (int j = 0; j < Constants.numDimensions; j++)
	            		result.append(m_project.m_vectors.getVertex(i).getEntries()[j] + " ");
	            	result.append("\n");
	            }
	            Utility.writeStringToFile(result.toString(), fileName);
	       }
		} else if (source == miExportStructs) {
			final JFileChooser fc = new JFileChooser();
			fc.setFileSelectionMode(JFileChooser.DIRECTORIES_ONLY);
	        int returnVal = fc.showOpenDialog(this);
	        if (returnVal == JFileChooser.APPROVE_OPTION) {
	            String dirname = fc.getSelectedFile().getAbsolutePath();
	            //This is where a real application would open the file.
	            for (IComputedEntry ent : ((PDData)pdd).getAllEntries()) {
	            	int sid = ((StructureOrg)(ent)).getID();
	            	(new VaspIn(ent.getCell(), null, null, null)).writePoscar(dirname+"/"+sid+".POSCAR", false);
	            }
	       }
		} else if (source == miExportVoltageCurve) {
			final JFileChooser fc = new JFileChooser();
	        int returnVal = fc.showOpenDialog(this);
	        if (returnVal == JFileChooser.APPROVE_OPTION) {
	        	Element anode = (Element)JOptionPane.showInputDialog(null, "Choose intercalation element:", "Anode", JOptionPane.PLAIN_MESSAGE, null, pdd.getElements().toArray(), pdd.getElements().get(0));
	        	java.util.LinkedList<Composition> cathodeCompChoices = new java.util.LinkedList<Composition>();
	        	for (IComputedEntry i : ((PDData)pdd).getStableEntries()) {
	        		cathodeCompChoices.add(i.getComposition());
	        	}
	        	Composition cathode = (Composition)JOptionPane.showInputDialog(null, "Choose other endpoint:", "Cathode", JOptionPane.PLAIN_MESSAGE, null, cathodeCompChoices.toArray(), cathodeCompChoices.get(0));
	        	
	        	// get anode energy
	        	double anodeEnergy = Double.NaN;
	        	Composition anodeComp = new Composition(anode);
	        	for (int basisIndx : ((PDData)pdd).getIndxStableEntries())
	        		if (((PDData)pdd).getAllEntries().get(basisIndx).getComposition().equals(anodeComp))
	        			anodeEnergy = ((PDData)pdd).getAllEntries().get(basisIndx).getEnergyPerAtom();
	        	if (Double.isNaN(anodeEnergy)) {
	        		System.out.println("ERROR: no anode comp found in pseudo pdd?");
	        		(new Exception()).printStackTrace();
	        	}
	        	
	        	StringBuilder result = new StringBuilder();
	        	for (List<Integer> facet : ((PDData)pdd).getIndxFacets()) {
	        		// we should always have a binary pseudo pdd here
	        		if (facet.size() != 2) {
	        			System.out.println("ERROR: facet.size() != 2 in TernPD3D.");
	        			(new Exception()).printStackTrace();
	        		}
	        		IComputedEntry endpoint0 = ((PDData)pdd).getAllEntries().get(facet.get(0));
	        		IComputedEntry endpoint1 = ((PDData)pdd).getAllEntries().get(facet.get(1));
	        		
	        		// if this facet has too much anode, continue
	        		if (endpoint0.getComposition().getFractionalCompo(anode) > cathode.getFractionalCompo(anode)
	        				|| endpoint1.getComposition().getFractionalCompo(anode) > cathode.getFractionalCompo(anode))
	        			continue;
	        		
	        		// make sure endpoint 1 has more anode (lithium)
	        		if (endpoint0.getComposition().getFractionalCompo(anode) > endpoint1.getComposition().getFractionalCompo(anode)) {
	        			IComputedEntry temp = endpoint0;
	        			endpoint0 = endpoint1;
	        			endpoint1 = temp;
	        		}
	        		
	        		// e.g. compound is Li_a Si_b  ==  Li_(a/b) Si_1
	        		double endpoint0_a = endpoint0.getComposition().getFractionalCompo(anode);
	        		double endpoint0_b = 1-endpoint0_a;
	        		double endpoint1_a = endpoint1.getComposition().getFractionalCompo(anode);
	        		double endpoint1_b = 1-endpoint1_a;
	        		double x0 = endpoint0_a / endpoint0_b;
	        		double x1 = endpoint1_a / endpoint1_b;
	        		double energy0 = endpoint0.getEnergyPerAtom() * (1 + x0); // total energy of Li_(x0)Si_1
	        		double energy1 = endpoint1.getEnergyPerAtom() * (1 + x1); // total energy of Li_(x1)Si_1

	        		double voltage = - (energy1 - energy0 - (x1-x0) * anodeEnergy) / (x1 - x0);
	        		
	        		result.append(x0 + " " + x1 + " " + voltage + "\n");
	        	}
	        	
	            String fileName = fc.getSelectedFile().getAbsolutePath();
	            Utility.writeStringToFile(result.toString(), fileName);

	       } /* else if (source == miExportVoltageCurve) {
				final JFileChooser fc = new JFileChooser();
		        int returnVal = fc.showOpenDialog(this);
		        if (returnVal == JFileChooser.APPROVE_OPTION) {
		        	Element anode = (Element)JOptionPane.showInputDialog(null, "Choose intercalation element:", "Anode", JOptionPane.PLAIN_MESSAGE, null, pdd.getElements().toArray(), pdd.getElements().get(0));
		        	java.util.LinkedList<Composition> cathodeCompChoices = new java.util.LinkedList<Composition>();
		        	for (IComputedEntry i : ((PDData)pdd).getStableEntries()) {
		        		cathodeCompChoices.add(i.getComposition());
		        	}
		        	Composition cathode = (Composition)JOptionPane.showInputDialog(null, "Choose other endpoint:", "Cathode", JOptionPane.PLAIN_MESSAGE, null, cathodeCompChoices.toArray(), cathodeCompChoices.get(0));
		        	
		        	java.util.ArrayList<Composition> ppdbComps = new java.util.ArrayList<Composition>();
		        	Composition anodeComp = new Composition(anode);
		        	ppdbComps.add(anodeComp);
		        	ppdbComps.add(cathode);
		        	PseudoPDData ppdb = (new PseudoPDBuilder((PDData)pdd)).getPseudoPDData(ppdbComps);
		        	
		        	// get anode energy
		        	double anodeEnergy = Double.NaN;
		        	for (int basisIndx : ppdb.getBasisEntryIndxs())
		        		if (ppdb.getAllPhases().get(basisIndx).getComposition().equals(anodeComp))
		        			anodeEnergy = ppdb.getAllPhases().get(basisIndx).getEnergy();
		        	if (Double.isNaN(anodeEnergy)) {
		        		System.out.println("ERROR: no anode comp found in pseudo pdd?");
		        		(new Exception()).printStackTrace();
		        	}
		        	
		        	StringBuilder result = new StringBuilder();
		        	for (List<Integer> facet : ppdb.getIndxFacets()) {
		        		// we should always have a binary pseudo pdd here
		        		if (facet.size() != 2) {
		        			System.out.println("ERROR: facet.size() != 2 in TernPD3D.");
		        			(new Exception()).printStackTrace();
		        		}
		        		MixedPhase endpoint0 = ppdb.getAllPhases().get(facet.get(0));
		        		MixedPhase endpoint1 = ppdb.getAllPhases().get(facet.get(1));
		        		
		        		// make sure endpoint 1 has more anode (lithium)
		        		if (endpoint0.getComposition().getFractionalCompo(anode) > endpoint1.getComposition().getFractionalCompo(anode)) {
		        			MixedPhase temp = endpoint0;
		        			endpoint0 = endpoint1;
		        			endpoint1 = temp;
		        		}
		        		
		        		// e.g. compound is Li_a Si_b  ==  Li_(a/b) Si_1
		        		double endpoint0_a = endpoint0.getComposition().getFractionalCompo(anode);
		        		double endpoint0_b = 1-endpoint0_a;
		        		double endpoint1_a = endpoint1.getComposition().getFractionalCompo(anode);
		        		double endpoint1_b = 1-endpoint1_a;
		        		double x0 = endpoint0_a / endpoint0_b;
		        		double x1 = endpoint1_a / endpoint1_b;
		        		double energy0 = endpoint0.getEnergy();
		        		double energy1 = endpoint1.getEnergy();

		        		double voltage = (energy1 - energy0 - (x1-x0) * anodeEnergy) / (x1 - x0);
		        		
		        		result.append(x0 + " " + x1 + " " + voltage + "\n");
		        	}
		        	
		            String fileName = fc.getSelectedFile().getAbsolutePath();
		            Utility.writeStringToFile(result.toString(), fileName);

		       } */
		}
	}
        
        
        /**
         * Shyue Ping Ong : Oct 24 2008, 1400
         * I added this routinue to allow me to set the max distance from hull
         * as a default setting for grand canonical PDs. 
         * This is because for the grand canonical diagram, having all the entries
         * above the hull in the PD is bad since there will be overlaps.
         * @param dist
         */
        public void setMaxDistanceFromHull(double dist) {
            maxDistFromHull = dist;
            miMaxDistFromHull.setLabel("Max dist from hull: " + maxDistFromHull);
            resetGraph();
	}
	
	// for the checkboxes
	public void itemStateChanged(ItemEvent e) {
	    Object source = e.getItemSelectable();
	    
	   // if (e.getStateChange() == ItemEvent.DESELECTED)

	    if (source == cbmiDisplayPosEnergyPoints) {
	    	// replot the figure
		    resetGraph();
	    } else if (source == cbmiEnergyType) {
	    	// replot the figure
		    resetGraph();
	    
		} else if (source == cbmiAutoCenter) {
			// nothin
		} else if (source == cbmiShowAxes) {
			// replot the figure
		    resetGraph();
		} else if (source == cbmiShowLabelsNone) {
			cbmiShowLabelsNone.setState(true);
			cbmiShowLabelsPseudo.setState(false);
			cbmiShowLabelsReal.setState(false);
			cbmiShowLabelsAll.setState(false);
			// replot the figure
		    resetGraph();
		} else if (source == cbmiShowLabelsPseudo) {
			cbmiShowLabelsNone.setState(false);
			cbmiShowLabelsPseudo.setState(true);
			cbmiShowLabelsReal.setState(false);
			cbmiShowLabelsAll.setState(false);
			// replot the figure
		    resetGraph();
		} else if (source == cbmiShowLabelsReal) {
			cbmiShowLabelsNone.setState(false);
			cbmiShowLabelsPseudo.setState(false);
			cbmiShowLabelsReal.setState(true);
			cbmiShowLabelsAll.setState(false);
			// replot the figure
		    resetGraph();
		} else if (source == cbmiShowLabelsAll) {
			cbmiShowLabelsNone.setState(false);
			cbmiShowLabelsPseudo.setState(false);
			cbmiShowLabelsReal.setState(false);
			cbmiShowLabelsAll.setState(true);
			// replot the figure
		    resetGraph();
		} else if (source == cbmiChempotNums) {
			// replot the figure
			resetGraph();
		} else if (colorByChempotChoices.contains(source)) {
			for (CheckboxMenuItem mi : colorByChempotChoices)
				if (source == mi) {
					// nothin
				} else {
					mi.setState(false);
				}
			resetGraph();
		} else if (source == cbmiColorBinary) {			
			cbmiColorNone.setState(false);
			cbmiColorBinary.setState(true);
			cbmiColorSpectrum.setState(false);
			cbmiColorType.setState(false);
			resetGraph();
		} else if (source == cbmiColorSpectrum) {			
			cbmiColorNone.setState(false);
			cbmiColorBinary.setState(false);
			cbmiColorSpectrum.setState(true);
			cbmiColorType.setState(false);
			resetGraph();
		} else if (source == cbmiColorNone) {
			cbmiColorNone.setState(true);
			cbmiColorBinary.setState(false);
			cbmiColorSpectrum.setState(false);
			cbmiColorType.setState(false);
			resetGraph();
		} else if (source == cbmiColorType) {
			cbmiColorNone.setState(false);
			cbmiColorBinary.setState(false);
			cbmiColorSpectrum.setState(false);
			cbmiColorType.setState(true);
			resetGraph();
		}else if (source == cbmiBigMode) {
			resetGraph();
		} else if (source == cbmiMakePDD) {
			// don't work with PseudoPDDatas
			if (pdd instanceof PseudoPDData)
				JOptionPane.showMessageDialog(this, "Making pseudo-pseudo phase diagrams not supported");
			
			if (cbmiMakePDD.getState()) {
				this.add(selectPtsInfoPanel, BorderLayout.SOUTH);
			} else {
				this.remove(selectPtsInfoPanel);
			}
			this.validate();
		}
	}
	
	public boolean getUseFormationEnergies() {
		return cbmiEnergyType.getState();
	}
	public boolean getDispPosEPoints() {
		return cbmiDisplayPosEnergyPoints.getState();
	}
	public boolean getShowAxes() {
		return cbmiShowAxes.getState();
	}
	public boolean getShowLabelsReal() {
		return cbmiShowLabelsReal.getState() || cbmiShowLabelsAll.getState();
	}
	public boolean getShowLabelsPseudo() {
		return cbmiShowLabelsPseudo.getState() || cbmiShowLabelsAll.getState();
	}
	public Double getMaxDistFromHull() {
		return maxDistFromHull;
	}
	public boolean getColorBinary() {
		return cbmiColorBinary.getState();
	}
	public boolean getColorSpectrum() {
		return cbmiColorSpectrum.getState();
	}
	public boolean getColorType() {
		return cbmiColorType.getState();
	}
	// returns an integer indicating the element whose chemical potential we
	// should use to color each facet or -1 if none
	public int getColorByChempots() {
		for (int i = 0; i < colorByChempotChoices.size(); i++)
			if (colorByChempotChoices.get(i).getState())
				return i - 1;
		
		// should never get here
		return -1;
	}
	public boolean getNumberChempots() {
		return cbmiChempotNums.getState();
	}
	public boolean getBigMode() {
		return cbmiBigMode.getState();
	}
	
	/**
	 * Standalone application support. The main() method acts as the applet's
	 * entry point when it is run as a standalone application. It is ignored
	 * if the applet is run from within an HTML page.
	 */
//	public static void main(String args[]) {
	//	new TernPD3D(null);

//	}
}
