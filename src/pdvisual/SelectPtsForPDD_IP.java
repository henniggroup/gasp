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

import java.awt.*;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.util.List;
import java.util.LinkedList;
import java.util.Map;
import java.util.HashMap;
import javax.swing.JOptionPane;

import chemistry.*;

import jv.object.PsObject;
import jv.object.PsUpdateIf;
import jv.project.PgGeometryIf;
import jv.project.PjProject_IP;
import jv.project.PvDisplayIf;
import jv.project.PvPickEvent;
import jv.project.PvPickListenerIf;
import jv.vecmath.PdVector;

public class SelectPtsForPDD_IP extends PjProject_IP implements ActionListener,PvPickListenerIf {
	static final long serialVersionUID = 1;
	
	protected	TernPD3D			m_ternpd3d;
	protected	Button				m_doneButton;
	protected	Button				m_resetButton;
	protected	Button				m_manualEntryButton;
	protected	Button				m_removeOneButton;
	protected	Button				m_selectorButton;
	protected	Choice				m_pointsChoice;
	
	// currently selected compositions
	List<Composition> comps;
	
	Map<Composition,Integer> compToIndxMap;

	public SelectPtsForPDD_IP(TernPD3D t) {
		super();
		
		m_ternpd3d = t;
		
		comps = new LinkedList<Composition>();
		compToIndxMap = new HashMap<Composition, Integer>();
		
		if (getClass() == SelectPtsForPDD_IP.class)
			init();
	}
	
	public void init() {
		super.init();
		addTitle("");

		Panel surfPanel = new Panel();
		{
			surfPanel.setLayout(new GridLayout(1,6));
			//surfPanel.add(new Label("Make a pseudo-phase diagram."));
			// add Selector button
			m_selectorButton = new Button("Selector");
			m_selectorButton.addActionListener(this);
			surfPanel.add(m_selectorButton);
			// add Done button
			m_doneButton = new Button("Do it!");
			m_doneButton.addActionListener(this);
			surfPanel.add(m_doneButton);
			// add Reset button
			m_resetButton = new Button("Reset");
			m_resetButton.addActionListener(this);
			surfPanel.add(m_resetButton);
			// add Manual Entry button
			m_manualEntryButton = new Button("Manual entry");
			m_manualEntryButton.addActionListener(this);
			surfPanel.add(m_manualEntryButton);
			// add Remove One button
			m_removeOneButton = new Button("Remove");
			m_removeOneButton.addActionListener(this);
			surfPanel.add(m_removeOneButton);
			// add points viewer choice
			m_pointsChoice = new Choice();
			surfPanel.add(m_pointsChoice);
		}
		add(surfPanel);
	}
	/**
	 * Set parent of panel which supplies the data inspected by the panel.
	 */
	public void setParent(PsUpdateIf parent) {
		super.setParent(parent);
		m_ternpd3d = (TernPD3D)parent;
	//	m_scalarChoice.removeAll();
	//	String [] listScalarFields = m_pjHeight.m_scalarFields;
	//	if (listScalarFields != null) {
	//		int len = listScalarFields.length;
	//		for (int i=0; i<len; i++)
	//			m_scalarChoice.add(listScalarFields[i]);
	//	}
	}
	/**
	 * Update the panel whenever the parent has changed somewhere else.
	 * Method is invoked from the parent or its superclasses.
	 */
	public boolean update(Object event) {
//		if (event == m_pjHeight) {
//			m_scalarChoice.select(m_pjHeight.getScalarName());
//			setTitle("Scalar Field "+m_pjHeight.getScalarName());
//		}
		return super.update(event);
	}
	/**
	 * Handle item state events invoked from choices, checkboxes, list items.
	 */
	public void actionPerformed(ActionEvent event) {
//		if (m_pjHeight==null)
//			return;
		Object source = event.getSource();
		if (source == m_doneButton) {
			finishMakingPseudoPDD();
		} else if (source == m_resetButton) {
			// clear comps list and the Choice and compToIndxMap 
			while (! comps.isEmpty())
				removePoint(0);
		} else if (source == m_manualEntryButton) {
			// make a popup to get it, then call addComposition
			CompoManualEntryDialog entryDialog = new CompoManualEntryDialog(m_ternpd3d, ((PDData)(m_ternpd3d.pdd)).getElements());
			Integer values[] = entryDialog.getValues();
			// if values is null, user pushed cancel
			if (values == null)
				return;
			Composition newC = new Composition(((PDData)m_ternpd3d.pdd).getElementsArray(), values);
			addComposition(newC);
		} else if (source == m_removeOneButton) {
			// remove the currently selected choice from Choice and comps
			int i = m_pointsChoice.getSelectedIndex();
			if (i != -1)
				removePoint(i);
		} else if (source == m_selectorButton) {
			m_ternpd3d.m_viewer.getDisplay().setMajorMode(PvDisplayIf.MODE_MARK);
		}
	}	
	
	private void removePoint(int i) {
		// if it was a point given by manual selection of a vertex, unmark it
		if (compToIndxMap.containsKey(comps.get(i))) {
			m_ternpd3d.m_project.m_vectors.getVertex(compToIndxMap.get(comps.get(i))).clearTag(PsObject.IS_SELECTED);
			compToIndxMap.remove(comps.get(i));
		}
		comps.remove(i);
		m_pointsChoice.remove(i);
		// redraw graph here
		m_ternpd3d.resetGraph();
	}
	
	private boolean addComposition(Composition newComp) {
		// check to make sure it's not a duplicate
		for (Composition c : comps) 
			if (c.equals(newComp)) 
				return false;

		// if all's ok, add it
		comps.add(newComp);
		// add it to the list
		m_pointsChoice.add(newComp.toString());

		return true;
	}
	
	private void finishMakingPseudoPDD() {
		
		// fail if we have too many or two few points selected
		if (comps.size() <= 1 || comps.size() >= m_ternpd3d.pdd.getDimension()) {
			JOptionPane.showMessageDialog(null, "Can't work with " + comps.size() + " points.");
			return;
		}

		// make the pseudopdd
		// catch and inform the user in the case of degenerate selected subspace
		PseudoPDBuilder pddb = null;
		PseudoPDData ppdd = null;
		try {
			if (m_ternpd3d.pdd instanceof PDData)
				pddb = new PseudoPDBuilder((PDData)m_ternpd3d.pdd);
			else if (m_ternpd3d.pdd instanceof PseudoPDData)
				pddb = new PseudoPDBuilder(((PseudoPDData)m_ternpd3d.pdd).getPDD());
			 ppdd = pddb.getPseudoPDData(comps);
		} catch (RuntimeException x) {
			JOptionPane.showMessageDialog(null, "Oops, looks like you gave a degenerate basis!");
			return;
		}
		new TernPD3D(ppdd);
		
	}
	
	public void dragDisplay(PvPickEvent arg0) { /* nothing */ }
	public void dragInitial(PvPickEvent arg0) { /* nothing */ }
	public void dragVertex(PgGeometryIf arg0, int arg1, PdVector arg2) { /* nothing */ }
	public void markVertices(PvPickEvent arg0) { /* nothing */ 
		if (this.isVisible()) {
			// get marked vertices
			int markedVertices[] = m_ternpd3d.m_project.m_vectors.getMarkedVertices();
			for (int v : markedVertices) {
				// add comp to Choice and to comps
				int pddIndx = m_ternpd3d.m_project.mapIndxsGraphToPDD(v);
				Composition c = m_ternpd3d.pdd.getComposition(pddIndx);
				if (addComposition(c)) {
					compToIndxMap.put(c, v);
				}
			}
		}
	}
	public void pickDisplay(PvPickEvent arg0) { /* nothing */ }
	public void pickInitial(PvPickEvent arg0) { /* nothing */ }
	public void pickVertex(PgGeometryIf arg0, int arg1, PdVector arg2) { /* nothing */	}
	public void selectGeometry(PgGeometryIf arg0) { /* nothing */ }
	public void unmarkVertices(PvPickEvent arg0) { /* nothing */ }
	
	
}

