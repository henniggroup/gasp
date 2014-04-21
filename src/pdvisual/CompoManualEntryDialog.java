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

import chemistry.*;
import java.util.List;
import java.util.LinkedList;
import java.awt.BorderLayout;
import java.awt.Frame;
import java.awt.FlowLayout;
import java.awt.event.ActionListener;
import javax.swing.*;

// thanks to: http://java.joycoding.com/java/544/194544-jdialog-multiple-inputs.html

public class CompoManualEntryDialog extends JDialog {
 static final long serialVersionUID = 1;

 private Integer values[];
 private Frame m_frame;
 List<JTextField> textFields;

 public CompoManualEntryDialog(Frame frame, List<Element> elements) {
		 super(frame, "Input composition", true);
		 m_frame = frame;
		 values = new Integer[elements.size()];
		 JPanel btnPanel = new JPanel();
		 setLayout(new FlowLayout());
		// getContentPane().setLayout(new GridLayout(1,2));
		 JButton okBtn = new JButton("OK");
		 JButton noBtn = new JButton("Cancel");
		 btnPanel.add(okBtn);
		 okBtn.addActionListener(new ActionListener() {
			 public void actionPerformed(java.awt.event.ActionEvent ae) {
				 okButton();
			 }
		 });
		 noBtn.addActionListener(new ActionListener() {
			 public void actionPerformed(java.awt.event.ActionEvent ae) {
				 noButton();
			 }
		 });
		 btnPanel.add(noBtn);
		 
		 textFields = new LinkedList<JTextField>();
		 for (Element e : elements) {
			 JLabel label = new JLabel(e.getSymbol() + ": ");
			 getContentPane().add(label, BorderLayout.CENTER);
			 JTextField input = new JTextField(4);
			 getContentPane().add(input, BorderLayout.CENTER);
			 textFields.add(input);
		 }
		 
		 getContentPane().add(btnPanel, BorderLayout.SOUTH);
		 pack();
		 
		 setVisible(true);
 	}

	 public Integer[] getValues() {
		 return values;
	 }
	
	 private void okButton() {
		 for (int i = 0; i < values.length; i++) {
			String s = textFields.get(i).getText();
			if ((s != null) && (s.length() > 0)) {
			    try {
			    	values[i] = Integer.parseInt(s);
			    } catch (NumberFormatException x) {
			    	JOptionPane.showMessageDialog(m_frame, "Ill formatted response: " + s);
			    	return;
			    }
			} else {
				values[i] = 0;
			}
		 }
		 setVisible(false);
	 }
	
	 private void noButton() {
		 values = null;
		 setVisible(false);
	 }
}