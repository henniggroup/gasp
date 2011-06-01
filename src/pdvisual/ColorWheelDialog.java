package pdvisual;

import java.awt.BorderLayout;
import java.awt.Frame;
import java.awt.event.ActionListener;
import java.awt.Color;
import java.util.Collections;
import java.util.List;
import java.awt.LinearGradientPaint;
import java.awt.Graphics2D;
import java.awt.geom.RoundRectangle2D;
import java.awt.image.BufferedImage;
import javax.swing.JButton;
import javax.swing.JDialog;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.ImageIcon;

import jv.number.PdColor;

public class ColorWheelDialog extends JDialog {
	 static final long serialVersionUID = 1;

	 private List<Double> values;
	 
	 private double maxValue;
	 private double minValue;
	 private int gradientWidth = 50;
	 private int gradientHeight = 200;
	 
	 // list of colors should be the same length as the list of fractions
	 Color colors[] = new Color[]{Color.blue,Color.green,Color.red};
	 // this should start with 0.0, end with 1.0, and have monotone-increasing fractions in between
	 float fractions[] = new float[]{(float)0.0, (float)0.5, (float)1.0};

	 public ColorWheelDialog(Frame frame, List<Double> vs) {
		 super(frame, "Color Wheel Dialog", true);
		 setModal(false);
		 values = vs;
		 
		 // calculate some things
		 maxValue = Collections.max(values);
		 minValue = Collections.min(values);

		 JPanel btnPanel = new JPanel();
		 setLayout(new BorderLayout());
		// getContentPane().setLayout(new GridLayout(1,2));
		 JButton doneBtn = new JButton("OK");
		 btnPanel.add(doneBtn);
		 doneBtn.addActionListener(new ActionListener() {
			 public void actionPerformed(java.awt.event.ActionEvent ae) {
				 doneButton();
			 }
		 });
		 
		 // make an image w/ some colors

		 JLabel scale = new JLabel(new ImageIcon(makeScale()));
		 getContentPane().add(scale, BorderLayout.EAST);
		 JLabel grad = new JLabel(new ImageIcon(makeGradient()));
		 getContentPane().add(grad, BorderLayout.WEST);
		 
		 getContentPane().add(btnPanel, BorderLayout.SOUTH);
		 pack();
		 
		 setVisible(true);
 	}
	 
	 private BufferedImage makeScale() {
		 BufferedImage buff = new BufferedImage(gradientWidth,gradientHeight,BufferedImage.TYPE_4BYTE_ABGR);
		 Graphics2D g2 = (Graphics2D)buff.getGraphics();
		 for (int i = 0; i < fractions.length; i++) {
			 String chempotStr = new Double(getUnnormedVal(fractions[i])).toString();
			 g2.setColor(getColor(getUnnormedVal(fractions[i])));
			 int yloc = (int) Math.max(fractions[i] * gradientHeight,12);
			 g2.drawString(chempotStr, 5, yloc);
		 }
		 
		 return buff;
	 }
	 
	 private BufferedImage makeGradient() {
		 BufferedImage buff = new BufferedImage(gradientWidth,gradientHeight,BufferedImage.TYPE_4BYTE_ABGR);
		 Graphics2D g2 = (Graphics2D)buff.getGraphics();
		 g2.setPaint(new LinearGradientPaint(0, 0, 0, gradientHeight,fractions,colors));
		 g2.fill(new RoundRectangle2D.Double(0, 0, gradientWidth, gradientHeight, 5, 5));
		 return buff;
	 }
	 
	 private double getFraction(double value) {
		 return (value - minValue) / (maxValue - minValue);
	 }
	 private double getUnnormedVal(double fraction) {
		 return fraction * (maxValue - minValue) + minValue;
	 }
	 
	 public Color getColor(double value) {
		 double normvalue = getFraction(value);
		 
		 int topIndx;
		 for (topIndx = 1; topIndx < fractions.length; topIndx++)
			 if (normvalue <= fractions[topIndx])
				 break;
		 
		 // so we want a linear gradient between:
		 Color bottomColor = colors[topIndx - 1];
		 Color topColor = colors[topIndx];
		 
		 double frac = (normvalue - fractions[topIndx-1]) / (fractions[topIndx] - fractions[topIndx-1]);
		 
		 Color resultC = PdColor.blend(1-frac, bottomColor, frac, topColor);
		 
		 return resultC;
	 }
	
	 private void doneButton() {
		 setVisible(false);
	 }
}