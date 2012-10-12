package pdvisual;

import jv.geom.*;
import jv.object.PsObject;
import jv.project.PgGeometryIf;
import jv.project.PvCameraIf;
import jv.project.PvDisplayIf;
import jv.project.PjProject;
import jv.vecmath.PdVector;
import jv.vecmath.PiVector;
import jv.number.PdColor;
import Jama.Matrix;
import chemistry.*;
import crystallography.*;
//import com.cmcweb.db.postgresql.PostgresConnection;
//import com.cmcweb.db.*;
//import com.cmcweb.db.postgresql.*;
//import java.sql.Connection;
import java.util.*;

import jv.thirdParty.ruler.PgAxes;
import javax.swing.*;

import utility.Utility;

import java.awt.Color;
import java.awt.MouseInfo;
import java.awt.PointerInfo;
import java.awt.event.*;
import jv.project.PvGeometryIf;
import java.awt.Component;
import java.text.DecimalFormat;
import java.text.NumberFormat;


public class TernPD3DProj extends PjProject implements ActionListener, MouseMotionListener {
	static final long serialVersionUID = 1;

	protected	PgPolygonSet			m_vectors; 
	protected   PgElementSet	m_face_colors;
	
	protected IPDPlottable pdd;
	
	// a reference to the TernPD3D lets us get menu options from the user
	private TernPD3D m_ternpd3d;
	
	private IPDAnalyzer pda;
	
	private JPopupMenu entrySelectorMenu;
	private Map<JMenuItem, IComputedEntry> entriesOnMenu;
	
	private List<Double> m_energies;
	
	ColorWheelDialog chempotColorDialog = null;
	
	// We need some way to map indices of vertices in the graph to indices of 
	// entries or phases in the original (pseudo)pddata.  They're almost the same
	// thing! except that we want to remove some undesirable (see displayVertex())
	// points from the graph.  So, we keep a mapping.
	private int indxMapping[];
 //   private List<EntryVisualizer> visualizers;
	

	/**
	 * Constructors creates two argument vectors and a result vector. Initializations
	 * of instance variables should be performed within the init() method
	 * to allow resetting the project later by simply calling the init() method.
	 */
	public TernPD3DProj(IPDPlottable _pdd, TernPD3D ternpd3d) {
		super("Interactive Phase diagram");					// Super class constructors assigns name of this project.
		pdd = _pdd;
		m_ternpd3d = ternpd3d;
		// Create a new polygon set in 2D space to store the two argument vectors.
		m_vectors	= new PgPolygonSet(3);
		m_vectors.setName("Phase diagram");							// Set unique name of polygon set.
		m_face_colors = new PgElementSet(3);
		m_face_colors.setName("Face colors");
		
		// set up the visualizer
		entrySelectorMenu = new JPopupMenu();
		entriesOnMenu = new HashMap<JMenuItem, IComputedEntry>();
		
		m_ternpd3d.m_viewer.getDisplay().addPickListener(m_ternpd3d.selectPtsInfoPanel);
		
		((Component)(m_ternpd3d.m_viewer.getDisplay())).addMouseMotionListener(this);
		
		
		// Create a new polygon in 3D space to store the result vector.
		//m_result		= new PgPolygon(3);
		//m_result.setName("Result");								// Set unique name of polygon
		if (getClass() == TernPD3DProj.class)						// Invoke the init() method. This construct refuses to call
		  init();														// the init method if this constructor is invoked by a subclass.
	}
	/**
	 * Initialize all instance variables of this project. This method
	 * is invoked during the construction of this project and whenever
	 * this project is resetted. Therefore, allocations of instance variables
	 * should be performed in the constructor of the project.
	 */
	public void init() {
		super.init();		
		
		m_vectors.init();		// Initialize the argument vectors.
		m_face_colors.init();
		List<Double[]> compFracs = pdd.getCompositionFractions();
		
		// integer indicating the element whose chemical potential we
		// should use to color each facet or -1 if none
		int colorFaces = m_ternpd3d.getColorByChempots();

		Element colorElement = null;
        if (colorFaces !=-1)
            colorElement = pdd.getElements().get(colorFaces);
		// initialize our index mapping
		indxMapping = new int[compFracs.size()];
		for (int i = 0; i < indxMapping.length; i++)
			indxMapping[i] = i;
		
		pda = PDAnalyzerFactory.getIPDAnalyzer(pdd);

		if (m_ternpd3d.getUseFormationEnergies())
			m_energies = pdd.getFormEnergiesPerAtom();
		else
			m_energies = pdd.getEnergiesPerAtom();

		List<String> names = pdd.getVertexLabels();
		m_vectors.setNumVertices(compFracs.size());
		if (colorFaces != -1) {
			addGeometry(m_face_colors);
			m_face_colors.setNumElements(compFracs.size());
			m_face_colors.assureElementColors();
			m_face_colors.showElementColors(true);	
			m_face_colors.showEdges(true);
			m_face_colors.showVertices(false);
		} else {
			removeGeometry(m_face_colors);
		}
		
		// this is just to initialize the font field in the m_vectors geometry.
		// the actual format of the string tha t setLabelFont expect seems to be nowhere documented,
		// but after initializing the field, we can use, e.g. setLabelSize()
		m_vectors.setLabelFont(PvGeometryIf.GEOM_ITEM_POINT, "");
		
		// Enable usage of individual polygon colors
		m_vectors.showVertexColors(true);
		// Assure that polygon colors are allocated, needs to be called only once.
		m_vectors.assureVertexColors();
		
		// if necessary, calculate maximum distance of displayed points above hull
		// to give us a reference for the coloring
		double maxDistAboveHull = 0.0;
		if (m_ternpd3d.getColorSpectrum()) {
			for (int j = 0; j < compFracs.size(); j++) {
				if (displayVertex(j, m_energies) && pda.getEnergyPerAtomAboveHull(j) > maxDistAboveHull)
					maxDistAboveHull = pda.getEnergyPerAtomAboveHull(j);
			}
			if (maxDistAboveHull <= 0) // don't divide by 0
				maxDistAboveHull = 1;
		}
		
		int dim = pdd.getDimension();
		
		// show element labels if we want to
		m_face_colors.showElementLabels(m_ternpd3d.getNumberChempots());
		m_face_colors.setTransparency(0.35);
		m_face_colors.showTransparency(true);
		
		// add a vertex at each for each entry
		for (int i = 0; i < compFracs.size(); i++) {
			double c[] = Utility.todouble(compFracs.get(i));
			double[] coord;
			// handle different dims
			if (dim == 2) {
				coord = c;
				m_vectors.setVertex(i, coord[0], m_energies.get(i), 0.0);
				if (colorFaces != -1)
					m_face_colors.setVertex(i, coord[0], m_energies.get(i), 0.0);
			} else if (dim == 3) {
				coord = getTernaryXYCoord(c);
				m_vectors.setVertex(i, coord[0], coord[1], m_energies.get(i));
				if (colorFaces != -1)
					m_face_colors.setVertex(i, coord[0], coord[1], m_energies.get(i));
			} else if (dim == 4) {
				coord = getQuaternaryXYZCoord(c);
				m_vectors.setVertex(i, coord[0], coord[1], coord[2]);
				if (colorFaces != -1)
					m_face_colors.setVertex(i, coord[0], coord[1], coord[2]);
			}
			// set the vertex labels 
			m_vectors.getVertex(i).setName("");
			if (m_ternpd3d.getShowLabelsReal() && pdd.isPurePhase(i))
				m_vectors.getVertex(i).setName(names.get(i));
			if (m_ternpd3d.getShowLabelsPseudo() && ! pdd.isPurePhase(i))
				m_vectors.getVertex(i).setName(names.get(i));
				

			// color vertices
			if (m_ternpd3d.getColorBinary()) {
				// make stable entries blue, the unstable red
				if (pdd.getIndxUnstableEntries().contains(i))
					m_vectors.setVertexColor(i, Color.BLUE);
				else
					m_vectors.setVertexColor(i, Color.RED);
			} else if (m_ternpd3d.getColorSpectrum()) {
				// make the gradient
				System.out.println(maxDistAboveHull);
				double blueness = (maxDistAboveHull - pda.getEnergyPerAtomAboveHull(i)) / maxDistAboveHull;
				if (blueness > 1) // numerical problem or undisplayed vertex
					blueness = 1;
				else if (blueness < 0)
					blueness = 0;
				m_vectors.setVertexColor(i, new Color((float)(1-blueness),0,(float)(blueness)));
			} else if (m_ternpd3d.getColorType()) {
				if ( ((PDData)pdd).getEntry(i).getClass().getName().contains("ManualComputedEntry")) {
					m_vectors.setVertexColor(i, Color.GREEN);
				} else { 
					m_vectors.setVertexColor(i, Color.BLUE);
				}
			} else { // set them all to the same color
				m_vectors.setVertexColor(i, Color.RED);
			}
		}
	
		// show vertex labels (if the option is set not to, we set the label to "")
		m_vectors.showVertexLabels(true);

		// show axes if we want to
		configureAxes();
		
		// add edges between the facets
		int counter = 0;
		m_vectors.setNumPolygons(pdd.getIndxFacets().size());
		List<List<Integer>> indxFacets = pdd.getIndxFacets();
		// color faces?
		if (chempotColorDialog != null) 
			chempotColorDialog.setVisible(false);
		if (colorFaces == -1) {
			chempotColorDialog = null;
		} else {
			List<Double> chempots = new ArrayList<Double>();
            for (List<Integer> facet : pdd.getIndxFacets())
				chempots.add(pda.getChemicalPotential(colorElement, facet));
			chempotColorDialog = new ColorWheelDialog(m_ternpd3d, chempots);
		}
		
		for (List<Integer> facet : indxFacets) {
			if (facet.size() == 2) {
				if (colorFaces != -1) {
                    m_face_colors.setElement(counter, new PiVector(facet.get(0),facet.get(1)));
					colorFace(counter, pda.getChemicalPotential(colorElement, facet));
				}
				m_vectors.setPolygon(counter++, new PiVector(facet.get(0),facet.get(1)));
			} else if (facet.size() == 3) {
				if (colorFaces != -1) {
                    m_face_colors.setElement(counter, new PiVector(facet.get(0),facet.get(1),facet.get(2)));
					colorFace(counter, pda.getChemicalPotential(colorElement, facet));
				}
				m_vectors.setPolygon(counter++, new PiVector(facet.get(0),facet.get(1),facet.get(2)));	
			} else if (facet.size() == 4) {
				// we actually need to draw 4 triangles here, to make a tetrahedron
				for (int i = 0; i < 4; i++) {
					List<Integer> f = new LinkedList<Integer>(facet);
					f.remove(i);
					if (colorFaces != -1) {
                    	m_face_colors.setElement(counter, new PiVector(f.get(0),f.get(1),f.get(2)));
						colorFace(counter, pda.getChemicalPotential(colorElement, facet));
					}
					m_vectors.setPolygon(counter++, new PiVector(f.get(0),f.get(1),f.get(2)));
				}
			}
		}

		// if we're in greater dimension than a binary, add lines between the basis phases
		if (dim > 2) {
			List<Integer> basisIndxs = pdd.getBasisEntryIndxs();
			for (Integer i : basisIndxs) {
				for (Integer j : basisIndxs) {
					if (i != j) {
						m_vectors.setPolygon(counter++, new PiVector(i, j));
					}
				}
			}
		}
		
		// set label look
		//m_vectors.setLabelAttribute(PvGeometryIf.GEOM_ITEM_POINT, 0, 0, PgGeometryIf.LABEL_CENTER, PgGeometryIf.LABEL_BASE, PsConfig.FONT_TEXT);
		if (m_ternpd3d.getBigMode()) {
			m_vectors.setGlobalPolygonSize(2.5);
			m_vectors.setGlobalVertexSize(6.);
			m_vectors.setLabelSize(PvGeometryIf.GEOM_ITEM_POINT,20.);
		} else {
			m_vectors.setGlobalPolygonSize(1);
			m_vectors.setGlobalVertexSize(4.);
			m_vectors.setLabelSize(PvGeometryIf.GEOM_ITEM_POINT,12.);
		}

		/*m_vectors.setVertex(0, 0., 0.);							// Base point (index = 0) of both vectors
		m_vectors.setVertex(1, .5, 1.);							// Endpoint (index = 1) of first vector
		m_vectors.setVertex(2, -2., 2.);							// Endpoint (index = 2) of second vector
		*/
		/*
		m_vectors.setNumPolygons(2);								// Allocate two lines being the two argument vectors
		// Determine the index of the basepoint (0) and endpoint (1) of the first vector
		m_vectors.setPolygon(0, new PiVector(0, 1));
		// Determine the index of the basepoint (0) and endpoint (2) of the second vector
		m_vectors.setPolygon(1, new PiVector(0, 2));
		m_vectors.getPolygon(0).setName("v");					// Set name of first argument vector
		m_vectors.getPolygon(1).setName("w");					// Set name of second argument vector
		m_vectors.showPolygonLabels(true);						// Show name of argument vectors at tip
		m_vectors.setGlobalPolygonColor(Color.blue);			// Set color of argument vectors
		m_vectors.setGlobalPolygonSize(2.);						// Set thickness of argument vectors
		m_vectors.showPolygonEndArrow(true); 
		*/
		
		// don't display points if user doesn't want to.  update the mapping accordingly.
		for (int i = m_vectors.getNumVertices() - 1; i >= 0; i--)
			if (!displayVertex(i,m_energies)) {
				m_vectors.removeVertex(i);
				updateIndxMapping(i);
			}
		
		
		m_vectors.setParent(this);									// Register the argument vectors as children
		if (colorFaces != -1)
			m_face_colors.update(m_face_colors);
		
	//	m_vectors.removeMarkedPolygons();
																			// allows to catch their update events in update(Object)
																			// and initiate recomputation of the result vector.
		//m_result.setNumVertices(2);								// Allocate place for both endpoints of result vector.
		//m_result.setGlobalEdgeColor(Color.red);					// Set color of result vector
		//m_result.setGlobalEdgeSize(2.);							// Set thickness of result vector
		//m_result.showPolygonEndArrow(true);						// Enable drawing of arrow of result vector
										// Set initial computation mode.

 //       visualizers = new ArrayList<EntryVisualizer>();
	}

//    public void addEntryVisualizer(EntryVisualizer visualizer)
 //   {
 //       visualizers.add(visualizer);
 //   }

	// update our indxMapping due to the removal from the graph of point i
	private void updateIndxMapping(int i) {
		for (int j = i; j < indxMapping.length - 1; j++)
			indxMapping[j] = indxMapping[j+1];
	}

	
	private void colorFace(int i, Double k) {

	/*	PdColor hsvCol = new PdColor();
		hsvCol.setColor(PdColor.hsv2rgb((int)(127+128*Math.acos(-k)/Math.PI),
				 (int)(255*(1.-Math.abs(k))), 255));
		hsvCol.enableAlpha(true);
		hsvCol.setAlpha(0.99);
		*/
		m_face_colors.setElementColor(i, chempotColorDialog.getColor(k));
		if(m_ternpd3d.getNumberChempots()) {
			m_face_colors.setLabelColor(i, Color.black);
			NumberFormat formatter = new DecimalFormat("0.###");
			m_face_colors.getElement(i).setName(formatter.format(k));
		}
	}
	
	private boolean displayVertex(int i, List<Double> energies) {
		// don't display positive energy points if user doesn't want to
		if ( (!m_ternpd3d.getDispPosEPoints() && energies.get(i) > PDBuilder.EFORM_TOL)
				// don't display things which are too unstable
			|| pda.getEnergyPerAtomAboveHull(i) > PDBuilder.EFORM_TOL + m_ternpd3d.getMaxDistFromHull()
			|| (m_ternpd3d.getMaxDistFromHull() == 0 && pdd.getIndxUnstableEntries().contains(i)))
			return false;
		return true;
	}
	
	public int mapIndxsGraphToPDD(int i) {
		return indxMapping[i];
	}
	
	private void configureAxes() {
		PvDisplayIf disp = getDisplay();							// Get the display from this project
		if (disp != null) {
			if (m_ternpd3d.getShowAxes()) {
				PgAxes axes = ((jv.viewer.PvDisplay)disp).getAxes();
				if (axes == null)
					axes = new PgAxes(3);
				axes.setMode(PgAxes.AXES_3DCENTRAL);
				axes.setEnabledAutoHashing(true);
				axes.setEnabledAutoLayout(true);
				// set labels to empty
				List<String> axisLabels = pdd.getAxisLabels();
				String labels[] = new String[axisLabels.size()];
				for (int i = 0; i < labels.length; i++)
					labels[i] = ""; // no labels for the composition axes
				labels[labels.length - 1] = "Energy (eV/atom)";
				axes.setNames(labels);
					
				((jv.viewer.PvDisplay)disp).setAxes(axes);
				disp.showAxes(true);
			} else {
				disp.showAxes(false);
			}
		}
	}
	
	public void start() {
		//addGeometry(m_result);									// Register result vector in display
		addGeometry(m_vectors);										// Register two argument vectors in display
		selectGeometry(m_vectors);									// Set two vectors as the active geometry in display
																		// to be able to pick their vertices.
		PvDisplayIf disp = getDisplay();							// Get the display from this project
		if (disp != null) {											// and set some display options.
			//disp.showGrid(true);									// Show grid in the background
			disp.setEnabledBoxRatio(true);
			// default to projection onto xy plane unless we're looking at a quaternary
			if (pdd.getDimension() < 4)
				disp.selectCamera(PvCameraIf.CAMERA_ORTHO_XY);	
			else
				disp.selectCamera(PvCameraIf.CAMERA_PERSPECTIVE);
			
			disp.setMajorMode(PvDisplayIf.MODE_ORBIT);				
			// set up some axes
			configureAxes();

		}
		super.start();
	}

	/**
	 * Perform recomputations whenever a child has changed. Here the two argument
	 * vectors have registered this class as parent, therefore, we are able to catch
	 * their update events and recompute the result vectors if they have changed.
	 * Method is usually invoked from the children.
	 */
	public boolean update(Object event) {
//		if (event == m_vectors) {									// If the two argument vector have changed...
//			computeResult(m_mode);									// Recompute the result vector based on the current mode.
		//	m_result.update(m_result);								// Redraw result vector in display.
//			return true;												// Event successfully handled, just return.
//		}
		// If we do not know about the event then just forward it to the superclass.
//		return super.update(event);
		return true;
	}
	
	// runs when we click on a choice in the entry selector menu
	public void actionPerformed(ActionEvent e) {
		Object source = e.getSource();
		
		if (source instanceof JMenuItem) {
			IComputedEntry ent = entriesOnMenu.get(source);
			
			//close menu after click
			entrySelectorMenu.setVisible(false);
			JOptionPane.showMessageDialog(null, ent.toString());
			new MyJmolViewer(ent);
//			for (EntryVisualizer visualizer : visualizers)
 //           {
 //               visualizer.setEntry(ent);
 //               visualizer.doVis();
 //           }
		} 
	}
	

	/**
	 * Method is called when a vertex is initial picked. We pop up an information dialog
	 * about the corresponding structure.
	 */
	public void pickVertex(PgGeometryIf geom, int indxGraph, PdVector vertex) {
			
	// return if we didn't click on a real point
	//	int indx = pos.getVertexInd();
		int indx = mapIndxsGraphToPDD(indxGraph);
		if (indx == -1) {
			return;
		}
		
		// otherwise, make the menu:    
		// first clear it:
		entrySelectorMenu.removeAll();
		entriesOnMenu.clear();
		// add composition as title
		JMenuItem title = new JMenuItem(pdd.getVertexLabels().get(indx));
		entrySelectorMenu.add(title);
		entrySelectorMenu.addSeparator();
		
		for (Integer i : pdd.getPhasesWSameCompAs(indx)) {
			// if pdd is a pseudopdd, we're dealing with mixed phases, and we need a submenu
			 if (pdd instanceof PseudoPDData) {  
				 MixedPhase m = ((PseudoPDData)pdd).getAllPhases().get(i);
				 JMenu subMenu = new JMenu("Phase: " + m.getEnergyPerAtom() + " eV");
				 entrySelectorMenu.add(subMenu);
				 for (IComputedEntry e : m.getPhases().keySet()) {
					 String percent = new Double(m.getPhases().get(e)).toString();
					 JMenuItem menuItem = new JMenuItem("  - " + percent + "%: " + e.getComposition().toString() + " : " + e.getEnergyPerAtom() + " eV");
					 menuItem.addActionListener(this);
					 entriesOnMenu.put(menuItem, e);
					 //subMenu.add(menuItem);
					 entrySelectorMenu.add(menuItem);
				 }
			 } else if (pdd instanceof PDData) {
			     String en = String.format("%9.6g", m_energies.get(i));
			     String eah = String.format("%8.5g", pda.getEnergyPerAtomAboveHull(i));
				 JMenuItem menuItem = new JMenuItem( en + " eV ( "+ eah + " )");
				 entriesOnMenu.put(menuItem, ((PDData)pdd).getEntry(i));
			     menuItem.addActionListener(this);
			     entrySelectorMenu.add(menuItem);
			 }
        }
        
        // show it in the right place
        PointerInfo pInfo = MouseInfo.getPointerInfo();
        entrySelectorMenu.show(null, (int)(pInfo.getLocation().getX()), (int)(pInfo.getLocation().getY()));	
	}
	
    /**
     * returns ternary XY coordinates for our composition with respect to
     * particular elements.
     * @param basis
     * @return
     */ /*
    public double[] getTernaryXYCoord(double[] rectangularCoords){
        if(rectangularCoords.length != 3) 
        	System.out.println("Error: in Composition.getTernaryXYCoord basis.length != 3");
        
        double sqrt3o2 = 0.5 * Math.sqrt(2);
        
        double[] abc2xy = new double[2];
        abc2xy[0] = (1.0-rectangularCoords[0])*Math.cos(Math.PI/6)/sqrt3o2 - 0.5;
        abc2xy[0] += 0.5-(1.0-rectangularCoords[1])*Math.cos(Math.PI/6)/sqrt3o2;
        abc2xy[0] *= 0.5;
        abc2xy[1] = rectangularCoords[2]*sqrt3o2;
        return abc2xy;
    } */
    
    public double[] getTernaryXYCoord(double[] rectangularCoords){
        if(rectangularCoords.length != 3) 
        	System.out.println("Error: in Composition.getTernaryXYCoord basis.length != 3");
        
        double sqrt3o2 = 0.5 * Math.sqrt(2);
        
        double[][] F = {{ 1.0, 0.5,     0.0}
        				,{0.0, sqrt3o2, 0.0}
        				,{0.0, 0.0,     0.0}}; 
        Matrix MF = new Matrix(F);
        
        double[][] B = new double[rectangularCoords.length][1];
        for (int i = 0; i < B.length; i++)
        	B[i][0] = rectangularCoords[i];
        Matrix MB = new Matrix(B);
        
        Matrix r = MF.times(MB);
        return r.getRowPackedCopy();
    }
    
    public double[] getQuaternaryXYZCoord(double[] rectangularCoords) {
        if(rectangularCoords.length != 4) 
        	System.out.println("Error: in Composition.getQuaternaryXYZCoord basis.length != 4");
        
        double sqrt3o2 = 0.5 * Math.sqrt(2);
        double sqrt3o4 = 0.25 * Math.sqrt(2);
        
        double[][] F = {{ 1.0, 0.5,     0.0, 0.5}
        				,{0.0, sqrt3o2, 0.0, sqrt3o4}
        				,{0.0, 0.0,     0.0, sqrt3o2}}; 
        Matrix MF = new Matrix(F);
        
        double[][] B = new double[rectangularCoords.length][1];
        for (int i = 0; i < B.length; i++)
        	B[i][0] = rectangularCoords[i];
        Matrix MB = new Matrix(B);
        
        Matrix r = MF.times(MB);
        return r.getRowPackedCopy();
    }
    
    public  void 	dragVertex(PgGeometryIf geom, int index, PdVector vertex) {
    	// don't let user move vertices!
    	init();
    }
    
    // implement the focus listener for the entry selector menu.  basically, we wanna
    // make it go away when it loses focus
    public void mouseDragged(MouseEvent me) { /*nothin*/ }
    public void mouseMoved(MouseEvent me) {	
    	// close the menu when we drag over the graph
		if (entrySelectorMenu.isVisible()) {
			entrySelectorMenu.setVisible(false);
			return;
		}
	}

}
