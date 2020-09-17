package name.mjw.jquante.config.impl;

import java.beans.PropertyVetoException;
import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.util.Enumeration;
import java.util.HashMap;
import java.util.Map.Entry;

import org.w3c.dom.Document;
import org.w3c.dom.NamedNodeMap;
import org.w3c.dom.Node;

import name.mjw.jquante.common.EventListenerList;
import name.mjw.jquante.common.Utility;
import name.mjw.jquante.common.resource.StringResource;
import name.mjw.jquante.config.Configuration;
import name.mjw.jquante.config.Parameter;
import name.mjw.jquante.config.event.AtomInfoChangeEvent;
import name.mjw.jquante.config.event.AtomInfoChangeListener;

/**
 * The default AtomProperty configuration.
 * 
 * .. follows a singleton pattern. .. and an observer pattern for notifying the
 * registered classes of the changes at runtime.
 * 
 * @author V.Ganesh
 * @version 2.0 (Part of MeTA v2.0)
 */
public final class AtomInfo implements Configuration {

	private static AtomInfo atomInfo;

	private static final int DEFAULT_TABLE_SIZE = 90;

	private String element = "element";
	private String symbol = "symbol";

	/** Holds value of property nameTable. */
	private HashMap<String, String> nameTable;
	private HashMap<String, String> originalNameTable;

	/** Holds value of property atomicNumberTable. */
	private HashMap<String, Integer> atomicNumberTable;
	private HashMap<String, Integer> originalAtomicNumberTable;

	/** Holds value of property atomicWeightTable. */
	private HashMap<String, Double> atomicWeightTable;
	private HashMap<String, Double> originalAtomicWeightTable;

	/** Utility field used by event firing mechanism. */
	private EventListenerList<AtomInfoChangeListener> listenerList = null;

	/** Creates a new instance of AtomInfo */
	public AtomInfo() {
		nameTable = new HashMap<>(DEFAULT_TABLE_SIZE);
		atomicNumberTable = new HashMap<>(DEFAULT_TABLE_SIZE);
		atomicWeightTable = new HashMap<>(DEFAULT_TABLE_SIZE);


		originalNameTable = new HashMap<>(DEFAULT_TABLE_SIZE);
		originalAtomicNumberTable = new HashMap<>(DEFAULT_TABLE_SIZE);
		originalAtomicWeightTable = new HashMap<>(DEFAULT_TABLE_SIZE);


		// the initial parameters
		try {
			setDefaultParams();
		} catch (PropertyVetoException ignored) {
			// because it never should happen in this context
			System.err.println(ignored.toString());
		}
	}

	/**
	 * Obtain an instance of this ...
	 * 
	 * @return atomInfo AtomProperty configuration
	 */
	public static AtomInfo getInstance() {
		if (atomInfo == null) {
			atomInfo = new AtomInfo();
		}

		return atomInfo;
	}

	/**
	 * method to reset the instance of AtomInfo .. so that default values are
	 * loaded. Note that the the values are not saved on to a persistent store.
	 * 
	 */
	public static void reset() {
		if (atomInfo == null)
			return;

		// re-initilize the parameters
		try {
			atomInfo.setDefaultParams();

			AtomInfoChangeEvent changeEvent = new AtomInfoChangeEvent(atomInfo);

			changeEvent.setChangeType(AtomInfoChangeEvent.ALL_CANGED);
			// fire the event!
			atomInfo.fireAtomInfoChangeListenerAtomInfoChanged(changeEvent);
		} catch (PropertyVetoException ignored) {
			// because it never should happen in this context
			System.err.println(ignored.toString());
		} // end of try catch
	}

	/**
	 * Method to reset the instance of AtomInfo .. so that default values are loaded
	 * for a particular atom. Note that the values are not saved on to a persistent
	 * store.
	 * 
	 * @param symbol Atomic symbol
	 */
	public void resetValues(String symbol) {
		nameTable.put(symbol, originalNameTable.get(symbol));
		atomicNumberTable.put(symbol, originalAtomicNumberTable.get(symbol));
		atomicWeightTable.put(symbol, originalAtomicWeightTable.get(symbol));

		AtomInfoChangeEvent changeEvent = new AtomInfoChangeEvent(this);

		changeEvent.setChangeType(AtomInfoChangeEvent.ALL_CANGED);
		// fire the event!
		fireAtomInfoChangeListenerAtomInfoChanged(changeEvent);
	}

	/**
	 * Save user level configuration file
	 * 
	 * @throws IOException Signals that an I/O exception of some sort has occurred.
	 */
	public void saveUserConfigFile() throws IOException {
		// TODO: we are being a bit paranoid by replicating every thing, though
		// that is not really necessary
		StringResource strings = StringResource.getInstance();
		FileOutputStream fos = new FileOutputStream(strings.getUserAtomInfo());

		fos.write(strings.getXmlHeader().getBytes());

		fos.write("<atominfo> \n".getBytes());
		for (String key : nameTable.keySet()) {
			fos.write("\t<atom> \n".getBytes());
			fos.write(("\t\t<physical symbol=\"" + key + "\" name=\"" + nameTable.get(key) + "\" atomicNumber=\""
					+ atomicNumberTable.get(key) + "\" atomicWeight=\"" + atomicWeightTable.get(key)
					+ "\"/> \n").getBytes());
			fos.write("\t\t<display> \n".getBytes());
			fos.write("\t\t</display>\n".getBytes());
			fos.write("\t</atom> \n".getBytes());
		} // end for
		fos.write("</atominfo> \n".getBytes());

		fos.close();
	}

	/**
	 * private method to set the default parameters
	 */
	private void setDefaultParams() throws PropertyVetoException {
		StringResource strings = StringResource.getInstance();
		try {
			// read the internal XML config file
			Document configDoc = Utility.parseXML(getClass().getResourceAsStream(strings.getDefaultAtomInfo()));

			// and save the info. properly
			saveOriginal(configDoc);

			// read in user XML config file
			// and override the settings
			saveUser();
		} catch (Exception e) {
			throw new PropertyVetoException("Exception in " + "AtomInfo.setParameter()" + e.toString(), null);
		} // end of try .. catch block
	}

	/**
	 * Recursive routine save DOM tree nodes
	 */
	private void saveOriginal(Node n) {
		int type = n.getNodeType(); // get node type

		switch (type) {
		case Node.ATTRIBUTE_NODE:
			String nodeName = n.getNodeName();

			if (nodeName.equals("name")) {
				originalNameTable.put(symbol, n.getNodeValue());
			} else if (nodeName.equals("atomicNumber")) {
				originalAtomicNumberTable.put(symbol, Integer.valueOf(n.getNodeValue()));
			} else if (nodeName.equals("atomicWeight")) {
				originalAtomicWeightTable.put(symbol, Double.valueOf(n.getNodeValue()));
			}
			break;
		case Node.ELEMENT_NODE:
			element = n.getNodeName();

			NamedNodeMap atts = n.getAttributes();
			if (element.equals("physical")) {
				symbol = atts.getNamedItem("symbol").getNodeValue();

				// save the others
				for (int i = 0; i < atts.getLength(); i++) {
					Node att = atts.item(i);
					saveOriginal(att);
				} // end for
			} else {
				if (atts == null)
					return;

				for (int i = 0; i < atts.getLength(); i++) {
					Node att = atts.item(i);
					saveOriginal(att);
				}
			}

			break;
		default:
			break;
		}

		// save children if any
		for (Node child = n.getFirstChild(); child != null; child = child.getNextSibling()) {
			saveOriginal(child);
		}
	}

	/**
	 * Recursive routine save DOM tree nodes
	 */
	private void saveUserNode(Node n) {
		int type = n.getNodeType(); // get node type

		switch (type) {
		case Node.ATTRIBUTE_NODE:
			String nodeName = n.getNodeName();

			if (nodeName.equals("name")) {
				nameTable.put(symbol, n.getNodeValue());
			} else if (nodeName.equals("atomicNumber")) {
				atomicNumberTable.put(symbol, Integer.valueOf(n.getNodeValue()));
			} else if (nodeName.equals("atomicWeight")) {
				atomicWeightTable.put(symbol, Double.valueOf(n.getNodeValue()));

			}
			break;
		case Node.ELEMENT_NODE:
			element = n.getNodeName();

			NamedNodeMap atts = n.getAttributes();
			if (element.equals("physical")) {
				symbol = atts.getNamedItem("symbol").getNodeValue();

				// save the others
				for (int i = 0; i < atts.getLength(); i++) {
					Node att = atts.item(i);
					saveUserNode(att);
				}

			} else {
				if (atts == null)
					return;

				for (int i = 0; i < atts.getLength(); i++) {
					Node att = atts.item(i);
					saveUserNode(att);
				}
			}

			break;
		default:
			break;
		}

		// save children if any
		for (Node child = n.getFirstChild(); child != null; child = child.getNextSibling()) {
			saveUserNode(child);
		}
	}

	/** make a copy of the original */
	private void copyOriginals() {
		// just make a copy of original
		copyTableS(originalNameTable, nameTable);
		copyTableI(originalAtomicNumberTable, atomicNumberTable);
		copyTableD(originalAtomicWeightTable, atomicWeightTable);
	}

	/**
	 * Routine save DOM tree nodes
	 */
	private void saveUser() throws Exception {
		StringResource strings = StringResource.getInstance();

		try {
			if ((new File(strings.getUserAtomInfo())).exists()) {
				// if a user configuration file exists, then we copy
				// its contents
				Document configDoc = Utility.parseXML(strings.getUserAtomInfo());

				saveUserNode(configDoc);
			} else {
				copyOriginals();
			}
		} catch (Exception e) {
			System.err.println("Restricted mode copying.");
			copyOriginals();
		}
	}

	/**
	 * Copy one table to another (of strings)
	 */
	private void copyTableS(HashMap<String, String> src, HashMap<String, String> dest) {
		for (Entry<String, String> entry : src.entrySet()) {
			dest.put(entry.getKey(), entry.getValue());
		}
	}

	/**
	 * Copy one table to another (of ints)
	 */
	private void copyTableI(HashMap<String, Integer> src, HashMap<String, Integer> dest) {
		for (Entry<String, Integer> entry : src.entrySet()) {
			dest.put(entry.getKey(), entry.getValue());
		}
	}

	/**
	 * Copy one table to another (of doubles)
	 */
	private void copyTableD(HashMap<String, Double> src, HashMap<String, Double> dest) {
		for (Entry<String, Double> entry : src.entrySet()) {
			dest.put(entry.getKey(), entry.getValue());
		}
	}


	/**
	 * Method returns the value of parameter pertaining to the key.
	 * 
	 * @return Parameter - the parameter value.
	 * @throws NullPointerException if the key is not found
	 */
	@Override
	public Parameter getParameter(String key) {
		return new AtomProperty((String) nameTable.get(key), (atomicNumberTable.get(key)).intValue(),
				(atomicWeightTable.get(key)).doubleValue());
	}

	/**
	 * method to set the new parameter value.
	 * 
	 * @param key       - the key whose value needs to be changed or added
	 * @param parameter - the new parameter value
	 * @throws PropertyVetoException - incase changing property value is not
	 *                               supported
	 */
	@Override
	public void setParameter(String key, Parameter parameter) throws PropertyVetoException {
		AtomProperty ap = (AtomProperty) parameter;

		AtomInfoChangeEvent changeEvent = new AtomInfoChangeEvent(this);

		changeEvent.setChangeType(AtomInfoChangeEvent.ALL_CANGED);

		changeEvent.setOldValue(new AtomProperty(Utility.getString(nameTable, key),
				Utility.getInteger(atomicNumberTable, key), Utility.getDouble(atomicWeightTable, key))); // old val
		changeEvent.setNewValue(ap); // new val
		changeEvent.setAtomSymbol(key); // the atomic symbol

		// change to new value
		nameTable.put(key, ap.getName());
		atomicNumberTable.put(key, Integer.valueOf(ap.getAtomicNumber()));
		atomicWeightTable.put(key, Double.valueOf(ap.getAtomicWeight()));		
		// fire the event!
		fireAtomInfoChangeListenerAtomInfoChanged(changeEvent);
	}

	/**
	 * is the configuration stored as file?
	 * 
	 * @return boolean - true / false
	 */
	@Override
	public boolean isStoredAsFile() {
		return true;
	}

	/**
	 * Method returns the absolute path of the configuration file if
	 * isStoredAsFile() call results in a true value, else a null is returned.
	 * 
	 * @return String - the absolute path of the configuration file
	 */
	@Override
	public String getConfigFile() {
		return StringResource.getInstance().getDefaultAtomInfo();
	}

	/**
	 * changes current value of file storage
	 * 
	 * @param storedAsFile new value of this storage
	 */
	@Override
	public void setStoredAsFile(boolean storedAsFile) {
		// igonored! ... we don't want to listen here
	}

	/**
	 * sets the value of a parameter
	 * 
	 * @param configFile - the new name of the config file
	 * @throws PropertyVetoException some problem can't change the stuff!
	 */
	@Override
	public void setConfigFile(String configFile) throws PropertyVetoException {
		throw new PropertyVetoException("Cannot change this property!", null);
	}

	/**
	 * Getter for property nameTable.
	 * 
	 * @return Value of property nameTable.
	 * 
	 */
	public HashMap<String, String> getNameTable() {
		return this.nameTable;
	}

	/**
	 * Setter for property nameTable.
	 * 
	 * @param nameTable New value of property nameTable.
	 * 
	 */
	public void setNameTable(HashMap<String, String> nameTable) {
		this.nameTable = nameTable;
	}

	/**
	 * Getter for property atomicNumberTable.
	 * 
	 * @return Value of property atomicNumberTable.
	 * 
	 */
	public HashMap<String, Integer> getAtomicNumberTable() {
		return this.atomicNumberTable;
	}

	/**
	 * Setter for property atomicNumberTable.
	 * 
	 * @param atomicNumberTable New value of property atomicNumberTable.
	 * 
	 */
	public void setAtomicNumberTable(HashMap<String, Integer> atomicNumberTable) {
		this.atomicNumberTable = atomicNumberTable;
	}

	/**
	 * Getter for property atomicWeightTable.
	 * 
	 * @return Value of property atomicWeightTable.
	 * 
	 */
	public HashMap<String, Double> getAtomicWeightTable() {
		return this.atomicWeightTable;
	}

	/**
	 * Setter for property atomicWeightTable.
	 * 
	 * @param atomicWeightTable New value of property atomicWeightTable.
	 * 
	 */
	public void setAtomicWeightTable(HashMap<String, Double> atomicWeightTable) {
		this.atomicWeightTable = atomicWeightTable;
	}





	/**
	 * Registers AtomInfoChangeListener to receive events.
	 * 
	 * @param listener The listener to register.
	 * 
	 */
	public synchronized void addAtomInfoChangeListener(AtomInfoChangeListener listener) {
		if (listenerList == null) {
			listenerList = new EventListenerList<>();
		}
		listenerList.add(AtomInfoChangeListener.class, listener);
	}

	/**
	 * Removes AtomInfoChangeListener from the list of listeners.
	 * 
	 * @param listener The listener to remove.
	 * 
	 */
	public synchronized void removeAtomInfoChangeListener(AtomInfoChangeListener listener) {
		listenerList.remove(AtomInfoChangeListener.class, listener);
	}

	/**
	 * Notifies all registered listeners about the event.
	 * 
	 * @param event The event to be fired
	 * 
	 */
	private void fireAtomInfoChangeListenerAtomInfoChanged(AtomInfoChangeEvent event) {
		if (listenerList == null)
			return;

		// the great hack ... bypass notification if there was no real
		// change!
		if (event.getOldValue() != null) {
			if (event.getOldValue().equals(event.getNewValue()))
				return;
		}

		for (Object listener : listenerList.getListenerList()) {
			((AtomInfoChangeListener) listener).atomInfoChanged(event);
		}

	}

	/**
	 * Getter for property atomicNumber.
	 * 
	 * @param symbol - the atom symbol, IUPAC name.
	 * @return Value of property atomicNumber for the specified symbol
	 */
	public int getAtomicNumber(String symbol) {
		try {
			return atomicNumberTable.get(symbol);
		} catch (Exception e) {
			return atomicNumberTable.get("X");
		}
	}

	/**
	 * Getter for property symbol.
	 * 
	 * @param atomicNumber for the required symbol
	 * @return symbol - the atom symbol, IUPAC name. If no such atomic number is
	 *         found then "X" is returned.
	 */
	public String getSymbol(int atomicNumber) {
		Enumeration<Integer> eles = (Enumeration<Integer>) atomicNumberTable.values();
		Enumeration<String> keys = (Enumeration<String>) atomicNumberTable.keySet();

		while (eles.hasMoreElements()) {
			if (eles.nextElement() == atomicNumber)
				return keys.nextElement();
			keys.nextElement();
		}

		return "X";
	}

	/**
	 * Setter for property atomicNumber.
	 * 
	 * @param symbol       - the atom symbol, IUPAC name!
	 * @param atomicNumber New value of property atomicNumber.
	 */
	public void setAtomicNumber(String symbol, int atomicNumber) {
		AtomInfoChangeEvent changeEvent = new AtomInfoChangeEvent(this);

		changeEvent.setChangeType(AtomInfoChangeEvent.ATOMIC_NUMBER);
		changeEvent.setAtomSymbol(symbol);
		changeEvent.setOldValue(atomicNumberTable.get(symbol));

		atomicNumberTable.put(symbol, Integer.valueOf(atomicNumber));

		changeEvent.setNewValue(Integer.valueOf(atomicNumber));

		fireAtomInfoChangeListenerAtomInfoChanged(changeEvent);
	}

	/**
	 * Getter for property atomicWeight.
	 * 
	 * @param symbol - the atom symbol, IUPAC name!
	 * @return Value of property atomicWeight for the specified symbol
	 */
	public double getAtomicWeight(String symbol) {
		try {
			return atomicWeightTable.get(symbol);
		} catch (Exception e) {
			return atomicWeightTable.get("X");
		}
	}

	/**
	 * Setter for property atomicWeight.
	 * 
	 * @param symbol       - the atom symbol, IUPAC name!
	 * @param atomicWeight New value of property atomicWeight.
	 */
	public void setAtomicWeight(String symbol, double atomicWeight) {
		AtomInfoChangeEvent changeEvent = new AtomInfoChangeEvent(this);

		changeEvent.setChangeType(AtomInfoChangeEvent.ATOMIC_WEIGHT);
		changeEvent.setAtomSymbol(symbol);
		changeEvent.setOldValue(atomicWeightTable.get(symbol));

		atomicWeightTable.put(symbol, Double.valueOf(atomicWeight));

		changeEvent.setNewValue(Double.valueOf(atomicWeight));

		fireAtomInfoChangeListenerAtomInfoChanged(changeEvent);
	}





	/**
	 * Getter for property name.
	 * 
	 * @param symbol - the atom symbol, IUPAC name!
	 * @return Value of property name for the specified symbol
	 */
	public String getName(String symbol) {
		try {
			return nameTable.get(symbol);
		} catch (Exception e) {
			return nameTable.get("X");
		}
	}





}
