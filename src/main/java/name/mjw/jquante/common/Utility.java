package name.mjw.jquante.common;

import org.w3c.dom.Document;
import org.xml.sax.ErrorHandler;
import org.xml.sax.SAXException;
import org.xml.sax.SAXParseException;

import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStream;
import java.util.HashMap;

import javax.xml.parsers.DocumentBuilder;
import javax.xml.parsers.DocumentBuilderFactory;
import javax.xml.parsers.ParserConfigurationException;

/**
 * Utility class providing basic functionality such as ID generation etc. The
 * class it self is declared as final, and all the methods exposed are static.
 * 
 * @author V.Ganesh
 * @version 2.0 (Part of MeTA v2.0)
 */
public final class Utility {

	/**
	 * 1 a.u. = 0.52917724924 angstroms When converting from a.u. to angstroms,
	 * multiply by this factor. When converting from angstroms to a.u., divide by
	 * this factor.
	 */
	public static final double AU_TO_ANGSTROM_FACTOR = 0.52917724924;

	/**
	 * 1 kcal/mol = 0.00159362 AU.
	 */
	public static final double KCAL_TO_AU_FACTOR = 0.00159362;

	/**
	 * private constructor ... avoid object creation
	 */
	private Utility() {
	}

	/**
	 * Method to capitalise a string.
	 * 
	 * @param theString the string to be capitalised
	 * @return the string which is capitalised
	 * @throws NullPointerException           If string is null.
	 * @throws ArrayIndexOutOfBoundsException If string is greater than two
	 *                                        characters.
	 */
	public static String capitalise(final String theString) {
		final String firstChar = theString.substring(0, 1);
		final String restString = theString.substring(1, theString.length());

		if (restString == null)
			return firstChar.toUpperCase();
		else
			return (firstChar.toUpperCase() + restString);
	}

	/**
	 * A simple utility method to parse XML using a default DOM parser and return
	 * the results as abstract instance of Document class. An helper utility for
	 * <code>workspace</code> subsystem.
	 * 
	 * @param xmlFile : fully qualified name of the XML file to be parsed.
	 * @return Document : parsed XML document if no error occurs.
	 * 
	 * @throws ParserConfigurationException / SAXException in case of error reading
	 *                                      xml file or creating the parser. An
	 *                                      IOException in case the XML file is not
	 *                                      present!
	 * @throws org.xml.sax.SAXException     org.xml.sax.SAXException
	 * @throws java.io.IOException          java.io.IOException
	 */
	public static Document parseXML(final String xmlFile)
			throws ParserConfigurationException, SAXException, IOException {
		return parseXML(new FileInputStream(xmlFile));
	}

	/**
	 * A simple utility method to parse XML using a default DOM parser and return
	 * the results as abstract instance of Document class. An helper utility for
	 * <code>workspace</code> subsystem.
	 * 
	 * @param xmlStream : a stream connecting to XML data.-
	 * @return Document : parsed XML document if no error occurs.
	 * 
	 * @throws ParserConfigurationException / SAXException in case of error reading
	 *                                      xml file or creating the parser. An
	 *                                      IOException in case the XML file is not
	 *                                      present!
	 * @throws org.xml.sax.SAXException     org.xml.sax.SAXException
	 * @throws java.io.IOException          java.io.IOException
	 */
	public static Document parseXML(final InputStream xmlStream)
			throws ParserConfigurationException, SAXException, IOException {
		// we use the JAXP DOM parser for the job!
		// get an instance of the parser
		final DocumentBuilderFactory dbf = DocumentBuilderFactory.newInstance();
		final DocumentBuilder db = dbf.newDocumentBuilder();

		// set the error handler .. to detect errors(?) in xml
		// in normal situations, error should never occur because
		// the xml files parsed here are automatically generated ones
		db.setErrorHandler(new ErrorHandler() {
			@Override
			public void warning(final SAXParseException e) throws SAXException {
				System.err.println("Warning : " + e);
				e.printStackTrace();
			}

			@Override
			public void error(final SAXParseException e) throws SAXException {
				System.err.println("Error : " + e);
				e.printStackTrace();
			}

			@Override
			public void fatalError(final SAXParseException e) throws SAXException {
				System.err.println("Error : " + e);
				e.printStackTrace();
			}
		}); // setErrorHandler

		// and now, parse the input file:
		return db.parse(xmlStream);
	}

	/**
	 * method to get an integer out of Hashtable. If the key is not present, then
	 * zero is returned
	 * 
	 * @param table - the hashtable containing Integer objects
	 * @param key   - the key
	 * @return int - the int value intended. if the key is not found in the table a
	 *         zero is returned instead.
	 */
	public static int getInteger(final HashMap<?, ?> table, final Object key) {
		int value = 0;

		if (table.containsKey(key)) {
			value = ((Integer) table.get(key)).intValue();
		}

		return value;
	}

	/**
	 * method to get an double out of Hashtable. If the key is not present, then
	 * zero is returned
	 * 
	 * @param table - the hashtable containing Double objects
	 * @param key   - the key
	 * @return double - the double value intended. if the key is not found in the
	 *         table a zero is returned instead.
	 */
	public static double getDouble(final HashMap<?, ?> table, final Object key) {
		double value = 0.0;

		if (table.containsKey(key)) {
			value = ((Double) table.get(key)).doubleValue();
		}
		return value;
	}

	/**
	 * method to get an String out of Hashtable. If the key is not present, then a
	 * blank string is returned
	 * 
	 * @param table - the hashtable containing Double objects
	 * @param key   - the key
	 * @return String - the String value intended. if the key is not found in the
	 *         table a blank String is returned instead.
	 */
	public static String getString(final HashMap<?, ?> table, final Object key) {
		String value = " ";

		if (table.containsKey(key)) {
			value = (String) table.get(key);
		}
		return value;
	}
}
