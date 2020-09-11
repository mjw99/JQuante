package name.mjw.jquante.math.qm.basis;

import name.mjw.jquante.common.Utility;

import org.w3c.dom.Document;
import org.w3c.dom.NamedNodeMap;
import org.w3c.dom.Node;
import org.xml.sax.SAXException;

import java.io.IOException;
import java.lang.ref.WeakReference;

import javax.xml.parsers.ParserConfigurationException;

/**
 * This class provides the means to read a basis set stored in XML format and
 * convert it into appropriate data structure representation.
 * 
 * Follows a singleton pattern.
 * 
 * @author V.Ganesh
 * @version 2.0 (Part of MeTA v2.0)
 */
public class BasisSetReader {

	private static WeakReference<BasisSetReader> _basisSetReader = null;

	private BasisSet basisSet;
	private AtomicBasis atomicBasis;
	private Orbital orbital;

	/** Creates a new instance of BasisSetReader */
	private BasisSetReader() {
	}

	/**
	 * Get an instance (and the only one) of BasisSetReader
	 * 
	 * @return BasisSetReader instance
	 */
	public static BasisSetReader getInstance() {
		if (_basisSetReader == null) {
			_basisSetReader = new WeakReference<>(new BasisSetReader());
		}

		BasisSetReader basisSetReader = _basisSetReader.get();

		if (basisSetReader == null) {
			basisSetReader = new BasisSetReader();
			_basisSetReader = new WeakReference<>(basisSetReader);
		}

		return basisSetReader;
	}

	/**
	 * Read a particular basis for the basis library (XML)
	 * 
	 * @param basisSetName
	 *            the name of basis.
	 * @return BasisSet object, representing the requested basis set.
	 * @throws IOException if stream from basisName cannot be read.
	 * @throws SAXException if there was a problem with parsing the XML.
	 * @throws ParserConfigurationException if there is a serious configuration error.
	 */
	public BasisSet readBasisSet(final String basisSetName)
			throws ParserConfigurationException, SAXException, IOException {
		final String basisLibraryPath = "/name/mjw/jquante/math/qm/basis/basis_";
		// read the XML config file
		final Document basisDoc = Utility
				.parseXML(getClass().getResourceAsStream(basisLibraryPath + basisSetName + ".xml"));

		// and save the basis info. properly
		saveIt(basisDoc);

		return basisSet;
	}

	/**
	 * Read a particular basis set from the basis set library (XML)
	 * 
	 * @param basisSetFileName the name of external basis set file name.
	 * @return BasisSet object, representing the requested basis set.
	 * @throws IOException                  if stream from basisFileName cannot be
	 *                                      read.
	 * @throws SAXException                 if there was a problem with parsing the
	 *                                      XML.
	 * @throws ParserConfigurationException if there is a serious configuration
	 *                                      error.
	 */
	public BasisSet readExternalBasisSet(final String basisSetFileName)
			throws ParserConfigurationException, SAXException, IOException {
		// read the XML basis set file
		final Document basisDoc = Utility.parseXML(getClass().getResourceAsStream(basisSetFileName));

		// and save the basis info. properly
		saveIt(basisDoc);

		return basisSet;
	}

	/**
	 * Recursive routine save DOM tree nodes ... and hoping that it reads the basis
	 * properly
	 */
	private void saveIt(final Node n) {
		final int type = n.getNodeType(); // get node type

		switch (type) {
			case Node.ATTRIBUTE_NODE:
				final String nodeName = n.getNodeName();

				if (nodeName.equals("name")) {
					// instance of a new basis set
					basisSet = new BasisSet(nodeName);
				}

				break;
			case Node.ELEMENT_NODE:
				final String element = n.getNodeName();

				final NamedNodeMap atts = n.getAttributes();

				if ("atom".equals(element)) {
					// a new atomic basis
					atomicBasis = new AtomicBasis(atts.getNamedItem("symbol").getNodeValue(),
							Integer.parseInt(atts.getNamedItem("atomicNumber").getNodeValue()));
					basisSet.addAtomicBasis(atomicBasis);
				} else if ("orbital".equals(element)) {
					// a orbital entry for atomic basis
					orbital = new Orbital(atts.getNamedItem("type").getNodeValue());
					atomicBasis.addOrbital(orbital);
				} else if ("entry".equals(element)) {
					// a orbital (coefficient, exponent) entry
					orbital.addEntry(Double.parseDouble(atts.getNamedItem("coeff").getNodeValue()),
							Double.parseDouble(atts.getNamedItem("exp").getNodeValue()));
				} else {
					if (atts == null)
						return;

					for (int i = 0; i < atts.getLength(); i++) {
						final Node att = atts.item(i);
					saveIt(att);
				}
			}

			break;
		default:
			break;
		}

		// save children if any
		for (Node child = n.getFirstChild(); child != null; child = child.getNextSibling()) {
			saveIt(child);
		} // end for
	}
}
