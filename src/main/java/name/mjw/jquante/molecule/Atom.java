package name.mjw.jquante.molecule;

import java.util.ArrayList;
import java.util.Enumeration;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Objects;

import org.hipparchus.geometry.euclidean.threed.Vector3D;

import name.mjw.jquante.common.Units;
import name.mjw.jquante.common.Utility;

/**
 * This class define the structure of an Atom (programmatically of course! ;) ).
 * 
 * @author V.Ganesh
 * @version 2.0 (Part of MeTA v2.0)
 */
public class Atom implements Cloneable {
	/** Holds value of property symbol. */
	private String symbol;

	/** Holds value of property atomCenter. */
	private Vector3D atomCenter;

	/**
	 * Holds information on all the connected atoms and the cardinalities, the
	 * atomIndices to which the present atom is connected serves as the
	 * <I>key</I> of the table and the value contains object of class BondType.
	 */
	private HashMap<Integer, BondType> connectedList;

	/**
	 * Holds information on the ZMatrix representation of this atom in the
	 * molecule. This can either contain zero, one, two or three entries
	 * indicating bond length, angle or dihedral connectivity with other atoms
	 * in the molecular graph.
	 */
	private ArrayList<ZMatrixItem> zMatrixElement;

	/** Holds value of property index. */
	private int index;

	/**
	 * Creates a new instance of Atom
	 * 
	 * @param symbol
	 *            The atom symbol
	 * @param atomCenter
	 *            The nuclear center of the atom in Cartesian coordinates
	 */
	public Atom(String symbol, Vector3D atomCenter) {
		this(symbol, atomCenter, new HashMap<>(1),
				null, 0);
	}

	/**
	 * Creates a new instance of Atom
	 * 
	 * @param symbol
	 *            The atom symbol
	 * @param atomCenter
	 *            The nuclear center of the atom in Cartesian coordinates
	 * @param atomIndex
	 *            The atom index of this atom
	 */
	public Atom(String symbol, Vector3D atomCenter, int atomIndex) {
		this(symbol, atomCenter, new HashMap<>(1),
				null, atomIndex);
	}

	/**
	 * Creates a new instance of Atom
	 * 
	 * @param symbol
	 *            The atom symbol
	 * @param atomCenter
	 *            The nuclear center of the atom in Cartesian coordinates
	 * @param connectedList
	 *            The connected list
	 * @param zMatrixElement
	 *            ZMatrix elements
	 * @param atomIndex
	 *            the atom index of this atom.
	 */
	private Atom(String symbol, Vector3D atomCenter,
			HashMap<Integer, BondType> connectedList,
			ArrayList<ZMatrixItem> zMatrixElement, int atomIndex) {
		this.symbol = Utility.capitalise(symbol.toLowerCase());
		this.index = atomIndex;
		this.atomCenter = atomCenter;

		this.connectedList = connectedList;
		this.zMatrixElement = zMatrixElement;

		// by default all atom centers are in angstrom units
		atomCenterUnits = Units.ANGSTROM;
	}

	/** initialize the ZMatrix Items */
	private void initZMatrixItems() {
		if (this.zMatrixElement == null) {
			this.zMatrixElement = new ArrayList<>(3);

			this.zMatrixElement.add(new ZMatrixItem(null, 0.0));
			this.zMatrixElement.add(new ZMatrixItem(null, 0.0));
			this.zMatrixElement.add(new ZMatrixItem(null, 0.0));
		}
	}

	/**
	 * Getter for property symbol.
	 * 
	 * @return Value of property symbol.
	 */
	public String getSymbol() {
		return this.symbol;
	}

	/**
	 * Setter for property symbol.
	 * 
	 * @param symbol
	 *            New value of property symbol.
	 */
	public void setSymbol(String symbol) {
		this.symbol = symbol;
	}

	/**
	 * Getter for property atomCenter.
	 * 
	 * @return Value of property atomCenter.
	 */
	public Vector3D getAtomCenter() {
		return this.atomCenter;
	}

	/**
	 * Getter for property atomCenterInAU.
	 * 
	 * @return Value of property atomCenterInAU.
	 */
	public Vector3D getAtomCenterInAU() {
		return this.atomCenter.scalarMultiply((atomCenterUnits == Units.AU)
				? 1.0
				: 1.0 / Utility.AU_TO_ANGSTROM_FACTOR);
	}

	protected Units atomCenterUnits;

	/**
	 * Get the value of atomCenterUnits
	 * 
	 * @return the value of atomCenterUnits
	 */
	public Units getAtomCenterUnits() {
		return atomCenterUnits;
	}

	/**
	 * Set the value of atomCenterUnits
	 * 
	 * @param atomCenterUnits
	 *            new value of atomCenterUnits
	 */
	public void setAtomCenterUnits(Units atomCenterUnits) {
		this.atomCenterUnits = atomCenterUnits;
	}

	/**
	 * Setter for property atomCenter.
	 * 
	 * @param atomCenter
	 *            New value of property atomCenter.
	 */
	public void setAtomCenter(Vector3D atomCenter) {
		this.atomCenter = atomCenter;
	}

	/**
	 * Getter for property x.
	 * 
	 * @return Value of property x.
	 */
	public double getX() {
		return atomCenter.getX();
	}

	/**
	 * Getter for property y.
	 * 
	 * @return Value of property y.
	 */
	public double getY() {
		return atomCenter.getY();
	}

	/**
	 * Getter for property z.
	 * 
	 * @return Value of property z.
	 */
	public double getZ() {
		return atomCenter.getZ();
	}

	/**
	 * overloaded toString() method.
	 * 
	 * @return A description of Atom object
	 */
	@Override
	public String toString() {
		return String.format("%s%2d % f % f % f", symbol, index, atomCenter.getX(), atomCenter.getY(),
				atomCenter.getZ());
	}

	/**
	 * for elaborate information
	 * 
	 * @return An extended description!
	 */
	public String toExtendedString() {
		return symbol + " " + atomCenter.toString() + " "
				+ connectedList.toString();
	}


	/**
	 * Simple method to find the distance between two atom centers.
	 * 
	 * @param atom
	 *            - the atom from which the distance is to be found
	 * @return the distance between the two atom centers
	 */
	public double distanceFrom(Atom atom) {
		return distanceFrom(atom.atomCenter);
	}

	/**
	 * Simple method to find the distance between atom center and an point.
	 * 
	 * @param point
	 *            - the point to which the distance is to be found
	 * @return the distance between the point and the atom center
	 */
	public double distanceFrom(Vector3D point) {
		return atomCenter.distance(point);
	}

	/**
	 * i do some cloning business ;)
	 * 
	 * @throws CloneNotSupportedException
	 *             If that isn't possible
	 * @return A copy of the present object
	 */
	@Override
	public Object clone() throws CloneNotSupportedException {
		HashMap<Integer, BondType> theClonedConnection = new HashMap<>();

		Enumeration<Integer> keys = (Enumeration<Integer>) connectedList.keySet();
		Integer key;

		while (keys.hasMoreElements()) {
			key = keys.nextElement();
			theClonedConnection.put(key, connectedList.get(key));
		}

		ArrayList<ZMatrixItem> theColonedZMatrixElement = new ArrayList<>();

		try {
			for (ZMatrixItem item : zMatrixElement) {
				theColonedZMatrixElement.add((ZMatrixItem) item.clone());
			} // end for
		} catch (NullPointerException npe) {
			theColonedZMatrixElement = null;
		} // end of try .. catch block

		return new Atom(this.symbol, new Vector3D(this.atomCenter.toArray()), theClonedConnection,
				theColonedZMatrixElement, index);
	}

	/**
	 * Method to add a connection between the present atom and the atom
	 * specified by atomIndex as specified by the implimentation of Molecule
	 * interface.
	 * 
	 * @param atomIndex
	 *            the integer index of atom to be connected
	 * @param bondType
	 *            the bond type
	 */
	public void addConnection(int atomIndex, BondType bondType) {
		Integer key = Integer.valueOf(atomIndex);

		if (bondType.equals(BondType.NO_BOND) && connectedList.containsKey(key)) {
			connectedList.remove(key);
			return;
		}
		
		connectedList.put(key, bondType);
	}

	/**
	 * Method to remove connection between the present atom and the atom
	 * specified by atomIndex as specified by the implementation of Molecule
	 * interface.
	 * 
	 * @param atomIndex
	 *            the index of atom to be disconnected
	 */
	public void removeConnection(int atomIndex) {
		connectedList.remove(Integer.valueOf(atomIndex));
	}

	/**
	 * remove all the connections to this atom
	 */
	public void removeAllConnections() {
		connectedList = new HashMap<>();
	}

	/**
	 * Getter for property connectedList.
	 * 
	 * @return Value of property connectedList.
	 * 
	 */
	public HashMap<Integer, BondType> getConnectedList() {
		return this.connectedList;
	}

	/**
	 * Setter for property connectedList.
	 * 
	 * @param connectedList
	 *            New value of property connectedList.
	 * 
	 */
	public void setConnectedList(HashMap<Integer, BondType> connectedList) {
		this.connectedList = connectedList;
	}

	/**
	 * check if an atom index is in the connected list.
	 * 
	 * @param atomIndex
	 *            the index to be checked
	 * @return boolean - true / false indicating the present connectivity
	 *         status.
	 */
	public boolean isConnected(int atomIndex) {
		return (connectedList.containsKey(Integer.valueOf(atomIndex)));
	}

	/**
	 * get the bond type of the current atom with the specified atom index.
	 * 
	 * @param atomIndex
	 *            the index to be checked
	 * @return BondType - type of bonding.
	 */
	public BondType getConnectivity(int atomIndex) {
		if (isConnected(atomIndex)) {
			return connectedList.get(Integer.valueOf(atomIndex));
		} else {
			return BondType.NO_BOND;
		}
	}

	/**
	 * @return the number of atoms that are strongly bonded to this atom
	 */
	public int getNumberOfStrongBonds() {
		int bonds = 0;
		Iterator<BondType> connectedListIter = connectedList.values()
				.iterator();
		while (connectedListIter.hasNext()) {
			BondType bond = connectedListIter.next();
			if (bond.isStrongBond()) {
				bonds++;
			}
		}

		return bonds;
	}

	/**
	 * @return the sum of the bond orders of this atom
	 */
	public double getSumOfBondOrders() {
		double sumBondOrders = 0;
		Iterator<BondType> connectedListIter = connectedList.values()
				.iterator();
		while (connectedListIter.hasNext()) {
			BondType bond = connectedListIter.next();
			sumBondOrders += bond.getBondOrder();
		}

		return sumBondOrders;
	}

	/**
	 * @return the number of atoms that are "double bonded" to this atom
	 */
	public int getNumberOfDoubleBonds() {
		int doubleBonds = 0;
		Iterator<BondType> connectedListIter = connectedList.values()
				.iterator();
		while (connectedListIter.hasNext()) {
			BondType bond = connectedListIter.next();
			if (bond == BondType.DOUBLE_BOND) {
				doubleBonds++;
			}
		}

		return doubleBonds;
	}

	/**
	 * @return true if this atom participates in a resonant bond structure
	 */
	public boolean hasResonantBonds() {
		Iterator<Integer> connectedAtomsIter = connectedList.keySet()
				.iterator();
		while (connectedAtomsIter.hasNext()) {
			Integer connectedIndex = connectedAtomsIter.next();
			BondType bond = connectedList.get(connectedIndex);
			if (bond == BondType.AROMATIC_BOND || bond == BondType.AMIDE_BOND) {
				return true;
			}
		}
		return false;
	}

	/**
	 * return the degree of connectivity of this atom.
	 * 
	 * @return the degree of connectivity of this atom
	 */
	public int getDegree() {
		return connectedList.size();
	}

	/**
	 * Getter for property index.
	 * 
	 * @return Value of property index.
	 * 
	 */
	public int getIndex() {
		return this.index;
	}

	/**
	 * Setter for property index.
	 * 
	 * @param index
	 *            New value of property index.
	 * 
	 */
	public void setIndex(int index) {
		this.index = index;
	}
	/** Position of length reference in the ZMatrixElement */
	private static final int LENGTH_REFERENCE_POS = 0;

	/**
	 * Getter for property lengthReference.
	 * 
	 * @return Value of property lengthReference.
	 */
	public ZMatrixItem getLengthReference() {
		initZMatrixItems();
		return zMatrixElement.get(LENGTH_REFERENCE_POS);
	}

	/**
	 * Setter for property lengthReference.
	 * 
	 * @param lengthReference
	 *            New value of property lengthReference.
	 */
	public void setLengthReference(ZMatrixItem lengthReference) {
		initZMatrixItems();
		zMatrixElement.set(LENGTH_REFERENCE_POS, lengthReference);
	}

	/** Position of angle reference in the ZMatrixElement */
	private static final int ANGLE_REFERENCE_POS = 1;

	/**
	 * Getter for property angleReference.
	 * 
	 * @return Value of property angleReference.
	 */
	public ZMatrixItem getAngleReference() {
		initZMatrixItems();
		return zMatrixElement.get(ANGLE_REFERENCE_POS);
	}

	/**
	 * Setter for property angleReference.
	 * 
	 * @param angleReference
	 *            New value of property angleReference.
	 */
	public void setAngleReference(ZMatrixItem angleReference) {
		initZMatrixItems();
		zMatrixElement.set(ANGLE_REFERENCE_POS, angleReference);
	}

	/** Position of dihedral reference in the ZMatrixElement */
	private static final int DIHEDRAL_REFERENCE_POS = 2;

	/**
	 * Getter for property dihedralReference.
	 * 
	 * @return Value of property dihedralReference.
	 */
	public ZMatrixItem getDihedralReference() {
		initZMatrixItems();
		return zMatrixElement.get(DIHEDRAL_REFERENCE_POS);
	}

	/**
	 * Setter for property dihedralReference.
	 * 
	 * @param dihedralReference
	 *            New value of property dihedralReference.
	 */
	public void setDihedralReference(ZMatrixItem dihedralReference) {
		initZMatrixItems();
		zMatrixElement.set(DIHEDRAL_REFERENCE_POS, dihedralReference);
	}

	protected ArrayList<UserDefinedAtomProperty> userProperties;

	/**
	 * Add a user defined atom property.
	 * 
	 * @param uProp
	 *            the new instance of property to be added
	 * @throws UnsupportedOperationException
	 *             if a property with the same name already exists
	 */
	public void addUserDefinedAtomProperty(UserDefinedAtomProperty uProp) {
		if (userProperties == null)
			userProperties = new ArrayList<>();

		// first check if this property already exists
		String name = uProp.getName();
		for (UserDefinedAtomProperty uprop : userProperties) {
			if (uprop.getName().equals(name)) {
				throw new UnsupportedOperationException("Property with name"
						+ " '" + name + "' is already defined for this object!");
			}
		}

		userProperties.add(uProp);
	}

	/**
	 * Remove a user defined property.
	 * 
	 * @param uProp
	 *            instance of UserDefinedAtomProperty to be removed
	 */
	public void removeUserDefinedAtomProperty(UserDefinedAtomProperty uProp) {
		if (userProperties == null)
			return;

		userProperties.remove(uProp);
	}

	/**
	 * Get a list of all the user defined atom properties, null if no user
	 * defined properties were ever added to this Atom object.
	 * 
	 * @return an Iterator object of list of instances of
	 *         UserDefinedAtomProperty
	 */
	public Iterator<UserDefinedAtomProperty> getUserDefinedAtomProperty() {
		if (userProperties == null)
			return null;

		return userProperties.iterator();
	}

	/**
	 * Get a user defined atom property with the specified name, null if no user
	 * defined property with specified name exists.
	 * 
	 * @param name
	 *            the name of the property want a reference to
	 * @return an Iterator object of list of instances of
	 *         UserDefinedAtomProperty
	 */
	public UserDefinedAtomProperty getUserDefinedAtomProperty(String name) {
		if (userProperties == null)
			return null;

		for (UserDefinedAtomProperty uprop : userProperties)
			if (uprop.getName().equals(name))
				return uprop;

		return null;
	}

	@Override
	public int hashCode() {
		return Objects.hash(atomCenter, index, symbol);
	}

	@Override
	public boolean equals(Object obj) {
		if (this == obj)
			return true;
		if (obj == null)
			return false;
		if (getClass() != obj.getClass())
			return false;
		Atom other = (Atom) obj;
		return Objects.equals(atomCenter, other.atomCenter)
				&& index == other.index
				&& Objects.equals(symbol, other.symbol);
	}
}
