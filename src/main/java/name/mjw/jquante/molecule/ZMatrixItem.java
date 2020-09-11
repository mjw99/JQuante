package name.mjw.jquante.molecule;

/**
 * Represents a single ZMatrix element: {@code <reference atom, value>} pair.
 * The pair may indicate distance, angle or dihedral with the reference atom.
 * 
 * @author V.Ganesh
 * @version 2.0 (Part of MeTA v2.0)
 */
public class ZMatrixItem implements Cloneable {

	/**
	 * Holds value of property referenceAtom.
	 */
	private Atom referenceAtom;
	/**
	 * Holds value of property value.
	 */
	private double value;

	/**
	 * Creates a new instance of ZMatrixItem
	 * 
	 * @param referenceAtom
	 *            Atom reference
	 * @param value
	 *            Value
	 */
	public ZMatrixItem(final Atom referenceAtom, final double value) {
		this.referenceAtom = referenceAtom;
		this.value = value;
	}

	/**
	 * to String!
	 * 
	 * @return the string representation of this item
	 */
	@Override
	public String toString() {
		return (referenceAtom.getSymbol() + referenceAtom.getIndex() + " " + value);
	}

	/**
	 * Getter for property referenceAtom.
	 * 
	 * @return Value of property referenceAtom.
	 */
	public Atom getReferenceAtom() {
		return this.referenceAtom;
	}

	/**
	 * Setter for property referenceAtom.
	 * 
	 * @param referenceAtom New value of property referenceAtom.
	 */
	public void setReferenceAtom(final Atom referenceAtom) {
		this.referenceAtom = referenceAtom;
	}

	/**
	 * Getter for property value.
	 * 
	 * @return Value of property value.
	 */
	public double getValue() {
		return this.value;
	}

	/**
	 * Setter for property value.
	 * 
	 * @param value New value of property value.
	 */
	public void setValue(final double value) {
		this.value = value;
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
		return new ZMatrixItem(this.referenceAtom, this.value);
	}
}