package name.mjw.jquante.molecule;

import java.io.Serializable;

/**
 * A simple framework for adding arbitrary user defined property to Atom
 * instance.
 * 
 * @author V.Ganesh
 * @version 2.0 (Part of MeTA v2.0)
 */
public class UserDefinedAtomProperty implements Serializable {

	/**
	 * Eclipse generated serialVersionUID
	 */
	private static final long serialVersionUID = -172897270553276800L;

	/** Creates new instance of UserDefinedAtomProperty */
	public UserDefinedAtomProperty() {
	}

	/**
	 * Creates new instance of UserDefinedAtomProperty
	 * 
	 * @param name
	 *            Property name.
	 * @param value
	 *            Property name.
	 */
	public UserDefinedAtomProperty(final String name, final Serializable value) {
		this.name = name;
		this.value = value;
	}

	protected String name;

	/**
	 * Get the value of name
	 * 
	 * @return the value of name
	 */
	public String getName() {
		return name;
	}

	/**
	 * Set the value of name
	 * 
	 * @param name new value of name
	 */
	public void setName(final String name) {
		this.name = name;
	}

	protected Serializable value;

	/**
	 * Get the value of value
	 * 
	 * @return the value of value
	 */
	public Serializable getValue() {
		return value;
	}

	/**
	 * Set the value of value
	 * 
	 * @param value new value of value
	 */
	public void setValue(final Serializable value) {
		this.value = value;
	}
}
