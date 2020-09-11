package name.mjw.jquante.molecule.property.electronic;

/**
 * Exception thrown when a Property is not defined (at a requested loaction, or
 * some criterion is not met).
 * 
 * @author V.Ganesh
 * @version 2.0 (Part of MeTA v2.0)
 */
public class PropertyNotDefinedException extends Exception {

	/**
	 * Eclipse generated serialVersionUID
	 */
	private static final long serialVersionUID = 302077363758923072L;

	/** Creates a new instance of PropertyNotDefinedException */
	public PropertyNotDefinedException() {
	}

	private String errMsg = "Undefined Property";

	/**
	 * Creates a new instance of PropertyNotDefinedException
	 * 
	 * @param errMsg
	 *            Error message.
	 */
	public PropertyNotDefinedException(final String errMsg) {
		this.errMsg = errMsg;
	}

	/**
	 * overloaded toString() method
	 */
	@Override
	public String toString() {
		return errMsg;
	}
}
