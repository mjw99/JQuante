package name.mjw.jquante.math.qm.basis;

/**
 * A runtime exception thrown is the requested basis is not found.
 * 
 * @author V.Ganesh
 * @version 2.0 (Part of MeTA v2.0)
 */
public final class BasisNotFoundException extends java.lang.RuntimeException {

	/**
	 * 
	 */
	private static final long serialVersionUID = 2348577368651609288L;

	/**
	 * Creates a new instance of <code>BasisNotFoundException</code> without
	 * detail message.
	 */
	public BasisNotFoundException() {
	}

	/**
	 * Constructs an instance of <code>BasisNotFoundException</code> with the
	 * specified detail message.
	 * 
	 * @param msg
	 *            the detail message.
	 */
	public BasisNotFoundException(String msg) {
		super(msg);
	}
}