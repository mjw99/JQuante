package name.mjw.jquante.config.event;

import java.util.EventObject;

/**
 * The EventObject indicating the AtomInfo change event to all the "listeners".
 * 
 * @author V.Ganesh
 * @version 2.0 (Part of MeTA v2.0)
 */
public final class AtomInfoChangeEvent extends EventObject {

	private static final long serialVersionUID = -5378842897630157145L;

	/** Holds value of property changeType. */
	private int changeType;

	/** Holds value of property newValue. */
	private Object newValue;

	/** Holds value of property oldValue. */
	private Object oldValue;

	/** the change types */
	public static final int ATOM_NAME = 0;
	public static final int ATOMIC_NUMBER = 1;
	public static final int ATOMIC_WEIGHT = 2;
	public static final int ALL_CANGED = 100;

	/** Holds value of property atomSymbol. */
	private String atomSymbol;

	/**
	 * Creates a new instance of AtomInfoChangeEvent *
	 * 
	 * @param source
	 *            The object on which the Event initially occurred.
	 */
	public AtomInfoChangeEvent(Object source) {
		super(source);
	}

	/**
	 * Getter for property changeType.
	 * 
	 * @return Value of property changeType.
	 * 
	 */
	public int getChangeType() {
		return this.changeType;
	}

	/**
	 * Setter for property changeType.
	 * 
	 * @param changeType
	 *            New value of property changeType.
	 * 
	 */
	public void setChangeType(int changeType) {
		this.changeType = changeType;
	}

	/**
	 * Getter for property newValue.
	 * 
	 * @return Value of property newValue.
	 * 
	 */
	public Object getNewValue() {
		return this.newValue;
	}

	/**
	 * Setter for property newValue.
	 * 
	 * @param newValue
	 *            New value of property newValue.
	 * 
	 */
	public void setNewValue(Object newValue) {
		this.newValue = newValue;
	}

	/**
	 * Getter for property oldValue.
	 * 
	 * @return Value of property oldValue.
	 * 
	 */
	public Object getOldValue() {
		return this.oldValue;
	}

	/**
	 * Setter for property oldValue.
	 * 
	 * @param oldValue
	 *            New value of property oldValue.
	 * 
	 */
	public void setOldValue(Object oldValue) {
		this.oldValue = oldValue;
	}

	/**
	 * Getter for property atomSymbol.
	 * 
	 * @return Value of property atomSymbol.
	 * 
	 */
	public String getAtomSymbol() {
		return this.atomSymbol;
	}

	/**
	 * Setter for property atomSymbol.
	 * 
	 * @param atomSymbol
	 *            New value of property atomSymbol.
	 * 
	 */
	public void setAtomSymbol(String atomSymbol) {
		this.atomSymbol = atomSymbol;
	}

}