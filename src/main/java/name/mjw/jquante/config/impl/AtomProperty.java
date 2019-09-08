package name.mjw.jquante.config.impl;

import name.mjw.jquante.config.Parameter;

/**
 * Represents an collection "wrapper" of all the properties of a single atom.
 * 
 * @author V.Ganesh
 * @version 2.0 (Part of MeTA v2.0)
 */
public final class AtomProperty implements Parameter {

	/** Holds value of property name. */
	private String name;

	/** Holds value of property atomicNumber. */
	private int atomicNumber;

	/** Holds value of property atomicWeight. */
	private double atomicWeight;

	/**
	 * Creates a new instance of AtomProperty
	 * 
	 * @param name           Name of atom.
	 * @param atomicNumber   Atomic nubmer of atom.
	 * @param atomicWeight   Atomic weight of atom.
	 */
	public AtomProperty(String name, int atomicNumber, double atomicWeight) {

		this.name = name;
		this.atomicNumber = atomicNumber;
		this.atomicWeight = atomicWeight;
	}

	/**
	 * Getter for property name.
	 * 
	 * @return Value of property name.
	 * 
	 */
	public String getName() {
		return this.name;
	}

	/**
	 * Setter for property name.
	 * 
	 * @param name New value of property name.
	 * 
	 */
	public void setName(String name) {
		this.name = name;
	}

	/**
	 * Getter for property atomicNumber.
	 * 
	 * @return Value of property atomicNumber.
	 * 
	 */
	public int getAtomicNumber() {
		return this.atomicNumber;
	}

	/**
	 * Setter for property atomicNumber.
	 * 
	 * @param atomicNumber New value of property atomicNumber.
	 * 
	 */
	public void setAtomicNumber(int atomicNumber) {
		this.atomicNumber = atomicNumber;
	}

	/**
	 * Getter for property atomicWeight.
	 * 
	 * @return Value of property atomicWeight.
	 * 
	 */
	public double getAtomicWeight() {
		return this.atomicWeight;
	}

	/**
	 * Setter for property atomicWeight.
	 * 
	 * @param atomicWeight New value of property atomicWeight.
	 * 
	 */
	public void setAtomicWeight(double atomicWeight) {
		this.atomicWeight = atomicWeight;
	}

	@Override
	public Object getValue() {
		return this;
	}

	@Override
	public void setValue(Object value) {
		AtomProperty ap = (AtomProperty) value;

		this.name = ap.name;
		this.atomicNumber = ap.atomicNumber;
		this.atomicWeight = ap.atomicWeight;
	}
}