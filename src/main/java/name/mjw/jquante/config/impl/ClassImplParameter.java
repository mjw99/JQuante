package name.mjw.jquante.config.impl;

import name.mjw.jquante.config.Parameter;

/**
 * This class defines the a simple implementation of Parameter interface for
 * storing the implementation class of an interface for MeTA Studio. The class
 * may be either a concrete Java class or a class generated from an embedded
 * script engine.
 * 
 * @author V.Ganesh
 * @version 2.0 (Part of MeTA v2.0)
 */
public final class ClassImplParameter implements Parameter {

	private String implClass;

	/**
	 * Creates a new instance of ClassImplParameter
	 * 
	 * @param implClass
	 *            Name of ClassImplParameter
	 */
	public ClassImplParameter(String implClass) {
		this.implClass = implClass;
	}

	/**
	 * Getter for property value.
	 * 
	 * @return Value of property value.
	 */
	@Override
	public Object getValue() {
		return implClass;
	}

	/**
	 * Setter for property value.
	 * 
	 * @param value
	 *            New value of property value.
	 * 
	 */
	@Override
	public void setValue(Object value) {
		this.implClass = (String) value;
	}

}