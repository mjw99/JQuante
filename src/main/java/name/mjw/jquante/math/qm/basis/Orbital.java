package name.mjw.jquante.math.qm.basis;

import java.util.ArrayList;
import java.util.List;

/**
 * Represents an orbital type and its coefficients and exponents in
 * <code>AtomicBasis</code>.
 * 
 * @author V.Ganesh
 * @version 2.0 (Part of MeTA v2.0)
 */
public final class Orbital {

	/**
	 * Holds value of property type.
	 */
	private final String type;

	/**
	 * Holds value of property coefficients.
	 */
	private final List<Double> coefficients;

	/**
	 * Holds value of property exponents.
	 */
	private final List<Double> exponents;

	/**
	 * Creates a new instance of Orbital
	 * 
	 * @param type the type of this orbital (e.g. 'S', 'P', 'D' etc...)
	 */
	public Orbital(final String type) {
		this.type = type;

		coefficients = new ArrayList<>(10);
		exponents = new ArrayList<>(10);
	}

	/**
	 * Getter for property type.
	 * 
	 * @return Value of property type.
	 */
	public String getType() {
		return this.type;
	}

	/**
	 * Getter for property coefficients.
	 * 
	 * @return Value of property coefficients.
	 */
	public List<Double> getCoefficients() {
		return this.coefficients;
	}

	/**
	 * Getter for property exponents.
	 * 
	 * @return Value of property exponents.
	 */
	public List<Double> getExponents() {
		return this.exponents;
	}

	/**
	 * adds a coefficient to this Orbital entry
	 * 
	 * @param coeff the coefficient to be added
	 */
	public void addCoefficient(final double coeff) {
		coefficients.add(Double.valueOf(coeff));
	}

	/**
	 * adds a exponent to this Orbital entry
	 * 
	 * @param exp the exponent to be added
	 */
	public void addExponent(final double exp) {
		exponents.add(Double.valueOf(exp));
	}

	/**
	 * adds a (coefficient, exponent) pair to this Orbital entry
	 * 
	 * @param coeff the coefficient to be added
	 * @param exp   the exponent to be added
	 */
	public void addEntry(final double coeff, final double exp) {
		addCoefficient(coeff);
		addExponent(exp);
	}
}
