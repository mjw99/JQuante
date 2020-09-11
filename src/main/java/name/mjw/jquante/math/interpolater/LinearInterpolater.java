package name.mjw.jquante.math.interpolater;

/**
 * A simple liner interpolater based on the following formula: <br>
 * <code>y1 * (1 - mu) + y2 * mu</code>
 * 
 * @author V.Ganesh
 * @version 2.0 (Part of MeTA v2.0)
 */
public class LinearInterpolater extends Interpolater {

	/** Creates a new instance of LinearInterpolater */
	public LinearInterpolater() {
	}

	/**
	 * Interpolate value at X depending upon value at Y Uses the formula: <br>
	 * <code>y1 * (1 - mu) + y2 * mu <br>  Where, <br>
	 * y1, y2 &ge; y[0] and y[1] respectively <br>
	 * mu &ge; x[0]
	 * </code>
	 * 
	 * @param y
	 *            the Y values (results of function evaluation)
	 * @param x
	 *            the X values at function evaluation is performed or is expected
	 * @return the interpolated value depending upon the interpolation formula
	 */
	@Override
	public double interpolate(final double[] y, final double[] x) {
		return (y[0] * (1.0 - x[0]) + y[1] * x[0]);
	}

	/**
	 * Return number of Y arguments required for this Interpolater
	 * 
	 * @return number of Y arguments
	 */
	@Override
	public int getNumberOfRequiredYArgs() {
		return 2;
	}

	/**
	 * Return number of X arguments required for this Interpolater
	 * 
	 * @return number of X arguments
	 */
	@Override
	public int getNumberOfRequiredXArgs() {
		return 1;
	}

	/**
	 * Overloaded toString()
	 * 
	 * @return string representation
	 */
	@Override
	public String toString() {
		return "Linear interpolater";
	}
}