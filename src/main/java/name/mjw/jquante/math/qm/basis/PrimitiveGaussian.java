package name.mjw.jquante.math.qm.basis;

import name.mjw.jquante.math.MathUtil;

import java.util.Objects;

import org.hipparchus.geometry.euclidean.threed.Vector3D;
import name.mjw.jquante.math.qm.integral.Integrals;
import net.jafama.FastMath;

/**
 * 
 * The class defines a primitive Gaussian (PG) and the operations on it.
 * 
 * Gaussian Type Orbitial based upon the following function:
 * <p>
 * g(x,y,z) = normalization*(x^i)*(y^j)*(z^k)*exp{-exponent*(r-ro)^2}
 * <p>
 * 
 * Refs: http://dx.doi.org/10.1143/JPSJ.21.2313
 * 
 * 
 * @author V.Ganesh
 * @author mw529
 * @version 2.0 (Part of MeTA v2.0)
 */
public final class PrimitiveGaussian implements Comparable<PrimitiveGaussian>{

	/**
	 * Holds value of property exponent.
	 */
	private final double exponent;

	/**
	 * Holds value of property origin.
	 */
	private final Vector3D origin;

	/**
	 * Holds value of property powers.
	 */
	private final Power powers;

	/**
	 * Holds value of property coefficient.
	 */
	private final double coefficient;

	/**
	 * normalization factor
	 */
	private final double normalization;

	private static final double PI_RAISE_TO_1DOT5 = FastMath.pow(Math.PI, 1.5);

	/**
	 * Creates a new instance of PrimitiveGaussian
	 * 
	 * @param origin      - the (x, y, z) on which this Gaussian is centered
	 * @param powers      - the powers of this Gaussian
	 * @param exponent    - the exponent for this PG
	 * @param coefficient - the coefficient for this PG
	 */
	public PrimitiveGaussian(Vector3D origin, Power powers, double exponent, double coefficient) {
		this.origin = origin;
		this.powers = powers;
		this.exponent = exponent;
		this.coefficient = coefficient;

		int l = powers.getL();
		int m = powers.getM();
		int n = powers.getN();

		/**
		 * Normalize this primitive Gaussian.
		 * 
		 * Equ 2.2 http://dx.doi.org/10.1143/JPSJ.21.2313 -
		 */
		normalization = FastMath.sqrt(Math.pow(2, 2 * (l + m + n) + 1.5) * FastMath.pow(exponent, l + m + n + 1.5)
				/ MathUtil.factorial2(2 * l - 1) / MathUtil.factorial2(2 * m - 1) / MathUtil.factorial2(2 * n - 1)
				/ PI_RAISE_TO_1DOT5);
	}

	/**
	 * Getter for property exponent.
	 * 
	 * @return Value of property exponent.
	 */
	public double getExponent() {
		return this.exponent;
	}

	/**
	 * Getter for property origin.
	 * 
	 * @return Value of property origin.
	 */
	public Vector3D getOrigin() {
		return this.origin;
	}

	/**
	 * Getter for property powers.
	 * 
	 * @return Value of property powers.
	 */
	public Power getPowers() {
		return this.powers;
	}

	/**
	 * Getter for property coefficient.
	 * 
	 * @return Value of property coefficient.
	 */
	public double getCoefficient() {
		return this.coefficient;
	}

	/**
	 * Gaussian product center. Return a new PG with product of this PG with a new
	 * PG. The powers of the resulting PG is set to 0,0,0 and the coefficient to 0.0
	 * 
	 * @param pg the PG with which to multiply
	 * @return the product of the PG
	 */
	public PrimitiveGaussian mul(PrimitiveGaussian pg) {
		double gamma = exponent + pg.exponent;
		Vector3D newOrigin = new Vector3D((exponent * origin.getX() + pg.exponent * pg.origin.getX()) / gamma,
				(exponent * origin.getY() + pg.exponent * pg.origin.getY()) / gamma,
				(exponent * origin.getZ() + pg.exponent * pg.origin.getZ()) / gamma);

		return new PrimitiveGaussian(newOrigin, new Power(0, 0, 0), gamma, 0.0);
	}

	/**
	 * Overlap matrix element with another PrimitiveGaussian
	 * 
	 * @param pg the PrimitiveGaussian with which the overlap is to be be
	 *           determined.
	 * @return the overlap value
	 */
	public double overlap(PrimitiveGaussian pg) {
		return (normalization * pg.normalization
				* Integrals.overlap(exponent, powers, origin, pg.exponent, pg.powers, pg.origin));
	}

	/**
	 * Kinetic Energy (KE) matrix element with another PrimitiveGaussian
	 * 
	 * @param pg the PrimitiveGaussian with which KE is to be determined.
	 * @return the KE value
	 */
	public double kinetic(PrimitiveGaussian pg) {
		return (normalization * pg.normalization
				* Integrals.kinetic(exponent, powers, origin, pg.exponent, pg.powers, pg.origin));
	}

	/**
	 * Nuclear matrix element with another PrimitiveGaussian
	 * 
	 * @param pg     the PrimitiveGaussian with which nuclear interaction is to be
	 *               determined.
	 * @param center the center at which nuclear energy is to be computed
	 * @return the nuclear value
	 */
	public double nuclear(PrimitiveGaussian pg, Vector3D center) {
		return (Integrals.nuclearAttraction(origin, normalization, powers, exponent, pg.origin, pg.normalization,
				pg.powers, pg.exponent, center));
	}

	/**
	 * Return the nuclear gradient of this PG w.r.t a center
	 * 
	 * @param pg     the other PG
	 * @param center the reference center
	 * @return partial derivatives w.r.t the center
	 */
	public Vector3D nuclearAttractionGradient(PrimitiveGaussian pg, Vector3D center) {
		Vector3D nder = new Vector3D(0, 0, 0);

		// TODO:

		return nder;
	}

	/**
	 * The amplitude of this primitive Gaussian at a given point.
	 * 
	 * @param point the reference point
	 * @return the amplitude of this PG at the specified point
	 */
	public double amplitude(Vector3D point) {
		final int l = powers.getL();
		final int m = powers.getM();
		final int n = powers.getN();

		final double x = point.getX() - origin.getX();
		final double y = point.getY() - origin.getY();
		final double z = point.getZ() - origin.getZ();

		final double d2 = FastMath.pow(x, l) * FastMath.pow(y, m) * FastMath.pow(z, n);

		return (normalization * coefficient * d2 * FastMath.exp(-exponent * (x * x + y * y + z * z)));
	}

	/**
	 * Calculate Laplacian at the specified point
	 * 
	 * @param point the point where Laplacian is to be computed
	 * @return the Laplacian at this point
	 */
	public double laplacian(Vector3D point) {
		double value = 0.0;
		double x = point.getX() - origin.getX(), y = point.getY() - origin.getY(), z = point.getZ() - origin.getZ();
		double x2 = x * x;
		double y2 = y * y;
		double z2 = z * z;
		int l = powers.getL();
		int m = powers.getM();
		int n = powers.getN();

		value = (l * (l - 1) / x2 + m * (m - 1) / y2 + n * (n - 1) / z2) + 4 * exponent * exponent * (x2 + y2 + z2)
				- 2 * exponent * (2 * (l + m + n) + 3);

		return value * normalization * coefficient * amplitude(point);
	}

	/**
	 * Evaluate gradient of this function at a given point
	 * 
	 * @param point the point where gradient is to be evaluated
	 * @return partial derivatives with respect to x, y, z
	 */
	public Vector3D gradient(Vector3D point) {
		int l = powers.getL();
		int m = powers.getM();
		int n = powers.getN();

		double x = point.getX() - origin.getX();
		double y = point.getY() - origin.getY();
		double z = point.getZ() - origin.getZ();

		double fx = FastMath.pow(x, l) * FastMath.exp(-exponent * FastMath.pow(x, 2));
		double fy = FastMath.pow(y, m) * FastMath.exp(-exponent * FastMath.pow(y, 2));
		double fz = FastMath.pow(z, n) * FastMath.exp(-exponent * FastMath.pow(z, 2));
		double gx = -2.0 * exponent * x * fx;
		double gy = -2.0 * exponent * y * fy;
		double gz = -2.0 * exponent * z * fz;

		if (l > 0)
			gx += FastMath.pow(x, l - 1) * FastMath.exp(-exponent * FastMath.pow(x, 2));
		if (m > 0)
			gy += FastMath.pow(y, m - 1) * FastMath.exp(-exponent * FastMath.pow(y, 2));
		if (n > 0)
			gz += FastMath.pow(z, n - 1) * FastMath.exp(-exponent * FastMath.pow(z, 2));

		double nc = normalization * coefficient;
		return new Vector3D(gx * fy * fz * nc, fx * gy * fz * nc, fx * fy * gz * nc);

	}

	/**
	 * Getter for property normalization.
	 * 
	 * @return Value of property normalization.
	 */
	public double getNormalization() {
		return this.normalization;
	}

	/**
	 * Return the maximum angular momentum for this primitive basis function
	 * 
	 * @return the maximum angular momentum of this primitive basis function
	 */
	public int getMaximumAngularMomentum() {
		return powers.getMaximumAngularMomentum();
	}

	/**
	 * Return the minimum angular momentum for this primitive basis function
	 * 
	 * @return the minimum angular momentum of this primitive basis function
	 */
	public int getMinimumAngularMomentum() {
		return powers.getMinimumAngularMomentum();
	}

	/**
	 * Return the total angular momentum of this primitive basis function
	 * 
	 * @return the maximum of the primitive basis function
	 */
	public int getTotalAngularMomentum() {
		return powers.getTotalAngularMomentum();
	}

	/**
	 * overloaded toString()
	 */
	@Override
	public String toString() {
		return "Origin : " + origin + " Powers : " + powers + " Exponent : " + exponent + " Coefficient : "
				+ coefficient + " Normalization : " + normalization + "\n";
	}

	@Override
	public int hashCode() {
		return Objects.hash(coefficient, exponent, normalization, origin, powers);
	}

	@Override
	public boolean equals(Object obj) {
		if (this == obj)
			return true;
		if (obj == null)
			return false;
		if (getClass() != obj.getClass())
			return false;
		PrimitiveGaussian other = (PrimitiveGaussian) obj;
		return Double.doubleToLongBits(coefficient) == Double.doubleToLongBits(other.coefficient)
				&& Double.doubleToLongBits(exponent) == Double.doubleToLongBits(other.exponent)
				&& Double.doubleToLongBits(normalization) == Double.doubleToLongBits(other.normalization)
				&& Objects.equals(origin, other.origin) && Objects.equals(powers, other.powers);
	}

	/**
	 * Sort exponents in decending order.
	 */
	@Override
	public int compareTo(PrimitiveGaussian other) {
		if (this.exponent > other.exponent) {
			return -1;
		} else {
			return 1;
		}
	}
}
