package name.mjw.jquante.math;

import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.ArrayRealVector;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.RealMatrixFormat;
import org.apache.commons.math3.linear.RealVector;
import org.apache.commons.math3.util.CombinatoricsUtils;

import net.jafama.FastMath;

import java.text.DecimalFormat;

import org.apache.commons.math3.geometry.euclidean.threed.Vector3D;

/**
 * A collection of few misc. utility math functions. All methods are static and
 * the class cannot be instantiated.
 * 
 * @author V.Ganesh
 * @version 2.0 (Part of MeTA v2.0)
 */
public final class MathUtil {

	/** Creates a new instance of MathUtil */
	private MathUtil() {
	}

	/**
	 * Method to find the angle in radians defined by three points v1-v2-v3.
	 * 
	 * @param v1
	 *            the first point
	 * @param v2
	 *            the second point (central angle)
	 * @param v3
	 *            the third point
	 * @return the angle defined
	 */
	public static final double findAngle(Vector3D v1, Vector3D v2, Vector3D v3) {
		Vector3D v12 = v2.subtract(v1);
		Vector3D v32 = v2.subtract(v3);

		return Vector3D.angle(v12, v32);
	}

	/**
	 * Method to find the dihedral angle defined by planes v1-v2-v3 and v2-v3-v4.
	 * 
	 * @param v1
	 *            first point
	 * @param v2
	 *            second point
	 * @param v3
	 *            third point (2nd and 3rd point define the angle)
	 * @param v4
	 *            the fourth angle
	 * @return the dihedral angle defined
	 */
	public static final double findDihedral(Vector3D v1, Vector3D v2, Vector3D v3, Vector3D v4) {
		// normal of plane 1
		Vector3D v12 = v2.subtract(v1);
		Vector3D v32 = v2.subtract(v3);
		Vector3D n123 = v12.crossProduct(v32).normalize();

		// normal of plane 2
		Vector3D v23 = v3.subtract(v2);
		Vector3D v43 = v3.subtract(v4);
		Vector3D n234 = v23.crossProduct(v43).normalize();

		// sign of the dihedral
		double sign = v32.dotProduct(n123.crossProduct(n234));

		if (sign >= 0.0) {
			sign = -1.0;
		} else {
			sign = 1.0;
		}

		// and find the angle between the two planes
		return Vector3D.angle(n123, n234) * sign;
	}

	/**
	 * compute double N! ... (1*3*5*...*n)
	 * 
	 * @param n
	 *            the n, whose factorial is to be found
	 * @return the factorial
	 */
	public static final int factorial2(int n) {
		int value = 1;

		while (n > 0) {
			value = value * n;
			n -= 2;
		}
		return value;
	}

	/**
	 * Does ( a! / b! / (a-2*b)! )
	 * 
	 * @param a
	 *            the first term
	 * @param b
	 *            the second term
	 * @return ( a! / b! / (a-2*b)! )
	 */
	public static final double factorialRatioSquared(int a, int b) {
		return CombinatoricsUtils.factorialDouble(a) / CombinatoricsUtils.factorialDouble(b)
				/ CombinatoricsUtils.factorialDouble(a - 2 * b);
	}

	/**
	 * Pre-factor of binomial expansion.
	 * 
	 * From Augspurger and Dykstra:
	 * <a href="http://dx.doi.org/10.1021/j100176a037">10.1021/j100176a037</a>
	 * 
	 * @param s
	 *            s
	 * @param ia
	 *            ia
	 * @param ib
	 *            ib
	 * @param xpa
	 *            xpa
	 * @param xpb
	 *            xpb
	 * @return Pre-factor of binomial expansion.
	 */
	public static final double binomialPrefactor(int s, int ia, int ib, double xpa, double xpb) {
		double sum = 0.0;

		for (int t = 0; t < (s + 1); t++) {
			if (((s - ia) <= t) && (t <= ib)) {
				sum += CombinatoricsUtils.binomialCoefficientDouble(ia, s - t)
						* CombinatoricsUtils.binomialCoefficientDouble(ib, t) * FastMath.pow(xpa, ia - s + t)
						* FastMath.pow(xpb, ib - t);
			}
		}

		return sum;
	}

	public static final RealVector realMatrixToRealVector(RealMatrix realMatrix) {
		RealVector vec = new ArrayRealVector(realMatrix.getRowDimension() * realMatrix.getColumnDimension());
		int ii = 0;

		for (int i = 0; i < realMatrix.getRowDimension(); i++) {
			for (int j = 0; j < realMatrix.getColumnDimension(); j++) {
				vec.setEntry(ii++, realMatrix.getEntry(i, j));
			}
		}
		return vec;

	}

	public static final String matrixToString(Array2DRowRealMatrix array2DRowRealMatrix) {
		DecimalFormat df = new DecimalFormat("+#,##0.000;-#");
		RealMatrixFormat mf = new RealMatrixFormat("\n", "", "", "", "\n", " ", df);
		return mf.format(array2DRowRealMatrix);
	}
}