package name.mjw.jquante.math.qm.integral;

import org.hipparchus.geometry.euclidean.threed.Vector3D;

import net.jafama.FastMath;

/**
 * Utility class for IntegralsPackage.
 * 
 * @author V.Ganesh
 * @version 2.0 (Part of MeTA v2.0)
 */
public final class IntegralsUtil {

	private IntegralsUtil() {
	}

	/**
	 * the Gaussian product theorem
	 * 
	 * @param alpha1
	 *            exponent of Gaussian a
	 * @param a
	 *            centre of Gaussian a
	 * @param alpha2
	 *            exponent of Gaussian b
	 * @param b
	 *            centre of Gaussian b
	 * 
	 * @return centre of resulting two-centre product Gaussian
	 */
	public static final Vector3D gaussianProductCenter(final double alpha1, final Vector3D a, final double alpha2,
			final Vector3D b) {
		final double gamma = alpha1 + alpha2;
		return new Vector3D((alpha1 * a.getX() + alpha2 * b.getX()) / gamma,
				(alpha1 * a.getY() + alpha2 * b.getY()) / gamma, (alpha1 * a.getZ() + alpha2 * b.getZ()) / gamma);
	}

	public static double[] gaussianProductCenter(double alpha1, double ax, double ay, double az, double alpha2,
			double bx, double by, double bz) {

		final double gamma = alpha1 + alpha2;

		return new double[] { (alpha1 * ax + alpha2 * bx) / gamma, (alpha1 * ay + alpha2 * by) / gamma,
				(alpha1 * az + alpha2 * bz) / gamma };
	}

	/**
	 *
	 * Index into the array using ⟨ij|kl⟩.
	 *
	 * The two particle integral ⟨ij|kl⟩ is:
	 *
	 * ⟨ij|kl⟩=∫ψi(x)ψj(x′)ψk(x)ψl(x′)|x−x′|d3xd3x′
	 *
	 * With the following 8 symmetries:
	 *
	 * ⟨ij|kl⟩=⟨ji|lk⟩=⟨kj|il⟩=⟨il|kj⟩=⟨kl|ij⟩=⟨lk|ji⟩=⟨il|kj⟩=⟨kj|il⟩
	 *
	 * This makes use of the XOR swap algorithm
	 * https://en.wikipedia.org/wiki/XOR_swap_algorithm
	 *
	 * @param i the index of Gaussian i
	 * @param j the index of Gaussian j
	 * @param k the index of Gaussian k
	 * @param l the index of Gaussian l
	 * @return the index into the array
	 */
	public static final int ijkl2intindex(int i, int j, int k, int l) {
		if (i < j) {
			i ^= j;
			j ^= i;
			i ^= j;
		}

		if (k < l) {
			k ^= l;
			l ^= k;
			k ^= l;
		}
		int ij = i * (i + 1) / 2 + j;
		int kl = k * (k + 1) / 2 + l;

		if (ij < kl) {
			ij ^= kl;
			kl ^= ij;
			ij ^= kl;
		}

		return (ij * (ij + 1) / 2 + kl);
	}

	/**
	 * Incomplete gamma function
	 * 
	 * @param m
	 *            The parameter of the integral
	 * @param x
	 *            The upper bound for the interval of integration
	 * @return &Gamma;
	 */
	public static final double computeFGamma(final int m, double x) {
		final double SMALL = 0.00000001;

		x = FastMath.max(FastMath.abs(x), SMALL);

		return (0.5 * FastMath.pow(x, -m - 0.5) * lowerIncompleteGamma(m + 0.5, x));
	}

	/**
	 * Computes the lower incomplete gamma function, &gamma;(a,x)= integral from
	 * zero to x of (exp(-t)t^(a-1))dt.
	 * <p>
	 * 
	 * @param a
	 *            The parameter of the integral
	 * @param x
	 *            The upper bound for the interval of integration
	 * @return the incomplete gamma function P(a,x)
	 */
	public static final double lowerIncompleteGamma(final double a, final double x) {
		final double EPS = 3.0e-7;
		final double FPMIN = 1.0e-30;
		final int MAX_ITERATION = 100;

		final double gln = logGamma(a);

		if (x < (a + 1.0)) {
			return FastMath.exp(gln) * seriesMethod(a, x, EPS, MAX_ITERATION, gln);
		} else {
			return FastMath.exp(gln) * continuedFractionMethod(a, x, EPS, FPMIN, MAX_ITERATION, gln);
		}

	}

	/**
	 * Computes the incomplete Gamma function by using the continued fraction
	 * representation. See NumRec sect 6.1
	 *
	 * @param a
	 *            The parameter of the integral
	 * @param x
	 *            The upper bound for the interval of integration
	 * @param EPS
	 * @param FPMIN
	 * @param MAX_ITERATION
	 * @param gln
	 * @return
	 */
	private static final double continuedFractionMethod(final double a, final double x, final double EPS, final double FPMIN,
			final int MAX_ITERATION, final double gln) {

		double b = (x + 1.0) - a;
		double c = 1.0 / FPMIN;
		double d = 1.0 / b;
		double h = d;
		double an;
		double delta;

		for (int i = 1; i < (MAX_ITERATION + 1); i++) {
			an = -i * (i - a);
			b += 2.0;
			d = an * d + b;

			if (FastMath.abs(d) < FPMIN) {
				d = FPMIN;
			}

			c = b + an / c;

			if (FastMath.abs(c) < FPMIN) {
				c = FPMIN;
			}

			d = 1.0 / d;
			delta = d * c;
			h *= delta;

			if (FastMath.abs(delta - 1.0) < EPS)
				break;
		}

		return 1.0 - (FastMath.exp(-x + a * FastMath.log(x) - gln) * h);
	}

	/**
	 * Computes the incomplete Gamma function by using the series representation.
	 * See NumRec sect 6.1.
	 *
	 * @param a
	 *            The parameter of the integral
	 * @param x
	 *            The upper bound for the interval of integration
	 * @param EPS
	 * @param MAX_ITERATION
	 * @param gln
	 * @return
	 */
	private static final double seriesMethod(final double a, final double x, final double EPS, final int MAX_ITERATION,
			final double gln) {
		if (x != 0.0) {
			double ap = a;
			double sum;
			double delta = sum = 1.0 / a;

			for (int i = 0; i < MAX_ITERATION; i++) {
				ap++;
				delta *= x / ap;
				sum += delta;
				if (FastMath.abs(delta) < FastMath.abs(sum) * EPS)
					break;
			}

			return sum * FastMath.exp(-x + a * FastMath.log(x) - gln);
		} else {
			return 0.0;
		}
	}

	/**
	 * Computes the value of ln(gamma(x)) Numerical recipes, section 6.1
	 * 
	 * @param x
	 *            Value
	 * @return the complete Gamma(x) function
	 */
	public static final double logGamma(final double x) {
		double y = x;
		double tmp = x + 5.5;

		tmp -= (x + 0.5) * FastMath.log(tmp);

		double ser = 1.000000000190015;

		y++;
		ser += 76.18009172947146 / y;
		y++;
		ser += -86.50532032941677 / y;
		y++;
		ser += 24.01409824083091 / y;
		y++;
		ser += -1.231739572450155 / y;
		y++;
		ser += 0.1208650973866179e-2 / y;
		y++;
		ser += -0.5395239384953e-5 / y;

		return (-tmp + FastMath.log((2.5066282746310005 * ser) / x));
	}
}
