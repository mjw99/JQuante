package name.mjw.jquante.math.qm.integral;

import name.mjw.jquante.math.geom.Point3D;

/**
 * Utility class for IntegralsPackage.
 * 
 * @author V.Ganesh
 * @version 2.0 (Part of MeTA v2.0)
 */
public final class IntegralsUtil {

	// for gammaIncomplete method
	private static final double SMALL = 0.00000001;
	private static final double EPS = 3.0e-7;
	private static final double FPMIN = 1.0e-30;

	private static final int MAX_ITERATION = 100;

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
	public static Point3D gaussianProductCenter(double alpha1, Point3D a, double alpha2, Point3D b) {
		double gamma = alpha1 + alpha2;
		return new Point3D((alpha1 * a.getX() + alpha2 * b.getX()) / gamma,
				(alpha1 * a.getY() + alpha2 * b.getY()) / gamma, (alpha1 * a.getZ() + alpha2 * b.getZ()) / gamma);
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
	 * @param i
	 *            the index of Gaussian i
	 * @param j
	 *            the index of Gaussian j
	 * @param k
	 *            the index of Gaussian k
	 * @param l
	 *            the index of Gaussian l
	 * @return the index into the array
	 */
	public static int ijkl2intindex(int i, int j, int k, int l) {
		int temp;

		if (i < j) {
			temp = i;
			i = j;
			j = temp;
		}
		if (k < l) {
			temp = k;
			k = l;
			l = temp;
		}

		int ij = i * (i + 1) / 2 + j;
		int kl = k * (k + 1) / 2 + l;

		if (ij < kl) {
			temp = ij;
			ij = kl;
			kl = temp;
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
	public static double computeFGamma(int m, double x) {
		x = Math.max(Math.abs(x), SMALL);

		return (0.5 * Math.pow(x, -m - 0.5) * gammaIncomplete(m + 0.5, x));
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
	public static double gammaIncomplete(double a, double x) {
		double gammap;
		double gln;

		gln = logGamma(a);

		if (x < (a + 1.0)) {
			// Computes the incomplete Gamma function by using the series representation.
			// See NumRec sect 6.1.
			if (x != 0.0) {
				double ap = a;
				double sum;
				double delta = sum = 1.0 / a;

				for (int i = 0; i < MAX_ITERATION; i++) {
					ap++;
					delta *= x / ap;
					sum += delta;
					if (Math.abs(delta) < Math.abs(sum) * EPS)
						break;
				}

				gammap = sum * Math.exp(-x + a * Math.log(x) - gln);
			} else {
				gammap = 0.0;
			}
		} else {
			// Computes the incomplete Gamma function by using the continued fraction
			// representation. See NumRec sect 6.1
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

				if (Math.abs(d) < FPMIN)
					d = FPMIN;

				c = b + an / c;

				if (Math.abs(c) < FPMIN)
					c = FPMIN;

				d = 1.0 / d;
				delta = d * c;
				h *= delta;

				if (Math.abs(delta - 1.0) < EPS)
					break;
			}

			gammap = 1.0 - (Math.exp(-x + a * Math.log(x) - gln) * h);
		}

		return (Math.exp(gln) * gammap);
	}

	/**
	 * Computes the value of ln(gamma(x)) Numerical recipes, section 6.1
	 * 
	 * @param x
	 *            Value
	 * @return the complete Gamma(x) function
	 */
	public static double logGamma(double x) {
		double y = x;
		double tmp = x + 5.5;

		tmp -= (x + 0.5) * Math.log(tmp);

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

		return (-tmp + Math.log((2.5066282746310005 * ser) / x));
	}
}
