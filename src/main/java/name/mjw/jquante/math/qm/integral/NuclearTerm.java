package name.mjw.jquante.math.qm.integral;

import org.hipparchus.geometry.euclidean.threed.Vector3D;
import org.hipparchus.util.CombinatoricsUtils;

import name.mjw.jquante.math.MathUtil;
import name.mjw.jquante.math.qm.basis.Power;
import net.jafama.FastMath;

/**
 * 
 * @author V.Ganesh
 * @version 2.0 (Part of MeTA v2.0)
 */
public final class NuclearTerm implements IntegralsPackage {
	/**
	 * The nuclear attraction term.
	 * 
	 * <i> Taken from http://dx.doi.org/10.1143/JPSJ.21.2313 eq. 2.15 </i>
	 * 
	 * @param a      the coefficient of primitive Gaussian a.
	 * @param norm1  the normalization factor of primitive Gaussian a.
	 * @param power1 the orbital powers of primitive Gaussian a.
	 * @param alpha1 the coefficient of primitive Gaussian a.
	 * @param b      the coefficient of primitive Gaussian b.
	 * @param norm2  the normalization factor of primitive Gaussian b.
	 * @param power2 the orbital powers of primitive Gaussian b.
	 * @param alpha2 the coefficient of primitive Gaussian b.
	 * @param c      the location of nuclear centre.
	 * @return the nuclear attraction integral.
	 */
	public final double nuclearAttraction(final Vector3D a, final double norm1, final Power power1, final double alpha1,
			final Vector3D b, final double norm2, final Power power2, final double alpha2, final Vector3D c) {

		final Vector3D product = IntegralsUtil.gaussianProductCenter(alpha1, a, alpha2, b);

		final double rABSquared = a.distanceSq(b);
		final double rCPSquared = c.distanceSq(product);

		final double gamma = alpha1 + alpha2;

		final double[] ax = constructAArray(power1.getL(), power2.getL(), product.getX() - a.getX(),
				product.getX() - b.getX(), product.getX() - c.getX(), gamma);

		final double[] ay = constructAArray(power1.getM(), power2.getM(), product.getY() - a.getY(),
				product.getY() - b.getY(), product.getY() - c.getY(), gamma);

		final double[] az = constructAArray(power1.getN(), power2.getN(), product.getZ() - a.getZ(),
				product.getZ() - b.getZ(), product.getZ() - c.getZ(), gamma);

		double sum = 0.0;
		for (int i = 0; i < ax.length; i++) {
			for (int j = 0; j < ay.length; j++) {
				for (int k = 0; k < az.length; k++) {
					sum += ax[i] * ay[j] * az[k] * IntegralsUtil.computeFGamma(i + j + k, rCPSquared * gamma);
				}
			}
		}

		return (-norm1 * norm2 * 2.0 * FastMath.PI / gamma * FastMath.exp(-alpha1 * alpha2 * rABSquared / gamma) * sum);
	}

	/**
	 * <i> http://dx.doi.org/10.1143/JPSJ.21.2313 eq. 2.18 and 3.1 </i>
	 * 
	 * @param l1    the angular momentum number of Gaussian 1.
	 * @param l2    the angular momentum number of Gaussian 2.
	 * @param pa    the distance of Gaussian 1 to the product centre.
	 * @param pb    the distance of Gaussian 2 to the product centre.
	 * @param pc    the distance of the nucleus to the product centre.
	 * @param gamma the sum of both Gaussians' exponent.
	 * @return the nuclear attraction integral.
	 */
	public final double[] constructAArray(final int l1, final int l2, final double pa, final double pb, final double pc,
			final double gamma) {
		final int iMax = l1 + l2 + 1;
		final double[] a = new double[iMax];

		int index;

		for (int i = 0; i < iMax; i++) {
			for (int r = 0; r < ((int) (FastMath.floor(i / 2.0) + 1.0)); r++) {
				for (int u = 0; u < ((int) (FastMath.floor((i - 2.0 * r) / 2.0) + 1.0)); u++) {
					index = i - 2 * r - u;

					a[index] += constructATerm(i, r, u, l1, l2, pa, pb, pc, gamma);
				}
			}
		}

		return a;
	}

	public final double[] constructGradAArray(final int l1, final int l2, final double pa, final double pb,
			final double pc, final double gamma) {
		final int iMax = l1 + l2 + 1;
		final double[] a = new double[iMax];

		int index;

		for (int i = 0; i < iMax; i++) {
			for (int r = 0; r < ((int) (FastMath.floor(i / 2.0) + 1.0)); r++) {
				for (int u = 0; u < ((int) (FastMath.floor((i - 2.0 * r) / 2.0) + 1.0)); u++) {
					index = i - 2 * r - u;

					a[index] += (i - 2.0 * r + 1.0) / (i - 2.0 * r - 2.0 * u + 1.0) * pc
							* constructATerm(i, r, u, l1, l2, pa, pb, pc, gamma);
				}
			}
		}

		return a;
	}

	/**
	 * The A term <br>
	 * <i> http://dx.doi.org/10.1143/JPSJ.21.2313 eq. 2.18 </i>
	 * 
	 * @param i     index i.
	 * @param r     index r.
	 * @param u     index u.
	 * @param l1    the angular momentum number of Gaussian 1.
	 * @param l2    the angular momentum number of Gaussian 2.
	 * @param pax   the distance of Gaussian 1 to the product centre.
	 * @param pbx   the distance of Gaussian 2 to the product centre.
	 * @param pcx   the distance of the nucleus to the product centre.
	 * @param gamma the sum of both Gaussians' exponent.
	 * @return the A term for index i,r,u.
	 */
	public final double constructATerm(final int i, final int r, final int u, final int l1, final int l2,
			final double pax, final double pbx, final double pcx, final double gamma) {
		return (FastMath.pow(-1.0, i) * MathUtil.binomialPrefactor(i, l1, l2, pax, pbx) * FastMath.pow(-1.0, u)
				* CombinatoricsUtils.factorialDouble(i) * FastMath.pow(pcx, i - 2 * r - 2 * u)
				* FastMath.pow(0.25 / gamma, r + u) / CombinatoricsUtils.factorialDouble(r)
				/ CombinatoricsUtils.factorialDouble(u) / CombinatoricsUtils.factorialDouble(i - 2 * r - 2 * u));
	}

	/**
	 * The nuclear attraction gradient term
	 * 
	 * @param a      the coefficient of primitive Gaussian a.
	 * @param power1 the orbital powers of primitive Gaussian a.
	 * @param alpha1 the coefficient of primitive Gaussian a.
	 * @param b      the coefficient of primitive Gaussian b.
	 * @param power2 the orbital powers of primitive Gaussian b.
	 * @param alpha2 the coefficient of primitive Gaussian b.
	 * @param c      the location of nuclear centre.
	 * @return the nuclear attraction gradient.
	 */
	public final Vector3D nuclearAttractionGradient(final Vector3D a, final Power power1, final double alpha1,
			final Vector3D b, final Power power2, final double alpha2, final Vector3D c) {
		final Vector3D product = IntegralsUtil.gaussianProductCenter(alpha1, a, alpha2, b);

		final double rABSquared = a.distanceSq(b);
		final double rCPSquared = c.distanceSq(product);

		final double gamma = alpha1 + alpha2;

		final double[] ax = constructAArray(power1.getL(), power2.getL(), product.getX() - a.getX(),
				product.getX() - b.getX(), product.getX() - c.getX(), gamma);

		final double[] ay = constructAArray(power1.getM(), power2.getM(), product.getY() - a.getY(),
				product.getY() - b.getY(), product.getY() - c.getY(), gamma);

		final double[] az = constructAArray(power1.getN(), power2.getN(), product.getZ() - a.getZ(),
				product.getZ() - b.getZ(), product.getZ() - c.getZ(), gamma);

		final double[] dax = constructGradAArray(power1.getL(), power2.getL(), product.getX() - a.getX(),
				product.getX() - b.getX(), product.getX() - c.getX(), gamma);

		final double[] day = constructGradAArray(power1.getM(), power2.getM(), product.getY() - a.getY(),
				product.getY() - b.getY(), product.getY() - c.getY(), gamma);

		final double[] daz = constructGradAArray(power1.getN(), power2.getN(), product.getZ() - a.getZ(),
				product.getZ() - b.getZ(), product.getZ() - c.getZ(), gamma);

		double gradX = 0.0;
		double gradY = 0.0;
		double gradZ = 0.0;

		for (int i = 0; i < dax.length; i++) {
			for (int j = 0; j < day.length; j++) {
				for (int k = 0; k < daz.length; k++) {
					final double fgamma = IntegralsUtil.computeFGamma(i + j + k + 1, rCPSquared * gamma);

					gradX += dax[i] * ay[j] * az[k] * fgamma;
					gradY += ax[i] * day[j] * az[k] * fgamma;
					gradZ += ax[i] * ay[j] * daz[k] * fgamma;
				}
			}
		}

		final double factor = -4 * FastMath.PI * FastMath.exp(-alpha1 * alpha2 * rABSquared / gamma);

		gradX = factor * gradX;
		gradY = factor * gradY;
		gradZ = factor * gradZ;

		return new Vector3D(gradX, gradY, gradZ);
	}
}
