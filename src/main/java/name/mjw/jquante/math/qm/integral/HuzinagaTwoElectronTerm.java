package name.mjw.jquante.math.qm.integral;

import java.util.ArrayList;

import org.hipparchus.geometry.euclidean.threed.Vector3D;
import org.hipparchus.linear.RealMatrix;

import name.mjw.jquante.math.MathUtil;
import name.mjw.jquante.math.qm.Density;
import name.mjw.jquante.math.qm.basis.ContractedGaussian;
import name.mjw.jquante.math.qm.basis.Power;
import net.jafama.FastMath;

/**
 * The Huzinaga integral package.
 * 
 * The equations herein are based upon: <br>
 * 'Gaussian Expansion Methods for Molecular Orbitals.'
 * <a href="http://dx.doi.org/10.1143/JPSJ.21.2313"> H. Taketa, S. Huzinaga, and
 * K. O-ohata. <i> H. Phys. Soc. Japan, </i> <b>21</b>, 2313, 1966.)</a> [THO
 * paper]. <br>
 * and PyQuante (<a href="http://pyquante.sf.net"> http://pyquante.sf.net </a>).
 * 
 * @author V.Ganesh
 * @version 2.0 (Part of MeTA v2.0)
 */
public final class HuzinagaTwoElectronTerm implements TwoElectronTerm {
	/**
	 * 2E coulomb interactions between four contracted Gaussians
	 */
	@Override
	public final double coulomb(final ContractedGaussian a, final ContractedGaussian b, final ContractedGaussian c,
			final ContractedGaussian d) {

		double jij = 0.0;

		double iaExp;
		double iaCoef;
		double iaNorm;
		double jbExp;
		double jbCoef;
		double jbNorm;
		double kcExp;
		double kcCoef;
		double kcNorm;
		double repulsionTerm;

		final ArrayList<Double> aExps = a.getExponents();
		final ArrayList<Double> aCoefs = a.getCoefficients();
		final ArrayList<Double> aNorms = a.getPrimNorms();
		final Vector3D aOrigin = a.getOrigin();
		final Power aPower = a.getPowers();

		final ArrayList<Double> bExps = b.getExponents();
		final ArrayList<Double> bCoefs = b.getCoefficients();
		final ArrayList<Double> bNorms = b.getPrimNorms();
		final Vector3D bOrigin = b.getOrigin();
		final Power bPower = b.getPowers();

		final ArrayList<Double> cExps = c.getExponents();
		final ArrayList<Double> cCoefs = c.getCoefficients();
		final ArrayList<Double> cNorms = c.getPrimNorms();
		final Vector3D cOrigin = c.getOrigin();
		final Power cPower = c.getPowers();

		final ArrayList<Double> dExps = d.getExponents();
		final ArrayList<Double> dCoefs = d.getCoefficients();
		final ArrayList<Double> dNorms = d.getPrimNorms();
		final Vector3D dOrigin = d.getOrigin();
		final Power dPower = d.getPowers();

		final int asz = aExps.size();
		final int bsz = bExps.size();
		final int csz = cExps.size();
		final int dsz = dExps.size();

		for (int i = 0; i < asz; i++) {
			iaCoef = aCoefs.get(i);
			iaExp = aExps.get(i);
			iaNorm = aNorms.get(i);

			for (int j = 0; j < bsz; j++) {
				jbCoef = bCoefs.get(j);
				jbExp = bExps.get(j);
				jbNorm = bNorms.get(j);

				for (int k = 0; k < csz; k++) {
					kcCoef = cCoefs.get(k);
					kcExp = cExps.get(k);
					kcNorm = cNorms.get(k);

					for (int l = 0; l < dsz; l++) {
						repulsionTerm = coulombRepulsion(aOrigin, iaNorm, aPower, iaExp, bOrigin, jbNorm, bPower, jbExp,
								cOrigin, kcNorm, cPower, kcExp, dOrigin, dNorms.get(l), dPower, dExps.get(l));

						jij += iaCoef * jbCoef * kcCoef * dCoefs.get(l) * repulsionTerm;
					}
				}
			}
		}

		return (a.getNormalization() * b.getNormalization() * c.getNormalization() * d.getNormalization() * jij);
	}

	/**
	 * coulomb repulsion term
	 */
	@Override
	public final double coulombRepulsion(final Vector3D a, final double aNorm, final Power aPower, final double aAlpha,
			final Vector3D b, final double bNorm, final Power bPower, final double bAlpha, final Vector3D c,
			final double cNorm, final Power cPower, final double cAlpha, final Vector3D d, final double dNorm,
			final Power dPower, final double dAlpha) {

		double sum = 0.0;

		final double radiusABSquared = a.distanceSq(b);
		final double radiusCDSquared = c.distanceSq(d);

		final Vector3D p = IntegralsUtil.gaussianProductCenter(aAlpha, a, bAlpha, b);
		final Vector3D q = IntegralsUtil.gaussianProductCenter(cAlpha, c, dAlpha, d);

		final double radiusPQSquared = p.distanceSq(q);

		final double gamma1 = aAlpha + bAlpha;
		final double gamma2 = cAlpha + dAlpha;
		final double delta = 0.25 * (1 / gamma1 + 1 / gamma2);

		final double quartRadiusPQSquaredOverDelta = 0.25 * radiusPQSquared / delta;

		final double[] bx = constructBArray(aPower.getL(), bPower.getL(), cPower.getL(), dPower.getL(), p.getX(),
				a.getX(), b.getX(), q.getX(), c.getX(), d.getX(), gamma1, gamma2, delta);

		final double[] by = constructBArray(aPower.getM(), bPower.getM(), cPower.getM(), dPower.getM(), p.getY(),
				a.getY(), b.getY(), q.getY(), c.getY(), d.getY(), gamma1, gamma2, delta);

		final double[] bz = constructBArray(aPower.getN(), bPower.getN(), cPower.getN(), dPower.getN(), p.getZ(),
				a.getZ(), b.getZ(), q.getZ(), c.getZ(), d.getZ(), gamma1, gamma2, delta);

		for (int i = 0; i < bx.length; i++) {
			for (int j = 0; j < by.length; j++) {
				for (int k = 0; k < bz.length; k++) {
					sum += bx[i] * by[j] * bz[k]
							* IntegralsUtil.computeFGamma(i + j + k, quartRadiusPQSquaredOverDelta);
				}
			}
		}

		return (2 * FastMath.pow(Math.PI, 2.5) / (gamma1 * gamma2 * FastMath.sqrt(gamma1 + gamma2))
				* FastMath.exp(-aAlpha * bAlpha * radiusABSquared / gamma1)
				* FastMath.exp(-cAlpha * dAlpha * radiusCDSquared / gamma2) * sum * aNorm * bNorm * cNorm * dNorm);
	}

	/**
	 * Construct B array.
	 * 
	 * <i> http://dx.doi.org/10.1143/JPSJ.21.2313 eq. 2.22 </i>
	 */
	private final double[] constructBArray(final int l1, final int l2, final int l3, final int l4, final double p,
			final double a, final double b, final double q, final double c, final double d, final double g1,
			final double g2, final double delta) {

		int i1;
		int i2;
		int r1;
		int r2;
		int u;
		int index;

		final int iMax = l1 + l2 + l3 + l4 + 1;
		final double[] bArr = new double[iMax];

		for (i1 = 0; i1 < (l1 + l2 + 1); i1++) {
			for (i2 = 0; i2 < (l3 + l4 + 1); i2++) {
				for (r1 = 0; r1 < (i1 / 2 + 1); r1++) {
					for (r2 = 0; r2 < (i2 / 2 + 1); r2++) {
						for (u = 0; u < ((i1 + i2) / 2 - r1 - r2 + 1); u++) {
							index = i1 + i2 - 2 * (r1 + r2) - u;

							bArr[index] += constructBTerm(i1, i2, r1, r2, u, l1, l2, l3, l4, p, a, b, q, c, d, g1, g2,
									delta);
						}
					}
				}
			}
		}

		return bArr;
	}

	/**
	 * Construct the B term
	 * 
	 * <i> http://dx.doi.org/10.1143/JPSJ.21.2313 eq. 2.22 </i>
	 */
	private final double constructBTerm(final int i1, final int i2, final int r1, final int r2, final int u,
			final int l1, final int l2, final int l3, final int l4, final double px, final double ax, final double bx,
			final double qx, final double cx, final double dx, final double gamma1, final double gamma2,
			final double delta) {

		return (functionB(i1, l1, l2, px, ax, bx, r1, gamma1) * FastMath.pow(-1, i2)
				* functionB(i2, l3, l4, qx, cx, dx, r2, gamma2) * FastMath.pow(-1, u)
				* MathUtil.factorialRatioSquared(i1 + i2 - 2 * (r1 + r2), u)
				* FastMath.pow(qx - px, i1 + i2 - 2 * (r1 + r2) - 2d * u)
				/ Math.pow(delta, i1 + i2 - 2d * (r1 + r2) - u));
	}

	/**
	 * the function B, taken from PyQuante
	 */
	private final double functionB(final int i, final int l1, final int l2, final double p, final double a,
			final double b, final int r, final double g) {
		return (MathUtil.binomialPrefactor(i, l1, l2, p - a, p - b) * functionB0(i, r, g));
	}

	/**
	 * the function B0, taken from PyQuante
	 */
	private final double functionB0(final int i, final int r, final double g) {
		return (MathUtil.factorialRatioSquared(i, r) * FastMath.pow(4 * g, r - i));
	}

	@Override
	public final double coulomb(final ContractedGaussian a, final ContractedGaussian b, final ContractedGaussian c,
			final ContractedGaussian d, final Density density, final RealMatrix jMat, final RealMatrix kMat) {
		throw new UnsupportedOperationException("Not supported yet.");
	}
}
