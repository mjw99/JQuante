package name.mjw.jquante.math.qm.integral;

import org.hipparchus.geometry.euclidean.threed.Vector3D;
import org.hipparchus.linear.RealMatrix;
import org.nd4j.linalg.api.ndarray.INDArray;
import org.nd4j.linalg.factory.Nd4j;

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
	public final double coulomb(ContractedGaussian a, ContractedGaussian b, ContractedGaussian c,
			ContractedGaussian d) {

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

		final INDArray aExps = Nd4j.create(a.getExponents());
		final INDArray aCoefs = Nd4j.create(a.getCoefficients());
		final INDArray aNorms = Nd4j.create(a.getPrimNorms());
		final INDArray aOrigin = Nd4j.create(a.getOrigin().toArray());
		final INDArray aPower = Nd4j.create(a.getPowers().toArray());

		final INDArray bExps = Nd4j.create(b.getExponents());
		final INDArray bCoefs = Nd4j.create(b.getCoefficients());
		final INDArray bNorms = Nd4j.create(b.getPrimNorms());
		final INDArray bOrigin = Nd4j.create(b.getOrigin().toArray());
		final INDArray bPower = Nd4j.create(b.getPowers().toArray());

		final INDArray cExps = Nd4j.create(c.getExponents());
		final INDArray cCoefs = Nd4j.create(c.getCoefficients());
		final INDArray cNorms = Nd4j.create(c.getPrimNorms());
		final INDArray cOrigin = Nd4j.create(c.getOrigin().toArray());
		final INDArray cPower = Nd4j.create(c.getPowers().toArray());

		final INDArray dExps = Nd4j.create(d.getExponents());
		final INDArray dCoefs = Nd4j.create(d.getCoefficients());
		final INDArray dNorms = Nd4j.create(d.getPrimNorms());
		final INDArray dOrigin = Nd4j.create(d.getOrigin().toArray());
		final INDArray dPower = Nd4j.create(d.getPowers().toArray());

		for (int i = 0; i < aExps.size(0); i++) {
			iaCoef = aCoefs.getDouble(i);
			iaExp = aExps.getDouble(i);
			iaNorm = aNorms.getDouble(i);

			for (int j = 0; j < bExps.size(0); j++) {
				jbCoef = bCoefs.getDouble(j);
				jbExp = bExps.getDouble(j);
				jbNorm = bNorms.getDouble(j);

				for (int k = 0; k < cExps.size(0); k++) {
					kcCoef = cCoefs.getDouble(k);
					kcExp = cExps.getDouble(k);
					kcNorm = cNorms.getDouble(k);

					for (int l = 0; l < dExps.size(0); l++) {
						repulsionTerm = coulombRepulsion(aOrigin.getDouble(0), aOrigin.getDouble(1),
								aOrigin.getDouble(2), iaNorm, aPower.getDouble(0), aPower.getDouble(1),
								aPower.getDouble(2), iaExp, bOrigin.getDouble(0), bOrigin.getDouble(1),
								bOrigin.getDouble(2), jbNorm, bPower.getDouble(0), bPower.getDouble(1),
								bPower.getDouble(2), jbExp, cOrigin.getDouble(0), cOrigin.getDouble(1),
								cOrigin.getDouble(2), kcNorm, cPower.getDouble(0), cPower.getDouble(1),
								cPower.getDouble(2), kcExp, dOrigin.getDouble(0), dOrigin.getDouble(1),
								dOrigin.getDouble(2), dNorms.getDouble(l), dPower.getDouble(0), dPower.getDouble(1),
								dPower.getDouble(2), dExps.getDouble(l));

						jij += iaCoef * jbCoef * kcCoef * dCoefs.getDouble(l) * repulsionTerm;
					}
				}
			}
		}

		return (a.getNormalization() * b.getNormalization() * c.getNormalization() * d.getNormalization() * jij);
	}

	public final double coulombRepulsion(double ax, double ay, double az, double aNorm, double al, double am, double an,
			double aAlpha, double bx, double by, double bz, double bNorm, double bl, double bm, double bn,
			double bAlpha, double cx, double cy, double cz, double cNorm, double cl, double cm, double cn,
			double cAlpha, double dx, double dy, double dz, double dNorm, double dl, double dm, double dn,
			double dAlpha) {

		double sum = 0.0;

		final double radiusABSquared = ((ax - bx) * (ax - bx) + (ay - by) * (ay - by) + (az - bz) * (az - bz));
		final double radiusCDSquared = ((cx - dx) * (cx - dx) + (cy - dy) * (cy - dy) + (cz - dz) * (cz - dz));

		final double[] p = IntegralsUtil.gaussianProductCenter(aAlpha, ax, ay, az, bAlpha, bx, by, bz);
		final double[] q = IntegralsUtil.gaussianProductCenter(cAlpha, cx, cy, cz, dAlpha, dx, dy, dz);

		final double radiusPQSquared = (p[0] - q[0]) * (p[0] - q[0]) + (p[1] - q[1]) * (p[1] - q[1])
				+ (p[2] - q[2]) * (p[2] - q[2]);

		final double gamma1 = aAlpha + bAlpha;
		final double gamma2 = cAlpha + dAlpha;
		final double delta = 0.25 * (1 / gamma1 + 1 / gamma2);

		final double quartRadiusPQSquaredOverDelta = 0.25 * radiusPQSquared / delta;

		final INDArray bxx = constructBArray((int) al, (int) bl, (int) cl, (int) dl, p[0], ax, bx, q[0], cx, dx, gamma1,
				gamma2, delta);

		final INDArray byy = constructBArray((int) am, (int) bm, (int) cm, (int) dm, p[1], ay, by, q[1], cy, dy, gamma1,
				gamma2, delta);

		final INDArray bzz = constructBArray((int) an, (int) bn, (int) cn, (int) dn, p[2], az, bz, q[2], cz, dz, gamma1,
				gamma2, delta);

		for (int i = 0; i < bxx.size(0); i++) {
			for (int j = 0; j < byy.size(0); j++) {
				for (int k = 0; k < bzz.size(0); k++) {
					sum += bxx.getDouble(i) * byy.getDouble(j) * bzz.getDouble(k)
							* IntegralsUtil.computeFGamma(i + j + k, quartRadiusPQSquaredOverDelta);
				}
			}
		}

		return (2 * FastMath.pow(Math.PI, 2.5) / (gamma1 * gamma2 * FastMath.sqrt(gamma1 + gamma2))
				* FastMath.exp(-aAlpha * bAlpha * radiusABSquared / gamma1)
				* FastMath.exp(-cAlpha * dAlpha * radiusCDSquared / gamma2) * sum * aNorm * bNorm * cNorm * dNorm);

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

		final INDArray bx = constructBArray(aPower.getL(), bPower.getL(), cPower.getL(), dPower.getL(), p.getX(), a.getX(),
				b.getX(), q.getX(), c.getX(), d.getX(), gamma1, gamma2, delta);

		final INDArray by = constructBArray(aPower.getM(), bPower.getM(), cPower.getM(), dPower.getM(), p.getY(), a.getY(),
				b.getY(), q.getY(), c.getY(), d.getY(), gamma1, gamma2, delta);

		final INDArray bz = constructBArray(aPower.getN(), bPower.getN(), cPower.getN(), dPower.getN(), p.getZ(), a.getZ(),
				b.getZ(), q.getZ(), c.getZ(), d.getZ(), gamma1, gamma2, delta);

		for (int i = 0; i < bx.size(0); i++) {
			for (int j = 0; j < by.size(0); j++) {
				for (int k = 0; k < bz.size(0); k++) {
					sum += bx.getDouble(i) * by.getDouble(j) * bz.getDouble(k)
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
	private final INDArray constructBArray(final int l1, final int l2, final int l3, final int l4, final double p,
			final double a, final double b, final double q, final double c, final double d, final double g1,
			final double g2, final double delta) {

		int i1;
		int i2;
		int r1;
		int r2;
		int u;
		int index;

		final int iMax = l1 + l2 + l3 + l4 + 1;
		INDArray bArr = Nd4j.zeros(iMax);

		for (i1 = 0; i1 < (l1 + l2 + 1); i1++) {
			for (i2 = 0; i2 < (l3 + l4 + 1); i2++) {
				for (r1 = 0; r1 < (i1 / 2 + 1); r1++) {
					for (r2 = 0; r2 < (i2 / 2 + 1); r2++) {
						for (u = 0; u < ((i1 + i2) / 2 - r1 - r2 + 1); u++) {
							index = i1 + i2 - 2 * (r1 + r2) - u;

							bArr.putScalar(index,
									constructBTerm(i1, i2, r1, r2, u, l1, l2, l3, l4, p, a, b, q, c, d, g1, g2, delta)
											+ bArr.getDouble(index));
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
	private final double functionB(final int i, final int l1, final int l2, final double p, final double a, final double b,
			final int r, final double g) {
		return (MathUtil.binomialPrefactor(i, l1, l2, p - a, p - b) * functionB0(i, r, g));
	}

	/**
	 * the function B0, taken from PyQuante
	 */
	private final double functionB0(final int i, final int r, final double g) {
		return (MathUtil.factorialRatioSquared(i, r) * FastMath.pow(4 * g, r - i));
	}

	@Override
	public final double coulomb(ContractedGaussian a, ContractedGaussian b, ContractedGaussian c, ContractedGaussian d,
			Density density, RealMatrix jMat, RealMatrix kMat) {
		throw new UnsupportedOperationException("Not supported yet.");
	}
}
