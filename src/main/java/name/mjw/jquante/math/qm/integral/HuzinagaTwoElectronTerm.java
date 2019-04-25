package name.mjw.jquante.math.qm.integral;

import org.apache.commons.math3.geometry.euclidean.threed.Vector3D;
import org.apache.commons.math3.linear.RealMatrix;
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

		int i;
		int j;
		int k;
		int l;
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

		INDArray aExps = Nd4j.create(a.getExponents());
		INDArray aCoefs = Nd4j.create(a.getCoefficients());
		INDArray aNorms = Nd4j.create(a.getPrimNorms());
		INDArray aOrigin = Nd4j.create(a.getOrigin().toArray());
		INDArray aPower = Nd4j.create(a.getPowers().toArray());

		INDArray bExps = Nd4j.create(b.getExponents());
		INDArray bCoefs = Nd4j.create(b.getCoefficients());
		INDArray bNorms = Nd4j.create(b.getPrimNorms());
		INDArray bOrigin = Nd4j.create(b.getOrigin().toArray());
		INDArray bPower = Nd4j.create(b.getPowers().toArray());

		INDArray cExps = Nd4j.create(c.getExponents());
		INDArray cCoefs = Nd4j.create(c.getCoefficients());
		INDArray cNorms = Nd4j.create(c.getPrimNorms());
		INDArray cOrigin = Nd4j.create(c.getOrigin().toArray());
		INDArray cPower = Nd4j.create(c.getPowers().toArray());

		INDArray dExps = Nd4j.create(d.getExponents());
		INDArray dCoefs = Nd4j.create(d.getCoefficients());
		INDArray dNorms = Nd4j.create(d.getPrimNorms());
		INDArray dOrigin = Nd4j.create(d.getOrigin().toArray());
		INDArray dPower = Nd4j.create(d.getPowers().toArray());

		for (i = 0; i < aExps.size(1); i++) {
			iaCoef = aCoefs.getDouble(i);
			iaExp = aExps.getDouble(i);
			iaNorm = aNorms.getDouble(i);

			for (j = 0; j < bExps.size(1); j++) {
				jbCoef = bCoefs.getDouble(j);
				jbExp = bExps.getDouble(j);
				jbNorm = bNorms.getDouble(j);

				for (k = 0; k < cExps.size(1); k++) {
					kcCoef = cCoefs.getDouble(k);
					kcExp = cExps.getDouble(k);
					kcNorm = cNorms.getDouble(k);

					for (l = 0; l < dExps.size(1); l++) {
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
		int i;
		int j;
		int k;

		double radiusABSquared = ((ax - bx) * (ax - bx) + (ay - by) * (ay - by) + (az - bz) * (az - bz));
		double radiusCDSquared = ((cx - dx) * (cx - dx) + (cy - dy) * (cy - dy) + (cz - dz) * (cz - dz));

		double[] p = IntegralsUtil.gaussianProductCenter(aAlpha, ax, ay, az, bAlpha, bx, by, bz);
		double[] q = IntegralsUtil.gaussianProductCenter(cAlpha, cx, cy, cz, dAlpha, dx, dy, dz);

		double radiusPQSquared = (p[0] - q[0]) * (p[0] - q[0]) + (p[1] - q[1]) * (p[1] - q[1])
				+ (p[2] - q[2]) * (p[2] - q[2]);

		double gamma1 = aAlpha + bAlpha;
		double gamma2 = cAlpha + dAlpha;
		double delta = 0.25 * (1 / gamma1 + 1 / gamma2);

		double quartRadiusPQSquaredOverDelta = 0.25 * radiusPQSquared / delta;

		INDArray bxx = constructBArray((int) al, (int) bl, (int) cl, (int) dl, p[0], ax, bx, q[0], cx, dx, gamma1,
				gamma2, delta);

		INDArray byy = constructBArray((int) am, (int) bm, (int) cm, (int) dm, p[1], ay, by, q[1], cy, dy, gamma1,
				gamma2, delta);

		INDArray bzz = constructBArray((int) an, (int) bn, (int) cn, (int) dn, p[2], az, bz, q[2], cz, dz, gamma1,
				gamma2, delta);

		for (i = 0; i < bxx.size(1); i++) {
			for (j = 0; j < byy.size(1); j++) {
				for (k = 0; k < bzz.size(1); k++) {
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
	public final double coulombRepulsion(Vector3D a, double aNorm, Power aPower, double aAlpha, Vector3D b,
			double bNorm, Power bPower, double bAlpha, Vector3D c, double cNorm, Power cPower, double cAlpha,
			Vector3D d, double dNorm, Power dPower, double dAlpha) {

		double sum = 0.0;
		int i;
		int j;
		int k;

		double radiusABSquared = a.distanceSq(b);
		double radiusCDSquared = c.distanceSq(d);

		Vector3D p = IntegralsUtil.gaussianProductCenter(aAlpha, a, bAlpha, b);
		Vector3D q = IntegralsUtil.gaussianProductCenter(cAlpha, c, dAlpha, d);

		double radiusPQSquared = p.distanceSq(q);

		double gamma1 = aAlpha + bAlpha;
		double gamma2 = cAlpha + dAlpha;
		double delta = 0.25 * (1 / gamma1 + 1 / gamma2);

		double quartRadiusPQSquaredOverDelta = 0.25 * radiusPQSquared / delta;

		INDArray bx = constructBArray(aPower.getL(), bPower.getL(), cPower.getL(), dPower.getL(), p.getX(), a.getX(),
				b.getX(), q.getX(), c.getX(), d.getX(), gamma1, gamma2, delta);

		INDArray by = constructBArray(aPower.getM(), bPower.getM(), cPower.getM(), dPower.getM(), p.getY(), a.getY(),
				b.getY(), q.getY(), c.getY(), d.getY(), gamma1, gamma2, delta);

		INDArray bz = constructBArray(aPower.getN(), bPower.getN(), cPower.getN(), dPower.getN(), p.getZ(), a.getZ(),
				b.getZ(), q.getZ(), c.getZ(), d.getZ(), gamma1, gamma2, delta);

		for (i = 0; i < bx.size(1); i++) {
			for (j = 0; j < by.size(1); j++) {
				for (k = 0; k < bz.size(1); k++) {
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
	private INDArray constructBArray(int l1, int l2, int l3, int l4, double p, double a, double b, double q, double c,
			double d, double g1, double g2, double delta) {

		int i1;
		int i2;
		int r1;
		int r2;
		int u;
		int index;

		int iMax = l1 + l2 + l3 + l4 + 1;
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
	private double constructBTerm(int i1, int i2, int r1, int r2, int u, int l1, int l2, int l3, int l4, double px,
			double ax, double bx, double qx, double cx, double dx, double gamma1, double gamma2, double delta) {

		return (functionB(i1, l1, l2, px, ax, bx, r1, gamma1) * FastMath.pow(-1, i2)
				* functionB(i2, l3, l4, qx, cx, dx, r2, gamma2) * FastMath.pow(-1, u)
				* MathUtil.factorialRatioSquared(i1 + i2 - 2 * (r1 + r2), u)
				* FastMath.pow(qx - px, i1 + i2 - 2 * (r1 + r2) - 2d * u)
				/ Math.pow(delta, i1 + i2 - 2d * (r1 + r2) - u));
	}

	/**
	 * the function B, taken from PyQuante
	 */
	private double functionB(int i, int l1, int l2, double p, double a, double b, int r, double g) {
		return (MathUtil.binomialPrefactor(i, l1, l2, p - a, p - b) * functionB0(i, r, g));
	}

	/**
	 * the function B0, taken from PyQuante
	 */
	private double functionB0(int i, int r, double g) {
		return (MathUtil.factorialRatioSquared(i, r) * FastMath.pow(4 * g, r - i));
	}

	@Override
	public final double coulomb(ContractedGaussian a, ContractedGaussian b, ContractedGaussian c, ContractedGaussian d,
			Density density, RealMatrix jMat, RealMatrix kMat) {
		throw new UnsupportedOperationException("Not supported yet.");
	}
}
