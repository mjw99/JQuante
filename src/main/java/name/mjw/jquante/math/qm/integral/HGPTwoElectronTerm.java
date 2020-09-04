package name.mjw.jquante.math.qm.integral;

import java.util.ArrayList;

import name.mjw.jquante.math.qm.Density;
import name.mjw.jquante.math.qm.basis.ContractedGaussian;
import name.mjw.jquante.math.qm.basis.Power;
import net.jafama.FastMath;

import org.hipparchus.geometry.euclidean.threed.Vector3D;
import org.hipparchus.linear.RealMatrix;

/**
 * Head-Gordon/Pople scheme of evaluating two-electron integrals.
 * 
 * The code is based upon PyQuante (<a href="http://pyquante.sf.net">
 * http://pyquante.sf.net </a>). See 'A method for two‐electron Gaussian
 * integral and integral derivative evaluation using recurrence relations'
 * <a href="http://dx.doi.org/10.1063/1.455553"> M. Head-Gordon and J. A. Pople,
 * J. Chem. Phys. <b>89</b>, 5777 (1988)</a> for more details.
 * 
 * @author V.Ganesh
 * @version 2.0 (Part of MeTA v2.0)
 */
public final class HGPTwoElectronTerm implements TwoElectronTerm {

	private final double sqrt2PI = FastMath.sqrt(2.0) * FastMath.pow(Math.PI, 1.25);

	/**
	 * 2E coulomb interactions between four contracted Gaussians.
	 */
	@Override
	public final double coulomb(final ContractedGaussian a, final ContractedGaussian b, final ContractedGaussian c,
			final ContractedGaussian d) {
		return (a.getNormalization() * b.getNormalization() * c.getNormalization() * d.getNormalization()
				* contractedHrr(a.getOrigin(), a.getPowers(), a.getCoefficients(), a.getExponents(), a.getPrimNorms(),
						b.getOrigin(), b.getPowers(), b.getCoefficients(), b.getExponents(), b.getPrimNorms(),
						c.getOrigin(), c.getPowers(), c.getCoefficients(), c.getExponents(), c.getPrimNorms(),
						d.getOrigin(), d.getPowers(), d.getCoefficients(), d.getExponents(), d.getPrimNorms()));
	}

	/**
	 * Coulomb repulsion term
	 */
	@Override
	public final double coulombRepulsion(final Vector3D a, final double aNorm, final Power aPower, final double aAlpha,
			final Vector3D b, final double bNorm, final Power bPower, final double bAlpha, final Vector3D c,
			final double cNorm, final Power cPower, final double cAlpha, final Vector3D d, final double dNorm,
			final Power dPower, final double dAlpha) {
		return vrr(a, aNorm, aPower, aAlpha, b, bNorm, bAlpha, c, cNorm, cPower, cAlpha, d, dNorm, dAlpha, 0);
	}

	/**
	 * Contracted HRR (Horizontal Recurrence Relation)
	 * 
	 * @param a      Center of contracted Gaussian function a.
	 * @param aPower Angular momentum of Gaussian function a.
	 * @param aCoeff Coefficients of primitives in contracted Gaussian function a.
	 * @param aExps  Orbital exponents of primitives in contracted Gaussian function
	 *               a.
	 * @param aNorms Primitive normalisation coefficients of primitives in
	 *               contracted Gaussian function a.
	 * 
	 * @param b      Center of contracted Gaussian function b.
	 * @param bPower Angular momentum of Gaussian function b.
	 * @param bCoeff Coefficients of primitives in contracted Gaussian function b.
	 * @param bExps  Orbital exponents of primitives in contracted Gaussian function
	 *               b.
	 * @param bNorms Primitive normalisation coefficients of primitives in
	 *               contracted Gaussian function b.
	 * @param c      Center of contracted Gaussian function c.
	 * @param cPower Angular momentum of Gaussian function c.
	 * @param cCoeff Coefficients of primitives in contracted Gaussian function c.
	 * @param cExps  Orbital exponents of primitives in contracted Gaussian function
	 *               c.
	 * @param cNorms Primitive normalisation coefficients of primitives in
	 *               contracted Gaussian function c.
	 * @param d      Center of contracted Gaussian function d.
	 * @param dPower Angular momentum of Gaussian function d.
	 * @param dCoeff Coefficients of primitives in contracted Gaussian function d.
	 * @param dExps  Orbital exponents of primitives in contracted Gaussian function
	 *               d.
	 * @param dNorms Primitive normalisation coefficients of primitives in
	 *               contracted Gaussian function d.
	 * 
	 * @return Contribution to Horizontal Recurrence Relation.
	 */
	protected final double contractedHrr(final Vector3D a, final Power aPower, final ArrayList<Double> aCoeff,
			final ArrayList<Double> aExps, final ArrayList<Double> aNorms, final Vector3D b, final Power bPower,
			final ArrayList<Double> bCoeff, final ArrayList<Double> bExps, final ArrayList<Double> bNorms,
			final Vector3D c, final Power cPower, final ArrayList<Double> cCoeff, final ArrayList<Double> cExps,
			final ArrayList<Double> cNorms, final Vector3D d, final Power dPower, final ArrayList<Double> dCoeff,
			final ArrayList<Double> dExps, final ArrayList<Double> dNorms) {

		final int la = aPower.getL();
		final int ma = aPower.getM();
		final int na = aPower.getN();

		final int lb = bPower.getL();
		final int mb = bPower.getM();
		final int nb = bPower.getN();

		final int lc = cPower.getL();
		final int mc = cPower.getM();
		final int nc = cPower.getN();

		final int ld = dPower.getL();
		final int md = dPower.getM();
		final int nd = dPower.getN();

		if (lb > 0) {
			final Power newBPower = new Power(lb - 1, mb, nb);
			return (contractedHrr(a, new Power(la + 1, ma, na), aCoeff, aExps, aNorms, b, newBPower, bCoeff, bExps,
					bNorms, c, cPower, cCoeff, cExps, cNorms, d, dPower, dCoeff, dExps, dNorms)
					+ (a.getX() - b.getX()) * contractedHrr(a, aPower, aCoeff, aExps, aNorms, b, newBPower, bCoeff,
							bExps, bNorms, c, cPower, cCoeff, cExps, cNorms, d, dPower, dCoeff, dExps, dNorms));
		} else if (mb > 0) {
			final Power newBPower = new Power(lb, mb - 1, nb);
			return (contractedHrr(a, new Power(la, ma + 1, na), aCoeff, aExps, aNorms, b, newBPower, bCoeff, bExps,
					bNorms, c, cPower, cCoeff, cExps, cNorms, d, dPower, dCoeff, dExps, dNorms)
					+ (a.getY() - b.getY()) * contractedHrr(a, aPower, aCoeff, aExps, aNorms, b, newBPower, bCoeff,
							bExps, bNorms, c, cPower, cCoeff, cExps, cNorms, d, dPower, dCoeff, dExps, dNorms));
		} else if (nb > 0) {
			final Power newBPower = new Power(lb, mb, nb - 1);
			return (contractedHrr(a, new Power(la, ma, na + 1), aCoeff, aExps, aNorms, b, newBPower, bCoeff, bExps,
					bNorms, c, cPower, cCoeff, cExps, cNorms, d, dPower, dCoeff, dExps, dNorms)
					+ (a.getZ() - b.getZ()) * contractedHrr(a, aPower, aCoeff, aExps, aNorms, b, newBPower, bCoeff,
							bExps, bNorms, c, cPower, cCoeff, cExps, cNorms, d, dPower, dCoeff, dExps, dNorms));
		} else if (ld > 0) {
			final Power newDPower = new Power(ld - 1, md, nd);
			return (contractedHrr(a, aPower, aCoeff, aExps, aNorms, b, bPower, bCoeff, bExps, bNorms, c,
					new Power(lc + 1, mc, nc), cCoeff, cExps, cNorms, d, newDPower, dCoeff, dExps, dNorms)
					+ (c.getX() - d.getX()) * contractedHrr(a, aPower, aCoeff, aExps, aNorms, b, bPower, bCoeff, bExps,
							bNorms, c, cPower, cCoeff, cExps, cNorms, d, newDPower, dCoeff, dExps, dNorms));
		} else if (md > 0) {
			final Power newDPower = new Power(ld, md - 1, nd);
			return (contractedHrr(a, aPower, aCoeff, aExps, aNorms, b, bPower, bCoeff, bExps, bNorms, c,
					new Power(lc, mc + 1, nc), cCoeff, cExps, cNorms, d, newDPower, dCoeff, dExps, dNorms)
					+ (c.getY() - d.getY()) * contractedHrr(a, aPower, aCoeff, aExps, aNorms, b, bPower, bCoeff, bExps,
							bNorms, c, cPower, cCoeff, cExps, cNorms, d, newDPower, dCoeff, dExps, dNorms));
		} else if (nd > 0) {
			final Power newDPower = new Power(ld, md, nd - 1);
			return (contractedHrr(a, aPower, aCoeff, aExps, aNorms, b, bPower, bCoeff, bExps, bNorms, c,
					new Power(lc, mc, nc + 1), cCoeff, cExps, cNorms, d, newDPower, dCoeff, dExps, dNorms)
					+ (c.getZ() - d.getZ()) * contractedHrr(a, aPower, aCoeff, aExps, aNorms, b, bPower, bCoeff, bExps,
							bNorms, c, cPower, cCoeff, cExps, cNorms, d, newDPower, dCoeff, dExps, dNorms));
		}

		return contractedVrr(a, aPower, aCoeff, aExps, aNorms, b, bCoeff, bExps, bNorms, c, cPower, cCoeff, cExps,
				cNorms, d, dCoeff, dExps, dNorms);
	}

	/**
	 * Contracted VRR (Vertical Recurrence Relation) contribution
	 * 
	 * @param a      Center of contracted Gaussian function a.
	 * @param aPower Angular momentum of Gaussian function a.
	 * @param aCoeff Coefficients of primitives in contracted Gaussian function a.
	 * @param aExps  Orbital exponents of primitives in contracted Gaussian function
	 *               a.
	 * @param aNorms Primitive normalisation coefficients of primitives in
	 *               contracted Gaussian function a.
	 * 
	 * @param b      Center of contracted Gaussian function b. *
	 * @param bCoeff Coefficients of primitives in contracted Gaussian function b.
	 * @param bExps  Orbital exponents of primitives in contracted Gaussian function
	 *               b.
	 * @param bNorms Primitive normalisation coefficients of primitives in
	 *               contracted Gaussian function b.
	 * 
	 * @param c      Center of contracted Gaussian function c.
	 * @param cPower Angular momentum of Gaussian function c.
	 * @param cCoeff Coefficients of primitives in contracted Gaussian function c.
	 * @param cExps  Orbital exponents of primitives in contracted Gaussian function
	 *               c.
	 * @param cNorms Primitive normalisation coefficients of primitives in
	 *               contracted Gaussian function c.
	 * 
	 * @param d      Center of contracted Gaussian function d. *
	 * @param dCoeff Coefficients of primitives in contracted Gaussian function d.
	 * @param dExps  Orbital exponents of primitives in contracted Gaussian function
	 *               d.
	 * @param dNorms Primitive normalisation coefficients of primitives in
	 *               contracted Gaussian function d.
	 * 
	 * @return Contribution to Vertical Recurrence Relation.
	 */
	protected final double contractedVrr(final Vector3D a, final Power aPower, final ArrayList<Double> aCoeff,
			final ArrayList<Double> aExps, final ArrayList<Double> aNorms, final Vector3D b,
			final ArrayList<Double> bCoeff, final ArrayList<Double> bExps, final ArrayList<Double> bNorms,
			final Vector3D c, final Power cPower, final ArrayList<Double> cCoeff, final ArrayList<Double> cExps,
			final ArrayList<Double> cNorms, final Vector3D d, final ArrayList<Double> dCoeff,
			final ArrayList<Double> dExps, final ArrayList<Double> dNorms) {

		double value = 0.0;

		double iaExp;
		double iaCoef;
		double iaNorm;

		double jbExp;
		double jbCoef;
		double jbNorm;

		double kcExp;
		double kcCoef;
		double kcNorm;

		for (int i = 0; i < aExps.size(); i++) {
			iaCoef = aCoeff.get(i);
			iaExp = aExps.get(i);
			iaNorm = aNorms.get(i);

			for (int j = 0; j < bExps.size(); j++) {
				jbCoef = bCoeff.get(j);
				jbExp = bExps.get(j);
				jbNorm = bNorms.get(j);

				for (int k = 0; k < cExps.size(); k++) {
					kcCoef = cCoeff.get(k);
					kcExp = cExps.get(k);
					kcNorm = cNorms.get(k);

					for (int l = 0; l < dExps.size(); l++) {
						value += iaCoef * jbCoef * kcCoef * dCoeff.get(l) * vrr(a, iaNorm, aPower, iaExp, b, jbNorm,
								jbExp, c, kcNorm, cPower, kcExp, d, dNorms.get(l), dExps.get(l), 0);
					}
				}
			}
		}
		return value;
	}

	/**
	 * VRR (Vertical Recurrence Relation) contribution
	 * 
	 * @param a      Center of primitive Gaussian function a.
	 * @param aNorm  Normalisation coefficients of primitive Gaussian function a.
	 * @param aPower Angular momentum of primitive Gaussian function a.
	 * @param aAlpha Orbital exponent of primitiveGaussian function a.
	 * 
	 * @param b      Center of primitive Gaussian function b.
	 * @param bNorm  Normalisation coefficients of primitive Gaussian function b.
	 * @param bAlpha Orbital exponent of primitiveGaussian function b.
	 * 
	 * @param c      Center of primitive Gaussian function c.
	 * @param cNorm  Normalisation coefficients of primitive Gaussian function c.
	 * @param cPower Angular momentum of primitive Gaussian function c.
	 * @param cAlpha Orbital exponent of primitiveGaussian function c.
	 * 
	 * @param d      Center of primitive Gaussian function d.
	 * @param dNorm  Normalisation coefficients of primitive Gaussian function d.
	 * @param dAlpha Orbital exponent of primitiveGaussian function d.
	 * 
	 * @param m      If equal to 0, ERI is true, else greater than 0, ERI is an
	 *               auxiliary integral.
	 * @return Contribution to Vertical Recurrence Relation.
	 */
	final double vrr(final Vector3D a, final double aNorm, final Power aPower, final double aAlpha, final Vector3D b,
			final double bNorm, final double bAlpha, final Vector3D c, final double cNorm, final Power cPower,
			final double cAlpha, final Vector3D d, final double dNorm, final double dAlpha, final int m) {
		double val;

		final Vector3D p = IntegralsUtil.gaussianProductCenter(aAlpha, a, bAlpha, b);
		final Vector3D q = IntegralsUtil.gaussianProductCenter(cAlpha, c, dAlpha, d);
		final double zeta = aAlpha + bAlpha;
		final double eta = cAlpha + dAlpha;
		final double zetaPlusEta = zeta + eta;
		final double zetaByZetaPlusEta = zeta / zetaPlusEta;
		final double etaByZetaPlusEta = eta / zetaPlusEta;
		final Vector3D w = IntegralsUtil.gaussianProductCenter(zeta, p, eta, q);

		final int la = aPower.getL();
		final int ma = aPower.getM();
		final int na = aPower.getN();
		final int lc = cPower.getL();
		final int mc = cPower.getM();
		final int nc = cPower.getN();

		if (nc > 0) {
			final Power newCPower = new Power(lc, mc, nc - 1);
			val = (q.getZ() - c.getZ())
					* vrr(a, aNorm, aPower, aAlpha, b, bNorm, bAlpha, c, cNorm, newCPower, cAlpha, d, dNorm, dAlpha, m)
					+ (w.getZ() - q.getZ()) * vrr(a, aNorm, aPower, aAlpha, b, bNorm, bAlpha, c, cNorm, newCPower,
							cAlpha, d, dNorm, dAlpha, m + 1);

			if (nc > 1) {
				final Power newCPower1 = new Power(lc, mc, nc - 2);
				val += 0.5 * (nc - 1) / eta
						* (vrr(a, aNorm, aPower, aAlpha, b, bNorm, bAlpha, c, cNorm, newCPower1, cAlpha, d, dNorm,
								dAlpha, m)
								- zetaByZetaPlusEta * vrr(a, aNorm, aPower, aAlpha, b, bNorm, bAlpha, c, cNorm,
										newCPower1, cAlpha, d, dNorm, dAlpha, m + 1));
			}

			if (na > 0) {
				val += 0.5 * na / zetaPlusEta * vrr(a, aNorm, new Power(la, ma, na - 1), aAlpha, b, bNorm, bAlpha, c,
						cNorm, newCPower, cAlpha, d, dNorm, dAlpha, m + 1);
			}
			return val;
		} else if (mc > 0) {
			final Power newCPower = new Power(lc, mc - 1, nc);
			val = (q.getY() - c.getY())
					* vrr(a, aNorm, aPower, aAlpha, b, bNorm, bAlpha, c, cNorm, newCPower, cAlpha, d, dNorm, dAlpha, m)
					+ (w.getY() - q.getY()) * vrr(a, aNorm, aPower, aAlpha, b, bNorm, bAlpha, c, cNorm, newCPower,
							cAlpha, d, dNorm, dAlpha, m + 1);

			if (mc > 1) {
				final Power newCPower1 = new Power(lc, mc - 2, nc);
				val += 0.5 * (mc - 1) / eta
						* (vrr(a, aNorm, aPower, aAlpha, b, bNorm, bAlpha, c, cNorm, newCPower1, cAlpha, d, dNorm,
								dAlpha, m)
								- zetaByZetaPlusEta * vrr(a, aNorm, aPower, aAlpha, b, bNorm, bAlpha, c, cNorm,
										newCPower1, cAlpha, d, dNorm, dAlpha, m + 1));
			}

			if (ma > 0) {
				val += 0.5 * ma / zetaPlusEta * vrr(a, aNorm, new Power(la, ma - 1, na), aAlpha, b, bNorm, bAlpha, c,
						cNorm, newCPower, cAlpha, d, dNorm, dAlpha, m + 1);
			}
			return val;
		} else if (lc > 0) {
			final Power newCPower = new Power(lc - 1, mc, nc);
			val = (q.getX() - c.getX())
					* vrr(a, aNorm, aPower, aAlpha, b, bNorm, bAlpha, c, cNorm, newCPower, cAlpha, d, dNorm, dAlpha, m)
					+ (w.getX() - q.getX()) * vrr(a, aNorm, aPower, aAlpha, b, bNorm, bAlpha, c, cNorm, newCPower,
							cAlpha, d, dNorm, dAlpha, m + 1);

			if (lc > 1) {
				final Power newCPower1 = new Power(lc - 2, mc, nc);
				val += 0.5 * (lc - 1) / eta
						* (vrr(a, aNorm, aPower, aAlpha, b, bNorm, bAlpha, c, cNorm, newCPower1, cAlpha, d, dNorm,
								dAlpha, m)
								- zetaByZetaPlusEta * vrr(a, aNorm, aPower, aAlpha, b, bNorm, bAlpha, c, cNorm,
										newCPower1, cAlpha, d, dNorm, dAlpha, m + 1));
			}

			if (la > 0) {
				val += 0.5 * la / zetaPlusEta * vrr(a, aNorm, new Power(la - 1, ma, na), aAlpha, b, bNorm, bAlpha, c,
						cNorm, newCPower, cAlpha, d, dNorm, dAlpha, m + 1);
			}
			return val;
		} else if (na > 0) {
			final Power newAPower = new Power(la, ma, na - 1);
			val = (p.getZ() - a.getZ())
					* vrr(a, aNorm, newAPower, aAlpha, b, bNorm, bAlpha, c, cNorm, cPower, cAlpha, d, dNorm, dAlpha, m)
					+ (w.getZ() - p.getZ()) * vrr(a, aNorm, newAPower, aAlpha, b, bNorm, bAlpha, c, cNorm, cPower,
							cAlpha, d, dNorm, dAlpha, m + 1);

			if (na > 1) {
				final Power newAPower1 = new Power(la, ma, na - 2);
				val += 0.5 * (na - 1) / zeta
						* (vrr(a, aNorm, newAPower1, aAlpha, b, bNorm, bAlpha, c, cNorm, cPower, cAlpha, d, dNorm,
								dAlpha, m)
								- etaByZetaPlusEta * vrr(a, aNorm, newAPower1, aAlpha, b, bNorm, bAlpha, c, cNorm,
										cPower, cAlpha, d, dNorm, dAlpha, m + 1));
			}

			return val;
		} else if (ma > 0) {
			final Power newAPower = new Power(la, ma - 1, na);
			val = (p.getY() - a.getY())
					* vrr(a, aNorm, newAPower, aAlpha, b, bNorm, bAlpha, c, cNorm, cPower, cAlpha, d, dNorm, dAlpha, m)
					+ (w.getY() - p.getY()) * vrr(a, aNorm, newAPower, aAlpha, b, bNorm, bAlpha, c, cNorm, cPower,
							cAlpha, d, dNorm, dAlpha, m + 1);

			if (ma > 1) {
				final Power newAPower1 = new Power(la, ma - 2, na);
				val += 0.5 * (ma - 1) / zeta
						* (vrr(a, aNorm, newAPower1, aAlpha, b, bNorm, bAlpha, c, cNorm, cPower, cAlpha, d, dNorm,
								dAlpha, m)
								- etaByZetaPlusEta * vrr(a, aNorm, newAPower1, aAlpha, b, bNorm, bAlpha, c, cNorm,
										cPower, cAlpha, d, dNorm, dAlpha, m + 1));
			}

			return val;
		} else if (la > 0) {
			final Power newAPower = new Power(la - 1, ma, na);
			val = (p.getX() - a.getX())
					* vrr(a, aNorm, newAPower, aAlpha, b, bNorm, bAlpha, c, cNorm, cPower, cAlpha, d, dNorm, dAlpha, m)
					+ (w.getX() - p.getX()) * vrr(a, aNorm, newAPower, aAlpha, b, bNorm, bAlpha, c, cNorm, cPower,
							cAlpha, d, dNorm, dAlpha, m + 1);

			if (la > 1) {
				final Power newAPower1 = new Power(la - 2, ma, na);
				val += 0.5 * (la - 1) / zeta
						* (vrr(a, aNorm, newAPower1, aAlpha, b, bNorm, bAlpha, c, cNorm, cPower, cAlpha, d, dNorm,
								dAlpha, m)
								- etaByZetaPlusEta * vrr(a, aNorm, newAPower1, aAlpha, b, bNorm, bAlpha, c, cNorm,
										cPower, cAlpha, d, dNorm, dAlpha, m + 1));
			}

			return val;
		}

		final double rab2 = a.distanceSq(b);
		final double Kab = sqrt2PI / zeta * FastMath.exp(-aAlpha * bAlpha / zeta * rab2);
		final double rcd2 = c.distanceSq(d);
		final double Kcd = sqrt2PI / eta * FastMath.exp(-cAlpha * dAlpha / eta * rcd2);
		final double rpq2 = p.distanceSq(q);
		final double T = zeta * eta / zetaPlusEta * rpq2;

		return aNorm * bNorm * cNorm * dNorm * Kab * Kcd / FastMath.sqrt(zetaPlusEta)
				* IntegralsUtil.computeFGamma(m, T);
	}

	@Override
	public final double coulomb(final ContractedGaussian a, final ContractedGaussian b, final ContractedGaussian c,
			final ContractedGaussian d, final Density density, final RealMatrix jMat, final RealMatrix kMat) {
		throw new UnsupportedOperationException("Not supported yet.");
	}
}
