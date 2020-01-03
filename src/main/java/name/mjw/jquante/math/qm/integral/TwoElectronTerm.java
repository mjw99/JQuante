package name.mjw.jquante.math.qm.integral;

import org.hipparchus.geometry.euclidean.threed.Vector3D;
import org.hipparchus.linear.RealMatrix;

import name.mjw.jquante.math.qm.Density;
import name.mjw.jquante.math.qm.basis.ContractedGaussian;
import name.mjw.jquante.math.qm.basis.Power;

/**
 * Top level interface for evaluating a 2E-integral term.
 * 
 * @author V.Ganesh
 * @version 2.0 (Part of MeTA v2.0)
 */
public interface TwoElectronTerm extends IntegralsPackage {

	/**
	 * 2E coulomb interactions between four contracted Gaussian functions.
	 * 
	 * @param a
	 *            Contracted Gaussian function a.
	 * @param b
	 *            Contracted Gaussian function b.
	 * @param c
	 *            Contracted Gaussian function c.
	 * @param d
	 *            Contracted Gaussian function d.
	 * @return Two-electron integral.
	 */
	public double coulomb(ContractedGaussian a, ContractedGaussian b,
			ContractedGaussian c, ContractedGaussian d);

	/**
	 * 2E coulomb interactions between four contracted Gaussians.
	 * 
	 * @param a
	 *            Contracted Gaussian function a.
	 * @param b
	 *            Contracted Gaussian function b.
	 * @param c
	 *            Contracted Gaussian function c.
	 * @param d
	 *            Contracted Gaussian function d.
	 * @param density
	 *            Density matrix.
	 * @param jMat
	 *            J matrix.
	 * @param kMat
	 *            K matrix.
	 * @return Two-electron integral.
	 */
	public double coulomb(ContractedGaussian a, ContractedGaussian b,
			ContractedGaussian c, ContractedGaussian d, Density density,
			RealMatrix jMat, RealMatrix kMat);

	/**
	 * The coulomb repulsion term between four centered Gaussians.
	 *
	 * @param a
	 *            Center of contracted Gaussian function a.
	 * @param aNorm
	 *            Normalisation factor of the contracted Gaussian function a.
	 * @param aPower
	 *            Angular momentum of Gaussian function a.
	 * @param aAlpha
	 *            Orbital exponent of Gaussian function a.
	 * @param b
	 *            Center of contracted Gaussian function b.
	 * @param bNorm
	 *            Normalisation factor of the contracted Gaussian function b.
	 * @param bPower
	 *            Angular momentum of Gaussian function b.
	 * @param bAlpha
	 *            Orbital exponent of Gaussian function b.
	 * @param c
	 *            Center of contracted Gaussian function c.
	 * @param cNorm
	 *            Normalisation factor of the contracted Gaussian function c.
	 * @param cPower
	 *            Angular momentum of Gaussian function c.
	 * @param cAlpha
	 *            Orbital exponent of Gaussian function c.
	 * @param d
	 *            Center of contracted Gaussian function d.
	 * @param dNorm
	 *            Normalisation factor of the contracted Gaussian function d.
	 * @param dPower
	 *            Angular momentum of Gaussian function d.
	 * @param dAlpha
	 *            Orbital exponent of Gaussian function d.
	 * @return Two-electron integral.
	 */
	public abstract double coulombRepulsion(Vector3D a, double aNorm,
			Power aPower, double aAlpha, Vector3D b, double bNorm, Power bPower,
			double bAlpha, Vector3D c, double cNorm, Power cPower,
			double cAlpha, Vector3D d, double dNorm, Power dPower, double dAlpha);

}
