package name.mjw.jquante.math.qm.integral;

import org.hipparchus.geometry.euclidean.threed.Vector3D;

import name.mjw.jquante.math.qm.basis.ContractedGaussian;
import name.mjw.jquante.math.qm.basis.Power;
import name.mjw.jquante.math.qm.basis.PrimitiveGaussian;

/**
 * Integral package dispatch class for MeTA Studio. <br>
 * 
 * The purpose of this class has substantially changed from the earlier code, in
 * that it no longer provides an implementation, but dispatches it to the
 * implementation classes that are actually handling the computation work. This
 * in turn gives an easy way to dynamically change the method used to actually
 * evaluate the integrals.
 * 
 * The static initialisation is done as a way to have good performance as well
 * as allow flexibility of using multiple integral packages.
 * 
 * @author V.Ganesh
 * @version 2.0 (Part of MeTA v2.0)
 */
public final class Integrals {

	private static NuclearTerm nuclearTerm;
	private static OneElectronTerm oneElectronTerm;
	private static TwoElectronTerm twoElectronTerm;

	/**
	 * No instantiation possible
	 */
	private Integrals() {
	}

	static {
		nuclearTerm = IntegralsPackageFactory.getInstance().getNuclearTerm();
		oneElectronTerm = IntegralsPackageFactory.getInstance()
				.getOneElectronTerm();
		twoElectronTerm = IntegralsPackageFactory.getInstance()
				.getTwoElectronTerm();
	}

	/**
	 * Re-initialize package references in case there are any changes..
	 */
	public static void reInitPackageReference() {
		nuclearTerm = IntegralsPackageFactory.getInstance().getNuclearTerm();
		oneElectronTerm = IntegralsPackageFactory.getInstance()
				.getOneElectronTerm();
		twoElectronTerm = IntegralsPackageFactory.getInstance()
				.getTwoElectronTerm();
	}

	/**
	 * Overlap integral taken form <i> Taken from THO eq. 2.12 </i> <br>
	 * 
	 * @param alpha1
	 *            the coefficient of primitive Gaussian a.
	 * @param power1
	 *            the orbital powers of primitive Gaussian a.
	 * @param a
	 *            the location of primitive Gaussian a.
	 * 
	 * @param alpha2
	 *            the coefficient of primitive Gaussian b.
	 * @param power2
	 *            the orbital powers of primitive Gaussian b.
	 * @param b
	 *            the location of primitive Gaussian b.
	 * @return the Overlap integral
	 */
	public static double overlap(double alpha1, Power power1, Vector3D a,
			double alpha2, Power power2, Vector3D b) {
		return oneElectronTerm.overlap(alpha1, power1, a, alpha2, power2, b);
	}

	/**
	 * The Kinetic Energy (KE) integral.
	 * 
	 * <i> Taken from THO eq. 2.13 </i>
	 * 
	 * @param alpha1
	 *            the coefficient of primitive Gaussian a.
	 * @param power1
	 *            the orbital powers of primitive Gaussian a.
	 * @param a
	 *            the location of primitive Gaussian a.
	 * 
	 * @param alpha2
	 *            the coefficient of primitive Gaussian b.
	 * @param power2
	 *            the orbital powers of primitive Gaussian b.
	 * @param b
	 *            the location of primitive Gaussian b.
	 * @return the Kinetic Energy integral
	 */
	public static double kinetic(double alpha1, Power power1, Vector3D a,
			double alpha2, Power power2, Vector3D b) {
		return oneElectronTerm.kinetic(alpha1, power1, a, alpha2, power2, b);
	}

	/**
	 * The nuclear attraction integral.
	 * 
	 * <i> Taken from THO eq. 2.15 </i>
	 * 
	 * @param a
	 *            the coefficient of primitive Gaussian a.
	 * @param norm1
	 *            the normalization factor of primitive Gaussian a.
	 * @param power1
	 *            the orbital powers of primitive Gaussian a.
	 * @param alpha1
	 *            the coefficient of primitive Gaussian a.
	 * @param b
	 *            the coefficient of primitive Gaussian b.
	 * @param norm2
	 *            the normalization factor of primitive Gaussian b.
	 * @param power2
	 *            the orbital powers of primitive Gaussian b.
	 * @param alpha2
	 *            the coefficient of primitive Gaussian b.
	 * @param c
	 *            the location of nuclear centre.
	 * @return the nuclear attraction integral.
	 */
	public static double nuclearAttraction(Vector3D a, double norm1,
			Power power1, double alpha1, Vector3D b, double norm2, Power power2,
			double alpha2, Vector3D c) {
		return nuclearTerm.nuclearAttraction(a, norm1, power1, alpha1, b,
				norm2, power2, alpha2, c);

	}

	/**
	 * The nuclear attraction gradient term.
	 * 
	 * @param a
	 *            the coefficient of primitive Gaussian a.
	 * @param power1
	 *            the orbital powers of primitive Gaussian a.
	 * @param alpha1
	 *            the coefficient of primitive Gaussian a.
	 * @param b
	 *            the coefficient of primitive Gaussian b.
	 * @param power2
	 *            the orbital powers of primitive Gaussian b.
	 * @param alpha2
	 *            the coefficient of primitive Gaussian b.
	 * @param c
	 *            the location of nuclear centre.
	 * @return the nuclear attraction gradient.
	 */
	public static Vector3D nuclearAttractionGradient(Vector3D a, Power power1,
			double alpha1, Vector3D b, Power power2, double alpha2, Vector3D c) {
		return nuclearTerm.nuclearAttractionGradient(a, power1, alpha1, b,
				power2, alpha2, c);

	}

	/**
	 * 2E coulomb interactions between 4 contracted Gaussians
	 * 
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
	public static double coulomb(ContractedGaussian a, ContractedGaussian b,
			ContractedGaussian c, ContractedGaussian d) {
		return twoElectronTerm.coulomb(a, b, c, d);
	}

	/**
	 * 2E coulomb interactions between 4 primitive Gaussians
	 * 
	 * @param a
	 *            Primitive Gaussian function a.
	 * @param b
	 *            Primitive Gaussian function b.
	 * @param c
	 *            Primitive Gaussian function c.
	 * @param d
	 *            Primitive Gaussian function d.
	 * @return Two-electron integral.
	 */
	public static double coulomb(PrimitiveGaussian a, PrimitiveGaussian b,
			PrimitiveGaussian c, PrimitiveGaussian d) {
		return twoElectronTerm.coulombRepulsion(a.getOrigin(),
				a.getNormalization(), a.getPowers(), a.getExponent(),
				b.getOrigin(), b.getNormalization(), b.getPowers(),
				b.getExponent(), c.getOrigin(), c.getNormalization(),
				c.getPowers(), c.getExponent(), d.getOrigin(),
				d.getNormalization(), d.getPowers(), d.getExponent());
	}
}