package name.mjw.jquante.math.qm;

import java.lang.ref.WeakReference;

import name.mjw.jquante.molecule.Molecule;

/**
 * Factory of SCF methods. <br>
 * Follows a singleton pattern.
 * 
 * @author V.Ganesh
 * @version 2.0 (Part of MeTA v2.0)
 */
public final class SCFMethodFactory {

	private static WeakReference<SCFMethodFactory> _scfMethodFactory = null;

	/** Creates a new instance of SCFMethodFactory */
	private SCFMethodFactory() {
	}

	/**
	 * Get an instance (and the only one) of SCFMethodFactory
	 * 
	 * @return SCFMethodFactory instance
	 */
	public static SCFMethodFactory getInstance() {
		if (_scfMethodFactory == null) {
			_scfMethodFactory = new WeakReference<>(new SCFMethodFactory());
		}

		SCFMethodFactory scfMethodFactory = _scfMethodFactory.get();

		if (scfMethodFactory == null) {
			scfMethodFactory = new SCFMethodFactory();
			_scfMethodFactory = new WeakReference<>(scfMethodFactory);
		}

		return scfMethodFactory;
	}

	/**
	 * Return an appropriate class appropriate <code>SCFType</code>
	 * 
	 * @param molecule
	 *            the molecule for which calculations are to be performed
	 * @param oneEI
	 *            1E integrals for this molecule at appropriate basis
	 * @param twoEI
	 *            2E integrals for this molecule at appropriate basis
	 * @param type
	 *            the SCFType instance
	 * 
	 * @return the the Self Consistent Field (SCF) method.
	 */
	public SCFMethod getSCFMethod(final Molecule molecule, final OneElectronIntegrals oneEI,
			final TwoElectronIntegrals twoEI, final SCFType type) {
		if (type.equals(SCFType.HARTREE_FOCK)) {
			return new RestrictedHartreeFockMethod(molecule, oneEI, twoEI);
		} else if (type.equals(SCFType.HARTREE_FOCK_DIRECT)) {
			return new RestrictedHartreeFockMethod(molecule, oneEI, twoEI, type);
		} else if (type.equals(SCFType.MOLLER_PLESSET)) {
			return new RestrictedMollerPlessetSCFMethod(molecule, oneEI, twoEI);
		} else if (type.equals(SCFType.UNRESTRICTED_HARTREE_FOCK_DIRECT)) {
			return new UnrestrictedHartreeFockMethod(molecule, oneEI, twoEI);
		} else {
			throw new UnsupportedOperationException("The type : " + type + " has no known implementing class!");
		}
	}
}