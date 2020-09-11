package name.mjw.jquante.math.qm;

import java.util.List;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

import name.mjw.jquante.config.impl.AtomInfo;
import name.mjw.jquante.math.qm.basis.BasisSetLibrary;
import name.mjw.jquante.math.qm.basis.ContractedGaussian;
import name.mjw.jquante.molecule.Molecule;

/**
 * The 1E integral (overlap S matrix) and 1E hCore matrices driver.
 * 
 * @author V.Ganesh
 * @version 2.0 (Part of MeTA v2.0)
 */
public final class OneElectronIntegrals {

	private static final Logger LOG = LogManager.getLogger(OneElectronIntegrals.class);

	/**
	 * The overlap S matrix.
	 */
	private Overlap overlap;

	/**
	 * The core Hamiltonian matrix contains integrals that represent the kinetic
	 * energy of an electron (T) and electron-nuclear potential energy (V).
	 */
	private HCore hCore;

	private final BasisSetLibrary basisSetLibrary;
	private final Molecule molecule;

	/**
	 * Creates a new instance of OneElectronIntegrals
	 * 
	 * @param basisSetLibrary the basis functions to be used
	 * @param mol             the Molecule object, of whose 1E integrals are to be
	 *                        evaluated
	 */
	public OneElectronIntegrals(final BasisSetLibrary basisSetLibrary, final Molecule mol) {
		this.basisSetLibrary = basisSetLibrary;
		this.molecule = mol;

		compute1E();
	}

	/**
	 * Get BasisFunctions associated with this 1E evaluation
	 * 
	 * @return instance of basisFunctions
	 */
	public BasisSetLibrary getBasisSetLibrary() {
		return basisSetLibrary;
	}

	/**
	 * compute the 1E integrals, form S matrix and hCore
	 */
	protected void compute1E() {
		computeOverlap();
		computeHcore();
	}

	private void computeOverlap() {
		final List<ContractedGaussian> bfs = basisSetLibrary.getBasisFunctions();
		final int noOfBasisFunctions = bfs.size();

		// Create the S matrix
		this.overlap = new Overlap(noOfBasisFunctions);

		// Populate the S matrix
		for (int i = 0; i < noOfBasisFunctions; i++) {
			final ContractedGaussian bfi = bfs.get(i);

			for (int j = 0; j < noOfBasisFunctions; j++) {
				final ContractedGaussian bfj = bfs.get(j);

				overlap.setEntry(i, j, bfi.overlap(bfj)); // the overlap matrix
			}
		}
	}

	private void computeHcore() {
		final List<ContractedGaussian> bfs = basisSetLibrary.getBasisFunctions();
		final int noOfBasisFunctions = bfs.size();

		this.hCore = new HCore(noOfBasisFunctions);

		// read in the atomic numbers
		final int[] atomicNumbers = new int[molecule.getNumberOfAtoms()];
		final AtomInfo ai = AtomInfo.getInstance();

		for (int atomIndex = 0; atomIndex < atomicNumbers.length; atomIndex++) {
			atomicNumbers[atomIndex] = ai.getAtomicNumber(molecule.getAtom(atomIndex).getSymbol());
		}

		// Populate the hCore matrix
		for (int i = 0; i < noOfBasisFunctions; i++) {
			final ContractedGaussian bfi = bfs.get(i);

			for (int j = 0; j < noOfBasisFunctions; j++) {
				final ContractedGaussian bfj = bfs.get(j);

				hCore.setEntry(i, j, bfi.kinetic(bfj)); // KE matrix elements

				for (int k = 0; k < atomicNumbers.length; k++) {
					hCore.setEntry(i, j, (hCore.getEntry(i, j)
							+ atomicNumbers[k] * bfi.nuclear(bfj, molecule.getAtom(k).getAtomCenterInAU())));
				}
			}

		}
	}

	/**
	 * Getter for property overlap.
	 * 
	 * @return Value of property overlap.
	 */
	public Overlap getOverlap() {
		return this.overlap;
	}

	/**
	 * Getter for property hCore.
	 * 
	 * @return Value of property hCore.
	 */
	public HCore getHCore() {
		return this.hCore;
	}

}
