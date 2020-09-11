package name.mjw.jquante.math.qm;

import org.hipparchus.linear.ArrayRealVector;
import org.hipparchus.linear.RealVector;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

import name.mjw.jquante.math.qm.event.SCFEvent;
import name.mjw.jquante.math.qm.integral.IntegralsUtil;
import name.mjw.jquante.molecule.Molecule;

/**
 * Moller Plesset Perturbation Theory as an extension of HF method. One of the
 * first and the most popular post HF methods. <br>
 * This class only implements second order perturbation (MP2), but may be
 * modified later to provide higher order perturbation.
 * 
 * @author V.Ganesh
 * @version 2.0 (Part of MeTA v2.0)
 */
public final class RestrictedMollerPlessetSCFMethod extends RestrictedHartreeFockMethod {

	private static final Logger LOG = LogManager.getLogger(RestrictedMollerPlessetSCFMethod.class);

	/**
	 * Holds value of property mpLevel.
	 */
	private int mpLevel;

	/**
	 * transformed MO integrals
	 */
	protected double[] moInts;

	/**
	 * Creates a new instance of MollerPlessetSCFMethod.
	 *
	 * @param molecule
	 *            the Molecule object.
	 * @param oneEI
	 *            the 1E integral driver.
	 * @param twoEI
	 *            the 2E integral driver.
	 */
	public RestrictedMollerPlessetSCFMethod(final Molecule molecule, final OneElectronIntegrals oneEI,
			final TwoElectronIntegrals twoEI) {
		super(molecule, oneEI, twoEI);

		mpLevel = 2; // second level correction (MP2) is default
	}

	/**
	 * Perform the SCF optimisation of the molecular wave function until the energy
	 * converges.
	 */
	@Override
	public void scf() {
		// first do the HF-SCF procedure
		super.scf();

		// transform AO -> MO ints
		transformAOIntsToMOInts();

		// now compute MP2 energy and then compute the total energy
		final int noOfBasisFunctions = mos.getRowDimension();
		final int noOfElectrons = molecule.getNumberOfElectrons();
		final int noOfOccupancies = noOfElectrons / 2;
		final int noOfVirtualOrbitals = noOfBasisFunctions - noOfOccupancies;
		final int noOfUnOccupied = noOfOccupancies + noOfVirtualOrbitals;

		double mp2Energy = 0.0;
		double arbs;
		double asbr;
		final double[] orbE = mos.getOrbitalEnergies();

		for (int a = 0; a < noOfOccupancies; a++) {
			for (int b = 0; b < noOfOccupancies; b++) {
				for (int r = noOfOccupancies; r < noOfUnOccupied; r++) {
					for (int s = noOfOccupancies; s < noOfUnOccupied; s++) {
						arbs = moInts[IntegralsUtil.ijkl2intindex(a, r, b, s)];
						asbr = moInts[IntegralsUtil.ijkl2intindex(a, s, b, r)];

						mp2Energy += arbs * (2.0 * arbs - asbr) / (orbE[a] + orbE[b] - orbE[r] - orbE[s]);
					}
				}
			}
		}

		LOG.debug(mp2Energy);

		energy += mp2Energy;

		scfEvent.setType(SCFEvent.CONVERGED_EVENT);
		scfEvent.setCurrentIteration(scfIteration);
		scfEvent.setCurrentEnergy(energy);
		fireSCFEventListenerScfEventOccured(scfEvent);
	}

	/**
	 * O(N^5) 4-index transformation of the two-electron integrals. Only transform
	 * the ones needed for MP2, which reduces the scaling to O(nN^4), where n are
	 * the occs (&lt;&lt;N).
	 */
	protected void transformAOIntsToMOInts() {
		// Start with (mu,nu|sigma,eta)
		// Unpack AOints and transform sigma -> b
		int mu;
		int nu;
		int sigma;
		int eta;
		int a;
		int b;
		int i;
		int j;

		final int noOfBasisFunctions = mos.getRowDimension();
		final int noOfMOS = mos.getColumnDimension();
		final int noOfElectrons = molecule.getNumberOfElectrons();
		final int noOfOccupancies = noOfElectrons / 2;

		final RealVector tempVector = new ArrayRealVector(noOfBasisFunctions);
		final double[] tempvec = tempVector.toArray();
		final double[] aoints = twoEI.getTwoEIntegrals();

		double[][][][] temp = new double[noOfBasisFunctions][noOfBasisFunctions][noOfOccupancies][noOfBasisFunctions];

		for (mu = 0; mu < noOfBasisFunctions; mu++) {
			for (nu = 0; nu < noOfBasisFunctions; nu++) {
				for (eta = 0; eta < noOfBasisFunctions; eta++) {
					for (b = 0; b < noOfOccupancies; b++) {
						for (sigma = 0; sigma < noOfBasisFunctions; sigma++) {
							tempvec[sigma] = aoints[IntegralsUtil.ijkl2intindex(mu, nu, sigma, eta)];
						} // end sigma

						temp[mu][nu][b][eta] = mos.getRowVector(b).dotProduct(tempVector);
					} // end b
				} // end eta
			} // end nu
		} // end mu

		final double[][][][] temp2 = new double[noOfOccupancies][noOfBasisFunctions][noOfOccupancies][noOfBasisFunctions];

		for (nu = 0; nu < noOfBasisFunctions; nu++) {
			for (eta = 0; eta < noOfBasisFunctions; eta++) {
				for (b = 0; b < noOfOccupancies; b++) {
					for (a = 0; a < noOfOccupancies; a++) {
						for (i = 0; i < noOfBasisFunctions; i++) {
							tempvec[i] = temp[i][nu][b][eta];
						} // end i

						temp2[a][nu][b][eta] = mos.getRowVector(a).dotProduct(tempVector);
					} // end a
				} // end b
			} // end eta
		} // end nu

		temp = new double[noOfBasisFunctions][noOfBasisFunctions][noOfOccupancies][noOfBasisFunctions];

		for (a = 0; a < noOfOccupancies; a++) {
			for (nu = 0; nu < noOfBasisFunctions; nu++) {
				for (b = 0; b < noOfOccupancies; b++) {
					for (j = 0; j < noOfMOS; j++) {
						for (i = 0; i < noOfBasisFunctions; i++) {
							tempvec[i] = temp2[a][nu][b][i];
						} // end i

						temp[a][nu][b][j] = mos.getRowVector(j).dotProduct(tempVector);
					} // end a
				} // end b
			} // end eta
		} // end nu

		moInts = new double[aoints.length];

		// Transform mu -> i and repack integrals
		for (a = 0; a < noOfOccupancies; a++) {
			for (j = 0; j < noOfMOS; j++) {
				for (b = 0; b < noOfOccupancies; b++) {
					for (i = 0; i < noOfMOS; i++) {
						for (sigma = 0; sigma < noOfBasisFunctions; sigma++) {
							tempvec[sigma] = temp[a][sigma][b][j];
						} // end sigma

						moInts[IntegralsUtil.ijkl2intindex(a, i, b, j)] = mos.getRowVector(i).dotProduct(tempVector);
					} // end i
				} // end b
			} // end j
		} // end a

	}

	/**
	 * Getter for property mpLevel.
	 * 
	 * @return Value of property mpLevel.
	 */
	public int getMpLevel() {
		return this.mpLevel;
	}

	/**
	 * Setter for property mpLevel.
	 * 
	 * @param mpLevel New value of property mpLevel.
	 */
	public void setMpLevel(final int mpLevel) {
		this.mpLevel = mpLevel;
	}

}
