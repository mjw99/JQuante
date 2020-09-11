package name.mjw.jquante.math.qm;

import org.hipparchus.linear.Array2DRowRealMatrix;
import org.hipparchus.linear.RealMatrix;

import name.mjw.jquante.math.MathUtil;

/**
 * Represents the Density P matrix or the charge order bond density matrix (DM).
 * 
 * @author V.Ganesh
 * @version 2.0 (Part of MeTA v2.0)
 */
public class Density extends Array2DRowRealMatrix {

	private static final long serialVersionUID = 4939105356956860039L;

	/**
	 * Creates a new instance of square (NxN) Matrix
	 * 
	 * @param n
	 *            the dimension
	 */
	public Density(final int n) {
		super(n, n);
	}

	public Density(final double[][] data) {
		super(data);
	}

	/**
	 * Compute the density matrix. It is computed using the number of occupied
	 * orbitals and the molecular orbitals coefficient matrix.
	 * 
	 * @param scfMethod       The SCF method for which the density is to be formed
	 * @param guessInitialDM  Do a guess for DM?
	 * @param densityGuesser  Is there any special way to make the guess?
	 * @param noOfOccupiedMOs Number of molecular orbitals occupied by electrons
	 * @param mos             the MolecularOrbitals, needed to compute the DM
	 */
	public void compute(final SCFMethod scfMethod, final boolean guessInitialDM, final DensityGuesser densityGuesser,
			final int noOfOccupiedMOs, final MolecularOrbitals mos) {
		if (guessInitialDM && densityGuesser != null) {
			this.setSubMatrix(densityGuesser.guessDM(scfMethod).getData(), 0, 0);
			return;
		}

		// else construct it from the MOs .. C*C'
		final RealMatrix dVector = new Array2DRowRealMatrix(noOfOccupiedMOs, mos.getRowDimension());

		final double[][] c = mos.getData();

		for (int i = 0; i < noOfOccupiedMOs; i++) {
			for (int j = 0; j < c.length; j++) {
				dVector.setEntry(i, j, mos.getEntry(i, j));
			}
		}

		this.setSubMatrix(dVector.transpose().multiply(dVector).getData(), 0, 0);

	}

	@Override
	public String toString() {
		return MathUtil.matrixToString(this);
	}
}
