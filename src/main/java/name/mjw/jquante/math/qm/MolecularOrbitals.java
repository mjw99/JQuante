package name.mjw.jquante.math.qm;

import org.hipparchus.linear.Array2DRowRealMatrix;
import org.hipparchus.linear.ArrayRealVector;
import org.hipparchus.linear.EigenDecomposition;
import org.hipparchus.linear.MatrixUtils;
import org.hipparchus.linear.RealMatrix;
import org.hipparchus.linear.RealVector;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

import name.mjw.jquante.math.MathUtil;

/**
 * Represents the Molecular orbitals as a coefficient matrix and the
 * corresponding eigenvalues representing the orbital energies.
 * 
 * @author V.Ganesh
 * @version 2.0 (Part of MeTA v2.0)
 */
public final class MolecularOrbitals extends Array2DRowRealMatrix {

	private static final long serialVersionUID = -3159550917810939744L;

	private static final Logger LOG = LogManager.getLogger(MolecularOrbitals.class);

	/**
	 * Creates a new instance of square (NxN) Matrix
	 * 
	 * @param n
	 *            the dimension
	 */
	public MolecularOrbitals(int n) {
		super(n, n);
	}

	/**
	 * Get the coefficient Matrix for this MolecularOrbitals
	 * 
	 * @return the coefficient matrix
	 */
	public double[][] getCoefficients() {
		return getData();
	}

	protected double[] orbitalEnergies;

	/**
	 * Get the value of orbitalEnergies
	 * 
	 * @return the value of orbitalEnergies
	 */
	public double[] getOrbitalEnergies() {
		return orbitalEnergies;
	}

	/**
	 * Compute the MO coefficients and the orbital energies
	 * 
	 * @param hCore
	 *            the HCore matrix
	 * @param overlap
	 *            the Overlap matrix
	 */
	public void compute(HCore hCore, Overlap overlap) {
		compute((RealMatrix) hCore, overlap);
	}

	/**
	 * Compute the MO coefficients and the orbital energies
	 * 
	 * @param fock
	 *            the Fock matrix
	 * @param overlap
	 *            the Overlap matrix
	 */
	public void compute(Fock fock, Overlap overlap) {
		compute((RealMatrix) fock, overlap);
	}

	/** The actual computation is irrelevant of the type of matrix */
	private void compute(RealMatrix theMat, Overlap overlap) {
		LOG.debug("");
		RealMatrix x = overlap.getSHalf();
		LOG.debug("x: " + x);
		RealMatrix a = x.multiply(theMat).multiply(x.transpose());
		LOG.debug("a: " + a);

		EigenDecomposition eig = new EigenDecomposition(a);

		SortedEigenDecomposition sortedEig = new SortedEigenDecomposition(eig);

		orbitalEnergies = sortedEig.getRealEigenvalues();

		this.setSubMatrix(sortedEig.getVT().multiply(x).getData(), 0, 0);

		LOG.debug("MO::values :" + this);
	}

	@Override
	public String toString() {
		return MathUtil.matrixToString(this);
	}

}

/**
 * Reorders commons-math's EigenDecomposition eigensystem, sorting by the
 * smallest eigenvalue first.
 *
 * @author mjw
 *
 */
class SortedEigenDecomposition {

	double[] realEigenvalues;
	RealVector[] eigenvectors;

	int n;

	public SortedEigenDecomposition(EigenDecomposition eig) {
		n = eig.getRealEigenvalues().length;

		realEigenvalues = eig.getRealEigenvalues();
		eigenvectors = new ArrayRealVector[n];

		for (int k = 0; k < n; ++k) {
			eigenvectors[k] = eig.getEigenvector(k);
		}

		RealVector tmpvector;
		for (int i = 0; i < n; i++) {

			int k = i;
			double p = realEigenvalues[i];

			for (int j = i + 1; j < n; j++) {
				if (realEigenvalues[j] < p) {
					k = j;
					p = realEigenvalues[j];
				}
			}

			if (k != i) {
				realEigenvalues[k] = realEigenvalues[i];
				realEigenvalues[i] = p;

				tmpvector = eigenvectors[i];
				eigenvectors[i] = eigenvectors[k];
				eigenvectors[k] = tmpvector;
			}
		}

	}

	RealMatrix getVT() {
		RealMatrix vt = MatrixUtils.createRealMatrix(n, n);
		for (int k = 0; k < n; ++k) {
			vt.setRowVector(k, eigenvectors[k]);
		}

		return vt;
	}

	double[] getRealEigenvalues() {
		return realEigenvalues;
	}

}
