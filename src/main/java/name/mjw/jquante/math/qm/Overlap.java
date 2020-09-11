package name.mjw.jquante.math.qm;

import java.util.ArrayList;
import java.util.List;

import org.hipparchus.geometry.euclidean.threed.Vector3D;
import org.hipparchus.linear.Array2DRowRealMatrix;
import org.hipparchus.linear.EigenDecomposition;
import org.hipparchus.linear.MatrixUtils;
import org.hipparchus.linear.RealMatrix;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

import name.mjw.jquante.math.MathUtil;
import name.mjw.jquante.math.qm.basis.BasisSetLibrary;
import name.mjw.jquante.math.qm.basis.ContractedGaussian;
import net.jafama.FastMath;

/**
 * Represents the Overlap matrix (S).
 * 
 * @author V.Ganesh
 * @version 2.0 (Part of MeTA v2.0)
 */
public final class Overlap extends Array2DRowRealMatrix {

	private static final Logger LOG = LogManager.getLogger(Overlap.class);

	private static final long serialVersionUID = 5209272241805533800L;

	/**
	 * Creates a new instance of square (NxN) Matrix
	 * 
	 * @param n
	 *            the dimension
	 */
	public Overlap(final int n) {
		super(n, n);
	}

	private transient RealMatrix sHalf = null;

	/**
	 * Get the S^1/2 matrix
	 * 
	 * Symmetric orthogonalization of the real symmetric matrix X (this). This is
	 * given by <code>U'(1/sqrt(lambda))U</code>, where lambda, U are the
	 * eigenvalues, vectors.
	 *
	 *
	 * @return return the symmetric orthogonalization matrix (S half)
	 */
	public RealMatrix getSHalf() {
		if (sHalf == null) {

			LOG.debug("Overlap::this " + this);
			final EigenDecomposition eig = new EigenDecomposition(this);

			final double[] eigenValues = eig.getRealEigenvalues();
			final RealMatrix eigenVectors = eig.getVT();

			LOG.trace("eigenVectors " + eigenVectors);

			this.sHalf = MatrixUtils.createRealIdentityMatrix(this.getRowDimension());

			for (int i = 0; i < this.getRowDimension(); i++) {
				sHalf.setEntry(i, i, (sHalf.getEntry(i, i) / FastMath.sqrt(eigenValues[i])));
			}

			this.sHalf = eigenVectors.transpose().multiply(sHalf).multiply(eigenVectors);
		}

		return this.sHalf;
	}

	private SCFMethod scfMethod;
	private int atomIndex;

	/**
	 * Compute Overlap partial derivative for an atom index.
	 * 
	 * @param atomIndex the atom index with respect to which the derivative are to
	 *                  be evaluated
	 * @param scfMethod the reference to the SCFMethod
	 * @return three element array of Overlap elements representing partial
	 *         derivatives with respect to x, y and z of atom position
	 */
	public List<Overlap> computeDerivative(final int atomIndex, final SCFMethod scfMethod) {
		this.scfMethod = scfMethod;
		this.atomIndex = atomIndex;

		final ArrayList<Overlap> dOverlap = new ArrayList<>(3);

		final int noOfBasisFunctions = this.getRowDimension();
		final Overlap dOverlapDx = new Overlap(noOfBasisFunctions);
		final Overlap dOverlapDy = new Overlap(noOfBasisFunctions);
		final Overlap dOverlapDz = new Overlap(noOfBasisFunctions);

		final double[][] hdx = dOverlapDx.getData();
		final double[][] hdy = dOverlapDy.getData();
		final double[][] hdz = dOverlapDz.getData();

		for (int i = 0; i < noOfBasisFunctions; i++) {
			for (int j = 0; j < noOfBasisFunctions; j++) {
				final Vector3D dOvrEle = computeOverlapDerElement(i, j);

				hdx[i][j] = dOvrEle.getX();
				hdy[i][j] = dOvrEle.getY();
				hdz[i][j] = dOvrEle.getZ();
			}
		}

		dOverlap.add(dOverlapDx);
		dOverlap.add(dOverlapDy);
		dOverlap.add(dOverlapDz);

		return dOverlap;
	}

	/**
	 * Compute one of the Overlap derivative elements, with respect to an atomIndex
	 */
	private Vector3D computeOverlapDerElement(final int i, final int j) {
		final BasisSetLibrary bsl = scfMethod.getOneEI().getBasisSetLibrary();
		final ContractedGaussian cgi = bsl.getBasisFunctions().get(i);
		final ContractedGaussian cgj = bsl.getBasisFunctions().get(j);

		return cgi.overlapDerivative(atomIndex, cgj);
	}

	@Override
	public String toString() {
		return MathUtil.matrixToString(this);
	}
}
