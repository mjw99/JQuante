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
	public Overlap(int n) {
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
			EigenDecomposition eig = new EigenDecomposition(this);

			double[] eigenValues = eig.getRealEigenvalues();
			RealMatrix eigenVectors = eig.getVT();

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
	 * @param atomIndex
	 *            the atom index with respect to which the derivative are to be
	 *            evaluated
	 * @param scfMethod
	 *            the reference to the SCFMethod
	 * @return three element array of Overlap elements representing partial
	 *         derivatives with respect to x, y and z of atom position
	 */
	public List<Overlap> computeDerivative(int atomIndex, SCFMethod scfMethod) {
		this.scfMethod = scfMethod;
		this.atomIndex = atomIndex;

		ArrayList<Overlap> dOverlap = new ArrayList<>(3);

		int noOfBasisFunctions = this.getRowDimension();
		Overlap dOverlapDx = new Overlap(noOfBasisFunctions);
		Overlap dOverlapDy = new Overlap(noOfBasisFunctions);
		Overlap dOverlapDz = new Overlap(noOfBasisFunctions);

		double[][] hdx = dOverlapDx.getData();
		double[][] hdy = dOverlapDy.getData();
		double[][] hdz = dOverlapDz.getData();

		for (int i = 0; i < noOfBasisFunctions; i++) {
			for (int j = 0; j < noOfBasisFunctions; j++) {
				Vector3D dOvrEle = computeOverlapDerElement(i, j);

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
	private Vector3D computeOverlapDerElement(int i, int j) {
		BasisSetLibrary bsl = scfMethod.getOneEI().getBasisSetLibrary();
		ContractedGaussian cgi = bsl.getBasisFunctions().get(i);
		ContractedGaussian cgj = bsl.getBasisFunctions().get(j);

		return cgi.overlapDerivative(atomIndex, cgj);
	}

	@Override
	public String toString() {
		return MathUtil.matrixToString(this);
	}
}
