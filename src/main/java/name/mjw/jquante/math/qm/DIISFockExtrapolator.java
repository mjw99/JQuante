package name.mjw.jquante.math.qm;

import java.text.DecimalFormat;
import java.util.ArrayList;

import name.mjw.jquante.math.MathUtil;
import org.hipparchus.linear.Array2DRowRealMatrix;
import org.hipparchus.linear.ArrayRealVector;
import org.hipparchus.linear.DecompositionSolver;
import org.hipparchus.linear.LUDecomposition;
import org.hipparchus.linear.RealMatrix;
import org.hipparchus.linear.RealMatrixFormat;
import org.hipparchus.linear.RealVector;
import org.hipparchus.linear.RealVectorFormat;
import org.hipparchus.exception.MathRuntimeException;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

/**
 * DIIS (Direct Inversion of iterative subspaces) proposed by Peter Pulay for
 * HF convergence acceleration.
 * 
 * @author V.Ganesh
 * @version 2.0 (Part of MeTA v2.0)
 */
public final class DIISFockExtrapolator implements FockExtrapolator {

	private static final Logger LOG = LogManager.getLogger(DIISFockExtrapolator.class);

	private ArrayList<Fock> fockMatrixList = new ArrayList<>();
	private ArrayList<RealVector> errorVectorList = new ArrayList<>();

	protected int diisStep = 0;

	private boolean isDiisStarted = false;

	private Fock oldFock = null;

	protected double errorThreshold = 0.00001;

	/**
	 * Get the next extrapolated fock matrix
	 * 
	 * @param currentFock
	 *            the current fock
	 * @param overlap
	 *            the overlap matrix
	 * @param density
	 *            the current density matrix
	 * @return extrapolated fock
	 */
	@Override
	public Fock next(Fock currentFock, Overlap overlap, Density density) {
		DecimalFormat df = new DecimalFormat("+#,##0.000;-#");
		df.setGroupingUsed(false);
		RealVectorFormat vf = new RealVectorFormat("{", "}", ",", df);
		RealMatrixFormat mf = new RealMatrixFormat("\n", "", "", "", "\n", " ", df);

		Fock newFock = new Fock(currentFock.getData());

		RealMatrix fPS = currentFock.multiply(density).multiply(overlap);
		RealMatrix sPF = overlap.multiply(density).multiply(currentFock);

		// The commutator of the Fock and density matrices (the orbital gradient)
		RealVector errorVector = MathUtil.realMatrixToRealVector(fPS.subtract(sPF));

		if (LOG.isDebugEnabled()) {
			LOG.debug("errorVector {} ", vf.format(errorVector));
		}
		double mxerr = errorVector.getNorm();
		LOG.debug("errorVector.getNorm() {}", mxerr);

		if (mxerr > errorThreshold && !isDiisStarted) {
			LOG.debug("starting DIIS");
			isDiisStarted = true;
		}

		// bootstrap
		if (!isDiisStarted) {
			if (oldFock == null) {
				oldFock = currentFock;
				return currentFock;
			} else {
				newFock.setSubMatrix(oldFock.scalarMultiply(0.5).add(currentFock.scalarMultiply(0.5)).getData(), 0, 0);
				oldFock = currentFock;
				return newFock;
			}
		}

		fockMatrixList.add(currentFock);
		errorVectorList.add(errorVector);

		int noOfIterations = errorVectorList.size();

		LOG.debug("noOfIterations: {}", noOfIterations);
		int noOfIterationsPlusOne = noOfIterations + 1;

		// B_{ij} matrix
		RealMatrix aMatrix = new Array2DRowRealMatrix(noOfIterationsPlusOne, noOfIterationsPlusOne);
		// 0,0...,-1 vector
		RealVector bVector = new ArrayRealVector(noOfIterationsPlusOne);

		// set up A x = b to be solved
		// Populate B_{ij} matrix
		for (int i = 0; i < noOfIterations; i++) {
			for (int j = 0; j < noOfIterations; j++) {
				aMatrix.setEntry(i, j, errorVectorList.get(i).dotProduct(errorVectorList.get(j)));
			}
		}

		for (int i = 0; i < noOfIterations; i++) {
			aMatrix.setEntry(noOfIterations, i, -1.0);
			aMatrix.setEntry(i, noOfIterations, -1.0);
			bVector.setEntry(i, 0.0);
		}

		aMatrix.setEntry(noOfIterations, noOfIterations, 0.0);
		bVector.setEntry(noOfIterations, -1.0);

		if (LOG.isDebugEnabled()) {
			LOG.debug("aMatrix {}", mf.format(aMatrix));
			LOG.debug("bVector {}", vf.format(bVector));
		}

		try {
			DecompositionSolver solver = new LUDecomposition(aMatrix).getSolver();
			RealVector solVec = solver.solve(bVector);
			if (LOG.isDebugEnabled()) {
				LOG.debug("solVec {}", vf.format(solVec));
			}

			newFock = new Fock(newFock.getRowDimension());
			Fock tmpFock;

			for (int i = 0; i < noOfIterations; i++) {
				tmpFock = new Fock(fockMatrixList.get(i).scalarMultiply(solVec.toArray()[i]).getData());
				newFock = new Fock(newFock.add(tmpFock).getData());
			}
			oldFock = currentFock;

		} catch (MathRuntimeException ignored) {
			if (LOG.isDebugEnabled()) {
				LOG.debug("No solution: {}", diisStep);
			}

			diisStep++;
			return currentFock;

		}
		diisStep++;
		if (LOG.isDebugEnabled()) {
			LOG.debug("final newFock {}", mf.format(newFock));
		}
		return newFock;
	}

	/**
	 * Get the value of errorThreshold
	 * 
	 * @return the value of errorThreshold
	 */
	public double getErrorThreshold() {
		return errorThreshold;
	}

	/**
	 * Set the value of errorThreshold
	 * 
	 * @param errorThreshold
	 *            new value of errorThreshold
	 */
	public void setErrorThreshold(double errorThreshold) {
		this.errorThreshold = errorThreshold;
	}
}