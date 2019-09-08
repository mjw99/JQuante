package name.mjw.jquante.math.qm;


import org.apache.commons.math3.linear.Array2DRowRealMatrix;

import name.mjw.jquante.math.MathUtil;

/**
 * Represents the Fock matrix.
 * 
 * @author V.Ganesh
 * @version 2.0 (Part of MeTA v2.0)
 */
public final class Fock extends Array2DRowRealMatrix {

	private static final long serialVersionUID = 8771239880594969731L;

	/**
	 * Creates a new instance of square (NxN) Matrix
	 * 
	 * @param n
	 *            the dimension
	 */
	public Fock(int n) {
		super(n, n);
	}

	/**
	 * Creates a new instance of Matrix, based on already allocated 2D array
	 * 
	 * @param a
	 *            the 2D array
	 */
	public Fock(double[][] a) {
		super(a);
	}

	/**
	 * Form the Fock matrix from the HCore and the GMatrix
	 * 
	 * @param hCore
	 *            the HCore matrix
	 * @param gMatrix
	 *            the GMatrix
	 */
	public void compute(HCore hCore, GMatrix gMatrix) {
		setSubMatrix((hCore.add(gMatrix)).getData(), 0, 0);
	}

	@Override
	public String toString() {
		return MathUtil.matrixToString(this);
	}
}
