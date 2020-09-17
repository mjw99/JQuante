package name.mjw.jquante.math.qm.basis;

/**
 * A shell is a set of basis function with the same centre and same contracted
 * exponent, for example, sp shell, d shell.
 * <p>
 * It also contains a mapping back to the original BasisFunctions that it is
 * composed of.
 * 
 */
public final class Shell extends ContractedGaussian {

	public Shell(ContractedGaussian cg) {
		super(cg);
	}

	private int firstBasisFunctionIndex;
	private int lastBasisFunctionIndex;

	public int getFirstBasisFunctionIndex() {
		return firstBasisFunctionIndex;
	}

	public void setFirstBasisFunctionIndex(int firstBasisFunctionIndex) {
		this.firstBasisFunctionIndex = firstBasisFunctionIndex;
	}

	public int getLastBasisFunctionIndex() {
		return lastBasisFunctionIndex;
	}

	public void setLastBasisFunctionIndex(int lastBasisFunctionIndex) {
		this.lastBasisFunctionIndex = lastBasisFunctionIndex;
	}

	@Override
	public String toString() {
		return Integer.toString(firstBasisFunctionIndex) + "-" + Integer.toString(lastBasisFunctionIndex) + "\n";
	}

	@Override
	public boolean equals(Object that) {
		if (this == that)
			return true;
		if (that == null)
			return false;
		if (getClass() != that.getClass())
			return false;
		Shell other = (Shell) that;
		return this.isSameShell(other);
	}

}
