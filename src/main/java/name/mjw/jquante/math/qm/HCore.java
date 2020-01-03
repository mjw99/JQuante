package name.mjw.jquante.math.qm;

import java.util.ArrayList;
import java.util.List;

import org.hipparchus.geometry.euclidean.threed.Vector3D;
import org.hipparchus.linear.Array2DRowRealMatrix;

import name.mjw.jquante.math.MathUtil;
import name.mjw.jquante.math.qm.basis.BasisSetLibrary;
import name.mjw.jquante.math.qm.basis.ContractedGaussian;

/**
 * The HCore matrix. The core Hamiltonian matrix contains integrals that
 * represent the kinetic energy of an electron (T) and electron-nuclear
 * potential energy (V).
 * 
 * 
 * @author V.Ganesh
 * @version 2.0 (Part of MeTA v2.0)
 */
public final class HCore extends Array2DRowRealMatrix {

	private static final long serialVersionUID = 290891895527849860L;
	protected int atomIndex;
	private transient SCFMethod scfMethod;
	private transient BasisSetLibrary bsl;
	private transient List<ContractedGaussian> cgs;

	/**
	 * Creates a new instance of square (NxN) Matrix
	 * 
	 * @param n
	 *            the dimension
	 */
	public HCore(int n) {
		super(n, n);
	}
	/**
	 * Compute HCore partial derivative for an atom index.
	 * 
	 * @param atomIndex
	 *            the atom index with respect to which the derivative are to be
	 *            evaluated
	 * @param scfMethod
	 *            the reference to the SCFMethod
	 * @return three element array of HCore elements representing partial
	 *         derivatives with respect to x, y and z of atom position
	 */
	public List<HCore> computeDerivative(int atomIndex, SCFMethod scfMethod) {
		this.atomIndex = atomIndex;
		this.scfMethod = scfMethod;
		this.bsl = scfMethod.getOneEI().getBasisSetLibrary();
		this.cgs = this.bsl.getBasisFunctions();

		ArrayList<HCore> dHCore = new ArrayList<>(3);

		int noOfBasisFunctions = this.getRowDimension();
		HCore dHCoreDx = new HCore(noOfBasisFunctions);
		HCore dHCoreDy = new HCore(noOfBasisFunctions);
		HCore dHCoreDz = new HCore(noOfBasisFunctions);

		double[][] hdx = dHCoreDx.getData();
		double[][] hdy = dHCoreDy.getData();
		double[][] hdz = dHCoreDz.getData();

		int i;
		int j;
		for (i = 0; i < noOfBasisFunctions; i++) {
			for (j = 0; j < noOfBasisFunctions; j++) {
				Vector3D dHCoreEle = computeHCoreDerElement(atomIndex, i, j);

				hdx[i][j] = dHCoreEle.getX();
				hdy[i][j] = dHCoreEle.getY();
				hdz[i][j] = dHCoreEle.getZ();
			}
		}

		dHCore.add(dHCoreDx);
		dHCore.add(dHCoreDy);
		dHCore.add(dHCoreDz);

		return dHCore;
	}

	/**
	 * Compute one of the HCore derivative elements, with respect to an
	 * atomIndex
	 */
	private Vector3D computeHCoreDerElement(int atomIndex, int i, int j) {		
		ContractedGaussian cgi = this.cgs.get(i);
		ContractedGaussian cgj = this.cgs.get(j);

		return cgi.kineticDerivative(atomIndex, cgj).add(
				cgi.nuclearAttractionDerivative(scfMethod.getMolecule(),
						atomIndex, cgj));

	}

	@Override
	public String toString() {
		return MathUtil.matrixToString(this);
	}
}
