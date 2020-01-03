package name.mjw.jquante.math.qm.property;

import java.util.List;

import org.hipparchus.geometry.euclidean.threed.Vector3D;
import org.hipparchus.linear.ArrayRealVector;
import org.hipparchus.linear.RealVector;

import name.mjw.jquante.math.qm.MolecularOrbitals;
import name.mjw.jquante.math.qm.SCFMethod;
import name.mjw.jquante.math.qm.basis.BasisSetLibrary;
import name.mjw.jquante.math.qm.basis.ContractedGaussian;
import name.mjw.jquante.molecule.Molecule;

/**
 * Computes MO density on specified points for a Molecule.
 * 
 * @author V.Ganesh
 * @version 2.0 (Part of MeTA v2.0)
 */
public class MODensity extends OneElectronProperty {

	private int monumber;

	private int nbf;
	private double[][] mos;

	private List<ContractedGaussian> bfs;

	/**
	 * Creates a new instance of MODensity.
	 *
	 * @param scfMethod
	 *            the Self Consistent Field (SCF) method.
	 */
	public MODensity(SCFMethod scfMethod) {
		super(scfMethod);

		bfs = scfMethod.getOneEI().getBasisSetLibrary().getBasisFunctions();
		nbf = bfs.size();

		mos = scfMethod.getMos().getCoefficients();

		monumber = scfMethod.getMolecule().getNumberOfElectrons() / 2;
	}

	/**
	 * Creates a new instance of MODensity.
	 *
	 * @param scfMethod
	 *            the Self Consistent Field (SCF) method.
	 * @param monumber
	 *            the number of molecular orbitals.
	 */
	public MODensity(SCFMethod scfMethod, int monumber) {
		this(scfMethod);

		this.monumber = monumber;
	}

	/**
	 * Creates a new instance of MODensity.
	 *
	 * @param molecule
	 *            the molecule object.
	 * @param bsl
	 *            the basis functions of a given molecule and a basis set.
	 * @param mos
	 *            the molecular orbitals as a coefficient matrix.
	 */
	public MODensity(Molecule molecule, BasisSetLibrary bsl, MolecularOrbitals mos) {

		this.bfs = bsl.getBasisFunctions();
		this.nbf = this.bfs.size();
		this.mos = mos.getCoefficients();
		this.monumber = molecule.getNumberOfElectrons() / 2;
	}

	/**
	 * Creates a new instance of MODensity.
	 *
	 * @param molecule
	 *            the molecule.
	 * @param bsl
	 *            the basis functions of a given molecule and a basis set.
	 * @param mos
	 *            the molecular orbitals as a coefficient matrix.
	 * @param monumber
	 *            the number of molecular orbitals.
	 *
	 */
	public MODensity(Molecule molecule, BasisSetLibrary bsl, MolecularOrbitals mos, int monumber) {

		this.bfs = bsl.getBasisFunctions();
		this.nbf = this.bfs.size();
		this.mos = mos.getCoefficients();
		this.monumber = monumber;
	}

	/**
	 * Computes the one electron property on the specified point and returns its
	 * value at the specified point. <br>
	 * Note that the unit of Point3D object must be a.u. No attempt is made to
	 * verify this.
	 *
	 * @param point
	 *            the point of interest
	 * @return the value of this property at this point
	 */
	@Override
	public double compute(Vector3D point) {
		RealVector amplitudes = new ArrayRealVector(nbf);
		double[] amp = amplitudes.toArray();
		double moden = 0.0;
		double amp_m;
		double mos_m;
		int m;
		int l;

		for (m = 0; m < nbf; m++)
			amp[m] = bfs.get(m).amplitude(point);

		moden = 0.0;

		for (m = 0; m < nbf; m++) {
			amp_m = amp[m];
			mos_m = mos[monumber][m];

			for (l = 0; l < nbf; l++) {
				moden += mos_m * mos[monumber][l] * amp[l] * amp_m;
			}
		}

		return moden;
	}
}
