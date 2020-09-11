package name.mjw.jquante.math.qm.property;

import java.util.List;

import org.hipparchus.geometry.euclidean.threed.Vector3D;
import org.hipparchus.linear.ArrayRealVector;
import org.hipparchus.linear.RealVector;

import name.mjw.jquante.math.qm.Density;
import name.mjw.jquante.math.qm.SCFMethod;
import name.mjw.jquante.math.qm.basis.BasisSetLibrary;
import name.mjw.jquante.math.qm.basis.ContractedGaussian;

/**
 * Computes electron density on specified points for a Molecule.
 * 
 * @author V.Ganesh
 * @version 2.0 (Part of MeTA v2.0)
 */
public class ElectronDensity extends OneElectronProperty {

	private final int nbf;

	private final List<ContractedGaussian> bfs;

	private final double[][] dm;

	/**
	 * Creates a new instance of ElectronDensity
	 * 
	 * @param scfMethod the Self Consistent Field (SCF) method
	 */
	public ElectronDensity(final SCFMethod scfMethod) {
		super(scfMethod);

		bfs = scfMethod.getOneEI().getBasisSetLibrary().getBasisFunctions();
		nbf = bfs.size();

		dm = scfMethod.getDensity().getData();
	}

	/**
	 * Creates a new instance of ElectronDensity
	 * 
	 * @param bsl the basis functions of a given molecule and a basis set
	 * @param den the density matrix
	 */
	public ElectronDensity(final BasisSetLibrary bsl, final Density den) {
		this.bfs = bsl.getBasisFunctions();
		this.nbf = this.bfs.size();

		this.dm = den.getData();
	}

	/**
	 * Computes the one electron property on the specified point and returns its
	 * value at the specified point. <br>
	 * Note that the unit of Point3D object must be a.u. No attempt is made to
	 * verify this.
	 * 
	 * @param point the point of interest
	 * @return the value of this property at this point
	 */
	@Override
	public double compute(final Vector3D point) {
		final RealVector amplitudes = new ArrayRealVector(nbf);
		final double[] amp = amplitudes.toArray();
		double density = 0.0;
		double amp_xyz;
		int m;
		int l;

		for (m = 0; m < nbf; m++)
			amp[m] = bfs.get(m).amplitude(point);

		for (m = 0; m < nbf; m++) {
			amp_xyz = 0.0;
			for (l = 0; l < nbf; l++) {
				amp_xyz += dm[m][l] * amp[l];
			}
			density += amp_xyz * amp[m];
		}

		// 2.0 represents double occupancy
		return 2.0 * density;
	}
}
