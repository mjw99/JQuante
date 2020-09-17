package name.mjw.jquante.math.qm.property;

import static org.junit.jupiter.api.Assertions.assertEquals;

import name.mjw.jquante.test.Fixtures;

import org.junit.jupiter.api.Test;

import name.mjw.jquante.math.qm.OneElectronIntegrals;
import name.mjw.jquante.math.qm.SCFMethod;
import name.mjw.jquante.math.qm.SCFMethodFactory;
import name.mjw.jquante.math.qm.SCFType;
import name.mjw.jquante.math.qm.TwoElectronIntegrals;
import name.mjw.jquante.math.qm.basis.BasisSetLibrary;
import name.mjw.jquante.molecule.Molecule;
import name.mjw.jquante.molecule.property.electronic.GridProperty;

class ElectronDensityTest {

	double diff = 0.0001;

	@Test
	void test() throws Exception {
		Molecule water = Fixtures.getWater();

		// Read Basis
		BasisSetLibrary basisFunctions = new BasisSetLibrary(water, "sto-3g");

		// Compute integrals
		OneElectronIntegrals oneElectronIntegrals = new OneElectronIntegrals(basisFunctions, water);
		TwoElectronIntegrals twoElectronIntegrals = new TwoElectronIntegrals(basisFunctions);

		// Set up HF method
		SCFMethod scfm = SCFMethodFactory.getInstance().getSCFMethod(water, oneElectronIntegrals, twoElectronIntegrals,
				SCFType.HARTREE_FOCK);

		// Actually do the SCF
		scfm.scf();

		// Set up the grid for where the electron density will be calculated:
		//
		// Set 5 grid points per x,y,z direction.
		// Extend the bounding box by 2 units in each direction.

		GridProperty gridProperty = new GridProperty(water.getBoundingBox().expand(2.0), 5);

		// Computes electron density on specified points for a Molecule
		ElectronDensity electronDensity = new ElectronDensity(scfm);
		electronDensity.compute(gridProperty);

		// Check the a few points
		assertEquals(0.08352817419051489, gridProperty.getFunctionValueAt(2, 2, 2), diff);

		assertEquals(0.1849781416342963, gridProperty.getFunctionValueAt(2, 2, 3), diff);

	}
}
