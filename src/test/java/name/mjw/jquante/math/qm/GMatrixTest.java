package name.mjw.jquante.math.qm;

import static org.junit.jupiter.api.Assertions.assertArrayEquals;

import org.junit.jupiter.api.BeforeAll;
import org.junit.jupiter.api.Test;

import name.mjw.jquante.math.qm.basis.BasisSetLibrary;
import name.mjw.jquante.molecule.Molecule;
import name.mjw.jquante.test.Fixtures;

class GMatrixTest {

	double diff = 0.00001;
	static Density density = null;
	static TwoElectronIntegrals twoEI = null;
	static BasisSetLibrary bf = null;

	@BeforeAll
	static void setUp() {
		// Create molecule
		Molecule hydrogen = Fixtures.getHydrogen();

		try {
			bf = new BasisSetLibrary(hydrogen, "sto-3g");
		} catch (Exception e) {
			e.printStackTrace();
		}
		twoEI = new TwoElectronIntegrals(bf);

		density = new Density(new double[][] { { 0.3012278366, 0.3012278366 }, { 0.3012278366, 0.3012278366 } });

	}

	@Test
	void testComputeInCore() {
		GMatrix gMatrix = new GMatrix(density.getRowDimension());

		gMatrix.compute(SCFType.HARTREE_FOCK, twoEI, density);

		double[][] actual = gMatrix.getData();
		double[][] expected = new double[][] { { 0.754933425, 0.3650758376 }, { 0.3650758376, 0.754933425 } };

		assertArrayEquals(expected[0], actual[0], diff);
		assertArrayEquals(expected[1], actual[1], diff);
	}

	@Test
	void testComputeDirect() {
		GMatrix gMatrix = new GMatrix(density.getRowDimension());

		gMatrix.compute(SCFType.HARTREE_FOCK_DIRECT, twoEI, density);

		double[][] actual = gMatrix.getData();
		double[][] expected = new double[][] { { 0.754933425, 0.3650758376 }, { 0.3650758376, 0.754933425 } };

		assertArrayEquals(expected[0], actual[0], diff);
		assertArrayEquals(expected[1], actual[1], diff);
	}
}
