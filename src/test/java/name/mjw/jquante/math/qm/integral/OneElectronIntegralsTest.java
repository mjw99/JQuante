package name.mjw.jquante.math.qm.integral;

import static org.junit.jupiter.api.Assertions.assertEquals;

import org.junit.jupiter.api.BeforeAll;
import org.junit.jupiter.api.Test;

import name.mjw.jquante.math.qm.OneElectronIntegrals;
import name.mjw.jquante.math.qm.basis.BasisSetLibrary;
import name.mjw.jquante.molecule.Molecule;
import name.mjw.jquante.test.Fixtures;

class OneElectronIntegralsTest {

	double diff = 0.0001;
	static OneElectronIntegrals e1 = null;

	@BeforeAll
	static void setUp() {

		// Create molecule
		Molecule water = Fixtures.getWater();

		// Read Basis
		BasisSetLibrary bf = null;

		try {
			bf = new BasisSetLibrary(water, "sto-3g");

		} catch (Exception e) {

			e.printStackTrace();
		}

		// compute integrals
		e1 = new OneElectronIntegrals(bf, water);

	}

	@Test
	void testgetHCore() {
		assertEquals(-32.70775590590571, e1.getHCore().getEntry(0, 0), diff);
		assertEquals(-1.7016641597, e1.getHCore().getEntry(6, 0), diff);
		assertEquals(-5.0309028634, e1.getHCore().getEntry(6, 6), diff);
	}

	@Test
	void testgetOverlap() {
		assertEquals(1.0, e1.getOverlap().getEntry(0, 0), diff);
		assertEquals(0.05251779176613285, e1.getOverlap().getEntry(6, 0), diff);
		assertEquals(1.0, e1.getOverlap().getEntry(6, 6), diff);
	}

	@Test
	void testgetBasisSetLibrary() {
		assertEquals("sto-3g", e1.getBasisSetLibrary().getBasisName());
	}

}
