package name.mjw.jquante.math.qm;

import name.mjw.jquante.math.qm.basis.BasisSetLibrary;
import name.mjw.jquante.molecule.Molecule;
import name.mjw.jquante.test.Fixtures;

import static org.junit.jupiter.api.Assertions.assertEquals;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.junit.jupiter.api.Test;

class SCFMethodFactoryTest {

	private final Logger LOG = LogManager.getLogger(SCFMethodFactoryTest.class);

	double diff = 0.0001;

	@Test
	void SinglePointHFHydrogenSTO3G() {

		// Create molecule
		Molecule hydrogen = Fixtures.getHydrogen();

		// Read Basis
		BasisSetLibrary bsl = null;

		try {
			bsl = new BasisSetLibrary(hydrogen, "sto-3g");

		} catch (Exception e) {

			e.printStackTrace();
		}

		// compute integrals
		OneElectronIntegrals e1 = new OneElectronIntegrals(bsl, hydrogen);
		TwoElectronIntegrals e2 = new TwoElectronIntegrals(bsl);

		// do SCF
		SCFMethod scfm = SCFMethodFactory.getInstance().getSCFMethod(hydrogen, e1, e2, SCFType.HARTREE_FOCK);
		scfm.scf();

		assertEquals(0.7151043908648649, scfm.nuclearEnergy(), diff);

		assertEquals(-1.1167593656305694, scfm.getEnergy(), diff);

		// orbital energies
		double[] ev = scfm.getOrbE();

		assertEquals(-0.578554081990778, ev[0], diff);
		assertEquals(0.6711435173832502, ev[1], diff);
	}

	@Test
	// Values checked with NWChem 6.6
	void SinglePointHFHydrogenFluorideSTO3G() {

		// Create molecule
		Molecule hydrogenFluoride = Fixtures.getHydrogenFluoride();

		// Read Basis
		BasisSetLibrary bf = null;

		try {
			bf = new BasisSetLibrary(hydrogenFluoride, "sto-3g");

		} catch (Exception e) {

			e.printStackTrace();
		}

		// compute integrals
		OneElectronIntegrals e1 = new OneElectronIntegrals(bf, hydrogenFluoride);
		TwoElectronIntegrals e2 = new TwoElectronIntegrals(bf);

		// do SCF
		SCFMethod scfm = SCFMethodFactory.getInstance().getSCFMethod(hydrogenFluoride, e1, e2, SCFType.HARTREE_FOCK);
		scfm.scf();

		assertEquals(5.1936698398691385, scfm.nuclearEnergy(), diff);
		assertEquals(-98.5707789400326, scfm.getEnergy(), diff);

		// orbital energies
		double[] ev = scfm.getOrbE();

		assertEquals(-25.89997382162182, ev[0], diff);
		assertEquals(-1.4711982487067965, ev[1], diff);
		assertEquals(-0.5851370412095496, ev[2], diff);
		assertEquals(-0.4641445641331199, ev[3], diff);
		assertEquals(-0.4641445641331199, ev[4], diff);
		assertEquals(0.6290476912112247, ev[5], diff);
	}

	@Test
	void SinglePointHFWaterSTO3G() {

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
		OneElectronIntegrals e1 = new OneElectronIntegrals(bf, water);
		TwoElectronIntegrals e2 = new TwoElectronIntegrals(bf);

		// do SCF
		SCFMethod scfm = SCFMethodFactory.getInstance().getSCFMethod(water, e1, e2, SCFType.HARTREE_FOCK);
		scfm.scf();

		assertEquals(9.087438510255588, scfm.nuclearEnergy(), diff);
		assertEquals(-74.964518362274, scfm.getEnergy(), diff);

		// orbital energies
		double[] ev = scfm.getOrbE();

		assertEquals(-20.24450742, ev[0], diff * 10);
		assertEquals(-1.26375686, ev[1], diff * 10);
		assertEquals(-0.61063305, ev[2], diff * 10);
		assertEquals(-0.45353394, ev[3], diff * 10);
		assertEquals(-0.39132131, ev[4], diff * 10);
		assertEquals(0.59589853, ev[5], diff * 10);
		assertEquals(0.72601218, ev[6], diff * 10);

	}

	@Test
	void SinglePointHFWaterSTO3GDirect() {

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
		OneElectronIntegrals e1 = new OneElectronIntegrals(bf, water);
		TwoElectronIntegrals e2 = new TwoElectronIntegrals(bf, water, true);

		// do SCF
		SCFMethod scfm = SCFMethodFactory.getInstance().getSCFMethod(water, e1, e2, SCFType.HARTREE_FOCK_DIRECT);
		scfm.scf();

		assertEquals(9.087438510255588, scfm.nuclearEnergy(), diff);

		assertEquals(-74.964518362274, scfm.getEnergy(), diff);

		// orbital energies
		double[] ev = scfm.getOrbE();

		assertEquals(-20.24450742, ev[0], diff * 10);
		assertEquals(-1.26375686, ev[1], diff * 10);
		assertEquals(-0.61063305, ev[2], diff * 10);
		assertEquals(-0.45353394, ev[3], diff * 10);
		assertEquals(-0.39132131, ev[4], diff * 10);
		assertEquals(0.59589853, ev[5], diff * 10);
		assertEquals(0.72601218, ev[6], diff * 10);

	}

	@Test
	void SinglePointHFWaterSTO3GShellPair() {

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
		OneElectronIntegrals e1 = new OneElectronIntegrals(bf, water);
		TwoElectronIntegrals e2 = new TwoElectronIntegrals(bf, water, false);

		// do SCF
		SCFMethod scfm = SCFMethodFactory.getInstance().getSCFMethod(water, e1, e2, SCFType.HARTREE_FOCK);
		scfm.scf();

		assertEquals(9.087438510255588, scfm.nuclearEnergy(), diff);
		assertEquals(-74.964518362274, scfm.getEnergy(), diff);

		// orbital energies
		double[] ev = scfm.getOrbE();

		assertEquals(-20.24450742, ev[0], diff * 10);
		assertEquals(-1.26375686, ev[1], diff * 10);
		assertEquals(-0.61063305, ev[2], diff * 10);
		assertEquals(-0.45353394, ev[3], diff * 10);
		assertEquals(-0.39132131, ev[4], diff * 10);
		assertEquals(0.59589853, ev[5], diff * 10);
		assertEquals(0.72601218, ev[6], diff * 10);

	}
}
