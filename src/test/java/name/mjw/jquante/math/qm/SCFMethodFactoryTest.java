package name.mjw.jquante.math.qm;


import name.mjw.jquante.math.qm.basis.BasisSetLibrary;
import name.mjw.jquante.molecule.Atom;
import name.mjw.jquante.molecule.Molecule;
import name.mjw.jquante.molecule.impl.MoleculeImpl;

import static org.junit.jupiter.api.Assertions.assertEquals;

import org.apache.commons.math3.geometry.euclidean.threed.Vector3D;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.junit.jupiter.api.Test;

public class SCFMethodFactoryTest {

	private final Logger LOG = LogManager.getLogger(SCFMethodFactoryTest.class);

	double diff = 0.0001;

	@Test
	public void SinglePointHFHydrogenSTO3G() {

		// Create molecule
		Atom H1 = new Atom("H", new Vector3D(0.00000000, 0.00000000, 0.00000000));
		Atom H2 = new Atom("H", new Vector3D(0.74000000, 0.00000000, 0.00000000));

		Molecule hydrogen = new MoleculeImpl("hydrogen");
		hydrogen.addAtom(H1);
		hydrogen.addAtom(H2);

		long t1 = System.currentTimeMillis();
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

		long t2 = System.currentTimeMillis();

		// do SCF
		SCFMethod scfm = SCFMethodFactory.getInstance().getSCFMethod(hydrogen, e1, e2, SCFType.HARTREE_FOCK);
		scfm.scf();

		long t3 = System.currentTimeMillis();

		LOG.debug("Time till 2E : " + (t2 - t1) + " ms");
		LOG.debug("Time for SCF : " + (t3 - t2) + " ms");

		assertEquals(0.7151043908648649, scfm.nuclearEnergy(), diff);

		assertEquals(-1.1167593656305694, scfm.getEnergy(), diff);

		// orbital energies
		double[] ev = scfm.getOrbE();

		assertEquals(-0.578554081990778, ev[0], diff);
		assertEquals(0.6711435173832502, ev[1], diff);
	}

	@Test
	// Values checked with NWChem 6.6
	public void SinglePointHFHydrogenFluorideSTO3G() {

		// Create molecule
		Atom H = new Atom("H", new Vector3D(0.00000000, 0.00000000, 0.00000000));
		Atom F = new Atom("F", new Vector3D(0.91700000, 0.00000000, 0.00000000));

		Molecule hydrogenFluoride = new MoleculeImpl("hydrogenFluoride");
		hydrogenFluoride.addAtom(H);
		hydrogenFluoride.addAtom(F);

		long t1 = System.currentTimeMillis();
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

		long t2 = System.currentTimeMillis();

		// do SCF
		SCFMethod scfm = SCFMethodFactory.getInstance().getSCFMethod(hydrogenFluoride, e1, e2, SCFType.HARTREE_FOCK);
		scfm.scf();

		long t3 = System.currentTimeMillis();

		LOG.debug("Time till 2E : " + (t2 - t1) + " ms");
		LOG.debug("Time for SCF : " + (t3 - t2) + " ms");

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
	public void SinglePointHFWaterSTO3G() {

		// Create molecule
		Atom O = new Atom("O", new Vector3D(0.00000000, 0.000000, 0.119748));
		Atom H1 = new Atom("H", new Vector3D(0.00000000, 0.761561, -0.478993));
		Atom H2 = new Atom("H", new Vector3D(0.00000000, -0.761561, -0.478993));

		Molecule water = new MoleculeImpl("water");
		water.addAtom(O);
		water.addAtom(H1);
		water.addAtom(H2);

		long t1 = System.currentTimeMillis();
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

		long t2 = System.currentTimeMillis();

		// do SCF
		SCFMethod scfm = SCFMethodFactory.getInstance().getSCFMethod(water, e1, e2, SCFType.HARTREE_FOCK);
		scfm.scf();

		long t3 = System.currentTimeMillis();

		LOG.debug("Time till 2E : " + (t2 - t1) + " ms");
		LOG.debug("Time for SCF : " + (t3 - t2) + " ms");

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
	public void SinglePointHFWaterSTO3GDirect() {

		// Create molecule
		Atom O = new Atom("O", new Vector3D(0.00000000, 0.000000, 0.119748));
		Atom H1 = new Atom("H", new Vector3D(0.00000000, 0.761561, -0.478993));
		Atom H2 = new Atom("H", new Vector3D(0.00000000, -0.761561, -0.478993));

		Molecule water = new MoleculeImpl("water");
		water.addAtom(O);
		water.addAtom(H1);
		water.addAtom(H2);

		long t1 = System.currentTimeMillis();
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

		long t2 = System.currentTimeMillis();

		// do SCF
		SCFMethod scfm = SCFMethodFactory.getInstance().getSCFMethod(water, e1, e2, SCFType.HARTREE_FOCK_DIRECT);
		scfm.scf();

		long t3 = System.currentTimeMillis();

		LOG.debug("Time till 2E : " + (t2 - t1) + " ms");
		LOG.debug("Time for SCF : " + (t3 - t2) + " ms");

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
	public void SinglePointHFWaterSTO3GShellPair() {

		// Create molecule
		Atom O = new Atom("O", new Vector3D(0.00000000, 0.000000, 0.119748));
		Atom H1 = new Atom("H", new Vector3D(0.00000000, 0.761561, -0.478993));
		Atom H2 = new Atom("H", new Vector3D(0.00000000, -0.761561, -0.478993));

		Molecule water = new MoleculeImpl("water");
		water.addAtom(O);
		water.addAtom(H1);
		water.addAtom(H2);

		long t1 = System.currentTimeMillis();
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

		long t2 = System.currentTimeMillis();

		// do SCF
		SCFMethod scfm = SCFMethodFactory.getInstance().getSCFMethod(water, e1, e2, SCFType.HARTREE_FOCK);
		scfm.scf();

		long t3 = System.currentTimeMillis();

		LOG.debug("Time till 2E : " + (t2 - t1) + " ms");
		LOG.debug("Time for SCF : " + (t3 - t2) + " ms");

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
