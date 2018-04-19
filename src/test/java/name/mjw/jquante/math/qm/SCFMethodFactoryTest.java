package name.mjw.jquante.math.qm;

import static org.junit.Assert.assertEquals;
import name.mjw.jquante.math.geom.Point3D;
import name.mjw.jquante.molecule.Atom;
import name.mjw.jquante.molecule.Molecule;
import name.mjw.jquante.molecule.impl.MoleculeImpl;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.junit.Test;

public class SCFMethodFactoryTest {

	private final Logger LOG = LogManager.getLogger(SCFMethodFactoryTest.class);

	double diff = 0.0001;

	@Test
	public void SinglePointHFHydrogenSTO3G() {

		// Create molecule
		Atom H1 = new Atom("H", 1.0, new Point3D(0.00000000, 0.00000000,
				0.00000000));
		Atom H2 = new Atom("H", 1.0, new Point3D(0.74000000, 0.00000000,
				0.00000000));

		Molecule hydrogen = new MoleculeImpl("hydrogen");
		hydrogen.addAtom(H1);
		hydrogen.addAtom(H2);

		long t1 = System.currentTimeMillis();
		// Read Basis
		BasisFunctions bf = null;

		try {
			bf = new BasisFunctions(hydrogen, "sto-3g");

		} catch (Exception e) {

			e.printStackTrace();
		}

		// compute integrals
		OneElectronIntegrals e1 = new OneElectronIntegrals(bf, hydrogen);
		TwoElectronIntegrals e2 = new TwoElectronIntegrals(bf);

		long t2 = System.currentTimeMillis();

		// do SCF
		SCFMethod scfm = SCFMethodFactory.getInstance().getSCFMethod(hydrogen,
				e1, e2, SCFType.HARTREE_FOCK);
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
		Atom H = new Atom("H", 1.0, new Point3D(0.00000000, 0.00000000,
				0.00000000));
		Atom F = new Atom("F", 7.0, new Point3D(0.91700000, 0.00000000,
				0.00000000));

		Molecule hydrogenFluoride = new MoleculeImpl("hydrogenFluoride");
		hydrogenFluoride.addAtom(H);
		hydrogenFluoride.addAtom(F);

		long t1 = System.currentTimeMillis();
		// Read Basis
		BasisFunctions bf = null;

		try {
			bf = new BasisFunctions(hydrogenFluoride, "sto-3g");

		} catch (Exception e) {

			e.printStackTrace();
		}

		// compute integrals
		OneElectronIntegrals e1 = new OneElectronIntegrals(bf, hydrogenFluoride);
		TwoElectronIntegrals e2 = new TwoElectronIntegrals(bf);

		long t2 = System.currentTimeMillis();

		// do SCF
		SCFMethod scfm = SCFMethodFactory.getInstance().getSCFMethod(
				hydrogenFluoride, e1, e2, SCFType.HARTREE_FOCK);
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
		Atom O = new Atom("O", 6.0, new Point3D(0.00000000, 0.000000, 0.119748));
		Atom H1 = new Atom("H", 1.0, new Point3D(0.00000000, 0.761561,
				-0.478993));
		Atom H2 = new Atom("H", 1.0, new Point3D(0.00000000, -0.761561,
				-0.478993));

		Molecule water = new MoleculeImpl("water");
		water.addAtom(O);
		water.addAtom(H1);
		water.addAtom(H2);

		long t1 = System.currentTimeMillis();
		// Read Basis
		BasisFunctions bf = null;

		try {
			bf = new BasisFunctions(water, "sto-3g");

		} catch (Exception e) {

			e.printStackTrace();
		}

		// compute integrals
		OneElectronIntegrals e1 = new OneElectronIntegrals(bf, water);
		TwoElectronIntegrals e2 = new TwoElectronIntegrals(bf);

		long t2 = System.currentTimeMillis();

		// do SCF
		SCFMethod scfm = SCFMethodFactory.getInstance().getSCFMethod(water, e1,
				e2, SCFType.HARTREE_FOCK);
		scfm.scf();

		long t3 = System.currentTimeMillis();

		LOG.debug("Time till 2E : " + (t2 - t1) + " ms");
		LOG.debug("Time for SCF : " + (t3 - t2) + " ms");

		assertEquals(9.087438510255588, scfm.nuclearEnergy(), diff);

		assertEquals(-74.9645174916, scfm.getEnergy(), diff);

		// orbital energies
		double[] ev = scfm.getOrbE();

		assertEquals(-20.244508003283794, ev[0], diff);
		assertEquals(-1.2637569211246862, ev[1], diff);
		assertEquals(-0.6106330109275634, ev[2], diff);
		assertEquals(-0.4535337774289125, ev[3], diff);
		assertEquals(-0.3913211137352477, ev[4], diff);
		assertEquals(0.5958982943824441, ev[5], diff);
		assertEquals(0.7260124013928394, ev[6], diff);

	}
}
