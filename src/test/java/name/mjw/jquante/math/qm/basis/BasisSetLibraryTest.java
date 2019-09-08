package name.mjw.jquante.math.qm.basis;

import static org.junit.jupiter.api.Assertions.assertEquals;

import org.apache.commons.math3.geometry.euclidean.threed.Vector3D;
import org.junit.jupiter.api.BeforeAll;
import org.junit.jupiter.api.Test;

import name.mjw.jquante.math.qm.basis.BasisSetLibrary;
import name.mjw.jquante.molecule.Atom;
import name.mjw.jquante.molecule.Molecule;
import name.mjw.jquante.molecule.impl.MoleculeImpl;

public class BasisSetLibraryTest {

	double diff = 0.00001;

	static Atom O;
	static Atom H1;
	static Atom H2;

	static Molecule water;
	BasisSetLibrary bsl = null;

	@BeforeAll
	public static void setup() {

		O = new Atom("O", new Vector3D(0.00000000, 0.000000, 0.119748));
		H1 = new Atom("H", new Vector3D(0.00000000, 0.761561, -0.478993));
		H2 = new Atom("H", new Vector3D(0.00000000, -0.761561, -0.478993));

		water = new MoleculeImpl("water");
		water.addAtom(O);
		water.addAtom(H1);
		water.addAtom(H2);
	}

	@Test
	public void testSTO3G() {

		try {
			bsl = new BasisSetLibrary(water, "sto-3g");

		} catch (Exception e) {

			e.printStackTrace();
		}

		// 7 functions
		assertEquals(7, bsl.getBasisFunctions().size());

		// 5 shells
		assertEquals(5, bsl.getShells().size());

		// 15 unique shell pairs
		assertEquals(15, bsl.getUniqueShellPairs().size());
	}

	@Test
	public void testCcPvtz() {

		try {
			bsl = new BasisSetLibrary(water, "cc-pvtz");

		} catch (Exception e) {

			e.printStackTrace();
		}
		// 65 functions
		assertEquals(65, bsl.getBasisFunctions().size());

		// 22 shells
		assertEquals(22, bsl.getShells().size());

		// 253 unique shell pairs
		assertEquals(253, bsl.getUniqueShellPairs().size());

	}

	@Test
	public void testPrimitiveOrdering() {

		BasisSetLibrary bsl = null;

		Atom H1 = new Atom("H", new Vector3D(0.752510, -0.454585, 0.000000));

		Molecule hydrogen = new MoleculeImpl("hydrogen");
		hydrogen.addAtom(H1);

		try {
			bsl = new BasisSetLibrary(hydrogen, "sto-3g");

		} catch (Exception e) {
			e.printStackTrace();
		}

		assertEquals(1, bsl.getBasisFunctions().size());

		assertEquals(3.42525091, bsl.getBasisFunctions().get(0).getExponents().get(0), diff);
		assertEquals(0.15432897, bsl.getBasisFunctions().get(0).getCoefficients().get(0), diff);

		assertEquals(0.62391373, bsl.getBasisFunctions().get(0).getExponents().get(1), diff);
		assertEquals(0.53532814, bsl.getBasisFunctions().get(0).getCoefficients().get(1), diff);

		assertEquals(0.1688554, bsl.getBasisFunctions().get(0).getExponents().get(2), diff);
		assertEquals(0.44463454, bsl.getBasisFunctions().get(0).getCoefficients().get(2), diff);

	}

	@Test
	public void testShellBasisFunctionIndexing() {

		try {
			bsl = new BasisSetLibrary(water, "sto-3g");

		} catch (Exception e) {

			e.printStackTrace();
		}

		Shell shell = bsl.getShells().get(1);

		assertEquals(1, shell.getFirstBasisFunctionIndex());
		assertEquals(1, shell.getLastBasisFunctionIndex());

		shell = bsl.getShells().get(2);

		assertEquals(2, shell.getFirstBasisFunctionIndex());
		assertEquals(4, shell.getLastBasisFunctionIndex());

	}

}
