package name.mjw.jquante.math.qm.basis;

import static org.junit.jupiter.api.Assertions.assertEquals;

import org.hipparchus.geometry.euclidean.threed.Vector3D;
import org.junit.jupiter.api.BeforeAll;
import org.junit.jupiter.api.Test;

class ShellTest {

	static ContractedGaussian cgtoS0;
	static Shell shell;

	@BeforeAll
	static void setUp() {
		cgtoS0 = new ContractedGaussian(new Vector3D(0, 0, 0), new Power(0, 0, 0));
		cgtoS0.addPrimitive(1.0, 1.0);
		cgtoS0.normalize();

		shell = new Shell(cgtoS0);

	}

	@Test
	void testFirstBasisFunctionIndex() {
		shell.setFirstBasisFunctionIndex(1);

		assertEquals(1, shell.getFirstBasisFunctionIndex());
	}

	@Test
	void testLastBasisFunctionIndex() {
		shell.setLastBasisFunctionIndex(2);

		assertEquals(2, shell.getLastBasisFunctionIndex());
	}

}
