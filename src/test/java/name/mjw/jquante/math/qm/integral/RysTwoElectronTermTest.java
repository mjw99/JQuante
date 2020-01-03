package name.mjw.jquante.math.qm.integral;


import static org.junit.jupiter.api.Assertions.assertEquals;

import org.hipparchus.geometry.euclidean.threed.Vector3D;
import org.junit.jupiter.api.BeforeAll;
import org.junit.jupiter.api.Test;

import name.mjw.jquante.math.qm.basis.ContractedGaussian;
import name.mjw.jquante.math.qm.basis.Power;

public class RysTwoElectronTermTest {

	private final double delta = 0.000001;

	RysTwoElectronTerm e2 = new RysTwoElectronTerm();

	static ContractedGaussian cgtoS0;
	static ContractedGaussian cgtoS1;
	static ContractedGaussian cgtoP0;
	static ContractedGaussian cgtoD0;
	static ContractedGaussian cgtoF0;

	@BeforeAll
	public static void setUp() {

		cgtoS0 = new ContractedGaussian(new Vector3D(0, 0, 0), new Power(0, 0, 0));
		cgtoS0.addPrimitive(1.0, 1.0);
		cgtoS0.normalize();

		cgtoS1 = new ContractedGaussian(new Vector3D(0, 0, 1), new Power(0, 0, 0));
		cgtoS1.addPrimitive(1.0, 1.0);
		cgtoS1.normalize();

		cgtoP0 = new ContractedGaussian(new Vector3D(0, 0, 0), new Power(1, 0, 0));
		cgtoP0.addPrimitive(1.0, 1.0);
		cgtoP0.normalize();

		cgtoD0 = new ContractedGaussian(new Vector3D(0, 0, 0), new Power(2, 0, 0));
		cgtoD0.addPrimitive(1.0, 1.0);
		cgtoD0.normalize();

		cgtoF0 = new ContractedGaussian(new Vector3D(0, 0, 0), new Power(3, 0, 0));
		cgtoF0.addPrimitive(1.0, 1.0);
		cgtoF0.normalize();
	}

	@Test
	public void ss00() {
		assertEquals(1.1283791633342477, e2.coulomb(cgtoS0, cgtoS0, cgtoS0, cgtoS0), delta);
	}

	@Test
	public void ss01() {
		assertEquals(0.8427007900292194, e2.coulomb(cgtoS0, cgtoS0, cgtoS1, cgtoS1), delta);
	}

	@Test
	public void sp00() {
		assertEquals(0.9403159699467084, e2.coulomb(cgtoS0, cgtoS0, cgtoP0, cgtoP0), delta);
	}

	@Test
	public void sd00() {
		assertEquals(0.8086717345109515, e2.coulomb(cgtoS0, cgtoS0, cgtoD0, cgtoD0), delta);
	}

	@Test
	public void sf00() {
		assertEquals(0.7132968291998493, e2.coulomb(cgtoS0, cgtoS0, cgtoF0, cgtoF0), delta);
	}

	@Test
	public void pd00() {
		assertEquals(0.8583741496984647, e2.coulomb(cgtoP0, cgtoP0, cgtoD0, cgtoD0), delta);
	}

	@Test
	public void dd00() {
		assertEquals(0.854418854766963, e2.coulomb(cgtoD0, cgtoD0, cgtoD0, cgtoD0), delta);
	}

	@Test
	public void ff00() {
		assertEquals(0.8186960564969021, e2.coulomb(cgtoF0, cgtoF0, cgtoF0, cgtoF0), delta);
	}

}
