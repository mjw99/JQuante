package name.mjw.jquante.math.qm.integral;

import static org.junit.jupiter.api.Assertions.assertEquals;

import org.junit.jupiter.api.Test;

import name.mjw.jquante.math.qm.basis.ContractedGaussian;
import name.mjw.jquante.test.Fixtures;

class HuzinagaTwoElectronTermTest {

	private final double delta = 0.000001;

	HuzinagaTwoElectronTerm e2 = new HuzinagaTwoElectronTerm();

	static ContractedGaussian cgtoS0 = Fixtures.getCgtoS0();
	static ContractedGaussian cgtoS1 = Fixtures.getCgtoS1();
	static ContractedGaussian cgtoP0 = Fixtures.getCgtoP0();
	static ContractedGaussian cgtoD0 = Fixtures.getCgtoD0();
	static ContractedGaussian cgtoF0 = Fixtures.getCgtoF0();

	@Test
	void ss00() {
		assertEquals(1.1283791633342477, e2.coulomb(cgtoS0, cgtoS0, cgtoS0, cgtoS0), delta);
	}

	@Test
	void ss01() {
		assertEquals(0.8427007900292194, e2.coulomb(cgtoS0, cgtoS0, cgtoS1, cgtoS1), delta);
	}

	@Test
	void sp00() {
		assertEquals(0.9403159699467084, e2.coulomb(cgtoS0, cgtoS0, cgtoP0, cgtoP0), delta);
	}

	@Test
	void sd00() {
		assertEquals(0.8086717345109515, e2.coulomb(cgtoS0, cgtoS0, cgtoD0, cgtoD0), delta);
	}

	@Test
	void sf00() {
		assertEquals(0.7132968291998493, e2.coulomb(cgtoS0, cgtoS0, cgtoF0, cgtoF0), delta);
	}

	@Test
	void pd00() {
		assertEquals(0.8583741496984647, e2.coulomb(cgtoP0, cgtoP0, cgtoD0, cgtoD0), delta);
	}

	@Test
	void dd00() {
		assertEquals(0.854418854766963, e2.coulomb(cgtoD0, cgtoD0, cgtoD0, cgtoD0), delta);
	}

	@Test
	void ff00() {
		assertEquals(0.8186960564969021, e2.coulomb(cgtoF0, cgtoF0, cgtoF0, cgtoF0), delta);
	}
}
