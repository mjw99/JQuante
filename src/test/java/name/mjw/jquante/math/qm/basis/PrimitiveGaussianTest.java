package name.mjw.jquante.math.qm.basis;


import static org.junit.jupiter.api.Assertions.assertEquals;

import org.hipparchus.geometry.euclidean.threed.Vector3D;
import org.junit.jupiter.api.Test;

class PrimitiveGaussianTest {

	private final double delta = 0.000001;

	@Test
	void s() {
		PrimitiveGaussian gto = new PrimitiveGaussian(new Vector3D(0, 0, 0),
				new Power(0, 0, 0), 1.0, 1.0);

		assertEquals(0.712705, gto.amplitude(new Vector3D(0, 0, 0)), delta);

	}

	@Test
	void d() {
		PrimitiveGaussian gto = new PrimitiveGaussian(new Vector3D(0, 0, 0),
				new Power(1, 0, 1), 1.0, 1.0);

		assertEquals(0.0, gto.amplitude(new Vector3D(0, 0, 0)), delta);

	}
	
	@Test
	void ordering() {
		PrimitiveGaussian a = new PrimitiveGaussian(new Vector3D(0, 0, 0), new Power(0, 0, 0), 1.0, 1.0);
		PrimitiveGaussian b = new PrimitiveGaussian(new Vector3D(0, 0, 0), new Power(0, 0, 0), 2.0, 1.0);

		assertEquals(-1, b.compareTo(a));
		assertEquals(1, a.compareTo(b));

	}
}
