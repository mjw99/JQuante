package name.mjw.jquante.math.qm.integral;


import static org.junit.jupiter.api.Assertions.assertEquals;

import org.hipparchus.geometry.euclidean.threed.Vector3D;
import org.junit.jupiter.api.Test;

import name.mjw.jquante.math.qm.basis.Power;

class OneElectronTermTest {

	private final double delta = 0.000001;

	OneElectronTerm oneElectronTerm = new OneElectronTerm();

	@Test
	void overlap1D() {
		assertEquals(1.0, oneElectronTerm.overlap1D(0, 0, 0, 0, 1), delta);
	}

	@Test
	void overlap() {
		assertEquals(1.968701, oneElectronTerm.overlap(1, new Power(0, 0, 0), new Vector3D(0, 0, 0), 1,
				new Power(0, 0, 0), new Vector3D(0, 0, 0)), delta);
	}

	@Test
	void kinetic() {
		assertEquals(2.953052, oneElectronTerm.kinetic(1, new Power(0, 0, 0), new Vector3D(0, 0, 0), 1,
				new Power(0, 0, 0), new Vector3D(0, 0, 0)), delta);
	}

}
