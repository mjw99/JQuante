package name.mjw.jquante.math.qm.integral;


import static org.junit.jupiter.api.Assertions.assertEquals;

import org.apache.commons.math3.geometry.euclidean.threed.Vector3D;
import org.junit.jupiter.api.Test;

import name.mjw.jquante.math.qm.basis.Power;

public class OneElectronTermTest {

	private final double delta = 0.000001;

	OneElectronTerm oneElectronTerm = new OneElectronTerm();

	@Test
	public void overlap1D() {
		assertEquals(1.0, oneElectronTerm.overlap1D(0, 0, 0, 0, 1), delta);
	}

	@Test
	public void overlap() {
		assertEquals(1.968701, oneElectronTerm.overlap(1, new Power(0, 0, 0), new Vector3D(0, 0, 0), 1,
				new Power(0, 0, 0), new Vector3D(0, 0, 0)), delta);
	}

	@Test
	public void kinetic() {
		assertEquals(2.953052, oneElectronTerm.kinetic(1, new Power(0, 0, 0), new Vector3D(0, 0, 0), 1,
				new Power(0, 0, 0), new Vector3D(0, 0, 0)), delta);
	}

}
