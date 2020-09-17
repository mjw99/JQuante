package name.mjw.jquante.math;

import static org.junit.jupiter.api.Assertions.assertEquals;

import org.hipparchus.geometry.euclidean.threed.Vector3D;
import org.junit.jupiter.api.Test;

class MathUtilTest {
	double diff = 0.00001;

	@Test
	void factorial2MinusOne() {

		assertEquals(1, MathUtil.factorial2(-1));

	}

	@Test
	void factorial2Zero() {

		assertEquals(1, MathUtil.factorial2(0));

	}

	@Test
	void factorial2Nine() {

		assertEquals(945, MathUtil.factorial2(9));

	}

	@Test
	void testfindDihedralAngle1() {

		Vector3D i = new Vector3D(60.614941, 67.371346, 12.713049);
		Vector3D j = new Vector3D(59.104156, 67.044731, 12.827555);
		Vector3D k = new Vector3D(58.827194, 65.911911, 13.835642);
		Vector3D l = new Vector3D(59.367950, 66.384277, 15.233553);

		assertEquals(59.66383, Math.toDegrees(MathUtil.findDihedral(i, j, k, l)), diff);

	}

	@Test
	void testfindDihedralAngle2() {

		Vector3D i = new Vector3D(30.958572, 73.709137, 38.031029);
		Vector3D j = new Vector3D(31.924915, 72.577698, 37.985279);
		Vector3D k = new Vector3D(31.643818, 71.563255, 36.873047);
		Vector3D l = new Vector3D(31.253489, 70.131218, 37.371773);

		assertEquals(-115.172063, Math.toDegrees(MathUtil.findDihedral(i, j, k, l)), diff);

	}

}
