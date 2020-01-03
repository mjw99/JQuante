package name.mjw.jquante.math;

import java.util.List;

import org.hipparchus.geometry.euclidean.threed.Vector3D;
import org.junit.jupiter.api.Assertions;
import org.junit.jupiter.api.Test;

public class LebedevRuleTest {

	private double hyperboloid(double x, double y, double z) {
		return (x * x) + (y * y) - (z * z);

	}

	@Test
	public void testBasicUsageExample() {
		final List<LebedevGridPoint> grid = LebedevRule.createGridByOrder(590);

		double result = 0.0;

		for (final LebedevGridPoint point : grid) {
			result += hyperboloid(point.getX(), point.getY(), point.getZ()) * point.getWeight();
		}

		Assertions.assertEquals(4.188790204786363, result * Math.PI * 4, 1E-14);

	}

	@Test
	public void testSphereIntegration() {

		for (int rule = 1; rule <= LebedevRule.MAX_RULE; rule++) {
			final List<LebedevGridPoint> grid = LebedevRule.createGridForRule(rule);

			if (grid == null) {
				Assertions.assertEquals(0, LebedevRule.isRuleAvailable(rule));
			} else {
				Assertions.assertEquals(1, LebedevRule.isRuleAvailable(rule));
				Assertions.assertEquals(LebedevRule.getRuleOrder(rule), grid.size());

				// integrate unitary function
				double result = 0.;
				Vector3D centroid = Vector3D.ZERO;
				for (final LebedevGridPoint point : grid) {
					result += point.getWeight();
					centroid = centroid.add(point.getVector3D());
				}
				// the integration of 1 should be 1 (the actual integral should be multiplied by
				// 4pi r^2)
				Assertions.assertEquals(1.0, result, 1E-13);
				// centroid is guaranteed to be zero due to symmetry of grid
				Assertions.assertEquals(0.0, centroid.getNorm(), 1E-14);
			}
		}

	}

	@Test
	public void testGetters() {
		final LebedevGridPoint point = LebedevRule.createGridForRule(1).get(0);
		Assertions.assertEquals(Math.PI / 2., point.getPhi(), 1E-15); // coelevation
		Assertions.assertEquals(0.0, point.getTheta(), 1E-15); // azimuth
		Assertions.assertEquals(1.0, point.getX(), 0.);
		Assertions.assertEquals(0.0, point.getY(), 0.);
		Assertions.assertEquals(0.0, point.getZ(), 0.);
	}

	@Test
	public void testLebedevBadOrder() {
		try {
			LebedevRule.createGridByOrder(0);
			Assertions.fail();
		} catch (final IllegalArgumentException e) {
			Assertions.assertTrue(true);
		}

		try {
			LebedevRule.getRuleOrder(0);
			Assertions.fail();
		} catch (final IllegalArgumentException e) {
			Assertions.assertTrue(true);
		}

		try {
			LebedevRule.getRuleOrder(1000);
			Assertions.fail();
		} catch (final IllegalArgumentException e) {
			Assertions.assertTrue(true);
		}

		try {
			LebedevRule.getRulePrecision(0);
			Assertions.fail();
		} catch (final IllegalArgumentException e) {
			Assertions.assertTrue(true);
		}

		try {
			LebedevRule.getRulePrecision(1000);
			Assertions.fail();
		} catch (final IllegalArgumentException e) {
			Assertions.assertTrue(true);
		}
	}

}
