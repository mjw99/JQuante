package name.mjw.jquante.math;

import org.apache.commons.math3.geometry.euclidean.threed.Vector3D;
import org.apache.commons.math3.linear.ArrayRealVector;
import org.apache.commons.math3.linear.RealVector;

import net.jafama.FastMath;

/**
 * Point information on a Lebedev grid.
 */
public final class LebedevGridPoint {

	/** x coordinate of the point in the unit sphere. */
	private final double x;
	/** y coordinate of the point in the unit sphere. */
	private final double y;
	/** z coordinate of the point in the unit sphere. */
	private final double z;
	/** Weight of the point. */
	private final double weight;

	/**
	 * Constructor (point must be in the unit sphere so x*x + y*x + z*z == 1).
	 * 
	 * @param xIn      x coordinate
	 * @param yIn      y coordinate
	 * @param zIn      z coordinate
	 * @param weightIn point weight
	 */
	public LebedevGridPoint(final double xIn, final double yIn, final double zIn, final double weightIn) {
		this.x = xIn;
		this.y = yIn;
		this.z = zIn;
		this.weight = weightIn;
	}

	/**
	 * Returns the grid point.
	 * 
	 * @return the grid point
	 */
	public final Vector3D getVector3D() {
		return new Vector3D(x, y, z);
	}

	/**
	 * Returns the grid point.
	 * 
	 * @return the grid point
	 */
	public final RealVector getXYZ() {
		return new ArrayRealVector(new double[] { x, y, z }, false);
	}

	/**
	 * Returns the grid point X.
	 * 
	 * @return the X grid point
	 */
	public final double getX() {
		return x;
	}

	/**
	 * Returns the grid point Y.
	 * 
	 * @return the Y grid point
	 */
	public final double getY() {
		return y;
	}

	/**
	 * Returns the grid point Z.
	 * 
	 * @return the Z grid point
	 */
	public final double getZ() {
		return z;
	}

	/**
	 * Returns the grid point's weight.
	 * 
	 * @return the grid point's weight
	 */
	public final double getWeight() {
		return weight;
	}

	/**
	 * Returns the grid point coelevation.
	 * 
	 * @return the grid point coelevation in [-PI,PI]
	 */
	public double getPhi() {
		return FastMath.acos(z);
	}

	/**
	 * Returns the grid point azimuth.
	 * 
	 * @return the grid point azimuth in [-PI,PI]
	 */
	public double getTheta() {
		return FastMath.atan2(y, x);
	}
}
