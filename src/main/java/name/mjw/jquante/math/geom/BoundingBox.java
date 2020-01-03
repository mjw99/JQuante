package name.mjw.jquante.math.geom;

import org.hipparchus.geometry.euclidean.threed.Vector3D;

import net.jafama.FastMath;

/**
 * A class representing a rectangular bonding box.
 * 
 * @author V.Ganesh
 * @version 2.0 (Part of MeTA v2.0)
 */
public class BoundingBox extends AbstractGeometricObject implements Cloneable {

	/**
	 * Eclipse generated serialVersionUID
	 */
	private static final long serialVersionUID = 2613593470614954917L;

	/** Holds value of property upperLeft. */
	private Vector3D upperLeft;

	/** Holds value of property bottomRight. */
	private Vector3D bottomRight;

	/** Creates a new instance of BoundingBox */
	public BoundingBox() {
		this(new Vector3D(0, 0, 0), new Vector3D(0, 0, 0));
	}

	public BoundingBox(Vector3D upperLeft, Vector3D bottomRight) {
		this.upperLeft = upperLeft;
		this.bottomRight = bottomRight;
	}

	/**
	 * Getter for property upperLeft.
	 * 
	 * @return Value of property upperLeft.
	 * 
	 */
	public Vector3D getUpperLeft() {
		return this.upperLeft;
	}

	/**
	 * Setter for property upperLeft.
	 * 
	 * @param upperLeft
	 *            New value of property upperLeft.
	 * 
	 */
	public void setUpperLeft(Vector3D upperLeft) {
		this.upperLeft = upperLeft;
	}

	/**
	 * Getter for property bottomRight.
	 * 
	 * @return Value of property bottomRight.
	 * 
	 */
	public Vector3D getBottomRight() {
		return this.bottomRight;
	}

	/**
	 * Setter for property bottomRight.
	 * 
	 * @param bottomRight
	 *            New value of property bottomRight.
	 * 
	 */
	public void setBottomRight(Vector3D bottomRight) {
		this.bottomRight = bottomRight;
	}

	/**
	 * return the width in X direction
	 * 
	 * @return width in X direction
	 */
	public double getXWidth() {
		return FastMath.abs(bottomRight.getX() - upperLeft.getX());
	}

	/**
	 * return the width in Y direction
	 * 
	 * @return width in Y direction
	 */
	public double getYWidth() {
		return FastMath.abs(bottomRight.getY() - upperLeft.getY());
	}

	/**
	 * return the width in Z direction
	 * 
	 * @return width in Z direction
	 */
	public double getZWidth() {
		return FastMath.abs(bottomRight.getZ() - upperLeft.getZ());
	}

	/**
	 * return the center of this box, defined as the midpoint of the line joining
	 * the points that define this BoundingBox
	 * 
	 * @return Point3D - the center
	 */
	public Vector3D center() {
		return new Vector3D((upperLeft.getX() + bottomRight.getX()) / 2.0,
				(upperLeft.getY() + bottomRight.getY()) / 2.0, (upperLeft.getZ() + bottomRight.getZ()) / 2.0);
	}

	/**
	 * Return all coordinates of all the 8 corner points associated with this
	 * bounding box. The order of points is: (upperLeft, the rest of the three
	 * points in clock wise order) (bottomRight, the rest of the three points in
	 * anti-clock wise order).
	 * 
	 * The size of the returned array is always 8, representing 8 corners of this
	 * bounding box.
	 * 
	 * @return the corner points
	 */
	public Vector3D[] corners() {
		Vector3D[] corners = new Vector3D[8];

		corners[0] = new Vector3D(upperLeft.toArray());

		corners[1] = new Vector3D(upperLeft.getX() + getXWidth(), upperLeft.getY(), upperLeft.getZ());
		corners[2] = new Vector3D(upperLeft.getX() + getXWidth(), upperLeft.getY() + getYWidth(), upperLeft.getZ());
		corners[3] = new Vector3D(upperLeft.getX(), upperLeft.getY() + getYWidth(), upperLeft.getZ());

		corners[4] = new Vector3D(bottomRight.toArray());

		corners[5] = new Vector3D(bottomRight.getX() - getXWidth(), bottomRight.getY(), bottomRight.getZ());
		corners[6] = new Vector3D(bottomRight.getX() - getXWidth(), bottomRight.getY() - getYWidth(),
				bottomRight.getZ());
		corners[7] = new Vector3D(bottomRight.getX(), bottomRight.getY() - getYWidth(), bottomRight.getZ());

		return corners;
	}

	/**
	 * Shift this bounding box to a new center
	 * 
	 * @param center
	 *            the new center
	 * @return a new instance of the shifted bounding box
	 */
	public BoundingBox shiftTo(Vector3D center) {
		Vector3D ul = new Vector3D((2.0 * center.getX() - getXWidth()) / 2.0, (2.0 * center.getY() - getYWidth()) / 2.0,
				(2.0 * center.getZ() - getZWidth() / 2.0));

		Vector3D br = new Vector3D(ul.getX() + getXWidth(), ul.getY() + getYWidth(), ul.getZ() + getZWidth());

		return new BoundingBox(ul, br);
	}

	/**
	 * volume of this box given by: <br>
	 * <code>V = a b c; where a, b and c are lengths along x, y and z 
	 *       respectively.</code>
	 * 
	 * @return the volume of this box
	 */
	@Override
	public double volume() {
		return (getXWidth() * getYWidth() * getZWidth());
	}

	/**
	 * total surface area of this box given by: <br>
	 * <code>S = 2(ab + bc + ca); where a, b and c are lengths along x, y and z 
	 *       respectively.</code>
	 * 
	 * @return the volume of this box
	 */
	@Override
	public double totalSurfaceArea() {
		return (2.0 * ((getXWidth() * getYWidth()) + (getYWidth() * getZWidth()) + (getZWidth() * getXWidth())));
	}

	/**
	 * combine 2 bounding boxes and return a third bigger one, enclosing both of
	 * them.
	 * 
	 * @param box
	 *            - the second box to be considered
	 * @return the bigger BoundingBox enclosing the two boxes
	 */
	public BoundingBox combine(BoundingBox box) {
		Vector3D ul = new Vector3D(FastMath.min(upperLeft.getX(), box.upperLeft.getX()),
				FastMath.min(upperLeft.getY(), box.upperLeft.getY()), FastMath.min(upperLeft.getZ(), box.upperLeft.getZ()));

		Vector3D br = new Vector3D(FastMath.max(bottomRight.getX(), box.bottomRight.getX()),
				FastMath.max(bottomRight.getY(), box.bottomRight.getY()),
				FastMath.max(bottomRight.getZ(), box.bottomRight.getZ()));

		return new BoundingBox(ul, br);

	}

	/**
	 * intersect 2 bounding boxes and return a third smaller one, enclosing the
	 * common part among both of them.
	 * 
	 * @param box
	 *            - the second box to be considered
	 * @return the bigger BoundingBox enclosing the two boxes
	 */
	public BoundingBox intersect(BoundingBox box) {
		BoundingBox bb = new BoundingBox();
		Vector3D bul = box.getUpperLeft();
		Vector3D bbr = box.getBottomRight();

		Vector3D ul = new Vector3D(FastMath.max(upperLeft.getX(), bul.getX()), FastMath.max(upperLeft.getY(), bul.getY()),
				FastMath.max(upperLeft.getZ(), bul.getZ()));
		bb.setUpperLeft(ul);

		Vector3D br = new Vector3D(FastMath.min(bottomRight.getX(), bbr.getX()), FastMath.min(bottomRight.getY(), bbr.getY()),
				FastMath.min(bottomRight.getZ(), bbr.getZ()));
		bb.setBottomRight(br);

		return bb;
	}

	/**
	 * Check whether a point is contained within this bounding box.
	 * 
	 * @param point
	 *            - the point to be checked
	 * @return boolen - returns true if the point is within the bounding box false
	 *         otherwise
	 */
	public boolean contains(Vector3D point) {
		if ((point.getX() >= upperLeft.getX()) && (point.getY() >= upperLeft.getY())
				&& (point.getZ() >= upperLeft.getZ()) && (point.getX() <= bottomRight.getX())
				&& (point.getY() <= bottomRight.getY()) && (point.getZ() <= bottomRight.getZ())) {
			return true;
		} else {
			return false;
		}
	}

	/**
	 * Expand this box along each side by a constant amount
	 * 
	 * @param e
	 *            the expansion parameter
	 * @return the new bounding box
	 */
	public BoundingBox expand(double e) {
		return expand(e, e, e);
	}

	/**
	 * Expand this bounding box with variable lengths in each direction
	 * 
	 * @param ex
	 *            expansion length in X direction
	 * @param ey
	 *            expansion length in Y direction
	 * @param ez
	 *            expansion length in Z direction
	 * @return the new bounding box
	 */
	public BoundingBox expand(double ex, double ey, double ez) {
		return new BoundingBox(upperLeft.add(new Vector3D(-ex, -ey, -ez)), bottomRight.add(new Vector3D(ex, ey, ez)));
	}

	/**
	 * overridden toString()
	 */
	@Override
	public String toString() {
		return ("UpperLeft: " + upperLeft + ", BottomRight: " + bottomRight);
	}

	/**
	 * i am clonable
	 */
	@Override
	public Object clone() throws CloneNotSupportedException {
		return (new BoundingBox(new Vector3D(upperLeft.toArray()), new Vector3D(bottomRight.toArray())));
	}

}