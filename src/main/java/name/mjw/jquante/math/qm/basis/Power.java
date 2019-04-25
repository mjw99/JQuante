package name.mjw.jquante.math.qm.basis;

import java.util.Objects;

import net.jafama.FastMath;

/**
 * Represents the powers on orbitals. <br>
 * They are also the magnetic quantum numbers.
 * 
 * @author V.Ganesh
 */
public class Power {

	/**
	 * Orbital power l.
	 */
	private final int l;

	/**
	 * Orbital power m.
	 */
	private final int m;

	/**
	 * Orbital power n.
	 */
	private final int n;

	/**
	 * Creates a new instance of Power.
	 * 
	 * @param l the orbital power l
	 * @param m the orbital power m
	 * @param n the orbital power n
	 */
	public Power(final int l, final int m, final int n) {
		this.l = l;
		this.m = m;
		this.n = n;
	}

	/**
	 * Getter for property l.
	 * 
	 * @return orbital power l
	 */
	public int getL() {
		return this.l;
	}

	/**
	 * Getter for property m.
	 * 
	 * @return orbital power m.
	 */
	public int getM() {
		return this.m;
	}

	/**
	 * Getter for property n.
	 * 
	 * @return orbital power n.
	 */
	public int getN() {
		return this.n;
	}

	/**
	 * Return the maximum angular momentum of this Power
	 * 
	 * @return the maximum of the powers
	 */
	public int getMaximumAngularMomentum() {
		return FastMath.max(FastMath.max(l, m), n);
	}

	/**
	 * Return the minimum angular momentum of this Power
	 * 
	 * @return the minimum of the powers
	 */
	public int getMinimumAngularMomentum() {
		return FastMath.min(FastMath.min(l, m), n);
	}

	/**
	 * Return the total angular momentum of this Power
	 * 
	 * @return the maximum of the powers
	 */
	public int getTotalAngularMomentum() {
		return l + m + n;
	}

	/**
	 * Get the power as a dimension 3 array.
	 *
	 * @return powers
	 */
	public double[] toArray() {
		return new double[] { l, m, n };
	}

	/**
	 * overloaded toString()
	 */
	@Override
	public String toString() {
		return "[" + l + ", " + m + ", " + n + "]";
	}

	@Override
	public int hashCode() {
		return Objects.hash(l, m, n);
	}

	@Override
	public boolean equals(Object obj) {
		if (this == obj)
			return true;
		if (obj == null)
			return false;
		if (getClass() != obj.getClass())
			return false;
		Power other = (Power) obj;
		return l == other.l && m == other.m && n == other.n;
	}
}