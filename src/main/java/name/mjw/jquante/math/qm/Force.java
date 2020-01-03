package name.mjw.jquante.math.qm;

import org.hipparchus.geometry.euclidean.threed.Vector3D;

/**
 * This interface defines the interfaces required for Force calculations on a
 * specific atomic position.
 * 
 * @author V.Ganesh
 * @version 2.0 (Part of MeTA v2.0)
 */
public interface Force {

	/**
	 * Compute the total force on the specified atom and return the results as a
	 * Vector3D object. Note that the SCF calculations should have been over
	 * before calling this method. This method in no way guarantees to check if
	 * the SCF was performed prior to calling this method.
	 * 
	 * @param atomIndex
	 *            the atom index for which the force is to be computed
	 * @param scfMethod
	 *            the instance of the SCF method
	 * @return the computed force
	 */
	public Vector3D computeForce(int atomIndex, SCFMethod scfMethod);
}
