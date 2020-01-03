package name.mjw.jquante.math.qm;

import java.util.List;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

import name.mjw.jquante.config.impl.AtomInfo;
import org.hipparchus.linear.RealMatrix;
import org.hipparchus.linear.Array2DRowRealMatrix;
import org.hipparchus.geometry.euclidean.threed.Vector3D;
import name.mjw.jquante.molecule.Atom;
import name.mjw.jquante.molecule.Molecule;
import net.jafama.FastMath;

/**
 * Force calculations for the Hartree-Fock method. This gradient (or Force)
 * calculation is based on Appendix C of Modern Quantum Chemistry by Szabo and
 * Ostland, which describes computing analytic gradients and geometry
 * optimisation.
 * 
 * @author V.Ganesh
 * @version 2.0 (Part of MeTA v2.0)
 */
public class HartreeFockForce implements Force {

	private static final Logger LOG = LogManager.getLogger(HartreeFockForce.class);

	private int atomIndex;
	private SCFMethod scfMethod;

	/**
	 * Compute the total force on the specified atom and return the results as a
	 * Vector3D object. Note that the SCF calculations should have been over before
	 * calling this method. This method in no way guarantees to check if the SCF was
	 * performed prior to calling this method.
	 * 
	 * @param atomIndex
	 *            the atom index for which the force is to be computed
	 * @param scfMethod
	 *            the instance of the SCF method
	 * @return the computed force
	 */
	@Override
	public Vector3D computeForce(int atomIndex, SCFMethod scfMethod) {
		this.atomIndex = atomIndex;
		this.scfMethod = scfMethod;

		Vector3D force = computeNuclearDerivative();
		LOG.debug("Nuclear /dr: " + force);
		force = force.add(computeDensityMatrixDerivative().scalarMultiply(2.0));
		LOG.debug("Density /dr: " + force);
		force = force.add(computeOneElectronDerivative().scalarMultiply(2.0));
		LOG.debug("1E /dr: " + force);
		force = force.add(computeTwoElectronDerivative());
		LOG.debug("2E /dr: " + force);

		return force.scalarMultiply(-1);
	}

	/** Compute the nuclear derivative contribution */
	private Vector3D computeNuclearDerivative() {
		Molecule mol = scfMethod.getMolecule();
		double nDer;

		AtomInfo ai = AtomInfo.getInstance();

		Atom a = mol.getAtom(atomIndex);
		Atom b;
		double aCharge = ai.getAtomicNumber(a.getSymbol());
		Vector3D aCenter = a.getAtomCenterInAU() ;
		double ndx = 0.0;
		double ndy = 0.0;
		double ndz = 0.0;

		for (int i = 0; i < mol.getNumberOfAtoms(); i++) {
			if (i != atomIndex) {
				b = mol.getAtom(i);
				Vector3D bCenter = b.getAtomCenterInAU();

				nDer = aCharge * ai.getAtomicNumber(b.getSymbol())
						/ FastMath.pow(bCenter.distanceSq(aCenter), 1.5);

				ndx += (nDer * (bCenter.getX() - aCenter.getX()));
				ndy += (nDer * (bCenter.getY() - aCenter.getY()));
				ndz += (nDer * (bCenter.getZ() - aCenter.getZ()));
			} // end if
		} // end for

		return new Vector3D(ndx, ndy, ndz);
	}

	/** Compute the one electron derivative contribution */
	private Vector3D computeOneElectronDerivative() {
		HCore hCore = scfMethod.getOneEI().getHCore();
		Density dens = scfMethod.getDensity();
		List<HCore> hCoreDer = hCore.computeDerivative(atomIndex, scfMethod);

		return new Vector3D(dens.multiply(hCoreDer.get(0)).getTrace(),
				dens.multiply(hCoreDer.get(1)).getTrace(),
				dens.multiply(hCoreDer.get(2)).getTrace());
	}

	/** Compute the derivative contribution from density matrix */
	private Vector3D computeDensityMatrixDerivative() {
		Overlap overlap = scfMethod.getOneEI().getOverlap();
		Density dens = scfMethod.getDensity();
		List<Overlap> overlapDer = overlap.computeDerivative(atomIndex, scfMethod);
		RealMatrix eMat = new Array2DRowRealMatrix(scfMethod.getOrbE());
		RealMatrix qMat = dens.multiply(eMat.multiply(dens.transpose()));

		return new Vector3D(qMat.multiply(overlapDer.get(0)).getTrace(),
				qMat.multiply(overlapDer.get(1)).getTrace(),
				qMat.multiply(overlapDer.get(2)).getTrace());
	}

	/** Compute the two electron derivative contribution */
	private Vector3D computeTwoElectronDerivative() {
		List<GMatrix> gDer = scfMethod.getGMatrix().computeDerivative(atomIndex, scfMethod);
		Density density = scfMethod.getDensity();

		return new Vector3D(density.multiply(gDer.get(0)).getTrace(),
				density.multiply(gDer.get(1)).getTrace(),
				density.multiply(gDer.get(2)).getTrace());
	}
}
