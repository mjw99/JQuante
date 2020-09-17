package name.mjw.jquante.math.qm;

import org.hipparchus.geometry.euclidean.threed.Vector3D;
import org.hipparchus.linear.RealMatrix;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

import name.mjw.jquante.math.optimizer.OptimizerFunction;
import name.mjw.jquante.math.qm.event.SCFEvent;
import name.mjw.jquante.molecule.Atom;
import name.mjw.jquante.molecule.Molecule;
import name.mjw.jquante.molecule.UserDefinedAtomProperty;
import net.jafama.FastMath;

/**
 * Implements the Hartree-Fock (HF) SCF method for single point energy
 * evaluation of a molecule.
 * 
 * @author V.Ganesh
 * @version 2.0 (Part of MeTA v2.0)
 * 
 *
 * @see <a href="https://en.wikipedia.org/wiki/Hartree%E2%80%93Fock_method">
 *      https://en.wikipedia.org/wiki/Hartree%E2%80%93Fock_method</a>
 * @see <a
 *      href="http://sirius.chem.vt.edu/~crawdad/programming/scf.pdf">http://sirius.chem.vt.edu/~crawdad/programming/scf.pdf</a>
 * @see <a href=
 *        "http://sirius.chem.vt.edu/wiki/doku.php?id=crawdad:programming:project3"
 *        >http://sirius.chem.vt.edu/wiki/doku.php?id=crawdad:programming:
 *        project3</a >
 * @see <a href="http://www.theoretical-physics.net/dev/quantum/hf.html">
 *      http://www.theoretical-physics.net/dev/quantum/hf.html</a>
 */
public class RestrictedHartreeFockMethod extends SCFMethod implements
		OptimizerFunction {

	/**
	 * Logger object
	 */
	private static final Logger LOG = LogManager
			.getLogger(RestrictedHartreeFockMethod.class);

	/**
	 * Represents an event in an SCF cycle
	 */
	protected SCFEvent scfEvent;

	/**
	 * SCF Type enumeration.
	 */
	private SCFType scfType;

	private boolean isDerivativeComputed = false;

	/**
	 * Creates a new instance of HartreeFockSCFMethod
	 * 
	 * @param molecule
	 *            The molecule under consideration.
	 * @param oneEI
	 *            The one electron integrals of the system
	 * @param twoEI
	 *            The two electron integrals of the system
	 */
	public RestrictedHartreeFockMethod(Molecule molecule,
			OneElectronIntegrals oneEI, TwoElectronIntegrals twoEI) {
		this(molecule, oneEI, twoEI, SCFType.HARTREE_FOCK);
	}

	/**
	 * Creates a new instance of HartreeFockSCFMethod
	 * 
	 * @param molecule
	 *            The molecule under consideration.
	 * @param oneEI
	 *            The one electron integrals of the system
	 * @param twoEI
	 *            The two electron integrals of the system
	 * @param scfType
	 *            Type of SCF Calculation.
	 */
	public RestrictedHartreeFockMethod(Molecule molecule,
			OneElectronIntegrals oneEI, TwoElectronIntegrals twoEI,
			SCFType scfType) {
		super(molecule, oneEI, twoEI);

		scfEvent = new SCFEvent(this);
		this.scfType = scfType;
	}

	/**
	 * Perform the SCF optimization of the molecular wave function until the
	 * energy converges.
	 */
	@Override
	public void scf() {
		// check first if closed shell run?
		int noOfElectrons = molecule.getNumberOfElectrons();
		int noOfOccupancies = noOfElectrons / 2;

		if (noOfElectrons % 2 != 0) {
			throw new UnsupportedOperationException("Open shell systems are"
					+ " not currently supported.");
		}

		Overlap overlap = oneEI.getOverlap();
		LOG.debug("Initial S matrix\n" + overlap);

		LOG.debug("S^-1/2 matrix\n" + overlap.getSHalf());

		HCore hCore = oneEI.getHCore();
		LOG.debug("Initial hCore\n" + hCore);

		boolean converged = false;
		double oldEnergy = 0.0;
		double nuclearEnergy = nuclearEnergy();
		double eOne;
		double eTwo;

		// init memory for the matrices
		gMatrix = new GMatrix(hCore.getRowDimension());
		mos = new MolecularOrbitals(hCore.getRowDimension());
		density = new Density(hCore.getRowDimension());
		fock = new Fock(hCore.getRowDimension());

		// compute initial MOs
		mos.compute(hCore, overlap);
		LOG.debug("Initial computed MO coefficient matrix as: \n" + mos);

		FockExtrapolator diis = new DIISFockExtrapolator();

		LOG.debug("Initial density matrix \n" + density);

		// start the SCF cycle
		for (scfIteration = 0; scfIteration < maxIteration; scfIteration++) {

			// make or guess density
			density.compute(this, guessInitialDM && (scfIteration == 0),
					densityGuesser, noOfOccupancies, mos);

			// make the G matrix
			gMatrix.compute(scfType, twoEI, density);

			// make fock matrix
			fock.compute(hCore, gMatrix);

			// apply DIIS
			fock = diis.next(fock, overlap, density);

			// compute the new MOs
			mos.compute(fock, overlap);

			// compute the total energy at this point
			eOne = density.multiply(hCore).getTrace();
			eTwo = density.multiply(fock).getTrace();

			energy = eOne + eTwo + nuclearEnergy;

			LOG.debug("SCF iteration: " + scfIteration + "\t Energy is : " + energy + "\tdelta_E: "
					+ (energy - oldEnergy));

			// fire the SCF event notification
			scfEvent.setType(SCFEvent.INFO_EVENT);
			scfEvent.setCurrentIteration(scfIteration);
			scfEvent.setCurrentEnergy(energy);
			fireSCFEventListenerScfEventOccured(scfEvent);

			// check for convergence
			if (FastMath.abs(energy - oldEnergy) < energyTolerance) {
				converged = true;
				scfEvent.setType(SCFEvent.CONVERGED_EVENT);
				scfEvent.setCurrentIteration(scfIteration);
				scfEvent.setCurrentEnergy(energy);
				fireSCFEventListenerScfEventOccured(scfEvent);
				break;
			}

			oldEnergy = energy;
		} // end of SCF iteration

		// not converged? then inform so...
		if (!converged) {
			scfEvent.setType(SCFEvent.FAILED_CONVERGENCE_EVENT);
			scfEvent.setCurrentIteration(scfIteration);
			scfEvent.setCurrentEnergy(energy);
			fireSCFEventListenerScfEventOccured(scfEvent);
		}
	}

	/**
	 * This gradient (or Force) calculation is based on Appendix C of Modern
	 * Quantum Chemistry by Szabo and Ostland, which describes computing
	 * analytic gradients and geometry optimization.
	 */
	private void computeForce() {
		Force hfForce = new HartreeFockForce();

		for (int i = 0; i < molecule.getNumberOfAtoms(); i++) {
			Atom atom = molecule.getAtom(i);
			Vector3D force = hfForce.computeForce(atom.getIndex(), this);

			UserDefinedAtomProperty atmForce = atom
					.getUserDefinedAtomProperty("force");
			if (atmForce == null) {
				atmForce = new UserDefinedAtomProperty("force", force);
				atom.addUserDefinedAtomProperty(atmForce);
			}

			atmForce.setValue(force);
		}

		isDerivativeComputed = true;
	}

	/**
	 * Evaluate the function with 'n' variables
	 * 
	 * @param variables
	 *            an array of variables
	 * @return a double value evaluating F(x1, x2, ...)
	 */
	@Override
	public double evaluate(double[] variables) {
		molecule.resetAtomCoordinates(variables, false);

		// perform scf, and return energy
		oneEI.compute1E();
		twoEI.compute2E();
		scf();
		isDerivativeComputed = false;
		return getEnergy();
	}

	/**
	 * For a method that encapsulates its own set of "base" variables (e.g. an
	 * initial set of atom positions), resets the base variables to new values.
	 * 
	 * @param variables
	 *            the new set of base variables
	 */
	@Override
	public void resetVariables(double[] variables) {
		throw new UnsupportedOperationException("Not supported yet.");
	}

	/**
	 * Getter for property derivativeAvailable.
	 * 
	 * @return Value of property derivativeAvailable.
	 */
	@Override
	public boolean isDerivativeAvailable() {
		throw new UnsupportedOperationException("Not supported yet.");
	}

	/**
	 * Getter for property hessianAvailable.
	 * 
	 * @return Value of property hessianAvailable.
	 */
	@Override
	public boolean isHessianAvailable() {
		throw new UnsupportedOperationException("Not supported yet.");
	}

	/**
	 * Getter for property derivatives.
	 * 
	 * @return Value of property derivatives.
	 */
	@Override
	public double[] getDerivatives() {
		if (!isDerivativeComputed) {
			computeForce();
		}

		// unpack and return the forces
		double[] forces = new double[molecule.getNumberOfAtoms() * 3];
		int ii = 0;
		for (int i = 0; i < molecule.getNumberOfAtoms(); i++) {
			UserDefinedAtomProperty atmForce = molecule.getAtom(i)
					.getUserDefinedAtomProperty("force");

			Vector3D force = (Vector3D) atmForce.getValue();
			forces[ii] = force.getX();
			forces[ii + 1] = force.getY();
			forces[ii + 2] = force.getZ();

			ii += 3;
		}

		return forces;
	}

	/**
	 * Return the max norm of the derivatives, if available
	 * 
	 * @return max norm of all the derivatives
	 */
	@Override
	public double getMaxNormOfDerivatives() {
		throw new UnsupportedOperationException("Not supported yet.");
	}

	/**
	 * Return the RMS of the derivatives, if available
	 * 
	 * @return max RMS all the derivatives
	 */
	@Override
	public double getRMSOfDerivatives() {
		throw new UnsupportedOperationException("Not supported yet.");
	}

	/**
	 * Getter for property Hessian.
	 * 
	 * @return Value of property Hessian.
	 */
	@Override
	public RealMatrix getHessian() {
		throw new UnsupportedOperationException("Not supported yet.");
	}
}
