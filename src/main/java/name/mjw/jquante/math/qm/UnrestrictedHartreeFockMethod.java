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
public class UnrestrictedHartreeFockMethod extends SCFMethod implements
		OptimizerFunction {

	/**
	 * Logger object
	 */
	private static final Logger LOG = LogManager
			.getLogger(UnrestrictedHartreeFockMethod.class);

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
	public UnrestrictedHartreeFockMethod(Molecule molecule,
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
	public UnrestrictedHartreeFockMethod(Molecule molecule,
			OneElectronIntegrals oneEI, TwoElectronIntegrals twoEI,
			SCFType scfType) {
		super(molecule, oneEI, twoEI, true);

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
		gMatrixList.set(0, new GMatrix(hCore.getRowDimension()));  // A
		gMatrixList.set(1, new GMatrix(hCore.getRowDimension()));  // B
		mosList.set(0, new MolecularOrbitals(hCore.getRowDimension()));
		mosList.set(1, new MolecularOrbitals(hCore.getRowDimension()));
		densityList.set(0, new Density(hCore.getRowDimension()));
		densityList.set(1, new Density(hCore.getRowDimension()));
		fockList.set(0, new Fock(hCore.getRowDimension()));
		fockList.set(1, new Fock(hCore.getRowDimension()));

		// TODO: compute initial MOs
		mosList.get(0).compute(hCore, overlap);
		LOG.debug("Initial computed MO coefficient matrix as: \n" + mosList.get(0));

		FockExtrapolator diis = new DIISFockExtrapolator();

		LOG.debug("Initial density matrix \n" + densityList.get(0));

		// start the SCF cycle
		for (scfIteration = 0; scfIteration < maxIteration; scfIteration++) {

			LOG.debug("");
			LOG.debug("SCF iteration: " + scfIteration);
            // make or guess density
            // TODO: will have two parts Da and Db
			densityList.get(0).compute(this, guessInitialDM && (scfIteration == 0),
					densityGuesser, noOfOccupancies, mosList.get(0));
			densityList.get(1).compute(this, guessInitialDM && (scfIteration == 0),
					densityGuesser, noOfOccupancies, mosList.get(1));

			LOG.debug("Density matrix:\n" + densityList);

            // make the G matrix
            // TODO: G will have two parts Ga = Ja+Jb-Ka and Gb = Ja+Jb-Kb 
			gMatrixList.get(0).compute(scfType, twoEI, densityList.get(0));
			gMatrixList.get(1).compute(scfType, twoEI, densityList.get(1));

            // make fock matrix
            // TODO: F will have two parts Fa = h+Ga and Fb = h+Gb
			fockList.get(0).compute(hCore, gMatrixList.get(0));
			fockList.get(1).compute(hCore, gMatrixList.get(1));

			// apply DIIS
			// TODO : diis needs to be rewritten to take care of Fa and Fb
			// fock = diis.next(fock, overlap, density);

            // compute the new MOs
            // TODO: two set of orbitals, mosA, mosB
			mosList.get(0).compute(fockList.get(0), overlap);
			mosList.get(1).compute(fockList.get(0), overlap);

            // compute the total energy at this point
            // TODO: eTwo (fock energy) will have two parts - this eqn needs to change
			eOne = densityList.get(0).multiply(hCore).getTrace();
			eTwo = densityList.get(0).multiply(fock).getTrace();

			energy = eOne + eTwo + nuclearEnergy;

			LOG.debug("Energy is : " + energy + "\tdelta_E: "
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
	 * Getter for property hessian.
	 * 
	 * @return Value of property hessian.
	 */
	@Override
	public RealMatrix getHessian() {
		throw new UnsupportedOperationException("Not supported yet.");
	}
}
