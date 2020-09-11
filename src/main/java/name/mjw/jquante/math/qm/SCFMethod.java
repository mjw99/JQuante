package name.mjw.jquante.math.qm;

import java.util.ArrayList;

import name.mjw.jquante.common.EventListenerList;
import name.mjw.jquante.config.impl.AtomInfo;
import name.mjw.jquante.math.optimizer.OptimizerFunction;
import name.mjw.jquante.math.qm.event.SCFEvent;
import name.mjw.jquante.math.qm.event.SCFEventListener;
import name.mjw.jquante.molecule.Atom;
import name.mjw.jquante.molecule.Molecule;

/**
 * An abstract class representing the Self Consistent Field (SCF) method like
 * Hartree-Fock, MP2 etc.
 * 
 * @author V.Ganesh
 * @version 2.0 (Part of MeTA v2.0)
 */
public abstract class SCFMethod implements OptimizerFunction {

	// default constants
	private static final int MAX_ITERATION = 20;
	private static final double ENERGY_TOLERANCE = 1.0e-4;
	private static final double DENSITY_TOLERANCE = 1.0e-4;

	/**
	 * The one electron integrals of the system. Provides Hcore and Overlap
	 * matrix.
	 */
	protected OneElectronIntegrals oneEI;

	/**
	 * The two electron integrals of the system.
	 */
	protected TwoElectronIntegrals twoEI;

	/**
	 * The Density matrix.
	 */
	protected Density density;
	protected ArrayList<Density> densityList;

	/**
	 * The Molecular Orbitals (coefficient) matrix.
	 */
	protected MolecularOrbitals mos;
	protected ArrayList<MolecularOrbitals> mosList;

	/**
	 * Preamble to Fock matrix construction; contains the contribution from the
	 * two electron integrals and the density matrix.
	 */
	protected GMatrix gMatrix;
	protected ArrayList<GMatrix> gMatrixList;

	/**
	 * Holds value of property fock - the Fock matrix.
	 */
	protected Fock fock;
	protected ArrayList<Fock> fockList;

	/**
	 * Holds value of property energyTolerance.
	 */
	protected double energyTolerance;

	/**
	 * Holds value of property densityTolerance.
	 */
	protected double densityTolerance;

	/**
	 * The molecule under consideration
	 */
	protected Molecule molecule;

	/**
	 * Holds value of property scfIteration.
	 */
	protected int scfIteration;

	/**
	 * Holds value of property maxIteration.
	 */
	protected int maxIteration;

	/**
	 * Holds value of property guessInitialDM.
	 */
	protected boolean guessInitialDM;

	/**
	 * Holds value of property densityGuesser.
	 */
	protected DensityGuesser densityGuesser;

	/**
	 * Holds value of property energy.
	 */
	protected double energy;

	/**
	 * Utility field used by event firing mechanism.
	 */
	private EventListenerList<SCFEventListener> listenerList = null;

	/**
	 * Is this an open shell SCF method? 
	 */
	private final boolean openShell;

	/**
	 * Creates a new instance of SCFMethod (closed shell).
	 *
	 * @param molecule the Molecule object.
	 * @param oneEI    the 1E integral driver.
	 * @param twoEI    the 2E integral driver.
	 */
	public SCFMethod(final Molecule molecule, final OneElectronIntegrals oneEI, final TwoElectronIntegrals twoEI) {
		this(molecule, oneEI, twoEI, false);
	}

	/**
	 * Creates a new instance of SCFMethod.
	 *
	 * @param molecule  the Molecule object.
	 * @param oneEI     the 1E integral driver.
	 * @param twoEI     the 2E integral driver.
	 * @param openShell is this a closed shell setup?
	 */
	public SCFMethod(final Molecule molecule, final OneElectronIntegrals oneEI, final TwoElectronIntegrals twoEI,
			final boolean openShell) {
		this.maxIteration = MAX_ITERATION;

		this.energyTolerance = ENERGY_TOLERANCE;
		this.densityTolerance = DENSITY_TOLERANCE;

		this.molecule = molecule;
		this.oneEI = oneEI;
		this.twoEI = twoEI;

		this.guessInitialDM = false;
		this.openShell = openShell;

		if (this.openShell) {
			this.densityList = new ArrayList<>(2);
			this.mosList = new ArrayList<>(2);
			this.gMatrixList = new ArrayList<>(2);
			this.fockList = new ArrayList<>(2);
		}
	}

	/**
	 * Perform the SCF
	 */
	public abstract void scf();

	/**
	 * compute nuclear repulsion energy
	 * 
	 * @return the nuclear repulsion energy
	 */
	public double nuclearEnergy() {
		double eNuke = 0.0;
		int i;
		int j;
		final int noOfAtoms = molecule.getNumberOfAtoms();

		Atom atomI;
		Atom atomJ;

		// read in the atomic numbers
		final int[] atomicNumbers = new int[noOfAtoms];
		final AtomInfo ai = AtomInfo.getInstance();

		for (i = 0; i < noOfAtoms; i++) {
			atomicNumbers[i] = ai.getAtomicNumber(molecule.getAtom(i).getSymbol());
		}

		// compute nuclear energy
		for (i = 0; i < noOfAtoms; i++) {
			atomI = molecule.getAtom(i);
			for (j = 0; j < i; j++) {
				atomJ = molecule.getAtom(j);

				eNuke += atomicNumbers[i] * atomicNumbers[j]
						/ atomI.getAtomCenterInAU().distance(atomJ.getAtomCenterInAU());
			}
		}

		return eNuke;
	}

	/**
	 * Getter for property energyTolerance.
	 * 
	 * @return Value of property energyTolerance.
	 */
	public double getEnergyTolerance() {
		return this.energyTolerance;
	}

	/**
	 * Setter for property energyTolerance.
	 * 
	 * @param energyTolerance New value of property energyTolerance.
	 */
	public void setEnergyTolerance(final double energyTolerance) {
		this.energyTolerance = energyTolerance;
	}

	/**
	 * Getter for property densityTolerance.
	 * 
	 * @return Value of property densityTolerance.
	 */
	public double getDensityTolerance() {
		return this.densityTolerance;
	}

	/**
	 * Setter for property densityTolerance.
	 * 
	 * @param densityTolerance New value of property densityTolerance.
	 */
	public void setDensityTolerance(final double densityTolerance) {
		this.densityTolerance = densityTolerance;
	}

	/**
	 * Getter for property maxIteration.
	 * 
	 * @return Value of property maxIteration.
	 */
	public int getMaxIteration() {
		return this.maxIteration;
	}

	/**
	 * Setter for property maxIteration.
	 * 
	 * @param maxIteration New value of property maxIteration.
	 */
	public void setMaxIteration(final int maxIteration) {
		this.maxIteration = maxIteration;
	}

	/**
	 * Getter for property density.
	 * 
	 * @return Value of property density.
	 */
	public Density getDensity() {
		return this.density;
	}

	/**
	 * Getter for property density.
	 * 
	 * @param denNo - 0 is A, 1 is B in case of open shell systems
	 * 
	 * @return Value of property density.
	 */
	public Density getDensity(final int denNo) {
		if (this.openShell) {
			return this.densityList.get(denNo);
		} else {
			return this.density;
		}
	}

	/**
	 * Getter for property mos.
	 * 
	 * @return Value of property mos.
	 */
	public MolecularOrbitals getMos() {
		return this.mos;
	}

	/**
	 * Getter for property mos.
	 * 
	 * @param mosNo - 0 is A, 1 is B in case of open shell systems
	 * 
	 * @return Value of property mos.
	 */
	public MolecularOrbitals getMos(final int mosNo) {
		if (this.openShell) {
			return this.mosList.get(mosNo);
		} else {
			return this.mos;
		}
	}

	/**
	 * Getter for property orbE.
	 * 
	 * @return Value of property orbE.
	 */
	public double[] getOrbE() {
		return this.mos.getOrbitalEnergies();
	}

	/**
	 * Getter for property orbE.
	 * 
	 * @param mosNo - 0 is A, 1 is B in case of open shell systems
	 * 
	 * @return Value of property orbE.
	 */
	public double[] getOrbE(final int mosNo) {
		if (this.openShell) {
			return this.mosList.get(mosNo).getOrbitalEnergies();
		} else {
			return this.mos.getOrbitalEnergies();
		}
	}

	/**
	 * Getter for property scfIteration.
	 * 
	 * @return Value of property scfIteration.
	 */
	public int getScfIteration() {
		return this.scfIteration;
	}

	/**
	 * Getter for property molecule.
	 * 
	 * @return Value of property molecule.
	 */
	public Molecule getMolecule() {
		return this.molecule;
	}

	/**
	 * Getter for property oneEI.
	 * 
	 * @return Value of property oneEI.
	 */
	public OneElectronIntegrals getOneEI() {
		return this.oneEI;
	}

	/**
	 * Getter for property twoEI.
	 * 
	 * @return Value of property twoEI.
	 */
	public TwoElectronIntegrals getTwoEI() {
		return this.twoEI;
	}

	/**
	 * Getter for property densityGuesser.
	 * 
	 * @return Value of property densityGuesser.
	 */
	public DensityGuesser getDensityGuesser() {
		return this.densityGuesser;
	}

	/**
	 * Setter for property densityGuesser.
	 * 
	 * @param densityGuesser New value of property densityGuesser.
	 */
	public void setDensityGuesser(final DensityGuesser densityGuesser) {
		this.densityGuesser = densityGuesser;
	}

	/**
	 * Getter for property guessInitialDM.
	 * 
	 * @return Value of property guessInitialDM.
	 */
	public boolean isGuessInitialDM() {
		return this.guessInitialDM;
	}

	/**
	 * Setter for property guessInitialDM.
	 * 
	 * @param guessInitialDM New value of property guessInitialDM.
	 */
	public void setGuessInitialDM(final boolean guessInitialDM) {
		this.guessInitialDM = guessInitialDM;
	}

	/**
	 * Getter for property fock.
	 * 
	 * @return Value of property fock.
	 */
	public Fock getFock() {
		return this.fock;
	}

	/**
	 * Getter for property fock.
	 * 
	 * @param fckNo - 0 is A, 1 is B in case of open shell systems
	 * 
	 * @return Value of property fock.
	 */
	public Fock getFock(final int fckNo) {
		if (this.openShell) {
			return this.fockList.get(fckNo);
		} else {
			return this.fock;
		}
	}

	/**
	 * Get the value of gMatrix
	 * 
	 * @return the value of gMatrix
	 */
	public GMatrix getGMatrix() {
		return gMatrix;
	}

	/**
	 * Get the value of gMatrix
	 * 
	 * @param gmatNo - 0 is A, 1 is B in case of open shell systems
	 * 
	 * @return the value of gMatrix
	 */
	public GMatrix getGMatrix(final int gmatNo) {
		if (this.openShell) {
			return this.gMatrixList.get(gmatNo);
		} else {
			return gMatrix;
		}
	}

	/**
	 * Registers SCFEventListener to receive events.
	 * 
	 * @param listener The listener to register.
	 */
	public synchronized void addSCFEventListener(final SCFEventListener listener) {
		if (listenerList == null) {
			listenerList = new EventListenerList<>();
		}
		listenerList.add(SCFEventListener.class, listener);
	}

	/**
	 * Removes SCFEventListener from the list of listeners.
	 * 
	 * @param listener The listener to remove.
	 */
	public synchronized void removeSCFEventListener(final SCFEventListener listener) {
		listenerList.remove(SCFEventListener.class, listener);
	}

	/**
	 * Notifies all registered listeners about the event.
	 * 
	 * @param event The event to be fired
	 */
	protected void fireSCFEventListenerScfEventOccured(final SCFEvent event) {
		if (listenerList == null)
			return;
		final Object[] listeners = listenerList.getListenerList();

		for (final Object listener : listeners) {
			((SCFEventListener) listener).scfEventOccured(event);
		} // end for
	}

	/**
	 * Getter for property energy.
	 * 
	 * @return Value of property energy.
	 */
	public double getEnergy() {
		return this.energy;
	}

	/**
	 * A string representation return the current computed energy
	 * 
	 * @return the String representation of this object
	 */
	@Override
	public String toString() {
		return Double.toString(getEnergy());
	}
}
