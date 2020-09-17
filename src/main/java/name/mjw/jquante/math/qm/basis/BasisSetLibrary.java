package name.mjw.jquante.math.qm.basis;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Iterator;
import java.util.List;
import java.util.Set;


import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

import com.google.common.collect.ImmutableSet;
import com.google.common.collect.Sets;

import name.mjw.jquante.molecule.Atom;
import name.mjw.jquante.molecule.Molecule;
import name.mjw.jquante.molecule.UserDefinedAtomProperty;

/**
 * Class to construct basis functions of a given molecule and a basis set
 * 
 * @author V.Ganesh
 * @author mw529
 */
public final class BasisSetLibrary {

	private static final Logger LOG = LogManager.getLogger(BasisSetLibrary.class);

	private String basisName;

	/**
	 * Holds value of property basisFunctions.
	 */
	private ArrayList<ContractedGaussian> basisFunctions;

	/**
	 * A shell is a set of basis functions sharing common orbital exponents and a
	 * common center, for example, sp shell or a d shell.
	 */
	private ArrayList<Shell> shells;

	private List<List<Shell>> uniqueShellPairs;


	/**
	 * Creates a new instance of BasisFunctions.
	 * 
	 * @param molecule  the Molecule whose basis function is requested.
	 * @param basisName the name of the basis set (like sto3g).
	 * @throws Exception the basisName was not found.
	 *
	 */
	public BasisSetLibrary(Molecule molecule, String basisName) throws Exception {
		// initialise the basis functions
		getBasisFunctions(molecule, basisName);
		this.basisName = basisName;

		// and initialise the shell list
		initShellList();

		// Shell Pair list
		initUniqueShellPairList();
	}

	/**
	 * Get the value of basisName
	 * 
	 * @return the value of basisName
	 */
	public String getBasisName() {
		return basisName;
	}

	/**
	 * Set the value of basisName
	 * 
	 * @param basisName new value of basisName
	 */
	public void setBasisName(String basisName) {
		this.basisName = basisName;
	}

	/**
	 * Getter for property basisFunctions.
	 * 
	 * @return Value of property basisFunctions.
	 */
	public List<ContractedGaussian> getBasisFunctions() {
		return this.basisFunctions;
	}

	/**
	 * Returns a set of basis function with the same centre and same contracted
	 * exponent.
	 *
	 * @return A shell Multimap for the current molecule and basis.
	 */
	public List<Shell> getShells() {
		return shells;
	}

	public List<List<Shell>> getUniqueShellPairs() {
		return uniqueShellPairs;
	}

	/**
	 * Getter for property basisFunctions.
	 * 
	 * @param molecule  the Molecule whose basis function is requested
	 * @param basisName the name of the basis set (like sto3g)
	 * @return Value of property basisFunctions.
	 */
	private ArrayList<ContractedGaussian> getBasisFunctions(Molecule molecule, String basisName) throws Exception {
		BasisSet basisSet = BasisSetReader.getInstance().readBasisSet(basisName);
		Iterator<Atom> atoms = molecule.getAtoms();

		basisFunctions = new ArrayList<>();

		Atom atom;
		AtomicBasis atomicBasis;
		while (atoms.hasNext()) { // loop over atoms
			atom = atoms.next();
			atomicBasis = basisSet.getAtomicBasis(atom.getSymbol());

			Iterator<Orbital> orbitals = atomicBasis.getOrbitals().iterator();
			Orbital orbital;
			ArrayList<ContractedGaussian> atomicFunctions = new ArrayList<>();

			while (orbitals.hasNext()) { // loop over atom orbitals
				orbital = orbitals.next();

				Iterator<Power> pList = PowerList.getInstance().getPowerList(orbital.getType());
				Power power;
				while (pList.hasNext()) { // and the power list, sp2 etc..
					power = pList.next();

					ContractedGaussian cg = new ContractedGaussian(atom, power);
					Iterator<Double> coeff = orbital.getCoefficients().iterator();
					Iterator<Double> exp = orbital.getExponents().iterator();

					while (coeff.hasNext()) { // build the CG from PGs
						cg.addPrimitive(exp.next().doubleValue(), coeff.next().doubleValue());
					}

					cg.normalize();
					cg.setBasisFunctionIndex(basisFunctions.size()); // set an index
					basisFunctions.add(cg); // add this CG to list
					atomicFunctions.add(cg); // add the reference to atom list
				}
			}

			// save a reference of the basis functions centered on
			// this atom as a user defined property of the atom
			try {
				atom.addUserDefinedAtomProperty(new UserDefinedAtomProperty("basisFunctions", atomicFunctions));
			} catch (UnsupportedOperationException e) {
				UserDefinedAtomProperty up = atom.getUserDefinedAtomProperty("basisFunctions");
				up.setValue(atomicFunctions);
			}
		}

		Collections.sort(basisFunctions);
		basisFunctions.trimToSize();
		return this.basisFunctions;
	}

	/**
	 * Initialise the shell list
	 */
	private void initShellList() {
		shells = new ArrayList<>();

		// First entry always goes in.
		Shell shell = new Shell(basisFunctions.get(0));
		shell.setFirstBasisFunctionIndex(0);
		shells.add(shell);


		for (int i = 1; i < basisFunctions.size(); i++) {
			ContractedGaussian contractedGaussian = basisFunctions.get(i);

			Shell tmpShell = new Shell(contractedGaussian);
			tmpShell.setFirstBasisFunctionIndex(i);
			tmpShell.setLastBasisFunctionIndex(i);

			// Increment lastBasisFunctionIndex if we have seen it before
			if (shells.contains(tmpShell)) {
				int index = shells.indexOf(tmpShell);
				shells.get(index).setLastBasisFunctionIndex(i);
			}

			if (!(shells.contains(tmpShell))) {
				shells.add(tmpShell);
			}

		}
		shells.trimToSize();

	}

	/**
	 * This is essentially itertools.combinations_with_replacement(shells, 2)
	 * <p>
	 * 
	 * https://docs.python.org/3/library/itertools.html#itertools.combinations_with_replacement
	 * 
	 * <p>
	 * itertools.combinations_with_replacement does not exist in Guava, but
	 * Sets.combinations can be augmented slightly to obtain the same result.
	 * <p>
	 * Overlall it is 
	 * {nShells}^C_{2} + nShells
	 * <p>
	 * where the {n}^C_{r} notation can be found at https://en.wikipedia.org/wiki/Combination
	 *
	 * <p>
	 * Also see: https://www.baeldung.com/java-combinations-algorithm
	 * 
	 */
	private void initUniqueShellPairList() {

		// Special case if shells only has one element
		if (shells.size() == 1) {
			uniqueShellPairs = new ArrayList<>();

			for (Shell shell : shells) {
				ArrayList<Shell> pair = new ArrayList<>();
				pair.add(shell);
				pair.add(shell);
				uniqueShellPairs.add(pair);
			}
			return;
		}

		// First, calculate {nShells}^C_{2} combinations
		Set<Set<Shell>> tmpShellPairs = Sets.combinations(ImmutableSet.copyOf(shells), 2);

		uniqueShellPairs = new ArrayList<>();
		for (Set<Shell> item : tmpShellPairs) {
			ArrayList<Shell> pair = new ArrayList<>(item);
			uniqueShellPairs.add(pair);
		}

		// Second, add in missing self pairs
		for (Shell shell : shells) {
			ArrayList<Shell> pair = new ArrayList<>();
			pair.add(shell);
			pair.add(shell);
			uniqueShellPairs.add(pair);
		}

	}

	public void printBasisFunctionList() {
		System.out.println("");
		System.out.println("Basis function list");
		System.out.println("===================");
		for (ContractedGaussian bfs : this.getBasisFunctions()) {
			System.out.println(String.format("%2d", bfs.getBasisFunctionIndex()) + " " + bfs.getCenteredAtom() + " "
					+ bfs.getPowers() + " " + bfs.getExponents() + " " + bfs.getCoefficients());

		}
	}

	public void printUniqueShellPairList() {
		System.out.println("");
		System.out.println("Unique shellpair list (" + this.uniqueShellPairs.size() + " pairs)");
		System.out.println("=================================");
		for (List<Shell> uniqueShellPair : this.uniqueShellPairs) {
			System.out.print(
					(shells.indexOf(uniqueShellPair.get(0)) + 1) + " " + (shells.indexOf(uniqueShellPair.get(1)) + 1)+ "\n"
					//+ uniqueShellPair.get(0)
					//+ uniqueShellPair.get(1) + "\n"
					);					
			
		}
	}

	public void printShellList() {
		System.out.println("");
		System.out.println("Shell list");
		System.out.println("==========");
		for (Shell shell : this.shells) {
			System.out.print(shell);
		}
	}

}
