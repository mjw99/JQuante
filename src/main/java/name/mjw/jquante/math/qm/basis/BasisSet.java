package name.mjw.jquante.math.qm.basis;

import java.util.HashMap;

/**
 * Represents an entire basis set (say sto-3g ... etc.)
 * 
 * @author V.Ganesh
 * @version 2.0 (Part of MeTA v2.0)
 */
public final class BasisSet {

	/**
	 * The name of the basis set.
	 */
	private String name;

	/**
	 * Holds the collection of atomic basis, with the symbol as the key
	 */
	private HashMap<String, AtomicBasis> atomicBasisSet;

	/**
	 * Creates a new instance of Basis
	 * 
	 * @param name
	 *            - the name of this basis set (say "sto-3g")
	 */
	public BasisSet(String name) {
		this.name = name;

		atomicBasisSet = new HashMap<>();
	}

	/**
	 * Getter for property name.
	 * 
	 * @return Value of property name.
	 */
	public String getName() {
		return this.name;
	}

	/**
	 * Add a relevant atomic basis to this basis set
	 * 
	 * @param atomicBasis
	 *            the instance of AtomicBasis to be added to this basis set
	 */
	public void addAtomicBasis(AtomicBasis atomicBasis) {
		atomicBasisSet.put(atomicBasis.getChemicalSymbol(), atomicBasis);
	}

	/**
	 * Returns appropriate basis for a given chemical symbol. Will throw
	 * <code>BasisNotFoundException</code> if the basis set does not contain
	 * atomic basis for the requested atomic symbol.
	 * 
	 * @param chemicalSymbol
	 *            for which the basis is requested
	 * @return instance of AtomicBasis
	 */
	public AtomicBasis getAtomicBasis(String chemicalSymbol) {
		AtomicBasis atomicBasis = atomicBasisSet.get(chemicalSymbol);

		if (atomicBasis == null) {
			throw new BasisNotFoundException("Basis for atom '" + chemicalSymbol
					+ "' is not defined in : " + name);
		}
		return atomicBasis;
	}
}