package name.mjw.jquante.math.qm.basis;

import java.util.ArrayList;
import java.util.List;

/**
 * Represents a single entity, i.e. an element, in a BasisSet
 * 
 * @author V.Ganesh
 * @version 2.0 (Part of MeTA v2.0)
 */
public class AtomicBasis {

	/**
	 * The chemical symbol of the atom.
	 */
	private String chemicalSymbol;

	/**
	 * The atomic number, Z, of the atom.
	 */
	private int atomicNumber;

	/**
	 * Holds value of property orbitals.
	 */
	private List<Orbital> orbitals;

	/**
	 * Creates a new instance of AtomicBasis
	 * 
	 * @param chemicalSymbol
	 *            the atomic symbol, of whose this is basis
	 * @param atomicNumber
	 *            its atomic number
	 */
	public AtomicBasis(String chemicalSymbol, int atomicNumber) {
		this.chemicalSymbol = chemicalSymbol;
		this.atomicNumber = atomicNumber;

		orbitals = new ArrayList<>();
	}

	/**
	 * Getter for property symbol.
	 * 
	 * @return Value of property symbol.
	 */
	public String getChemicalSymbol() {
		return this.chemicalSymbol;
	}

	/**
	 * Setter for property symbol.
	 * 
	 * @param symbol
	 *            New value of property symbol.
	 */
	public void setChemicalSymbol(String symbol) {
		this.chemicalSymbol = symbol;
	}

	/**
	 * Getter for property atomicNumber.
	 * 
	 * @return Value of property atomicNumber.
	 */
	public int getAtomicNumber() {
		return this.atomicNumber;
	}

	/**
	 * Setter for property atomicNumber.
	 * 
	 * @param atomicNumber
	 *            New value of property atomicNumber.
	 */
	public void setAtomicNumber(int atomicNumber) {
		this.atomicNumber = atomicNumber;
	}

	/**
	 * Getter for property orbitals.
	 * 
	 * @return Value of property orbitals.
	 */
	public List<Orbital> getOrbitals() {
		return this.orbitals;
	}

	/**
	 * Setter for property orbitals.
	 * 
	 * @param orbitals
	 *            New value of property orbitals.
	 */
	public void setOrbitals(List<Orbital> orbitals) {
		this.orbitals = orbitals;
	}

	/**
	 * Add an orbital object to this atomic basis
	 * 
	 * @param orbital
	 *            the Orbital object to be added to this atomic basis
	 */
	public void addOrbital(Orbital orbital) {
		orbitals.add(orbital);
	}
}