package name.mjw.jquante.molecule;

import java.util.Iterator;
import java.util.LinkedHashMap;

import name.mjw.jquante.molecule.event.MoleculeStateChangeEvent;
import name.mjw.jquante.molecule.event.MoleculeStateChangeListener;

/**
 * Represents a molecular formula.
 * 
 * @author V.Ganesh
 * @version 2.0 (Part of MeTA v2.0)
 */
public class MolecularFormula {

	// a private reference to the molecule object
	private Molecule molecule;

	// the formula string
	private String formulaString;

	// make a hashtable hashed by symbols
	private LinkedHashMap<String, Integer> uniqueSymbols;

	/** Creates a new instance of MolecularFormula */
	private MolecularFormula() {
	}

	/**
	 * Creates a new instance of MolecularFormula
	 * 
	 * @param molecule
	 *            Molecule object.
	 */
	public MolecularFormula(final Molecule molecule) {
		this.molecule = molecule;

		// add a change listener
		molecule.addMoleculeStateChangeListener(new MoleculeStateChangeListener() {
			@Override
			public void moleculeChanged(final MoleculeStateChangeEvent msce) {
				makeFormulaString();
			}
		});

		// in the beginning make the string
		makeFormulaString();
	}

	/**
	 * make the molecule string
	 */
	private void makeFormulaString() {
		// make a hashtable hashed by symbols
		uniqueSymbols = new LinkedHashMap<>(5);

		// iterate through each atoms and collect the unique ones and their
		// counts...
		final Iterator<Atom> atoms = molecule.getAtoms();
		Atom atom;

		// TODO: a bit expensive (?) operation, may need optimization
		while (atoms.hasNext()) {
			atom = atoms.next();

			if (uniqueSymbols.containsKey(atom.getSymbol())) {
				uniqueSymbols.put(atom.getSymbol(),
						Integer.valueOf((uniqueSymbols.get(atom.getSymbol())).intValue() + 1));
			} else {
				uniqueSymbols.put(atom.getSymbol(), Integer.valueOf(1));
			} // end if
		} // end while

		// now sort the stuff based on atomic symbols, alphabetically
		final Object[] keys = uniqueSymbols.keySet().toArray();
		final Object[] values = uniqueSymbols.values().toArray();
		sortKeys(keys, values);

		// now make the real molecular formula string
		makeFormulaString(keys, values);
	}

	/**
	 * Return atom count of a particular symbol. If none is found, zero is returned.
	 * 
	 * @param atomSymbol the atom symbol to query for
	 * @return the count of the atoms for the queried symbol
	 */
	public int getAtomCount(final String atomSymbol) {
		final Integer atmCnt = uniqueSymbols.get(atomSymbol);

		return (atmCnt == null ? 0 : atmCnt);
	}

	/**
	 * sort the keys
	 */
	private void sortKeys(final Object[] keys, final Object[] values) {
		Object temp;
		int i;
		int j;

		for (i = 0; i < keys.length; i++) {
			for (j = 0; j < i; j++) {
				if (((String) keys[i]).compareTo((String) keys[j]) < 0) {
					// exchange keys ...
					temp = keys[i];
					keys[i] = keys[j];
					keys[j] = temp;
					// and the values ...
					temp = values[i];
					values[i] = values[j];
					values[j] = temp;
				}
			}
		}
	}

	/**
	 * make the real string
	 */
	private void makeFormulaString(final Object[] keys, final Object[] values) {
		formulaString = "";

		for (int i = 0; i < keys.length; i++) {
			formulaString += keys[i].toString() + values[i].toString() + " ";

		}
	}

	/**
	 * Returns the formula as a string object
	 * 
	 * @return The formula as a string.
	 */
	public String getFormula() {
		return formulaString;
	}

	/**
	 * overridden toString()
	 */
	@Override
	public String toString() {
		return getFormula();
	}

	/**
	 * do add / sub depending on sign
	 */
	private MolecularFormula doAddSub(final MolecularFormula formula, final int sign) {
		final MolecularFormula newFormula = new MolecularFormula();

		// make a hashtable hashed by symbols
		newFormula.uniqueSymbols = new LinkedHashMap<>(5);

		// run through what is common and add/sub them
		Object[] keys = uniqueSymbols.keySet().toArray();

		for (int i = 0; i < keys.length; i++) {
			if (formula.uniqueSymbols.containsKey(keys[i])) {
				newFormula.uniqueSymbols.put(keys[i].toString(),
						uniqueSymbols.get(keys[i]) + (sign * formula.uniqueSymbols.get(keys[i])));
			} else {
				newFormula.uniqueSymbols.put(keys[i].toString(), uniqueSymbols.get(keys[i]));
			} // end if
		} // end for

		// and then add what ever is found in the formula to be added, but has
		// not alreay been added
		keys = formula.uniqueSymbols.keySet().toArray();
		for (int i = 0; i < keys.length; i++) {
			if (!newFormula.uniqueSymbols.containsKey(keys[i])) {
				newFormula.uniqueSymbols.put(keys[i].toString(), formula.uniqueSymbols.get(keys[i]));
			} // end if
		} // end for

		// then make the formula strings
		keys = newFormula.uniqueSymbols.keySet().toArray();
		final Object[] values = newFormula.uniqueSymbols.values().toArray();

		newFormula.sortKeys(keys, values);
		newFormula.makeFormulaString(keys, values);

		return newFormula;
	}

	/**
	 * Add a molecular formula and return the result as a new instance.
	 * 
	 * @param formula - is the formula to be added
	 * @return the new instance of the molecular formula
	 */
	public MolecularFormula add(final MolecularFormula formula) {
		return doAddSub(formula, 1);
	}

	/**
	 * Subtract a molecular formula and return the result as a new instance.
	 * 
	 * @param formula - is the formula to be subtracted
	 * @return the new instance of the molecular formula
	 */
	public MolecularFormula sub(final MolecularFormula formula) {
		return doAddSub(formula, -1);
	}

	/**
	 * overloaded equals()
	 * 
	 * @param obj the object to be compared
	 * @return true if they are equal, false otherwise
	 */
	@Override
	public boolean equals(final Object obj) {
		if (this == obj)
			return true;

		if (!(obj instanceof MolecularFormula)) {
			return false;
		} else {
			// then do a deep compare
			final MolecularFormula o = (MolecularFormula) obj;

			final Object[] keys = uniqueSymbols.keySet().toArray();

			for (int i = 0; i < keys.length; i++) {
				if (o.uniqueSymbols.containsKey(keys[i])) {
					if (!o.uniqueSymbols.get(keys[i]).equals(
							uniqueSymbols.get(keys[i]))) {
						return false;
					} // end if
				} else {
					return false;
				} // end if
			} // end for

			return true;
		} // end if
	}

	@Override
	public int hashCode() {
		int hash = 7;
		hash = 47
				* hash
				+ (this.uniqueSymbols != null ? this.uniqueSymbols.hashCode()
						: 0);
		return hash;
	}
}
