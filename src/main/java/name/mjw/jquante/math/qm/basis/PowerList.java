package name.mjw.jquante.math.qm.basis;

import java.lang.ref.WeakReference;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;

/**
 * Orbital symbol to power list map. Follows a singleton pattern.
 * 
 * @author V.Ganesh
 * @version 2.0 (Part of MeTA v2.0)
 */
public final class PowerList {

	private static WeakReference<PowerList> _powerList = null;

	private HashMap<String, ArrayList<Power>> thePowerList;

	/** Creates a new instance of PowerList */
	private PowerList() {
		// put in the standard powers for the orbital symbols

		thePowerList = new HashMap<>(6);

		thePowerList.put("S", generatePowerList(0));
		thePowerList.put("P", generatePowerList(1));
		thePowerList.put("D", generatePowerList(2));
		thePowerList.put("F", generatePowerList(3));
		thePowerList.put("G", generatePowerList(4));
		thePowerList.put("H", generatePowerList(5));
	}

	/**
	 * Generate a power list for the given maximum angular momentum value
	 * 
	 * @param maxAngularMomentum
	 *            the maximum angular momentum
	 * @return the Power list in order for this maximum angular momentum
	 */
	public ArrayList<Power> generatePowerList(int maxAngularMomentum) {
		ArrayList<Power> pList = new ArrayList<>();

		for (int i = maxAngularMomentum; i >= 0; i--)
			for (int j = maxAngularMomentum - i; j >= 0; j--)
				pList.add(new Power(i, j, maxAngularMomentum - i - j));

		return pList;
	}

	/**
	 * Get an instance (and the only one) of PowerList
	 * 
	 * @return PowerList instance
	 */
	public static PowerList getInstance() {
		if (_powerList == null) {
			_powerList = new WeakReference<>(new PowerList());
		}

		PowerList powerList = _powerList.get();

		if (powerList == null) {
			powerList = new PowerList();
			_powerList = new WeakReference<>(powerList);
		}

		return powerList;
	}

	/**
	 * get the power list for the specified orbital symbol ('S', 'P', 'D' or 'F'
	 * .. no explicit error checking done, but will throw a RuntimeException if
	 * the arguments are incorrect)
	 * 
	 * @param orbital
	 *            - 'S', 'P', 'D' or 'F'
	 * @return Iterator of Power object representing the powers
	 */
	public Iterator<Power> getPowerList(String orbital) {
		return thePowerList.get(orbital).iterator();
	}
}