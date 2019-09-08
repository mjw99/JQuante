package name.mjw.jquante.math.qm.integral;

/**
 * Factory interface for instantiating appropriate integral packages. Follows a
 * singleton pattern.
 * 
 * @author V.Ganesh
 * @version 2.0 (Part of MeTA v2.0)
 */
public final class IntegralsPackageFactory {

	private static IntegralsPackageFactory theInstance;
	protected IntegralPackageType defaultTwoElectronIntegralPackage;
	private OneElectronTerm oneElectronTerm;
	private NuclearTerm nuclearTerm;

	private IntegralsPackageFactory() {
		defaultTwoElectronIntegralPackage = IntegralPackageType.TWO_ELECTRON_HUZINAGA;
	}

	/**
	 * Get the instance of IntegralsPackageFactory
	 * 
	 * @return instance of IntegralsPackageFactory
	 */
	public static IntegralsPackageFactory getInstance() {
		if (theInstance == null) {
			theInstance = new IntegralsPackageFactory();
		}

		return theInstance;
	}

	/**
	 * Request a nuclear term integral package
	 * 
	 * @return instance of NuclearTerm
	 */
	public NuclearTerm getNuclearTerm() {
		if (nuclearTerm == null)
			nuclearTerm = new NuclearTerm();

		return nuclearTerm;
	}

	/**
	 * Request a one electron integral package
	 * 
	 * @return instance of OneElectronTerm
	 */
	public OneElectronTerm getOneElectronTerm() {
		if (oneElectronTerm == null)
			oneElectronTerm = new OneElectronTerm();

		return oneElectronTerm;
	}

	/**
	 * Request a 2E integral package.
	 * 
	 * @param type
	 *            requests a 2-E integral package
	 * @return the instance of 2-E integral package, if no suitable one
	 *         available, then UnsupportedOperationException is thrown
	 */
	public TwoElectronTerm getTwoElectronTerm(IntegralPackageType type) {
		TwoElectronTerm twoElectronTerm;

		switch (type) {
		case TWO_ELECTRON_HUZINAGA:
			twoElectronTerm = new HuzinagaTwoElectronTerm();
			break;
		case TWO_ELECTRON_RYS:
			twoElectronTerm = new RysTwoElectronTerm();
			break;
		case TWO_ELECTRON_HGP:
			twoElectronTerm = new HGPTwoElectronTerm();
			break;
		default:
			throw new UnsupportedOperationException("No 2E integral "
					+ "package yet for: " + type);
		}

		return twoElectronTerm;
	}

	/**
	 * Get the default 2E integral package. Currently this is:
	 * TWO_ELECTRON_HUZINAGA
	 * 
	 * @return the instance of 2-E integral package, if no suitable one
	 *         available, then UnsupportedOperationException is thrown
	 */
	public TwoElectronTerm getTwoElectronTerm() {
		return getTwoElectronTerm(defaultTwoElectronIntegralPackage);
	}

	/**
	 * Get the value of defaultTwoElectronIntegralPackage
	 * 
	 * @return the value of defaultTwoElectronIntegralPackage
	 */
	public IntegralPackageType getDefaultTwoElectronIntegralPackage() {
		return defaultTwoElectronIntegralPackage;
	}

	/**
	 * Set the value of defaultTwoElectronIntegralPackage
	 * 
	 * @param defaultTwoElectronIntegralPackage
	 *            new value of defaultTwoElectronIntegralPackage
	 */
	public void setDefaultTwoElectronIntegralPackage(
			IntegralPackageType defaultTwoElectronIntegralPackage) {
		this.defaultTwoElectronIntegralPackage = defaultTwoElectronIntegralPackage;

		// change here affects the statically initialised Integral package
		// so update reference there
		Integrals.reInitPackageReference();
	}

	/**
	 * Get a requested integral package.
	 * 
	 * @param type
	 *            the requested IntegralPackageType
	 * @return the instance of the requested IntegralPackageType, if not
	 *         available the an UnsupportedOperationException is thrown
	 */
	public IntegralsPackage getPackage(IntegralPackageType type) {
		switch (type) {
		case NUCLEAR_TERM:
			return getNuclearTerm();
		case ONE_ELECTRON_TERM:
			return getOneElectronTerm();
		case TWO_ELECTRON_HUZINAGA:
		case TWO_ELECTRON_HGP:
		case TWO_ELECTRON_RYS:
			return getTwoElectronTerm(type);
		default:
			throw new UnsupportedOperationException("No integral "
					+ "package yet for: " + type);
		}
	}
}
