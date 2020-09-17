package name.mjw.jquante.math.qm;

import java.util.ArrayList;
import java.util.List;
import java.util.stream.IntStream;

import org.hipparchus.geometry.euclidean.threed.Vector3D;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

import name.mjw.jquante.math.qm.basis.BasisSetLibrary;
import name.mjw.jquante.math.qm.basis.ContractedGaussian;
import name.mjw.jquante.math.qm.basis.Power;
import name.mjw.jquante.math.qm.basis.PrimitiveGaussian;
import name.mjw.jquante.math.qm.integral.Integrals;
import name.mjw.jquante.math.qm.integral.IntegralsUtil;
import name.mjw.jquante.molecule.Molecule;

import net.jafama.FastMath;

/**
 * The 2E integral driver.
 * 
 * @author V.Ganesh
 * @version 2.0 (Part of MeTA v2.0)
 */
public final class TwoElectronIntegrals {
	private static final Logger LOG = LogManager.getLogger(TwoElectronIntegrals.class);

	private BasisSetLibrary basisSetLibrary;

	private Molecule molecule;

	/**
	 * Holds value of property twoEIntegrals.
	 */
	private double[] twoEIntegrals;

	private int atomIndex;
	protected SCFMethod scfMethod;

	protected ArrayList<double[]> twoEDer;

	/** Calculate integrals on the fly instead of storing them in an array */
	protected boolean onTheFly;

	/**
	 * Creates a new instance of TwoElectronIntegrals
	 * 
	 * @param basisSetLibrary the basis functions to be used
	 */
	public TwoElectronIntegrals(BasisSetLibrary basisSetLibrary) {
		this(basisSetLibrary, false);
	}

	/**
	 * Creates a new instance of TwoElectronIntegrals
	 * 
	 * @param basisSetLibrary the basis functions to be used
	 * @param onTheFly        if true, the 2E integrals are not calculated and
	 *                        stored, they must be calculated individually by
	 *                        calling compute2E(i,j,k,l)
	 */
	public TwoElectronIntegrals(BasisSetLibrary basisSetLibrary, boolean onTheFly) {
		this.basisSetLibrary = basisSetLibrary;

		this.onTheFly = onTheFly;

		// compute the 2E integrals
		if (!onTheFly) {
			try {
				compute2E(); // try to do compute 2E incore
			} catch (OutOfMemoryError e) {
				// if no memory, resort to direct SCF
				LOG.error("No memory for in-core integral evaluation" + ". Switching to direct integral evaluation.");
				this.onTheFly = true;
			}
		}
	}

	/**
	 * Creates a new instance of TwoElectronIntegrals, and computes two electron
	 * integrals using shell-pair method rather than the conventional loop over all
	 * basis functions approach. Note that, this requires the Molecule object to be
	 * passed as an argument for it to function correctly. Note that the Atom
	 * objects which are a part of Molecule object must be specially configured to
	 * have a user defined property called "basisFunctions" that is essentially a
	 * list of ContractedGaussians that are centered on the atom. If this breaks,
	 * then the code will silently default to using the conventional way of
	 * computing the two electron integrals.
	 * 
	 * @param basisSetLibrary Basis functions for the molecule.
	 * @param molecule        Molecule for the two electron integrals.
	 * @param onTheFly        if true, the 2E integrals are not calculated and
	 *                        stored, they must be calculated individually by
	 *                        calling compute2EShellPair(i,j,k,l)
	 */
	public TwoElectronIntegrals(BasisSetLibrary basisSetLibrary, Molecule molecule, boolean onTheFly) {
		this.basisSetLibrary = basisSetLibrary;
		this.molecule = molecule;

		this.onTheFly = onTheFly;

		// compute the 2E integrals
		if (!onTheFly) {
			try {
				// first try in - core method
				try {
					compute2EShellPair(); // try to use shell pair algorithm
				} catch (Exception ignored) {
					LOG.error("WARNING: Using a slower algorithm" + " to compute two electron integrals. Error is: "
							+ ignored.toString());
					LOG.error(ignored);

					compute2E(); // if it doesn't succeed then fall back to
									// normal
				}
			} catch (OutOfMemoryError e) {
				// if no memory, resort to direct SCF
				LOG.error("No memory for in-core integral evaluation" + ". Switching to direct integral evaluation.");
				this.onTheFly = true;
			}
		}
	}

	/**
	 * compute the 2E integrals, and store it in a single 1D array, in the form
	 * [ijkl].
	 * 
	 * This method has been modified to take advantage of multi core systems where
	 * available.
	 */
	protected void compute2E() {
		LOG.debug("compute2E() called");
		final List<ContractedGaussian> bfs = basisSetLibrary.getBasisFunctions();

		// allocate required memory
		final int noOfBasisFunctions = bfs.size();
		final int noOfIntegrals = noOfBasisFunctions * (noOfBasisFunctions + 1)
				* (noOfBasisFunctions * noOfBasisFunctions + noOfBasisFunctions + 2) / 8;

		LOG.debug("noOfIntegrals is {}", noOfIntegrals);

		twoEIntegrals = new double[noOfIntegrals];

		// we only need i <= j, k <= l, and ij >= kl
		IntStream.range(0, noOfBasisFunctions).parallel().forEach(i -> {
			IntStream.range(0, i + 1).parallel().forEach(j -> {
				int ij = i * (i + 1) / 2 + j;
				IntStream.range(0, noOfBasisFunctions).parallel().forEach(k -> {
					IntStream.range(0, k + 1).parallel().forEach(l -> {
						int kl = k * (k + 1) / 2 + l;

						if (ij >= kl) {
							int ijkl = IntegralsUtil.ijkl2intindex(i, j, k, l);
							// record the 2E integrals
							twoEIntegrals[ijkl] = Integrals.coulomb(bfs.get(i), bfs.get(j), bfs.get(k), bfs.get(l));
						}
					});
				});
			});
		});

	}

	/** Compute a single 2E derivative element */
	private Vector3D compute2EDerivativeElement(ContractedGaussian bfi, ContractedGaussian bfj, ContractedGaussian bfk,
			ContractedGaussian bfl) {

		Vector3D twoEDerEle = new Vector3D(0, 0, 0);
		int[] paramIdx;
		Power currentPower;
		Vector3D currentOrigin;
		double currentAlpha;

		if (bfi.getCenteredAtom().getIndex() == atomIndex) {
			paramIdx = new int[] { 4, 1, 2, 3 };
			for (PrimitiveGaussian iPG : bfi.getPrimitives()) {
				currentOrigin = iPG.getOrigin();
				currentPower = iPG.getPowers();
				currentAlpha = iPG.getExponent();
				for (PrimitiveGaussian jPG : bfj.getPrimitives()) {
					for (PrimitiveGaussian kPG : bfk.getPrimitives()) {
						for (PrimitiveGaussian lPG : bfl.getPrimitives()) {
							compute2EDerivativeElementHelper(iPG, jPG, kPG, lPG, currentOrigin, currentPower,
									currentAlpha, paramIdx, twoEDerEle);
						}
					}
				}
			}
		}

		if (bfj.getCenteredAtom().getIndex() == atomIndex) {
			paramIdx = new int[] { 0, 4, 2, 3 };
			for (PrimitiveGaussian iPG : bfi.getPrimitives()) {
				for (PrimitiveGaussian jPG : bfj.getPrimitives()) {
					currentOrigin = jPG.getOrigin();
					currentPower = jPG.getPowers();
					currentAlpha = jPG.getExponent();
					for (PrimitiveGaussian kPG : bfk.getPrimitives()) {
						for (PrimitiveGaussian lPG : bfl.getPrimitives()) {
							compute2EDerivativeElementHelper(iPG, jPG, kPG, lPG, currentOrigin, currentPower,
									currentAlpha, paramIdx, twoEDerEle);
						} // end for
					} // end for
				} // end for
			} // end for
		} // end if

		if (bfk.getCenteredAtom().getIndex() == atomIndex) {
			paramIdx = new int[] { 0, 1, 4, 3 };
			for (PrimitiveGaussian iPG : bfi.getPrimitives()) {
				for (PrimitiveGaussian jPG : bfj.getPrimitives()) {
					for (PrimitiveGaussian kPG : bfk.getPrimitives()) {
						currentOrigin = kPG.getOrigin();
						currentPower = kPG.getPowers();
						currentAlpha = kPG.getExponent();
						for (PrimitiveGaussian lPG : bfl.getPrimitives()) {
							compute2EDerivativeElementHelper(iPG, jPG, kPG, lPG, currentOrigin, currentPower,
									currentAlpha, paramIdx, twoEDerEle);
						} // end for
					} // end for
				} // end for
			} // end for
		} // end if

		if (bfl.getCenteredAtom().getIndex() == atomIndex) {
			paramIdx = new int[] { 0, 1, 2, 4 };
			for (PrimitiveGaussian iPG : bfi.getPrimitives()) {
				for (PrimitiveGaussian jPG : bfj.getPrimitives()) {
					for (PrimitiveGaussian kPG : bfk.getPrimitives()) {
						for (PrimitiveGaussian lPG : bfl.getPrimitives()) {
							currentOrigin = lPG.getOrigin();
							currentPower = lPG.getPowers();
							currentAlpha = lPG.getExponent();

							compute2EDerivativeElementHelper(iPG, jPG, kPG, lPG, currentOrigin, currentPower,
									currentAlpha, paramIdx, twoEDerEle);
						} // end for
					} // end for
				} // end for
			} // end for
		} // end if

		return twoEDerEle;
	}

	/** helper for computing 2E derivative */
	private void compute2EDerivativeElementHelper(PrimitiveGaussian iPG, PrimitiveGaussian jPG, PrimitiveGaussian kPG,
			PrimitiveGaussian lPG, Vector3D currentOrigin, Power currentPower, double currentAlpha, int[] paramIdx,
			Vector3D derEle) {

		int l = currentPower.getL();
		int m = currentPower.getM();
		int n = currentPower.getN();

		double coeff = iPG.getCoefficient() * jPG.getCoefficient() * kPG.getCoefficient() * lPG.getCoefficient();

		PrimitiveGaussian xPG = new PrimitiveGaussian(currentOrigin, new Power(l + 1, m, n), currentAlpha, coeff);

		PrimitiveGaussian[] pgs = new PrimitiveGaussian[] { iPG, jPG, kPG, lPG, xPG };

		double terma = FastMath.sqrt(currentAlpha * (2.0 * l + 1.0)) * coeff
				* Integrals.coulomb(pgs[paramIdx[0]], pgs[paramIdx[1]], pgs[paramIdx[2]], pgs[paramIdx[3]]);
		double termbx = 0.0;
		double termby = 0.0;
		double termbz = 0.0;

		if (l > 0) {
			xPG = new PrimitiveGaussian(currentOrigin, new Power(xPG.getPowers().getL() - 1, xPG.getPowers().getM(), xPG.getPowers().getN()), coeff, currentAlpha);
			termbx = -2.0 * l * FastMath.sqrt(currentAlpha / (2. * l - 1)) * coeff
					* Integrals.coulomb(pgs[paramIdx[0]], pgs[paramIdx[1]], pgs[paramIdx[2]], pgs[paramIdx[3]]);
		}

		if (m > 0) {
			xPG = new PrimitiveGaussian(currentOrigin, new Power(xPG.getPowers().getL(), xPG.getPowers().getM() - 1, xPG.getPowers().getN()), coeff, currentAlpha);
			termby = -2.0 * m * FastMath.sqrt(currentAlpha / (2. * m - 1)) * coeff
					* Integrals.coulomb(pgs[paramIdx[0]], pgs[paramIdx[1]], pgs[paramIdx[2]], pgs[paramIdx[3]]);
		}

		if (n > 0) {
			xPG = new PrimitiveGaussian(currentOrigin, new Power(xPG.getPowers().getL(), xPG.getPowers().getM(), xPG.getPowers().getN() - 1), coeff, currentAlpha);
			termbz = -2.0 * n * FastMath.sqrt(currentAlpha / (2. * n - 1)) * coeff
					* Integrals.coulomb(pgs[paramIdx[0]], pgs[paramIdx[1]], pgs[paramIdx[2]], pgs[paramIdx[3]]);
		}

		derEle = new Vector3D(derEle.getX() + terma + termbx, derEle.getY() + terma + termby,
				derEle.getZ() + terma + termbz);
	}

	/**
	 * Compute the 2E integrals using shell pair based method, and store it in a
	 * single 1D array, in the form [ijkl].
	 * 
	 * TODO
	 * 1) This does not appear to be the shell pair method (mjw).
	 * 2) Add parallel support.
	 *
	 */
	protected void compute2EShellPair() {
		List<ContractedGaussian> bfs = basisSetLibrary.getBasisFunctions();

		// allocate required memory
		int noOfBasisFunctions = bfs.size();
		int noOfIntegrals = noOfBasisFunctions * (noOfBasisFunctions + 1)
				* (noOfBasisFunctions * noOfBasisFunctions + noOfBasisFunctions + 2) / 8;

		twoEIntegrals = new double[noOfIntegrals];

		int noOfAtoms = molecule.getNumberOfAtoms();
		int naFunc;
		int nbFunc;
		int ncFunc;
		int ndFunc;
		int twoEIndx;
		ArrayList<ContractedGaussian> aFunc;
		ArrayList<ContractedGaussian> bFunc;
		ArrayList<ContractedGaussian> cFunc;
		ArrayList<ContractedGaussian> dFunc;
		ContractedGaussian iaFunc;
		ContractedGaussian jbFunc;
		ContractedGaussian kcFunc;
		ContractedGaussian ldFunc;

		// center a
		for (int a = 0; a < noOfAtoms; a++) {
			aFunc = (ArrayList<ContractedGaussian>) molecule.getAtom(a).getUserDefinedAtomProperty("basisFunctions")
					.getValue();
			naFunc = aFunc.size();
			// basis functions on a
			for (int i = 0; i < naFunc; i++) {
				iaFunc = aFunc.get(i);

				// center b
				for (int b = 0; b <= a; b++) {
					bFunc = (ArrayList<ContractedGaussian>) molecule.getAtom(b)
							.getUserDefinedAtomProperty("basisFunctions").getValue();
					nbFunc = (b < a) ? bFunc.size() : i + 1;
					// basis functions on b
					for (int j = 0; j < nbFunc; j++) {
						jbFunc = bFunc.get(j);

						// center c
						for (int c = 0; c < noOfAtoms; c++) {
							cFunc = (ArrayList<ContractedGaussian>) molecule.getAtom(c)
									.getUserDefinedAtomProperty("basisFunctions").getValue();
							ncFunc = cFunc.size();
							// basis functions on c
							for (int k = 0; k < ncFunc; k++) {
								kcFunc = cFunc.get(k);

								// center d
								for (int d = 0; d <= c; d++) {
									dFunc = (ArrayList<ContractedGaussian>) molecule.getAtom(d)
											.getUserDefinedAtomProperty("basisFunctions").getValue();
									ndFunc = (d < c) ? dFunc.size() : k + 1;
									// basis functions on d
									for (int l = 0; l < ndFunc; l++) {
										ldFunc = dFunc.get(l);

										twoEIndx = IntegralsUtil.ijkl2intindex(iaFunc.getBasisFunctionIndex(), jbFunc.getBasisFunctionIndex(),
												kcFunc.getBasisFunctionIndex(), ldFunc.getBasisFunctionIndex());

										twoEIntegrals[twoEIndx] = compute2E(iaFunc, jbFunc, kcFunc, ldFunc);
									} // end for (l)
								} // end for (d)
							} // end for (k)
						} // end for (c)
					} // end for (j)
				} // end for (b)
			} // end for (i)
		} // end for (a)
	}

	/**
	 * 
	 * Compute an integral centered at &lt;ij&#x7C;kl&gt;
	 * 
	 * @param i Index of contracted Gaussian function i.
	 * @param j Index of contracted Gaussian function j.
	 * @param k Index of contracted Gaussian function k.
	 * @param l Index of contracted Gaussian function l.
	 * @return the value of two integral electron
	 */
	public double compute2E(int i, int j, int k, int l) {
		// if we really need to compute, go ahead!
		List<ContractedGaussian> bfs = basisSetLibrary.getBasisFunctions();

		return Integrals.coulomb(bfs.get(i), bfs.get(j), bfs.get(k), bfs.get(l));
	}

	/**
	 * Compute an integral centered at &lt;ij&#x7C;kl&gt;
	 * 
	 * @param cgi Contracted Gaussian function i.
	 * @param cgj Contracted Gaussian function j.
	 * @param cgk Contracted Gaussian function k.
	 * @param cgl Contracted Gaussian function l.
	 * @return the value of two integral electron
	 */
	public double compute2E(ContractedGaussian cgi, ContractedGaussian cgj, ContractedGaussian cgk,
			ContractedGaussian cgl) {
		return Integrals.coulomb(cgi, cgj, cgk, cgl);
	}

	/**
	 * Getter for property twoEIntegrals.
	 * 
	 * @return Value of property twoEIntegrals.
	 */
	public double[] getTwoEIntegrals() {
		return this.twoEIntegrals;
	}

	/**
	 * Setter for property twoEIntegrals.
	 * 
	 * @param twoEIntegrals New value of property twoEIntegrals.
	 */
	public void setTwoEIntegrals(double[] twoEIntegrals) {
		this.twoEIntegrals = twoEIntegrals;
	}

	/**
	 * Return the derivatives of the two electron integrals.
	 * 
	 * @param atomIndex the reference atom index
	 * @param scfMethod the reference SCF method
	 */
	public void compute2EDerivatives(int atomIndex, SCFMethod scfMethod) {
		final List<ContractedGaussian> bfs = basisSetLibrary.getBasisFunctions();
		final int noOfBasisFunctions = bfs.size();

		this.atomIndex = atomIndex;
		this.scfMethod = scfMethod;

		twoEDer = new ArrayList<>();

		double[] dxTwoE = new double[twoEIntegrals.length];
		double[] dyTwoE = new double[twoEIntegrals.length];
		double[] dzTwoE = new double[twoEIntegrals.length];

		twoEDer.add(dxTwoE);
		twoEDer.add(dyTwoE);
		twoEDer.add(dzTwoE);

		// we only need i <= j, k <= l, and ij >= kl
		IntStream.range(0, noOfBasisFunctions).parallel().forEach(i -> {
			ContractedGaussian bfi = bfs.get(i);

			IntStream.range(0, i + 1).parallel().forEach(j -> {
				ContractedGaussian bfj = bfs.get(j);
				int ij = i * (i + 1) / 2 + j;
				IntStream.range(0, noOfBasisFunctions).parallel().forEach(k -> {
					ContractedGaussian bfk = bfs.get(k);
					IntStream.range(0, k + 1).parallel().forEach(l -> {
						ContractedGaussian bfl = bfs.get(l);
						int kl = k * (k + 1) / 2 + l;

						if (ij >= kl) {
							int ijkl = IntegralsUtil.ijkl2intindex(i, j, k, l);

							// record derivative of the 2E integrals
							Vector3D twoEDerEle = compute2EDerivativeElement(bfi, bfj, bfk, bfl);

							dxTwoE[ijkl] = twoEDerEle.getX();
							dyTwoE[ijkl] = twoEDerEle.getY();
							dzTwoE[ijkl] = twoEDerEle.getZ();

						}
					});
				});
			});
		});
	}

	/**
	 * The currently computed list of two electron derivatives
	 * 
	 * @return partial derivatives of 2E integrals, computed in previous call to
	 *         compute2EDerivatives()
	 */
	public ArrayList<double[]> getTwoEDer() {
		return twoEDer;
	}

	/**
	 * Get the value of onTheFly
	 * 
	 * @return the value of onTheFly
	 */
	public boolean isOnTheFly() {
		return onTheFly;
	}

	/**
	 * Set the value of onTheFly
	 * 
	 * @param onTheFly new value of onTheFly
	 */
	public void setOnTheFly(boolean onTheFly) {
		this.onTheFly = onTheFly;
	}

	public void printCompleteIntegralList() {

		int noOfBasisFunctions = basisSetLibrary.getBasisFunctions().size();

		System.out.println("");
		System.out.println("Complete integral list");
		System.out.println("======================");
		for (int i = 0; i < noOfBasisFunctions; i++) {
			for (int j = 0; j < noOfBasisFunctions; j++) {
				for (int k = 0; k < noOfBasisFunctions; k++) {
					for (int l = 0; l < noOfBasisFunctions; l++) {
						int ijkl = IntegralsUtil.ijkl2intindex(i, j, k, l);
						System.out.println(i + " " + j + " " + k + " " + l + "\t\t" + twoEIntegrals[ijkl]);

					}

				}

			}

		}
	}
}
