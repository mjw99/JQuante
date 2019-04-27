package name.mjw.jquante.math.qm.basis;

import java.util.ArrayList;
import java.util.List;
import java.util.Objects;

import name.mjw.jquante.config.impl.AtomInfo;
import org.apache.commons.math3.geometry.euclidean.threed.Vector3D;
import name.mjw.jquante.molecule.Atom;
import name.mjw.jquante.molecule.Molecule;
import net.jafama.FastMath;

/**
 * The class defines a contracted Gaussian and the operations on it.
 * 
 * @author V.Ganesh
 * @version 2.0 (Part of MeTA v2.0)
 */
public final class ContractedGaussian {

	/**
	 * Holds value of property origin.
	 */
	private final Vector3D origin;

	/**
	 * Holds value of property powers, whole array size is 3.
	 */
	private final Power powers;

	/**
	 * The Primitive Gaussians (PGs) that make up this Contracted Gaussian.
	 */
	private final ArrayList<PrimitiveGaussian> primitives;

	/**
	 * The list of exponents for this Gaussian.
	 */
	private final ArrayList<Double> exponents;

	/**
	 * Holds value of property coefficients.
	 */
	private final ArrayList<Double> coefficients;

	/**
	 * Holds value of property normalization.
	 */
	private double normalization;

	/**
	 * Holds value of property primNorms.
	 * 
	 * normalization factors for PGs
	 */
	private ArrayList<Double> primNorms;

	protected Atom centeredAtom;

	protected int index;

	/**
	 * Creates a new instance of ContractedGaussian
	 * 
	 * @param origin - the (x, y, z) on which this gaussian is centered
	 * @param powers - the powers of this gaussian
	 */
	public ContractedGaussian(Vector3D origin, Power powers) {
		this.origin = origin;
		this.powers = powers;

		primitives = new ArrayList<>(10);
		exponents = new ArrayList<>(10);
		coefficients = new ArrayList<>(10);
		primNorms = new ArrayList<>(10);

		this.normalization = 1;
	}

	/**
	 * Creates a new instance of ContractedGaussian
	 * 
	 * @param centeredAtom - CG centered on an atom
	 * @param powers       - the powers of this gaussian
	 */
	public ContractedGaussian(Atom centeredAtom, Power powers) {
		this(centeredAtom.getAtomCenterInAU(), powers);

		this.centeredAtom = centeredAtom;
	}

	/**
	 * Get the value of centeredAtom
	 * 
	 * @return the value of centeredAtom
	 */
	public Atom getCenteredAtom() {
		return centeredAtom;
	}

	/**
	 * Adds a primitive gaussian (PG) to this contracted gaussian list
	 * 
	 * @param exponent    the exponent for this PG
	 * @param coefficient the coefficient of this PG
	 */
	public void addPrimitive(double exponent, double coefficient) {
		primitives.add(new PrimitiveGaussian(origin, powers, exponent, coefficient));

		exponents.add(exponent);
		coefficients.add(coefficient);
	}

	/**
	 * Getter for property origin.
	 * 
	 * @return Value of property origin.
	 */
	public Vector3D getOrigin() {
		return this.origin;
	}

	/**
	 * Getter for property powers.
	 * 
	 * @return Value of property powers.
	 */
	public Power getPowers() {
		return this.powers;
	}

	/**
	 * Getter for property primitives.
	 * 
	 * @return Value of property primitives.
	 */
	public ArrayList<PrimitiveGaussian> getPrimitives() {
		return this.primitives;
	}

	/**
	 * Getter for property exponents.
	 * 
	 * @return Value of property exponents.
	 */
	public ArrayList<Double> getExponents() {
		return this.exponents;
	}

	/**
	 * Getter for property coefficients.
	 * 
	 * @return Value of property coefficients.
	 */
	public ArrayList<Double> getCoefficients() {
		return this.coefficients;
	}

	/**
	 * Getter for property normalisation.
	 * 
	 * @return Value of property normalisation.
	 */
	public double getNormalization() {
		return this.normalization;
	}

	/**
	 * Normalise this basis function.
	 * 
	 */
	public void normalize() {
		normalization = 1.0 / FastMath.sqrt(this.overlap(this));

		for (int i = 0; i < primitives.size(); i++) {
			primNorms.add(primitives.get(i).getNormalization());
		}
	}

	/**
	 * Overlap matrix element with another ContractedGaussian
	 * 
	 * @param cg the ContractedGaussian with which the overlap is to be be
	 *           determined.
	 * @return the overlap value
	 */
	public double overlap(ContractedGaussian cg) {
		double sij = 0.0;

		ArrayList<PrimitiveGaussian> cgPrimitives = cg.getPrimitives();

		for (int i = 0; i < primitives.size(); i++) {
			PrimitiveGaussian iPG = primitives.get(i);
			for (int j = 0; j < cgPrimitives.size(); j++) {
				PrimitiveGaussian jPG = cgPrimitives.get(j);

				sij += iPG.getCoefficient() * jPG.getCoefficient() * iPG.overlap(jPG);
			}
		}

		return normalization * cg.normalization * sij;
	}

	/**
	 * Compute the overlap energy derivative term w.r.t the specified atom index
	 * 
	 * @param atomIndex the reference atomIndex
	 * @param cg        the other cg
	 * @return the derivative term as an instance of Vector3D
	 */
	public Vector3D overlapDerivative(int atomIndex, ContractedGaussian cg) {
		Vector3D ovrDer = new Vector3D(0, 0, 0);

		// case 1: atomIndex is centered on this CG
		if (this.centeredAtom.getIndex() == atomIndex) {
			for (PrimitiveGaussian jPG : cg.primitives) {
				for (PrimitiveGaussian iPG : primitives) {
					overlapDerivativeHelper(iPG, jPG, origin, ovrDer);
				}
			}
		}

		// case 2: atomIndex is centered on the other CG
		if (cg.centeredAtom.getIndex() == atomIndex) {
			for (PrimitiveGaussian iPG : primitives) {
				for (PrimitiveGaussian jPG : cg.primitives) {
					overlapDerivativeHelper(jPG, iPG, cg.origin, ovrDer);
				}
			}
		} // end if

		// case 3: is atomIndex is not centered on any of the CGs, then
		// then the result is zero.

		return ovrDer;
	}

	/** helper method for overlap derivative */
	private void overlapDerivativeHelper(PrimitiveGaussian iPG,
			PrimitiveGaussian jPG, Vector3D pOrigin, Vector3D ovrDer) {
		int l = iPG.getPowers().getL();
		int m = iPG.getPowers().getM();
		int n = iPG.getPowers().getN();
		double coeff = iPG.getCoefficient() * jPG.getCoefficient();
		double alpha = iPG.getExponent();

		PrimitiveGaussian xPG = new PrimitiveGaussian(pOrigin, new Power(l + 1,
				m, n), coeff, alpha);

		double terma = FastMath.sqrt(alpha * (2.0 * l + 1.0)) * coeff
				* xPG.overlap(jPG);
		double termbx = 0;
		double termby;
		double termbz;

		if (l > 0) {
			xPG = new PrimitiveGaussian(pOrigin, new Power(xPG.getPowers().getL() - 1, xPG.getPowers().getM(), xPG.getPowers().getN()), coeff, alpha);
			termbx = -2.0 * l * FastMath.sqrt(alpha / (2.0 * l - 1.0)) * coeff
					* xPG.overlap(jPG);
		}

		xPG = new PrimitiveGaussian(pOrigin, new Power(xPG.getPowers().getL(), xPG.getPowers().getM() + 1, xPG.getPowers().getN()), coeff, alpha);
		terma = FastMath.sqrt(alpha * (2.0 * m + 1.0)) * coeff * xPG.overlap(jPG);

		if (m > 0) {
			xPG = new PrimitiveGaussian(pOrigin, new Power(xPG.getPowers().getL(), xPG.getPowers().getM() - 1, xPG.getPowers().getN()), coeff, alpha);
			termby = -2 * m * FastMath.sqrt(alpha / (2.0 * m - 1.0)) * coeff
					* xPG.overlap(jPG);
		} else {
			termby = 0.0;
		}

		xPG = new PrimitiveGaussian(pOrigin, new Power(xPG.getPowers().getL(), xPG.getPowers().getM(), xPG.getPowers().getN() + 1), coeff, alpha);
		terma = FastMath.sqrt(alpha * (2.0 * n + 1.0)) * coeff * xPG.overlap(jPG);

		if (n > 0) {
			xPG = new PrimitiveGaussian(pOrigin, new Power(xPG.getPowers().getL(), xPG.getPowers().getM(), xPG.getPowers().getN() -1), coeff, alpha);
			termbz = -2.0 * n * FastMath.sqrt(alpha / (2.0 * n - 1.0)) * coeff
					* xPG.overlap(jPG);
		} else {
			termbz = 0.0;
		}

		ovrDer = new Vector3D(ovrDer.getX() + terma + termbx,
				ovrDer.getY() + terma + termby,
				ovrDer.getZ() + terma + termbz);
	}

	/**
	 * Kinetic Energy (KE) matrix element with another ContractedGaussian
	 * 
	 * @param cg the ContractedGaussian with which KE is to be determined.
	 * @return the KE value
	 */
	public double kinetic(ContractedGaussian cg) {
		double tij = 0.0;
		List<PrimitiveGaussian> cgPrimitives = cg.getPrimitives();

		for (int i = 0; i < primitives.size(); i++) {
			PrimitiveGaussian iPG = primitives.get(i);
			for (int j = 0; j < cgPrimitives.size(); j++) {
				PrimitiveGaussian jPG = cgPrimitives.get(j);

				tij += iPG.getCoefficient() * jPG.getCoefficient() * iPG.kinetic(jPG);
			}
		}

		return normalization * cg.normalization * tij;
	}

	/**
	 * Compute the kinetic energy derivative term w.r.t the specified atom index
	 * 
	 * @param atomIndex the reference atomIndex
	 * @param cg        the other cg
	 * @return the derivative term as an instance of Vector3D
	 */
	public Vector3D kineticDerivative(int atomIndex, ContractedGaussian cg) {
		Vector3D kder = new Vector3D(0, 0, 0);

		// case 1: atomIndex is centered on this CG
		if (this.centeredAtom.getIndex() == atomIndex) {
			for (PrimitiveGaussian jPG : cg.primitives) {
				for (PrimitiveGaussian iPG : primitives) {
					kineticDerivativeHelper(iPG, jPG, origin, kder);
				} // end for
			} // end for
		} // end if

		// case 2: atomIndex is centered on the other CG
		if (cg.centeredAtom.getIndex() == atomIndex) {
			for (PrimitiveGaussian iPG : primitives) {
				for (PrimitiveGaussian jPG : cg.primitives) {
					kineticDerivativeHelper(jPG, iPG, cg.origin, kder);
				} // end for
			} // end for
		} // end if

		// case 3: is atomIndex is not centered on any of the CGs, then
		// then the result is zero.

		return kder;
	}

	/** helper method for kinetic energy derivative */
	private void kineticDerivativeHelper(PrimitiveGaussian iPG, PrimitiveGaussian jPG, Vector3D pOrigin,
			Vector3D kder) {
		int l = iPG.getPowers().getL();
		int m = iPG.getPowers().getM();
		int n = iPG.getPowers().getN();
		double coeff = iPG.getCoefficient() * jPG.getCoefficient();
		double alpha = iPG.getExponent();

		PrimitiveGaussian xPG = new PrimitiveGaussian(pOrigin, new Power(l + 1, m, n), coeff, alpha);

		double terma = FastMath.sqrt(alpha * (2.0 * l + 1.0)) * coeff * xPG.kinetic(jPG);
		double termbx = 0.0;
		double termby = 0.0;
		double termbz = 0.0;

		if (l > 0) {
			xPG = new PrimitiveGaussian(pOrigin, new Power(xPG.getPowers().getL() - 1, xPG.getPowers().getM(), xPG.getPowers().getN()), coeff, alpha);
			termbx = -2.0 * l * FastMath.sqrt(alpha / (2.0 * l - 1.0)) * coeff * xPG.kinetic(jPG);
		}

		xPG = new PrimitiveGaussian(pOrigin, new Power(xPG.getPowers().getL(), xPG.getPowers().getM() + 1, xPG.getPowers().getN()), coeff, alpha);
		terma = FastMath.sqrt(alpha * (2.0 * m + 1.0)) * coeff * xPG.kinetic(jPG);

		if (m > 0) {
			xPG = new PrimitiveGaussian(pOrigin, new Power(xPG.getPowers().getL(), xPG.getPowers().getM() - 1, xPG.getPowers().getN()), coeff, alpha);
			termby = -2 * m * FastMath.sqrt(alpha / (2.0 * m - 1.0)) * coeff * xPG.kinetic(jPG);
		} else {
			termby = 0.0;
		}

		
		xPG = new PrimitiveGaussian(pOrigin, new Power(xPG.getPowers().getL(), xPG.getPowers().getM(), xPG.getPowers().getN() + 1), coeff, alpha);
		terma = FastMath.sqrt(alpha * (2.0 * n + 1.0)) * coeff * xPG.kinetic(jPG);

		if (n > 0) {
			xPG = new PrimitiveGaussian(pOrigin, new Power(xPG.getPowers().getL(), xPG.getPowers().getM(), xPG.getPowers().getN() - 1), coeff, alpha);
			termbz = -2.0 * n * FastMath.sqrt(alpha / (2.0 * n - 1.0)) * coeff * xPG.kinetic(jPG);
		} else {
			termbz = 0.0;
		}

		kder = new Vector3D(kder.getX() + terma + termbx, kder.getY() + terma + termby, kder.getZ() + terma + termbz);
	}

	/**
	 * Nuclear matrix element with another ContractedGaussian and a center
	 * 
	 * @param cg     the ContractedGaussian with which nuclear interaction is to
	 *               determined.
	 * @param center the center at which nuclear energy is to be computed
	 * @return the nuclear value
	 */
	public double nuclear(ContractedGaussian cg, Vector3D center) {
		double vij = 0.0;
		int i;
		int j;

		List<PrimitiveGaussian> cgPrimitives = cg.getPrimitives();

		for (i = 0; i < primitives.size(); i++) {
			PrimitiveGaussian iPG = primitives.get(i);
			for (j = 0; j < cgPrimitives.size(); j++) {
				PrimitiveGaussian jPG = cgPrimitives.get(j);

				vij += iPG.getCoefficient() * jPG.getCoefficient() * iPG.nuclear(jPG, center);
			} // end for
		} // end for

		return normalization * cg.normalization * vij;
	}

	/**
	 * Derivative of nuclear attraction term, with respect to an atom index
	 * 
	 * @param mol       the reference Molecule object
	 * @param atomIndex the reference atomIndex
	 * @param cg        the other ContractedGaussian
	 * @return partial derivative of nuclear attraction integral term
	 */
	public Vector3D nuclearAttractionDerivative(Molecule mol, int atomIndex, ContractedGaussian cg) {
		Vector3D nder = new Vector3D(0, 0, 0);

		// case 1: atom index is centered on this CG
		if (this.centeredAtom.getIndex() == atomIndex) {
			for (PrimitiveGaussian jPG : cg.primitives) {
				for (PrimitiveGaussian iPG : primitives) {
					nuclearDerivativeHelper(mol, iPG, jPG, origin, nder);
				}
			}
		}

		// case 2: atomIndex is centered on the other CG
		if (cg.centeredAtom.getIndex() == atomIndex) {
			for (PrimitiveGaussian iPG : primitives) {
				for (PrimitiveGaussian jPG : cg.primitives) {
					nuclearDerivativeHelper(mol, jPG, iPG, cg.origin, nder);
				}
			}
		}

		double factor;
		double atno = AtomInfo.getInstance().getAtomicNumber(centeredAtom.getSymbol());
		for (PrimitiveGaussian iPG : primitives) {
			for (PrimitiveGaussian jPG : cg.primitives) {
				factor = iPG.getCoefficient() * jPG.getCoefficient() * atno;
				nder = nder.add(jPG.nuclearAttractionGradient(iPG, centeredAtom.getAtomCenter()).scalarMultiply(factor));
			}
		}

		return nder;
	}

	/** helper method for kinetic energy derivative */
	private void nuclearDerivativeHelper(Molecule mol, PrimitiveGaussian iPG, PrimitiveGaussian jPG, Vector3D pOrigin,
			Vector3D nder) {
		int l = iPG.getPowers().getL();
		int m = iPG.getPowers().getM();
		int n = iPG.getPowers().getN();
		double coeff = iPG.getCoefficient() * jPG.getCoefficient();
		double alpha = iPG.getExponent();

		PrimitiveGaussian xPG = new PrimitiveGaussian(pOrigin, new Power(l + 1, m, n), coeff, alpha);

		double terma = 0.0;

		AtomInfo ai = AtomInfo.getInstance();
		for (int i = 0; i < mol.getNumberOfAtoms(); i++) {
			Atom atom = mol.getAtom(i);
			terma += ai.getAtomicNumber(atom.getSymbol()) * FastMath.sqrt(alpha * (2.0 * l + 1.0)) * coeff
					* xPG.nuclear(jPG, atom.getAtomCenter());
		}

		double termbx = 0.0;
		double termby = 0.0;
		double termbz = 0.0;

		if (l > 0) {
			xPG = new PrimitiveGaussian(pOrigin, new Power(xPG.getPowers().getL() - 1, xPG.getPowers().getM(), xPG.getPowers().getN()), coeff, alpha);

			for (int i = 0; i < mol.getNumberOfAtoms(); i++) {
				Atom atom = mol.getAtom(i);
				termbx += -2.0 * l * ai.getAtomicNumber(atom.getSymbol()) * FastMath.sqrt(alpha / (2.0 * l - 1.0))
						* coeff * xPG.nuclear(jPG, atom.getAtomCenter());
			}
		}

		xPG = new PrimitiveGaussian(pOrigin, new Power(xPG.getPowers().getL(), xPG.getPowers().getM() + 1, xPG.getPowers().getN()), coeff, alpha);
		terma = 0.0;
		for (int i = 0; i < mol.getNumberOfAtoms(); i++) {
			Atom atom = mol.getAtom(i);
			terma += ai.getAtomicNumber(atom.getSymbol()) * FastMath.sqrt(alpha * (2.0 * m + 1.0)) * coeff
					* xPG.nuclear(jPG, atom.getAtomCenter());
		}

		if (m > 0) {
			xPG = new PrimitiveGaussian(pOrigin, new Power(xPG.getPowers().getL(), xPG.getPowers().getM() - 1, xPG.getPowers().getN()), coeff, alpha);
			for (int i = 0; i < mol.getNumberOfAtoms(); i++) {
				Atom atom = mol.getAtom(i);
				termby += -2.0 * m * ai.getAtomicNumber(atom.getSymbol()) * FastMath.sqrt(alpha / (2.0 * m - 1.0))
						* coeff * xPG.nuclear(jPG, atom.getAtomCenter());
			}
		}

		xPG = new PrimitiveGaussian(pOrigin, new Power(xPG.getPowers().getL(), xPG.getPowers().getM(), xPG.getPowers().getN() + 1), coeff, alpha);
		terma = 0.0;
		for (int i = 0; i < mol.getNumberOfAtoms(); i++) {
			Atom atom = mol.getAtom(i);
			terma += ai.getAtomicNumber(atom.getSymbol()) * FastMath.sqrt(alpha * (2.0 * n + 1.0)) * coeff
					* xPG.nuclear(jPG, atom.getAtomCenter());
		}

		if (n > 0) {
			xPG = new PrimitiveGaussian(pOrigin, new Power(xPG.getPowers().getL(), xPG.getPowers().getM(), xPG.getPowers().getN() - 1), coeff, alpha);
			
			for (int i = 0; i < mol.getNumberOfAtoms(); i++) {
				Atom atom = mol.getAtom(i);
				termbz += -2.0 * n * ai.getAtomicNumber(atom.getSymbol()) * FastMath.sqrt(alpha / (2.0 * n - 1.0))
						* coeff * xPG.nuclear(jPG, atom.getAtomCenter());
			}
		} else
			termbz = 0.0;

		nder = new Vector3D(nder.getX() + terma + termbx, nder.getY() + terma + termby, nder.getZ() + terma + termbz);
	}

	/**
	 * The amplitude of this primitive Gaussian at a given point.
	 * 
	 * @param point the reference point
	 * @return the amplitude of this PG at the specified point
	 */
	public double amplitude(Vector3D point) {
		double value = 0.0;

		for (int i = 0; i < primitives.size(); i++) {
			value += primitives.get(i).amplitude(point);
		}
		return (normalization * value);
	}

	/**
	 * Calculate laplacian at the specified point
	 * 
	 * @param point the point where Laplacian is to be computed
	 * @return the laplacian at this point
	 */
	public double laplacian(Vector3D point) {
		double value = 0.0;

		for (PrimitiveGaussian pg : primitives)
			value += pg.laplacian(point);

		return normalization * value;
	}

	/**
	 * Evaluate gradient of this function at a given point
	 * 
	 * @param point the point where gradient is to be evaluated
	 * @return partial derivatives with respect to x, y, z
	 */
	public Vector3D gradient(Vector3D point) {
		Vector3D grad = new Vector3D(0, 0, 0);

		for (PrimitiveGaussian pg : primitives)
			grad.add(pg.gradient(point));

		return grad.scalarMultiply(normalization);
	}

	/**
	 * Getter for property primNorms.
	 * 
	 * @return Value of property primNorms.
	 */
	public ArrayList<Double> getPrimNorms() {
		return this.primNorms;
	}

	/**
	 * Return the maximum angular momentum for this contracted basis function
	 * 
	 * @return the maximum angular momentum of this contracted basis function
	 */
	public int getMaximumAngularMomentum() {
		return powers.getMaximumAngularMomentum();
	}

	/**
	 * Return the minimum angular momentum for this contracted basis function
	 * 
	 * @return the minimum angular momentum of this contracted basis function
	 */
	public int getMinimumAngularMomentum() {
		return powers.getMinimumAngularMomentum();
	}

	/**
	 * Return the total angular momentum of this contracted basis function
	 * 
	 * @return the maximum of the contracted basis function
	 */
	public int getTotalAngularMomentum() {
		return powers.getTotalAngularMomentum();
	}

	/**
	 * Get the value of index
	 * 
	 * @return the value of index
	 */
	public int getIndex() {
		return index;
	}

	/**
	 * Set the value of index
	 * 
	 * @param index new value of index
	 */
	public void setIndex(int index) {
		this.index = index;
	}

	/**
	 * isSameShell determines that two ContractedGaussian belong to the same shell.
	 * A shell is a set of basis function with the same centre and same contracted
	 * exponent, for example, sp shell, d shell.
	 *
	 * @param otherCg The other contracted Gaussian.
	 * @return True if otherCg is in the same shell.
	 */
	public boolean isSameShell(ContractedGaussian otherCg) {

		if (!this.origin.equals(otherCg.origin)) {
			return false;
		}

		if (this.getPrimitives().size() != otherCg.getPrimitives().size()) {
			return false;
		}

		if (this.getTotalAngularMomentum() != otherCg.getTotalAngularMomentum()) {
			return false;
		}

		for (int i = 0; i < this.getExponents().size(); i++) {
			if (!this.getExponents().get(i).equals(otherCg.getExponents().get(i))) {

				return false;
			}

		}

		for (int i = 0; i < this.getCoefficients().size(); i++) {
			if (!this.getCoefficients().get(i).equals(otherCg.getCoefficients().get(i))) {

				return false;
			}

		}

		return true;
	}

	/**
	 * overloaded toString()
	 */
	@Override
	public String toString() {
		return "Origin : " + origin + " Powers : " + powers + " " + exponents + " " + coefficients + "\n ";
	}

	@Override
	public int hashCode() {
		return Objects.hash(centeredAtom, coefficients, exponents, index, normalization, origin, powers, primNorms,
				primitives);
	}

	@Override
	public boolean equals(Object obj) {
		if (this == obj)
			return true;
		if (obj == null)
			return false;
		if (getClass() != obj.getClass())
			return false;
		ContractedGaussian other = (ContractedGaussian) obj;
		return Objects.equals(centeredAtom, other.centeredAtom) && Objects.equals(coefficients, other.coefficients)
				&& Objects.equals(exponents, other.exponents) && index == other.index
				&& Double.doubleToLongBits(normalization) == Double.doubleToLongBits(other.normalization)
				&& Objects.equals(origin, other.origin) && Objects.equals(powers, other.powers)
				&& Objects.equals(primNorms, other.primNorms) && Objects.equals(primitives, other.primitives);
	}
}
