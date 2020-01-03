package name.mjw.jquante.molecule.impl;

import java.util.ArrayList;
import java.util.Enumeration;
import java.util.HashMap;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.PriorityQueue;

import org.hipparchus.geometry.euclidean.threed.Vector3D;

import name.mjw.jquante.config.impl.AtomInfo;
import name.mjw.jquante.math.geom.BoundingBox;
import name.mjw.jquante.molecule.Atom;
import name.mjw.jquante.molecule.BondType;
import name.mjw.jquante.molecule.MolecularFormula;
import name.mjw.jquante.molecule.Molecule;
import name.mjw.jquante.molecule.UserDefinedAtomProperty;
import name.mjw.jquante.molecule.event.MoleculeStateChangeEvent;

/**
 * The default implementation of Molecule interface.
 * 
 * @author V.Ganesh
 * @version 2.0 (Part of MeTA v2.0)
 */
public class MoleculeImpl extends Molecule {

	private String title;

	private ArrayList<Atom> atomList;

	private BoundingBox boundingBox;

	private boolean stateChanged;

	private MoleculeStateChangeEvent msce;

	private AtomInfo atomInfo;

	private MolecularFormula molecularFormula;

	/** Creates a new instance of MoleculeImpl */
	public MoleculeImpl() {
		this("Molecule");
	}

	public MoleculeImpl(String title) {
		this.zMatrixComputed = false;
		this.additionalInformationAvailable = false;
		this.title = title.trim();

		// init values
		atomList = new ArrayList<>();


		stateChanged = true;

		numberOfElectrons = 0;
		molecularMass = 0.0;

		// event listener
		msce = new MoleculeStateChangeEvent(this);

		// atom info..
		atomInfo = AtomInfo.getInstance();

		// default value for center of mass
		centerOfMass = new Vector3D(0.0, 0.0, 0.0);
	}

	/**
	 * Returns the title.
	 * 
	 * @return String identifying the title of the Molecule object
	 */
	@Override
	public String getTitle() {
		return title;
	}

	/**
	 * Overloaded addAtom() method.
	 * 
	 * @param atom
	 *            the instance of atom class
	 */
	@Override
	public void addAtom(Atom atom) {
		atom.setIndex(atomList.size());
		atomList.add(atom);

		// adjust the molecular weight, and center of mass
		double mass = atomInfo.getAtomicWeight(atom.getSymbol());
		molecularMass += mass;
	
		centerOfMass = new Vector3D(centerOfMass.getX() + (mass * atom.getX()),
				centerOfMass.getY() + (mass * atom.getY()),
				centerOfMass.getZ() + (mass * atom.getZ()));

		// and numberOfElectrons
		numberOfElectrons += atomInfo.getAtomicNumber(atom.getSymbol());

		// fire the change event
		stateChanged = true;
		msce.setEventType(MoleculeStateChangeEvent.ATOM_ADDED);
		msce.setAtom1(atom);
		fireMoleculeStateChangeListenerMoleculeChanged(msce);
	}

	/**
	 * remove an atom from this Molecule at a given index May throw runtime
	 * exception, if invalid index
	 * 
	 * @param atomIndex
	 *            the atom index to be removed
	 */
	@Override
	public void removeAtomAt(int atomIndex) {
		Atom atom = atomList.get(atomIndex);

		// adjust the molecular weight and COM
		double mass = atomInfo.getAtomicWeight(atom.getSymbol());
		molecularMass -= mass;
		centerOfMass = new Vector3D(centerOfMass.getX() - (mass * atom.getX()),
				centerOfMass.getY() - (mass * atom.getY()),
				centerOfMass.getZ() - (mass * atom.getZ()));

		// and numberOfElectrons
		numberOfElectrons -= atomInfo.getAtomicNumber(atom.getSymbol());

		// remove the atom
		atomList.remove(atomIndex);

		// then reindex the atoms
		Iterator<Atom> atoms = atomList.iterator();
		int index = 0;
		Atom atom1;

		while (atoms.hasNext()) {
			atom1 = atoms.next();
			atom1.setIndex(index++);
			atom1.removeAllConnections();
		} // emd while

		// and fire the change event
		stateChanged = true;
		msce.setEventType(MoleculeStateChangeEvent.ATOM_REMOVED);
		msce.setAtom1(atom);
		fireMoleculeStateChangeListenerMoleculeChanged(msce);
	}

	/**
	 * remove an atom from this Molecule
	 * 
	 * @param atom
	 *            the instance of atom class
	 */
	@Override
	public void removeAtom(Atom atom) {
		removeAtomAt(atom.getIndex());
	}

	/**
	 * Overloaded addAtom() method.
	 * 
	 * @param symbol
	 *            the atom symbol
	 * @param atomCenter
	 *            the Cartesian coordinates of the atom stored as Point3D object
	 */
	@Override
	public void addAtom(String symbol, Vector3D atomCenter) {
		addAtom(new Atom(symbol, atomCenter));
	}

	/**
	 * Adds an atom to this molecule object.
	 * 
	 * @param symbol
	 *            the atom symbol
	 * @param x
	 *            X coordinate of the atom
	 * @param y
	 *            Y coordinate of the atom
	 * @param z
	 *            Z coordinate of the atom
	 */
	@Override
	public void addAtom(String symbol, double x, double y, double z) {
		addAtom(new Atom(symbol, new Vector3D(x, y, z)));
	}

	/**
	 * Adds an atom to this molecule object.
	 * 
	 * @param symbol
	 *            the atom symbol
	 * @param x
	 *            X coordinate of the atom
	 * @param y
	 *            Y coordinate of the atom
	 * @param z
	 *            Z coordinate of the atom
	 * @param xi
	 *            X coordinate of the atom base (vector)
	 * @param yj
	 *            Y coordinate of the atom base (vector)
	 * @param zk
	 *            Z coordinate of the atom base (vector)
	 * @param index
	 *            the atom index
	 */
	public void addAtom(String symbol, double x, double y, double z, double xi,
			double yj, double zk, int index) {
		Atom atm = new Atom(symbol, new Vector3D(x, y, z));
		atm.addUserDefinedAtomProperty(new UserDefinedAtomProperty(
				"baseCenter", new Vector3D(xi, yj, zk)));
		addAtom(atm);
	}

	/**
	 * Adds an atom to this molecule object.
	 * 
	 * @param symbol
	 *            the atom symbol
	 * @param x
	 *            X coordinate of the atom
	 * @param y
	 *            Y coordinate of the atom
	 * @param z
	 *            Z coordinate of the atom
	 * @param index
	 *            the atom index of this atom
	 */
	@Override
	public void addAtom(String symbol, double x, double y, double z, int index) {
		addAtom(new Atom(symbol, new Vector3D(x, y, z), index));
	}


	/**
	 * Method to get a particular atom from the lists of atoms.
	 * 
	 * @return Atom the instance of atom class
	 * @throws IndexOutOfBoundsException
	 *             If not a valid atomIndex
	 */
	@Override
	public Atom getAtom(int atomIndex) {
		return atomList.get(atomIndex);
	}

	/**
	 * Method to set a particular atom in the lists of atoms.
	 * 
	 * @param atomIndex
	 *            the index in the list, which is to be changed
	 * @param atom
	 *            the instance of atom class
	 * @throws IndexOutOfBoundsException
	 *             If not a valid atomIndex
	 */
	@Override
	public void setAtom(int atomIndex, Atom atom) {
		atomList.set(atomIndex, atom);
	}

	/**
	 * Method to get a particular atom from the lists of atoms.
	 * 
	 * @return an iterator object containing a linear list of atoms in the
	 *         Molecule!
	 */
	@Override
	public Iterator<Atom> getAtoms() {
		return atomList.iterator();
	}

	/**
	 * A way to indicate the implementation that all Atom objects are added and
	 * now the Molecule object could be "cleansed" to minimize memory usage.
	 */
	@Override
	public void pack() {
		atomList.trimToSize();
	}

	/**
	 * Return an atom index that seems similar to the one object passed.
	 * 
	 * @param atom
	 *            the Atom object to be searched for
	 * @return the atoms index which is similar to the Atom object passed (-1 is
	 *         no match found).
	 */
	@Override
	public int soundsLike(Atom atom) {
		int index = -1;

		int i = 0;
		for (Atom atm : atomList) {
			if (atm.distanceFrom(atom) <= 0.0001) {
				if (atm.getSymbol().equals(atom.getSymbol())) {
					index = i;
					break;
				}
			}
			i++;
		}

		return index;
	}

	/**
	 * This method returns true or false indicating whether an atom pair is
	 * bonded.
	 * 
	 * @param atomIndex1
	 *            the reference atom
	 * @param atomIndex2
	 *            the atom index for which the bonding information is queried.
	 * @return true if there is a bond between atomIndex1 and atomIndex2, false
	 *         otherwise.
	 * @throws IndexOutOfBoundsException
	 *             If not a valid atomIndex1.
	 */
	@Override
	public boolean isBonded(int atomIndex1, int atomIndex2) {
		return (atomList.get(atomIndex1).isConnected(atomIndex2));
	}

	/**
	 * This method returns true or false indicating whether an atom pair is
	 * having a bond between them.
	 * 
	 * @param atomIndex1
	 *            the reference atom
	 * @param atomIndex2
	 *            the atom index for which the bonding information is queried.
	 * @return true if there is a strong bond between atomIndex1 and atomIndex2,
	 *         false otherwise.
	 * @throws IndexOutOfBoundsException
	 *             If not a valid atomIndex1.
	 */
	@Override
	public boolean isStronglyBonded(int atomIndex1, int atomIndex2) {
		BondType bt = atomList.get(atomIndex1)
				.getConnectivity(atomIndex2);

		return bt.isStrongBond();
	}

	/**
	 * Method to obtain the bonding type between an atom pair.
	 * 
	 * @param atomIndex1
	 *            the reference atom
	 * @param atomIndex2
	 *            the atom index for which the bonding information is queried.
	 * @return instance of BondType indicating the type of bond.
	 */
	@Override
	public BondType getBondType(int atomIndex1, int atomIndex2) {
		return ((atomList.get(atomIndex1)).getConnectivity(atomIndex2));
	}

	/**
	 * Method to define the bonding type of an atom pair.
	 * 
	 * @param atomIndex1
	 *            the reference atom
	 * @param atomIndex2
	 *            the atom index between which the bonding is to be defined
	 * @param bondType
	 *            An instance of BondType defining the type of the bond
	 * @throws IndexOutOfBoundsException
	 *             If not a valid atomIndex1 or atomIndex2.
	 */
	@Override
	public void setBondType(int atomIndex1, int atomIndex2, BondType bondType) {
		atomList.get(atomIndex1).addConnection(atomIndex2, bondType);
		atomList.get(atomIndex2).addConnection(atomIndex1, bondType);

		// no boolean value set as this doesn't affect the molecules bound
		if (!bondType.equals(BondType.NO_BOND)) {
			msce.setEventType(MoleculeStateChangeEvent.BOND_MODIFIED);
		} else {
			msce.setEventType(MoleculeStateChangeEvent.BOND_REMOVED);
		} // end if

		msce.setAtom1((Atom) atomList.get(atomIndex1));
		msce.setAtom2((Atom) atomList.get(atomIndex2));
		fireMoleculeStateChangeListenerMoleculeChanged(msce);
	}

	/**
	 * Method to remove the bond between an atom pair
	 * 
	 * @param atomIndex1
	 *            the reference atom
	 * @param atomIndex2
	 *            the atom index between which the bonding is to be removed
	 */
	@Override
	public void removeBondBetween(int atomIndex1, int atomIndex2) {
		setBondType(atomIndex1, atomIndex2, BondType.NO_BOND);
	}

	/**
	 * sets the title.
	 * 
	 * @param title
	 *            string identifying the title of the Molecule object
	 */
	@Override
	public void setTitle(String title) {
		this.title = title.trim();

		msce.setEventType(MoleculeStateChangeEvent.TITLE_MODIFIED);
		fireMoleculeStateChangeListenerMoleculeChanged(msce);
	}

	/**
	 * Method to get total number of atoms in the molecule.
	 * 
	 * @return number of atoms in the Molecule
	 */
	@Override
	public int getNumberOfAtoms() {
		return atomList.size();
	}

	/**
	 * Method to return the bonding box enclosing the molecule.
	 * 
	 * @return BoundingBox inclosing the molecule.
	 */
	@Override
	public BoundingBox getBoundingBox() {
		Vector3D ul = null;
		Vector3D br = null;
		if (boundingBox == null || stateChanged) {



			if (atomList.isEmpty()) {
				return new BoundingBox();
			}

			Atom atom = atomList.get(0);

			double xmin = atom.getX();
			double xmax = xmin;
			double ymin = atom.getY();
			double ymax = ymin;
			double zmin = atom.getZ();
			double zmax = zmin;

			double x;
			double y;
			double z;

			for (int i = 1; i < atomList.size(); i++) {
				atom = atomList.get(i);

				x = atom.getX();

				if (x < xmin)
					xmin = x;
				if (x > xmax)
					xmax = x;

				y = atom.getY();

				if (y < ymin)
					ymin = y;
				if (y > ymax)
					ymax = y;

				z = atom.getZ();

				if (z < zmin)
					zmin = z;
				if (z > zmax)
					zmax = z;
			} // end for

			ul = new Vector3D(xmin,ymin,zmin);
			br = new Vector3D(xmax,ymax,zmax);
			
			
		} // end if

		stateChanged = false;

		return new BoundingBox(ul, br);
	}

	/**
	 * Method to return the bonding box center on COM and enclosing the molecule
	 * 
	 * @return BoundingBox inclosing the molecule and centered on COM.
	 */
	@Override
	public BoundingBox getMassCenteredBoundingBox() {
		BoundingBox bb = getBoundingBox();

		if (atomList.isEmpty())
			return bb;

		Vector3D com = getCenterOfMass();

		BoundingBox bborg = getBoundingBox();

		double dist = bborg.center().distance(com);

		Vector3D ul = new Vector3D(bb.getUpperLeft().getX() - dist,
				bb.getUpperLeft().getY() - dist,
				bb.getUpperLeft().getZ() - dist);
		
		Vector3D br = new Vector3D(bb.getBottomRight().getX() + dist,
					bb.getBottomRight().getY() + dist,
					bb.getBottomRight().getZ() + dist);

		return new BoundingBox(ul, br);
	}



	/**
	 * Return an object of molecular formula, pertaining to this molecule
	 * 
	 * @return MolecularFormula the molecular formula of this molecule
	 */
	@Override
	public MolecularFormula getFormula() {
		if (molecularFormula == null) {
			molecularFormula = new MolecularFormula(this);
		} // end if

		return molecularFormula;
	}

	/**
	 * The implementation of this method should update the coordinates of
	 * exixting Atom objects which for this molecule. The new coordinates are
	 * provided in a linear array, of size 3N, with the first three indices
	 * representing coordinates for first atom as (x,y,z) and so on. In case
	 * <code>updateConnectivity</code> is set to true, an appropriate event
	 * should be fired, with the listener taking care of updating the
	 * connectivity information. The implementation of this method should also
	 * update fields dependent on the atomic center such as center of mass and
	 * bounding box etc.
	 * 
	 * @param coords
	 *            the coordinates with size 3N
	 * @param updateConnectivity
	 *            is the connectivity to be updated?
	 * @throws IllegalArgumentException
	 *             If the array size does not match 3N.
	 */
	@Override
	public void resetAtomCoordinates(double[] coords, boolean updateConnectivity)
			throws IllegalArgumentException {
		// first ckeck if the size of coords is correct
		if (coords.length != getNumberOfAtoms() * 3)
			throw new IllegalArgumentException("Size of coords array should "
					+ "be exactly: " + getNumberOfAtoms() * 3
					+ ". You provided an array with size: " + coords.length);

		// note that the state is changed, before its actually changed!
		stateChanged = true;

		// reset center of mass
		centerOfMass = new Vector3D(0,0,0);

		// if its ok, we assume all is well
		double x;
		double y; 
		double z;
		int idxCoord = 0;
		for (Atom atom : atomList) {
			x = coords[idxCoord];
			y = coords[idxCoord + 1];
			z = coords[idxCoord + 2];

			atom.setAtomCenter(new Vector3D(x, y, z));

			double mass = atomInfo.getAtomicWeight(atom.getSymbol());
	
			centerOfMass = new Vector3D(centerOfMass.getX() + (mass * x),
					centerOfMass.getY() + (mass * y),
					centerOfMass.getZ() + (mass * z));

			idxCoord += 3;
		} // end for

		if (updateConnectivity)
			msce.setEventType(MoleculeStateChangeEvent.MAJOR_MODIFICATION);
		else
			msce.setEventType(MoleculeStateChangeEvent.REFRESH_EVENT);

		fireMoleculeStateChangeListenerMoleculeChanged(msce);
	}

	/**
	 * Do a breadth first traversal of the the molecular graph starting with the
	 * specified atomIndex and return the traversed path as an Iterator of atom
	 * indices.
	 * 
	 * @param atomIndex
	 *            the index from where to begin traversal
	 * @return Iterator object containing atom indices visited in that order
	 */
	@Override
	public java.util.Iterator<Integer> traversePath(int atomIndex) {
		ArrayList<Integer> visited = new ArrayList<>(atomList.size());

		// push the current index
		visited.add(Integer.valueOf(atomIndex));

		// then traverse others and record
		traverseAndRecordBFS(atomIndex, visited);

		for (int i = 0; i < atomList.size(); i++) {
			if (visited.size() == atomList.size())
				break;

			traverseAndRecordBFS(i, visited);
		} // end for

		return visited.iterator();
	}

	/**
	 * Do a breadth first traversal of the the molecular graph starting with the
	 * specified atomIndex and return the traversed path as an Iterator of atom
	 * indices.
	 * 
	 * @param atomIndex
	 *            the index from where to begin traversal
	 * @param noDepth
	 *            do not go to higher depths
	 * @return Iterator object containing atom indices visited in that order
	 */
	@Override
	public java.util.Iterator<Integer> traversePath(int atomIndex,
			boolean noDepth) {
		if (!noDepth)
			return traversePath(atomIndex);

		ArrayList<Integer> visited = new ArrayList<>(atomList.size());

		// then traverse others and record
		traverseAndRecordDFSNoWeak(atomIndex, visited);

		return visited.iterator();
	}

	/**
	 * traverse and record the path ... DFS
	 */
	private void traverseAndRecordDFS(int atomIndex,
			ArrayList<Integer> visited, int recursionLevel) {

		if (recursionLevel == 4)
			return;

		Enumeration<Integer> connectList = (Enumeration<Integer>) atomList.get(atomIndex)
				.getConnectedList().keySet();

		while (connectList.hasMoreElements()) {
			Integer toVisit = connectList.nextElement();

			// visited already? then cycle
			if (visited.contains(toVisit))
				continue;

			// push the current index
			visited.add(toVisit);

			// avoid visit of "dead-end" atom
			if (atomList.get(toVisit).getDegree() < 2)
				continue;

			// recursive call
			traverseAndRecordDFS(toVisit.intValue(), visited,
					recursionLevel + 1);
		} // end while
	}

	/**
	 * traverse and record the path ... DFS ... without traversing paths
	 * represented by a weak bond
	 */
	private void traverseAndRecordDFSNoWeak(int atomIndex,
			ArrayList<Integer> visited) {
		PriorityQueue<Integer> toVisitList = new PriorityQueue<>();
		toVisitList.add(atomIndex);

		while (!toVisitList.isEmpty()) {
			Integer toVisit = toVisitList.remove();

			// visited already? then cycle
			if (visited.contains(toVisit))
				continue;

			visited.add(toVisit);

			HashMap<Integer, BondType> connedAtoms = atomList.get(toVisit)
					.getConnectedList();
			Enumeration<Integer> connectList = (Enumeration<Integer>) connedAtoms.keySet();

			while (connectList.hasMoreElements()) {
				Integer visitThis = connectList.nextElement();

				// if already visited dont worry much
				if (visited.contains(visitThis))
					continue;

				// also do no traverse weak bond path
				if (connedAtoms.get(visitThis).equals(BondType.WEAK_BOND))
					continue;

				// mark it to be visted
				toVisitList.add(visitThis);
			}
		}
	}

	/**
	 * traverse and record the path .. BFS
	 */
	private void traverseAndRecordBFS(int atomIndex, ArrayList<Integer> visited) {
		LinkedList<Integer> queue = new LinkedList<>();
		Enumeration<Integer> connectList;

		// do a BFS
		while (true) {
			// traverse through its branches
			connectList = (Enumeration<Integer>) atomList.get(atomIndex).getConnectedList().keySet();

			while (connectList.hasMoreElements()) {
				Integer toVisit = connectList.nextElement();

				// visited already? then cycle
				if (visited.contains(toVisit))
					continue;

				// push the visited index
				visited.add(toVisit);

				// avoid visit of "dead-end" atom
				if (atomList.get(toVisit).getDegree() < 2)
					continue;

				// go to the depth of 4 so as to track the dihedrals
				traverseAndRecordDFS(toVisit, visited, 0);

				// and add to the queue
				queue.add(toVisit);
			} // end while

			// if queue is empty, get out
			if (queue.isEmpty())
				break;

			// else get the next unexplored node
			atomIndex = queue.remove();
		} // end while

		// finally visit atoms with degree 1
		int noOfAtoms = atomList.size();

		for (int i = 0; i < noOfAtoms; i++) {
			connectList = (Enumeration<Integer>) atomList.get(i).getConnectedList().keySet();

			while (connectList.hasMoreElements()) {
				Integer toVisit = connectList.nextElement();

				// visited already? then cycle
				if (visited.contains(toVisit))
					continue;

				// do avoid visit of an atom with degree 1
				if (atomList.get(toVisit).getDegree() > 1)
					continue;

				// push the visited index
				visited.add(toVisit);
			}
		}

	}

	/**
	 * Get a canonical ordering of atoms as a new Molecule object. It is up to
	 * the implementation as to how this ordering is done. Note: for the new
	 * Molecule object, the connectivity must be explicitely build by the
	 * caller.
	 * 
	 * @return a new Molecule object with different ordering
	 */
	@Override
	public Molecule getCanonicalOrdering() {
		Molecule mol = new MoleculeImpl(this.getTitle());

		// make a copy first
		try {
			for (Atom atm : atomList) {
				Atom atmCpy = (Atom) atm.clone();

				mol.addAtom(atmCpy);
				atmCpy.removeAllConnections();
			} // end for
		} catch (Exception e) {
			System.err.println("Exception in getCanonicalOrdering(): "
					+ e.toString());
			e.printStackTrace();
		}

		// get center of mass
		Vector3D com = mol.getCenterOfMass();
		for (Atom atm : atomList) {
			atm.setIndex(mol.closestTo(com));
		}

		return mol;
	}

	/**
	 * Return an atom index that is closest to a given point
	 * 
	 * @param point
	 *            the reference point
	 * @return the atom index that is closest to the point
	 */
	@Override
	public int closestTo(Vector3D point) {
		int index = -1;

		int i = 0;
		double dist = Double.MAX_VALUE;
		for (Atom atm : atomList) {
			if (atm.distanceFrom(point) <= dist) {
				index = i;
				break;
			}
			i++;
		}

		return index;
	}

	/**
	 * Do a simple re-indexing of the atoms
	 */
	@Override
	public void reIndexAtom() {
		int noOfAtoms = getNumberOfAtoms();

		for (int i = 0; i < noOfAtoms; i++) {
			getAtom(i).setIndex(i);
		}
	}
}
