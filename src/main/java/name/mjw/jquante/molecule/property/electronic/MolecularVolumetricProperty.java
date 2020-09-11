package name.mjw.jquante.molecule.property.electronic;

import java.util.ArrayList;

import name.mjw.jquante.molecule.Molecule;

/**
 * Represents one or many volumetric properties connected to a molecule object.
 * 
 * @author V.Ganesh
 * @version 2.0 (Part of MeTA v2.0)
 */
public class MolecularVolumetricProperty {

	protected Molecule molecule;

	protected ArrayList<VolumeItem> volumeItemList;

	/**
	 * Creates a new instance of MolecularVolumetricProperty
	 * 
	 * @param molecule
	 *            Molecule object.
	 */
	public MolecularVolumetricProperty(final Molecule molecule) {
		this.molecule = molecule;

		this.volumeItemList = new ArrayList<>();
	}

	/**
	 * Get the value of molecule
	 * 
	 * @return the value of molecule
	 */
	public Molecule getMolecule() {
		return molecule;
	}

	/**
	 * Add a new VolumeItem object to this property
	 * 
	 * @param volumeItem instance of new VolumeItem
	 */
	public void addVolumeItem(final VolumeItem volumeItem) {
		volumeItem.setReferenceMolecule(molecule);
		volumeItemList.add(volumeItem);
	}

	/**
	 * Remove an existing VolumeItem from this property
	 * 
	 * @param volumeItem instance of VolumeItem to be removed
	 */
	public void removeVolumeItem(final VolumeItem volumeItem) {
		volumeItem.setReferenceMolecule(null);
		volumeItemList.remove(volumeItem);
	}
}
