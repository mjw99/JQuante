package name.mjw.jquante.common.resource;

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Paths;

import java.security.DigestInputStream;
import java.security.MessageDigest;
import java.security.NoSuchAlgorithmException;
import java.text.DecimalFormat;
import java.text.DecimalFormatSymbols;
import java.text.ParseException;
import java.util.ArrayList;
import java.util.Locale;

import java.util.List;

import javax.xml.bind.DatatypeConverter;
import javax.xml.bind.annotation.XmlAttribute;
import javax.xml.bind.annotation.XmlElement;
import javax.xml.bind.annotation.XmlRootElement;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

import name.mjw.fortranformat.FortranFormat;
import name.mjw.jquante.config.impl.AtomInfo;

@XmlRootElement(name = "basis")
public class NWChemBasisSetFile {
	private static final Logger LOG = LogManager
			.getLogger(NWChemBasisSetFile.class);

	private static final int MARK_LENGTH = 5000;

	@XmlAttribute(name = "name")
	private String basisSetName;

	@XmlAttribute(name = "md5sum")
	private String md5sumOfInputBasisFile;

	@XmlElement(name = "atom")
	private List<Library> atoms;

	/**
	 * A representation of a NWChem basis file.
	 * 
	 * 
	 * @see <a
	 *      href="https://github.com/jeffhammond/nwchem/blob/HEAD/src/basis/basis_dox.c">src/basis/basis_dox.c</a>
	 * @see <a
	 *      href="http://www.ccl.net/cca/documents/basis-sets/basis.html">General
	 *      basis set naming information</a>
	 * 
	 */
	public NWChemBasisSetFile() {
		atoms = new ArrayList<>();
	}

	/**
	 * Reads a NWChem basis file.
	 * 
	 * @param fileName
	 *            Filename of NWChem basis file.
	 * 
	 * @see <a
	 *      href="https://github.com/jeffhammond/nwchem/blob/HEAD/src/basis/bas_input.F">src/basis/bas_input.F</a>
	 * 
	 * 
	 */
	public void read(String fileName) {

		File f = new File(fileName);
		basisSetName = f.getName();
		MessageDigest md = null;

		try {
			md = MessageDigest.getInstance("MD5");
		} catch (NoSuchAlgorithmException e) {
			LOG.warn(e);
		}

		// https://stackoverflow.com/questions/304268/getting-a-files-md5-checksum-in-java
		try (InputStream is = Files.newInputStream(Paths.get(fileName));
				DigestInputStream dis = new DigestInputStream(is, md);
				BufferedReader in = new BufferedReader(new InputStreamReader(
						dis, StandardCharsets.UTF_8))) {

			String line;

			while ((line = in.readLine()) != null) {

				if (line.matches("^basis.*")) {

					parseBasisSection(line, in);

				}

			}

			byte[] digest = md.digest();
			md5sumOfInputBasisFile = DatatypeConverter.printHexBinary(digest)
					.toLowerCase();

		} catch (IOException | ParseException e) {
			LOG.warn(e);
		}

	}

	/**
	 * Parse the section:
	 * 
	 * <pre>
	 * basis "He_aug-cc-pVTZ" SPHERICAL 
	 * ... 
	 * end
	 * </pre>
	 * 
	 * 
	 * @param inputLine
	 * @param br
	 * @throws ParseException
	 * @throws IOException
	 */
	void parseBasisSection(String inputLine, BufferedReader br)
			throws ParseException, IOException {

		String elementName;
		String []words;

		words = inputLine.split("\\s+");

		// basis "Na_cc-pVTZ" SPHERICAL
		LOG.debug(words[0] + " " + words[1] + " " + words[2]);

		if (!(words[0].contentEquals("basis"))) {
			throw new ParseException("Expecting basis term", 1);
		}

		elementName = words[1].replace("\"", "").split("_")[0];

		Library atom = new Library(elementName);

		List<Shell> shells = new ArrayList<>();

		br.mark(MARK_LENGTH);
		while ((inputLine = br.readLine()) != null) {

			if (inputLine.matches("^end.*")) {

				br.reset();
				break;

			} else {

				LOG.debug("dispatching to parseShellSection() with line ");
				LOG.debug(inputLine);

				br.reset();
				shells.addAll(parseShellSection(br, elementName));

			}

			br.mark(MARK_LENGTH);

		}

		atom.addShell(shells);
		atoms.add(atom);
	}

	/**
	 * Parse the section:
	 * 
	 * <pre>
	 * H    S
	 *      33.8700000              0.0060680
	 *       5.0950000              0.0453080
	 *       1.1590000              0.2028220
	 * </pre>
	 * 
	 * 
	 * @param line
	 * @param br
	 * @throws ParseException
	 * @throws IOException
	 */
	static List<Shell> parseShellSection(BufferedReader br, String elementName)
			throws ParseException, IOException {

		int numberOfContractedShellsToParse = 1;
		Double exp;
		Double coeff;
		String line;
		ArrayList<Shell> shells = null;
		ArrayList<Object> lineObjects;
		final String shellFormat = "(1X, E18.7, 9X, 6(E15.7, 8X))";

		DecimalFormat df = new DecimalFormat("0",
				DecimalFormatSymbols.getInstance(Locale.ENGLISH));
		df.setMaximumFractionDigits(340); // 340 =
		// DecimalFormat.DOUBLE_FRACTION_DIGITS

		line = br.readLine();
		LOG.debug("");
		LOG.debug("parseShellSection() is starting on line: ");
		LOG.debug(line);

		lineObjects = FortranFormat.read(line, "(2A5)");

		String shellType = (String) lineObjects.get(1);

		// Peek ahead here to work out how many contractedShells to return
		br.mark(MARK_LENGTH);
		line = br.readLine();

		// Hack since FortranFormat does not
		// understand Double-Precision Exponent
		line = line.replace("D", "E");
		lineObjects = FortranFormat.read(line, shellFormat);
		numberOfContractedShellsToParse = getNonNullLength(lineObjects) - 1;
		br.reset();

		shells = new ArrayList<>(numberOfContractedShellsToParse);

		// For the basis sets in which s- and p-type functions share the same
		// exponents, the term SP-shell is used.
		if (shellType.equals("SP")) {
			LOG.debug("SP Found");

			for (int i = 0; i < numberOfContractedShellsToParse; i++) {
				LOG.debug(ShellType.values()[i]);
				shells.add(new Shell(ShellType.values()[i].toString()));

			}
		} else {
			for (int i = 0; i < numberOfContractedShellsToParse; i++) {
				shells.add(new Shell(shellType));
			}
		}

		while (true) {
			line = br.readLine();

			if (line.matches("^" + elementName + ".*") || line.matches("^end.*")) {

				LOG.debug("Seen line: " + line);
				LOG.debug("New shell time..");
				LOG.debug("");

				br.reset();
				return shells;
			}

			// Hack since FortranFormat does not understand Double-Precision
			// Exponent
			line = line.replace("D", "E");

			lineObjects = FortranFormat.read(line, shellFormat);
			LOG.debug(lineObjects);

			exp = (Double) lineObjects.get(0);

			for (int i = 0; i < numberOfContractedShellsToParse; i++) {

				coeff = (Double) lineObjects.get(i + 1);

				shells.get(i).addExponentCoefficientPair(df.format(exp),
						df.format(coeff));
			}

			br.mark(MARK_LENGTH);
		}

	}

	static class Library {
		@XmlAttribute(name = "symbol")
		String element;

		@XmlAttribute
		int atomicNumber;

		@XmlElement(name = "orbital")
		List<Shell> shells;

		Library(String element) {
			this.element = element;
			atomicNumber = AtomInfo.getInstance().getAtomicNumber(element);
			shells = new ArrayList<>();
		}

		void addShell(List<Shell> orbitals) {
			this.shells.addAll(orbitals);
		}

	}

	static class Shell {

		@XmlAttribute(name = "type")
		String type;

		@XmlElement(name = "entry")
		List<ExponentCoefficientPair> exponentCoefficientPairs = null;

		Shell(String type) {
			this.type = type;
			exponentCoefficientPairs = new ArrayList<>();
		}

		void addExponentCoefficientPair(String exponent, String coefficient) {

			ExponentCoefficientPair pair = new ExponentCoefficientPair(
					exponent, coefficient);

			exponentCoefficientPairs.add(pair);

		}

	}

	static class ExponentCoefficientPair {

		@XmlAttribute(name = "exp")
		String exponent;

		@XmlAttribute(name = "coeff")
		String coefficient;

		ExponentCoefficientPair(String exponent, String coefficient) {
			this.exponent = exponent;
			this.coefficient = coefficient;
		}

	}

	/**
	 * Count the number of non-null elements in an array
	 * 
	 * @param lineObjects
	 * @return number of non-null elements
	 */
	private static int getNonNullLength(ArrayList<Object> lineObjects) {
		int count = 0;
		for (Object el : lineObjects)
			if (el != null)
				++count;
		return count;
	}

	private enum ShellType {
		S(0), P(1), D(2), F(3), G(4), H(5);

		private final int value;

		private ShellType(int value) {
			this.value = value;
		}
	}

	public String getBasisSetName() {
		return basisSetName;
	}
}
