import org.openscience.cdk.DefaultChemObjectBuilder;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.io.MDLV2000Reader;
import org.openscience.cdk.tools.manipulator.AtomContainerManipulator;

import org.openscience.cdk.interfaces.IAtom;
import java.util.Locale;
import org.openscience.nmrshiftdb.PredictionTool;
import org.openscience.cdk.aromaticity.CDKHueckelAromaticityDetector;
import org.openscience.nmrshiftdb.util.AtomUtils;
import org.openscience.cdk.DefaultChemObjectBuilder;
import org.openscience.cdk.interfaces.IChemObject;
import org.openscience.cdk.interfaces.IAtomContainer;
import java.io.Reader;
import org.openscience.cdk.io.MDLReader;
import java.io.FileReader;

import java.io.StringReader;

public class MolFileParser {
    public static void main(String[] args) {
        if (args.length < 1) {
            System.out.println("Please provide a Molfile string as a command-line argument.");
            return;
        }

        String molfileString = args[0];

        try {
            // Create a reader for the Molfile string
            StringReader stringReader = new StringReader(molfileString);
            MDLV2000Reader reader = new MDLV2000Reader(stringReader);

            // Read the molecule
            IAtomContainer molecule = reader.read(DefaultChemObjectBuilder.getInstance().newAtomContainer());

            // Optionally, you can do further manipulation or analysis of the molecule
            AtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(molecule);

            // Output some basic information about the molecule
            //System.out.println("Number of atoms: " + molecule.getAtomCount());
            //System.out.println("Number of bonds: " + molecule.getBondCount());

            AtomUtils.addAndPlaceHydrogens(molecule);
            CDKHueckelAromaticityDetector.detectAromaticity(molecule);
            final PredictionTool predictor = new PredictionTool();
            System.out.println("C,min,mean,max");
            String solvent = "Unreported";
            boolean use3d = true;
            for (int i = 0; i < molecule.getAtomCount(); ++i) {
                final IAtom curAtom = molecule.getAtom(i);
                if (curAtom.getAtomicNumber() == 6) {
                    final float[] result = predictor.predict(molecule, curAtom, use3d, solvent);
                    System.out.format(Locale.US, "%3d,%8.2f,%8.2f,%8.2f\n", i + 1, result[0], result[1], result[2]);
                }
            }

        } catch (Exception e) {
            e.printStackTrace();
        }
    }
}
