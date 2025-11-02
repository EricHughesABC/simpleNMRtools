import org.openscience.cdk.DefaultChemObjectBuilder;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.io.MDLV2000Reader;
import org.openscience.cdk.tools.manipulator.AtomContainerManipulator;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.nmrshiftdb.PredictionTool;
import org.openscience.cdk.aromaticity.CDKHueckelAromaticityDetector;
import org.openscience.nmrshiftdb.util.AtomUtils;

import java.io.StringReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.Locale;

/**
 * NMR Shift Predictor that can be called from Python via JPype
 */
public class MolFileParser {
    
    private PredictionTool predictor;
    private String defaultSolvent;
    private boolean use3d;
    private boolean isInitialized;
    
    /**
     * Constructor - initializes the predictor
     */
    public MolFileParser() {
        this.defaultSolvent = "Unreported";
        this.use3d = false;
        this.isInitialized = false;
        
        try {
            this.predictor = new PredictionTool();
            this.isInitialized = true;
        } catch (IOException e) {
            System.err.println("Failed to initialize PredictionTool: " + e.getMessage());
            e.printStackTrace();
            // Don't throw exception from constructor - set flag instead
        }
    }
    
    /**
     * Constructor with custom settings
     */
    public MolFileParser(String solvent, boolean use3d) {
        this.defaultSolvent = solvent;
        this.use3d = use3d;
        this.isInitialized = false;
        
        try {
            this.predictor = new PredictionTool();
            this.isInitialized = true;
        } catch (IOException e) {
            System.err.println("Failed to initialize PredictionTool: " + e.getMessage());
            e.printStackTrace();
        }
    }
    
    /**
     * Check if predictor is ready to use
     */
    public boolean isInitialized() {
        return this.isInitialized;
    }
    
    /**
     * Prediction result for a single atom
     */
    public static class PredictionResult {
        private int atomIndex;
        private float min;
        private float mean;
        private float max;
        
        public PredictionResult(int atomIndex, float min, float mean, float max) {
            this.atomIndex = atomIndex;
            this.min = min;
            this.mean = mean;
            this.max = max;
        }
        
        public int getAtomIndex() { return atomIndex; }
        public float getMin() { return min; }
        public float getMean() { return mean; }
        public float getMax() { return max; }
        
        @Override
        public String toString() {
            return String.format(Locale.US, "%d,%.2f,%.2f,%.2f", 
                               atomIndex, min, mean, max);
        }
    }
    
    /**
     * Parse molfile and predict NMR shifts
     * Returns a list of prediction results
     */
    public List<PredictionResult> predictFromMolfile(String molfileString) throws Exception {
        if (!isInitialized) {
            throw new IllegalStateException("PredictionTool failed to initialize. Check if required data files are present.");
        }
        
        List<PredictionResult> results = new ArrayList<>();
        
        if (molfileString == null || molfileString.trim().isEmpty()) {
            throw new IllegalArgumentException("Molfile string cannot be null or empty");
        }
        
        StringReader stringReader = null;
        MDLV2000Reader reader = null;
        
        try {
            // Create a reader for the Molfile string
            stringReader = new StringReader(molfileString);
            reader = new MDLV2000Reader(stringReader);
            
            // Read the molecule
            IAtomContainer molecule = reader.read(
                DefaultChemObjectBuilder.getInstance().newAtomContainer()
            );
            
            // Configure atoms
            AtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(molecule);
            
            // Add hydrogens and detect aromaticity
            AtomUtils.addAndPlaceHydrogens(molecule);
            CDKHueckelAromaticityDetector.detectAromaticity(molecule);
            
            // Predict for each carbon atom
            for (int i = 0; i < molecule.getAtomCount(); ++i) {
                final IAtom curAtom = molecule.getAtom(i);
                if (curAtom.getAtomicNumber() == 6) {  // Carbon atoms only
                    final float[] result = predictor.predict(
                        molecule, curAtom, use3d, defaultSolvent
                    );
                    results.add(new PredictionResult(
                        i + 1,      // 1-indexed atom number
                        result[0],  // min
                        result[1],  // mean
                        result[2]   // max
                    ));
                }
            }
            
        } catch (Exception e) {
            throw new Exception("Error parsing molfile or predicting shifts: " + e.getMessage(), e);
        } finally {
            // Clean up resources
            if (reader != null) {
                try {
                    reader.close();
                } catch (IOException e) {
                    // Ignore close errors
                }
            }
        }
        
        return results;
    }
    
    /**
     * Parse molfile and return results as CSV string
     * Useful for simple Python integration
     */
    public String predictFromMolfileAsCSV(String molfileString) throws Exception {
        List<PredictionResult> results = predictFromMolfile(molfileString);
        
        StringBuilder csv = new StringBuilder();
        csv.append("AtomIndex,Min,Mean,Max\n");
        
        for (PredictionResult result : results) {
            csv.append(result.toString()).append("\n");
        }
        
        return csv.toString();
    }
    
    /**
     * Get the number of carbon atoms that would be predicted
     * Useful for quick checks without full prediction
     */
    public int getCarbonCount(String molfileString) throws Exception {
        if (molfileString == null || molfileString.trim().isEmpty()) {
            throw new IllegalArgumentException("Molfile string cannot be null or empty");
        }
        
        StringReader stringReader = null;
        MDLV2000Reader reader = null;
        
        try {
            stringReader = new StringReader(molfileString);
            reader = new MDLV2000Reader(stringReader);
            
            IAtomContainer molecule = reader.read(
                DefaultChemObjectBuilder.getInstance().newAtomContainer()
            );
            
            int carbonCount = 0;
            for (int i = 0; i < molecule.getAtomCount(); ++i) {
                if (molecule.getAtom(i).getAtomicNumber() == 6) {
                    carbonCount++;
                }
            }
            
            return carbonCount;
            
        } catch (Exception e) {
            throw new Exception("Error parsing molfile: " + e.getMessage(), e);
        } finally {
            if (reader != null) {
                try {
                    reader.close();
                } catch (IOException e) {
                    // Ignore close errors
                }
            }
        }
    }
    
    /**
     * Set the solvent for predictions
     */
    public void setSolvent(String solvent) {
        this.defaultSolvent = solvent;
    }
    
    /**
     * Get current solvent setting
     */
    public String getSolvent() {
        return this.defaultSolvent;
    }
    
    /**
     * Set whether to use 3D coordinates
     */
    public void setUse3D(boolean use3d) {
        this.use3d = use3d;
    }
    
    /**
     * Get current 3D setting
     */
    public boolean getUse3D() {
        return this.use3d;
    }
    
    /**
     * Main method for command-line compatibility (optional)
     */
    public static void main(String[] args) {
        if (args.length < 1) {
            System.out.println("Please provide a Molfile string as a command-line argument.");
            return;
        }

        String molfileString = args[0];
        
        try {
            MolFileParser parser = new MolFileParser();
            
            if (!parser.isInitialized()) {
                System.err.println("Failed to initialize parser. Check if NMRShiftDB data files are present.");
                System.exit(1);
            }
            
            String results = parser.predictFromMolfileAsCSV(molfileString);
            System.out.println(results);
        } catch (Exception e) {
            System.err.println("Error: " + e.getMessage());
            e.printStackTrace();
            System.exit(1);
        }
    }
}