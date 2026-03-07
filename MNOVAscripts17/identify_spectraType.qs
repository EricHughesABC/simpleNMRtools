// Function to check if one string is a substring of another
// function isSubstring(sub, str) {
//     return str.indexOf(sub) !== -1;
// }

var validHSQCPulseSequenceStrings =  ["hsqcedetgpsisp2.3.ptg",
                                      "hsqcedetgpsisp2.3",
                                      "gHSQCAD", 
                                      "hsqcedetgpsp.3",
                                      "gHSQC",
                                      "inv4gp.wu",
                                      "hsqcetgp",
                                      "gns_noah3-BSScc.eeh",
                                      "hsqcedetgpsisp2.4"];

var validHMBCPulseSequenceStrings = ["ghmbc.wu", 
                                     "gHMBC", 
                                     "hmbcetgpl3nd", 
                                     "hmbcetgpl3nd.ptg",
                                     "gHMBCAD",
                                     "hmbcgpndqf",
                                     "gns_noah3-BSScc.eeh",
                                     "shmbcctetgpl2nd",
                                     "hmbcgplpndqf"];

var validCOSYCPulseSequenceStrings = ["cosygpqf", 
                                      "cosygp", 
                                      "gcosy", 
                                      "cosygpmfppqf", 
                                      "cosygpmfqf", 
                                      "gCOSY",
                                      "cosygpppqf_ptype",
                                      "cosyqf45",
                                      "cosygpmfphpp"];

var validDEPT135PulseSequenceStrings = ["dept135.wu", "DEPT", "deptsp135"];

var validDOUBLEDEPTCH3PulseSequenceStrings = ["hcdeptedetgpzf"];

var validHSQCCLIPCOSYPulseSequenceStrings = ["hsqc_clip_cosy_mc_notation.eeh",
                                             "gns_noah3-BSScc.eeh"
];

var valid1DH1PulseSequenceStrings = ["zg30",  
                                     "s2pul", 
                                     "zg", 
                                     "zgcppr",
                                     ];

var valid1DC13PulseSequenceStrings = ["zgdc30", "s2pul", "zgpg30",
                                      "zgzrse", "zg0dc.fr"
];

var validPURESHIFT1DH1PulseSequenceStrings = ["ja_PSYCHE_pr_03b", 
                                              "reset_psychetse.ptg"];

var validNOESYPulseSequenceStrings = ["noesygpphppzs"];



// Module-level constants - created once when script loads
const VALID_1D_PULSE_SEQUENCES = {
    "H1_1D": valid1DH1PulseSequenceStrings,
    "C13_1D": valid1DC13PulseSequenceStrings,
    "H1_pureshift": validPURESHIFT1DH1PulseSequenceStrings,
    "DEPT135": validDEPT135PulseSequenceStrings,
};

const VALID_2D_PULSE_SEQUENCES = {
    "HSQC": validHSQCPulseSequenceStrings,
    "HMBC": validHMBCPulseSequenceStrings,
    "COSY": validCOSYCPulseSequenceStrings,
    "DDEPT_CH3_ONLY": validDOUBLEDEPTCH3PulseSequenceStrings,
    "HSQC_CLIPCOSY": validHSQCCLIPCOSYPulseSequenceStrings,
    "NOESY": validNOESYPulseSequenceStrings
};

function identify_spectraType(spectra){

    print("identify_spectraType called");
    print("spectra type: " + Object.prototype.toString.call(spectra));

    // Function to handle 2D spectrum identification
    function process2DSpectrum(spectrum) {
        const pulseSequence = spectrum.getParam("Pulse Sequence");
        
        for( var ky in VALID_2D_PULSE_SEQUENCES){
            if(comparePulseSequence(VALID_2D_PULSE_SEQUENCES[ky], pulseSequence)){
                print( pulseSequence + " " + ky + " " + comparePulseSequence( VALID_2D_PULSE_SEQUENCES[ky], pulseSequence) );
                spectrum.experimentEEH = ky;
                spectrum.pulseSeqIdentified = true;
                return true; // Found a match
            }
        }
        return false; // No match found
    }

    // Function to handle 1D spectrum identification
    function process1DSpectrum(spectrum) {
        const pulseSequence = spectrum.getParam("Pulse Sequence");
        const subtype = spectrum.subtype;
        
        for( var ky in VALID_1D_PULSE_SEQUENCES){
            if(comparePulseSequence(VALID_1D_PULSE_SEQUENCES[ky], pulseSequence)){
                var pulseSeqIdentified = false;
                
                if(ky == "H1_1D" && subtype == "1H"){
                    pulseSeqIdentified = true;
                }
                else if(ky == "C13_1D" && subtype == "13C"){
                    pulseSeqIdentified = true;
                }
                else if(ky == "H1_pureshift"){
                    pulseSeqIdentified = true;
                }
                else if(ky == "DEPT135" && subtype == "13C"){
                    pulseSeqIdentified = true;
                }
                else{
                    print("Pulse sequence " + pulseSequence + " does not match subtype " + subtype);
                    continue;
                }
                
                print( pulseSequence + " " + ky + " " + comparePulseSequence( VALID_1D_PULSE_SEQUENCES[ky], pulseSequence) );
                spectrum.experimentEEH = ky;
                spectrum.pulseSeqIdentified = pulseSeqIdentified;
                return true; // Found a match
            }
        }
        return false; // No match found
    }

    // Function to initialize spectrum properties
    function initializeSpectrumProperties(spectrum) {
        spectrum.experimentEEH = "unidentified";
        spectrum.pulseSeqIdentified = false;
        spectrum.hasPeaks = spectrum.peaks().count > 0;
        spectrum.hasMultiplets = spectrum.multiplets().count > 0;
        spectrum.hasIntegrals = spectrum.integrals().count > 0;
    }

    // var doc = Application.mainWindow.activeDocument;
    // var simpleutils = new simpleUtils();
    
    // Handle both single spectrum and array of spectra (ES5 compatible)
    var isArray = Object.prototype.toString.call(spectra) === '[object Array]';
    var spectra_lst = isArray ? spectra : [spectra];

    print(Object.keys(VALID_2D_PULSE_SEQUENCES).length);

    for( var i=0; i < spectra_lst.length; i++){
        var spectrum = spectra_lst[i];

        // Initialize spectrum properties
        initializeSpectrumProperties(spectrum);

        // Process spectrum based on dimensions
        if(spectrum.dimCount == 2){
            process2DSpectrum(spectrum);
        }
        else if(spectrum.dimCount == 1){
            process1DSpectrum(spectrum);
        }
    }

    for( var i=0; i<spectra_lst.length; i++){
        var spectrum = spectra_lst[i];
        print("spectrum.experimentEEH " + spectrum.experimentEEH + " " + spectrum.getParam("Pulse Sequence"));
    }
    
    // Return single spectrum or array based on input type
    return isArray ? spectra_lst : spectra_lst[0];
}

function comparePulseSequence(validPulseSequences, pulseSequence){
    var pulseSeqIdentified = false;
    for (var i=0; i<validPulseSequences.length; i++){
        var pulseSeq = validPulseSequences[i];
        // check if the pulse sequence string contains the valid pulse sequence string
        if(isSubstring(pulseSeq, pulseSequence)){
            pulseSeqIdentified = true;
            break;
        }
    }
    return pulseSeqIdentified;
}

function isHSQC(spectrum) {

    var simpleutils = new simpleUtils();
    var pulseSequence = spectrum.getParam("Pulse Sequence");
    var isHSQC = false;

    // if not 2d return false
    print("spectrum dimCount " + spectrum.dimCount);
    if (spectrum.dimCount != 2) {
        return false;
    }

    for (var i = 0; i < validHSQCPulseSequenceStrings.length; i++) {
        var hsqcString = validHSQCPulseSequenceStrings[i];
        // check if the pulse sequence string contains the valid pulse sequence string
        if (isSubstring(hsqcString, pulseSequence)) {
            isHSQC = true;
            break;
        }
    }

    return isHSQC;
}




