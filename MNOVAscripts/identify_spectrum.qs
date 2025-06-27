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
                                    "shmbcctetgpl2nd"];

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
                                      "zgzrse"
];

var validPURESHIFT1DH1PulseSequenceStrings = ["ja_PSYCHE_pr_03b", 
                                              "reset_psychetse.ptg"];

var validNOESYPulseSequenceStrings = ["noesygpphppzs"];

// Function to check if one string is a substring of another
// function isSubstring(sub, str) {
//     return str.indexOf(sub) !== -1;
// }

function identify_spectrum(){

    const valid1DPulseSequences = {
        "H1_1D": valid1DH1PulseSequenceStrings,
        "C13_1D": valid1DC13PulseSequenceStrings,
        "PURESHIFT": validPURESHIFT1DH1PulseSequenceStrings,
        "DEPT135": validDEPT135PulseSequenceStrings,
    };

    const valid2DPulseSequences = {
        "HSQC": validHSQCPulseSequenceStrings,
        "HMBC": validHMBCPulseSequenceStrings,
        "COSY": validCOSYCPulseSequenceStrings,
        "DOUBLEDEPTCH3": validDOUBLEDEPTCH3PulseSequenceStrings,
        "HSQCCLIPCOSY": validHSQCCLIPCOSYPulseSequenceStrings,
        "NOESY": validNOESYPulseSequenceStrings
    };

    var doc = Application.mainWindow.activeDocument;
    var simpleutils = new simpleUtils();
    var spectra_lst = simpleutils.identifySpectra(doc);

    print(Object.keys(valid2DPulseSequences).length);

    for( var i=0; i < spectra_lst.length; i++){
        var spectrum = spectra_lst[i];
        var ndims = spectrum.dimCount;
        var pulseSequence = spectrum.getParam("Pulse Sequence");
        var subtype = spectrum.subtype;
        var spectrumFound = false;
        var exptIdentifiedName = "";

        if(ndims == 2){
            var pulseSeqIdentified = false;
            for( var ky in valid2DPulseSequences){
                if(comparePulseSequence(valid2DPulseSequences[ky], pulseSequence)){
                    print( pulseSequence + " " + ky + " " + comparePulseSequence( valid2DPulseSequences[ky], pulseSequence) );
                    pulseSeqIdentified = true;
                    spectrumFound = true;
                    exptIdentifiedName = ky;
                    spectrum["experimentEEH"] = ky;
                    break;
                }
            }
        }
        else if(ndims == 1){
            var pulseSeqIdentified = false;
            for( var ky in valid1DPulseSequences){
                if(comparePulseSequence(valid1DPulseSequences[ky], pulseSequence)){
                    if(ky == "H1_1D" && subtype == "1H"){
                        pulseSeqIdentified = true;
                    }
                    else if(ky == "C13_1D" && subtype == "13C"){
                        pulseSeqIdentified = true;
                    }
                    else if(ky == "PURESHIFT"){
                        pulseSeqIdentified = true;
                    }
                    else if(ky == "DEPT135" && subtype == "13C"){
                        pulseSeqIdentified = true;
                    }
                    else{
                        print("Pulse sequence " + pulseSequence + " does not match subtype " + subtype);
                        continue;
                    }
                    print( pulseSequence + " " + ky + " " + comparePulseSequence( valid1DPulseSequences[ky], pulseSequence) );
                    pulseSeqIdentified = true;
                    spectrumFound = true;
                    exptIdentifiedName = ky;
                    spectrum["experimentEEH"] = ky;
                    break;
                }
            }
        }
        
    }

    for( var i=0; i<spectra_lst.length; i++){
        var spectrum = spectra_lst[i];
        print("spectrum.experimentEEH " + spectrum.experimentEEH + " " + spectrum.getParam("Pulse Sequence"));

    }
    return spectra_lst;
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




