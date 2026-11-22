// Function to check if one string is a substring of another
// function isSubstring(sub, str) {
//     return str.indexOf(sub) !== -1;
// }

// Pulse-sequence lists, experiment names and the subtype a 1D spectrum must have all come
// from SIMPLENMR_EXPERIMENTS in simpleUtils.qs.

function identify_spectraType(spectra){

    print("identify_spectraType called");
    print("spectra type: " + Object.prototype.toString.call(spectra));

    // Try the experiment types of this spectrum's dimension in table order; first match wins.
    function processSpectrum(spectrum, ndims) {
        const pulseSequence = spectrum.getParam("Pulse Sequence");
        const subtype = spectrum.subtype;
        for (var e = 0; e < SIMPLENMR_EXPERIMENTS.length; e++) {
            var expt = SIMPLENMR_EXPERIMENTS[e];
            if (expt.dims != ndims || !comparePulseSequence(expt.pulseSequences, pulseSequence)) {
                continue;
            }
            if (expt.subtype !== null && expt.subtype != subtype) {
                print("Pulse sequence " + pulseSequence + " does not match subtype " + subtype);
                continue;
            }
            print(pulseSequence + " " + expt.name);
            spectrum.experimentEEH = expt.name;
            spectrum.pulseSeqIdentified = true;
            return true; // Found a match
        }
        return false; // No match found
    }

    // Function to initialize spectrum properties
    function initializeSpectrumProperties(spectrum) {
        spectrum.experimentEEH = EXPT.UNIDENTIFIED;
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


    for( var i=0; i < spectra_lst.length; i++){
        var spectrum = spectra_lst[i];

        // Initialize spectrum properties
        initializeSpectrumProperties(spectrum);

        // Process spectrum based on dimensions
        if(spectrum.dimCount == 1 || spectrum.dimCount == 2){
            processSpectrum(spectrum, spectrum.dimCount);
        }
    }

    for( var i=0; i<spectra_lst.length; i++){
        var spectrum = spectra_lst[i];
        print("spectrum.experimentEEH " + spectrum.experimentEEH + " " + spectrum.getParam("Pulse Sequence"));
    }
    
    // Return single spectrum or array based on input type
    return isArray ? spectra_lst : spectra_lst[0];
}

// comparePulseSequence() and isHSQC() are in simpleUtils.qs




