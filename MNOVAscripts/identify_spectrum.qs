// Pulse-sequence lists, experiment names and the subtype a 1D spectrum must have all come
// from SIMPLENMR_EXPERIMENTS in simpleUtils.qs.

function identify_spectrum(){

    const validSpectraKeys = spectra_keys;

    var doc = Application.mainWindow.activeDocument;
    var simpleutils = new simpleUtils();
    var spectra_lst = simpleutils.identifySpectra(doc);


    for( var i=0; i < spectra_lst.length; i++){
        var spectrum = spectra_lst[i];
        var ndims = spectrum.dimCount;
        var pulseSequence = spectrum.getParam("Pulse Sequence");
        var subtype = spectrum.subtype;
        var spectrumFound = false;
        var exptIdentifiedName = "";

        spectrum["experimentIdentified"] = false;
        spectrum["experimentEEH"] = "undefined";

        // try the experiment types in table order; the first match wins
        for (var e = 0; e < SIMPLENMR_EXPERIMENTS.length; e++) {
            var expt = SIMPLENMR_EXPERIMENTS[e];
            if (expt.dims != ndims || !comparePulseSequence(expt.pulseSequences, pulseSequence)) {
                continue;
            }
            if (expt.subtype !== null && expt.subtype != subtype) {
                print("Pulse sequence " + pulseSequence + " does not match subtype " + subtype);
                continue;
            }
            spectrumFound = true;
            exptIdentifiedName = expt.name;
            spectrum["experimentEEH"] = expt.name;
            spectrum["experimentIdentified"] = true;
            break;
        }
        
    }

    // check if any unidentified spectra remain
    // check if class has been used to identify the spectrum
    for (var i=0; i<spectra_lst.length; i++){
        var spectrum = spectra_lst[i];
        if (!spectrum["experimentIdentified"]) {
            var class_str =spectrum.getParam("Class");
            var note_str = spectrum['page']['notes'];
            print("Checking spectrum " + spectrum.getParam("Pulse Sequence") + " class: " + class_str + " note: " + note_str);
            for (var j=0; j < validSpectraKeys.length; j++){
                var ky = validSpectraKeys[j];
                print("Comparing to " + ky);
                if (class_str === ky) {
                    spectrum["experimentIdentified"] = true;
                    spectrum["experimentEEH"] = ky;
                    break;
                }
                if (note_str === ky) {
                    spectrum["experimentIdentified"] = true;
                    spectrum["experimentEEH"] = ky;
                    break;
                }
            }
        }
    }

    // check if there are any unidentifed spectra and if there are ask the user to identify them
    var unidentified_spectra = [];
    print();
    for (var i=0; i<spectra_lst.length; i++){
        var spectrum = spectra_lst[i];
        print("Pulse sequence: " + spectrum.getParam("Pulse Sequence") + ", experimentEEH: " + spectrum["experimentEEH"] + ", experimentIdentified: " + spectrum["experimentIdentified"]);
        if (!spectrum["experimentIdentified"]) {
            unidentified_spectra.push(spectrum);
        }
    }
    print();
    print("unidentified_spectra.length " + unidentified_spectra.length);
    print();
    if (unidentified_spectra.length > 0) {
        idspectra_dialog( spectra_lst );
    }

    for( var i=0; i<spectra_lst.length; i++){
        var spectrum = spectra_lst[i];
        print("spectrum.experimentEEH " + spectrum.experimentEEH + " " + spectrum.getParam("Pulse Sequence"));
    }
    return spectra_lst;
}


// comparePulseSequence() and isHSQC() are in simpleUtils.qs