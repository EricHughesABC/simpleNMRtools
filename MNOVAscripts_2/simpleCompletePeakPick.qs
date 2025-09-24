function simpleCompletePeakPick() {

    // get simpleutils functions
    var simpleutils = new simpleUtils();



    // //  get the active 2D spectrum
    // var nmr2Dspec = new NMRSpectrum(nmr.activeSpectrum());

    var doc = Application.mainWindow.activeDocument;
    var simpleutils = new simpleUtils();
    var spectra_lst = simpleutils.identifySpectra(doc);
    spectra_lst = identify_spectraType(spectra_lst);

    var nmr2Dspec = new NMRSpectrum(nmr.activeSpectrum());



    // if not 2D spectrum, exit with message
    if(!simpleutils.is2DSpectrum(nmr2Dspec)){
        return;
    }
    
    // find all the spectra in the document
    // identify the spectra in the document
    var currentSpectrum = identify_spectraType(nmr2Dspec);

    // print out the active spectrum type
    print("\nActive spectrum type: " + currentSpectrum["experimentEEH"]);
    print();

    const expt_types = ["H1_1D", "C13_1D", "PURESHIFT", "DEPT135", "HSQC", "HMBC", "COSY", "NOESY", "DOUBLEDEPTCH3", "HSQCCLIPCOSY"];

    var spectra_dict = {};

    for( var i=0; i<expt_types.length; i++){
        var spectrum = find_NMR_spectrum(spectra_lst, expt_types[i]);
        if (spectrum) {
          spectra_dict[expt_types[i]] = spectrum;
        }
    }

    // print(spectra_dict);

    print("\nSpectra summary:\n");
    for ( var key in spectra_dict){
      print("\t" + key + ": has peaks? " + spectra_dict[key]["hasPeaks"] + ", has multiplets? " + spectra_dict[key]["hasMultiplets"] + ", has integrals? " + spectra_dict[key]["hasIntegrals"]);
    }

    // dictionary of spectra with peaks
    var spectra_with_peaks = {};
    for (var key in spectra_dict) {
        if (spectra_dict[key]["hasPeaks"]) {
            spectra_with_peaks[key] = spectra_dict[key];
        }
    }

    // print(spectra_with_peaks);
    print("\nSpectra with peaks:\n");
    for( var key in spectra_with_peaks){
      print("\t" + key);
    }

    print( spectra_with_peaks.C13_1D.hasPeaks );

    var peakPickingPath = 0;

    if (spectra_with_peaks.C13_1D.hasPeaks) peakPickingPath |= 1;
    if (spectra_with_peaks.PURESHIFT.hasPeaks) peakPickingPath |= 2;

    if (spectra_with_peaks.HSQC && spectra_with_peaks.HSQC.hasPeaks) peakPickingPath |= 4;

    print( "peakPickingPath: " + peakPickingPath );

    switch(peakPickingPath){

        case 1:
            print("Using 1D 13C spectrum only for peak picking.");
            break;
        case 2:
            print("Using 1D PureShift spectrum only for peak picking.");
            break;
        case 3:
            print("Using both 1D 13C and PureShift spectra for peak picking.");
            break;

        case 4:
            print("Using 2D HSQC spectrum only for peak picking.");
            break;

        case 5:
            print("Using 1D 13C and 2D HSQC spectra for peak picking.");
            break;
    // default:
        default:
            print("No case matched for peakPickingPath: " + peakPickingPath);
            break;
    }


    // get active spectrum
    var activeSpectrum = new NMRSpectrum(nmr.activeSpectrum());

    // check if active spectrum is 2D
    if(!simpleutils.is2DSpectrum(activeSpectrum)){
        print("Active spectrum is not 2D. Please select a 2D spectrum and try again.");
        return;
    }

    
}


/**
 * Finds and returns the first 1D 13C spectrum object from a list of spectra.
 *
 * @param {Array<Object>} spectra_lst - List of spectrum objects to search.
 * @param {string} experiment - The experiment type to search for.
 * @returns {Object|undefined} The first spectrum object with "experimentEEH" equal to the specified experiment, or undefined if not found.
 */
function find_NMR_spectrum( spectra_lst, experiment ) {
    var nmr_spectrum = undefined;

    // search field "experimentEEH" for the specified experiment
    for (var i=0; i<spectra_lst.length; i++){
        if (spectra_lst[i]["experimentEEH"] == experiment) {
            nmr_spectrum = spectra_lst[i];
            break;
        }
    }
    return nmr_spectrum;
}


