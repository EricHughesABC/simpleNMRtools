


/**
 * Checks if a given substring exists within a string.
 *
 * @param {string} sub - The substring to search for.
 * @param {string} str - The string to search within.
 * @returns {boolean} True if the substring is found in the string, otherwise false.
 */
function isSubstring(sub, str) {
    return str.indexOf(sub) !== -1;
}

spectra_keys = ["HSQC", 
                    "HMBC", 
                    "COSY", 
                    "NOESY", 
                    "H1_1D", 
                    "C13_1D", 
                    "DEPT135", 
                    "PureShift", 
                    "DDEPT_CH3_ONLY", 
                    "SKIP", 
                    "HSQC_CLIPCOSY"];
function simpleUtils() {

    var doc = Application.mainWindow.activeDocument;

    // simpleUtils.spectra_keys = ["HSQC", 
    //                     "HMBC", 
    //                     "COSY", 
    //                     "NOESY", 
    //                     "H1_1D", 
    //                     "C13_1D", 
    //                     "DEPT135", 
    //                     "PureShift", 
    //                     "DDEPT_CH3_ONLY", 
    //                     "SKIP", 
    //                     "HSQC_CLIPCOSY"];

}

simpleUtils.prototype.isObjectEmpty = function(obj) {
    return Object.keys(obj).length === 0;
}

simpleUtils.prototype.getActiveMolecule = function (aDocWin, aMolPlugin) {
    var molec = aMolPlugin.activeMolecule();

    if (molec.isValid()) {
            return molec;
    }
    if (aDocWin.itemCount("Molecule") === 1) {
            molec = new Molecule(aDocWin.item(0, "Molecule"));
        return molec;
    }
    return undefined;
}	

simpleUtils.prototype.splitPathAndFilename = function (path) {
    var separatorIndex = path.lastIndexOf('/');
    var directoryPath = path.substring(0, separatorIndex);
    var filename = path.substring(separatorIndex + 1);

    var dotIndex = filename.lastIndexOf('.');
    var name = filename.substring(0, dotIndex);
    var extension = filename.substring(dotIndex + 1);

    return {
        directoryPath: directoryPath,
        filename: filename,
        name: name,
        extension: extension
    };
}

        // var json_spectra_string = JSON.stringify(mnova_to_iupac, null, 4);
        // var fout1 = new File(jsonfilename);
        // if (fout1.open(File.WriteOnly)) {
        //     sout1 = new TextStream(fout1, 'UTF-8');
        //     sout1.writeln(json_spectra_string);
        // }
        // fout1.close();   

simpleUtils.prototype.writeAsJsonToFile = function (data, filename) {

    print("writing to file ", filename);
    var json_spectra_string = JSON.stringify(data, null, 4);
    var fout1 = new File(filename);
    if( fout1.open(File.WriteOnly)){
        var sout1 = new TextStream(fout1, 'UTF-8');
        sout1.writeln(json_spectra_string);
    }
    else{
        print("Error opening file for writing\n", filename);
        messageBox.warning("Error opening file for writing\n" +  filename);
    }
    fout1.close();
}

simpleUtils.prototype.readJsonFromFile = function (filename) {
    var fin = new File(filename);
    var rtn = undefined;
    if( fin.open(File.ReadOnly)){
        var sin = new TextStream(fin, 'UTF-8');
        var json_string = sin.readAll();
        fin.close();
        rtn = JSON.parse(json_string);
    }
    else{
        print("Error opening file for reading\n", filename);
        messageBox.warning("Error opening file for reading\n" + filename);
    }
    return rtn;
}

simpleUtils.prototype.identifySpectra = function(doc) {

    var spectra = [];
    for (i = 0, pageCount = doc.pageCount(); i < pageCount; i++) {
        page = doc.page(i);
        for (j = 0, itemCount = page.itemCount(); j < itemCount; j++) {
            spec = new NMRSpectrum(page.item(j));
      
            if( spec.isValid() ){
                spectra.push(spec);
            }
        }
    }
    return spectra;
}

simpleUtils.prototype.findC13spectrum = function( spectra_lst ) {
    // find the first 1D 13C spectrum in the list of spectra
    // parameters
    //    spectra_lst:  list of spectra
    // returns
    //    c13spectrum: the spectrum object

    var c13spectrum = undefined;
    for (var i=0; i<spectra_lst.length; i++){
        if (spectra_lst[i].subtype == "predicted, 13C"){
            // skip and c13 spectrum that is predicted
            continue
        }
        if (spectra_lst[i].dimCount == 1 && spectra_lst[i].getParam("Nucleus") == "13C") {
            c13spectrum = spectra_lst[i];
            // print("C13spectrum index = ", i)
            // print("spectrum dimCount & subtype", spectra_lst[i].dimCount, spectra_lst[i].subtype, spectra_lst[i].getParam("Nucleus"));
            break;
        }
    }
    return c13spectrum;
}

simpleUtils.prototype.findPureShiftSpectrum = function( spectra_lst ) {
    // find the first 1D PureShift spectrum in the list of spectra
    // parameters
    //    spectra_lst:  list of spectra
    // returns
    //    psSpectrum: the spectrum object

    var psSpectrum = undefined;
    for (var i=0; i<spectra_lst.length; i++){
        var pulseSequence = spectra_lst[i].getParam("Pulse Sequence");
        // change string to lowercase
        pulseSequence = pulseSequence.toLowerCase();

        if (spectra_lst[i].dimCount == 1 &&  isSubstring( "psyc", pulseSequence)){
            psSpectrum = spectra_lst[i];
            break;
        }
    }
    return psSpectrum;
}   


simpleUtils.prototype.hasPeaks = function(spectrum) {
    // check if the spectrum has peaks
    // parameters
    //    spectrum:  spectrum object
    // returns
    //    boolean: true if the spectrum has peaks, false otherwise
    if(spectrum.peaks().count > 0){
        return true;
    }
    else{
        return false;
    }
}

simpleUtils.prototype.isCOSY = function(spectrum) {
    // check if the spectrum is a COSY spectrum
    // parameters
    //    spectrum:  spectrum object
    // returns
    //    boolean: true if the spectrum is a COSY spectrum, false otherwise
    var pulseSequence = spectrum.getParam("Pulse Sequence");
    // change string to lowercase
    pulseSequence = pulseSequence.toLowerCase();
    if(spectrum.dimCount == 2 && isSubstring("cosy", pulseSequence)){
        return true;
    }
    else{
        return false;
    }
}

simpleUtils.prototype.isHSQCCLIPCOSY = function(spectrum) {
    // check if the spectrum is a HSQC-CLIP-COSY spectrum
    // parameters
    //    spectrum:  spectrum object
    // returns
    //    boolean: true if the spectrum is a HSQC-CLIP-COSY spectrum, false otherwise
    var pulseSequence = spectrum.getParam("Pulse Sequence");
    // change string to lowercase
    pulseSequence = pulseSequence.toLowerCase();
    if(spectrum.dimCount == 2 && isSubstring("clip", pulseSequence)){
        return true;
    }
    else{
        return false;
    }
}

// check if the spectrum is a 2D spectrum
simpleUtils.prototype.is2DSpectrum = function(spectrum) {
    // check if the spectrum is a 2D spectrum
    // parameters
    //    spectrum:  spectrum object
    // returns
    //    boolean: true if the spectrum is a 2D spectrum, false otherwise
    
    print("spectrum type", spectrum.type);
    if(spectrum.dimCount == 2){
        return true;
    }
    else{
        MessageBox.warning("Please select a 2D spectrum\n Spectrum is: " + spectrum.type);
        return false;
    }
}

// is HSQC spectrum
simpleUtils.prototype.isHSQC = function(spectrum) {
    // check if the spectrum is a HSQC spectrum
    // parameters
    //    spectrum:  spectrum object
    // returns
    //    boolean: true if the spectrum is a HSQC spectrum, false otherwise
    var pulseSequence = spectrum.getParam("Pulse Sequence");
    // change string to lowercase
    pulseSequence = pulseSequence.toLowerCase();
    print("pulse sequence", pulseSequence);
    if(spectrum.dimCount == 2 && isSubstring("hsqc", pulseSequence) && !isSubstring("clip", pulseSequence)){
        return true;
    }
    else{
        return false;
    }
}

// find HSQC spectrum

simpleUtils.prototype.findHSQC = function(spectra_lst) {
    var hsqcSpectrum = undefined;
    for (var i=0; i<spectra_lst.length; i++){
        if (this.isHSQC(spectra_lst[i])){
            hsqcSpectrum = spectra_lst[i];
            break;
        }
    }
    return hsqcSpectrum;
}


// spectitle = spectitle.replace(/\s/g, "_");
// spectra[spectitle] = {};
// spectra[spectitle]["datatype"] = "nmrspectrum";
// spectra[spectitle]["origin"] = spec.originalFormat;
// spectra[spectitle]["type"] = spec.type;
// print("spec type", spec.type);
// spectra[spectitle]["subtype"] = spec.subtype;
// spectra[spectitle]["experimenttype"] = spec.experimentType;
// spectra[spectitle]["experiment"] = spec.getParam("Experiment");
// spectra[spectitle]["pulsesequence"] = spec.getParam("Pulse Sequence");
// spectra[spectitle]["intrument"] = spec.getParam("Instrument");
// spectra[spectitle]["probe"] = spec.getParam("Probe");
// spectra[spectitle]["datafilename"] = spec.getParam("Data File Name");
// // spectra[spectitle]["comment"] = spec.getParam("Comment");
// spectra[spectitle]["nucleus"] = spec.getParam("Nucleus");
// spectra[spectitle]["specfrequency"] = spec.getParam("Spectrometer Frequency");
// test identifySpectra
// var doc = Application.mainWindow.activeDocument;
// var simpleutils = new simpleUtils();
// var spectra = simpleutils.identifySpectra(doc);



