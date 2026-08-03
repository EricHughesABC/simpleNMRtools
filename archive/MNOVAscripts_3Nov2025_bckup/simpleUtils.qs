// Function to check if one string is a substring of another
function isSubstring(sub, str) {
    return str.indexOf(sub) !== -1;
}

// spectra_keys = ["HSQC", 
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

spectra_keys = ["HSQC", 
                    "HMBC", 
                    "COSY", 
                    "NOESY", 
                    "H1_1D", 
                    "C13_1D", 
                    "DEPT135", 
                    "PureShift", 
                    "HSQC_CLIPCOSY", 
                    "DDEPT_CH3_ONLY",
                    "SKIP", 
                ];

function simpleUtils() {

    var doc = Application.mainWindow.activeDocument;

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

simpleUtils.prototype.findC13spectrum = function( spectra ) {
    var c13spectrum = undefined;
    for (var i=0; i<spectra.length; i++){
        if (spectra[i].subtype == "predicted, 13C"){
            continue
        }
        if (spectra[i].dimCount == 1 && spectra[i].getParam("Nucleus") == "13C") {
            c13spectrum = spectra[i];
            print("C13spectrum index = ", i)
            print("spectrum dimCount & subtype", spectra[i].dimCount, spectra[i].subtype, spectra[i].getParam("Nucleus"));
            break;
        }
    }
    return c13spectrum;
}

simpleUtils.prototype.findPureShiftSpectrum = function( spectra ) {
    var psSpectrum = undefined;
    for (var i=0; i<spectra.length; i++){
        var pulseSequence = spectra[i].getParam("Pulse Sequence");
        // change string to lowercase
        pulseSequence = pulseSequence.toLowerCase();

        if (spectra[i].dimCount == 1 &&  isSubstring( "psyc", pulseSequence)){
            psSpectrum = spectra[i];
            break;
        }
    }
    return psSpectrum;
}   


simpleUtils.prototype.hasPeaks = function(spectrum) {
    if(spectrum.peaks().count > 0){
        return true;
    }
    else{
        return false;
    }
}

simpleUtils.prototype.isCOSY = function(spectrum) {
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
simpleUtils.prototype.findHSQC = function(spectra) {
    var hsqcSpectrum = undefined;
    for (var i=0; i<spectra.length; i++){
        if (this.isHSQC(spectra[i])){
            hsqcSpectrum = spectra[i];
            break;
        }
    }
    return hsqcSpectrum;
}

simpleUtils.prototype.exptUniqueIdentifier = function(spectrum) {
    var spec = spectrum;
    const type = spec.type;
    const subtype = spec.subtype.split(",")[0];
    const pulsesequence = spec.getParam("Pulse Sequence");
    const nucleus = spec.getParam("Nucleus");
    const experiment = spec.getParam("Experiment");
    const spectrum = spec.getParam("Title");
    return nucleus + " " + experiment + " " + pulsesequence + " " + spectrum;
}

// var doc = Application.mainWindow.activeDocument;
// var simpleutils = new simpleUtils();
// var spectra_lst = simpleutils.identifySpectra(doc);

// print("Number of spectra identified: " + spectra_lst.length + "\n");

// for( var i=0; i < spectra_lst.length; i++){
//     var spec = spectra_lst[i];
//     print("spectrum: " + i);
//     print("\tPulse Sequence: " + spec.getParam("Pulse Sequence"));
//     print("\ttype: " + spec.type);
//     print("\tsubtype: " + spec.subtype);
//     print("\texperimentType: " + spec.experimentType);
//     print("\tdimCount: " + spec.dimCount);
//     print("\tnotes: " + spectra_lst[i]['page']['notes']);
//     print("\tnotes: " + spec.getParam("Notes"));
//     print("\tComment: " + spec.getParam("Comment"));
//     print("\tExperiment: " + spec.getParam("Experiment"));
//     print("\tClass: " + spec.getParam("Class"));
//     print("\tNucleus: " + spec.getParam("Nucleus"));
//     print("\tSpectrometer Frequency: " + spec.getParam("Spectrometer Frequency"));
//     print("");
// }

// for( var i=0; i<spectra_lst.length; i++){
//     var spectrum = spectra_lst[i];
//     if(simpleutils.hasPeaks(spectrum)){
//         print("spectrum with peaks found: " + spectrum.getParam("Experiment"));
//     }
//     else{
//         print("spectrum with NO peaks: " + spectrum.getParam("Experiment"));
//     }
// }
// // naming all spectra
// for ( var i=0; i < spectra_lst.length; i++ ) {

//     var unique_ID = simpleutils.exptUniqueIdentifier( spectra_lst[i] );
//     print("Unique ID: " + unique_ID );
// }
    
// // print out exxperiment_eeh property
// for ( var i=0; i < spectra_lst.length; i++ ) {

//     var expt_eeh = spectra_lst[i]["experimentEEH"];
//     print("spectrum experimentEEH: " + expt_eeh );
// }





