// Read in simpleNMRglobals.json file
// and parse the JSON content to get experiments_dict

// read in the json file and return the json object
function readJsonFile(jsonfilename){

    var jsonobj = {};
    var fin = new File(jsonfilename);
    if (fin.open(File.ReadOnly)) {
        var sin = new TextStream(fin, 'UTF-8');
        var jsonstr = sin.readAll();
        jsonobj = JSON.parse(jsonstr);
    }
    else{
        MessageBox.warning("Could not open json file\n" + jsonfilename);
    }
    fin.close();
    return jsonobj;
}

// var jsonFilePath = "simpleNMRglobals.json";
// var globals = readJsonFile(jsonFilePath);
// var experiments_dict = globals.experiments_dict;
// var spectra_keys = Object.keys(experiments_dict);

// Function to check if one string is a substring of another
function isSubstring(sub, str) {
    return str.indexOf(sub) !== -1;
}

const spectra_keys = ["unidentified",
                      "Predicted",    
                      "SKIP", 
                      "HSQC", 
                      "HMBC", 
                      "COSY", 
                      "H1_1D", 
                      "C13_1D", 
                      "DEPT135", 
                      "H1_pureshift", 
                      "DDEPT_CH3_ONLY", 
                      "HSQC_CLIPCOSY",
                      "NOESY", 
                ];




// create a one to one spectra_dict

const experiments_dict = {
                            "unidentified": "unidentified",
                            "Predicted": "Predicted",
                            "SKIP": "SKIP",
                            "HSQC": "HSQC",
                            "HMBC": "HMBC",
                            "COSY": "COSY",
                            "H1_1D": "H1_1D",
                            "C13_1D": "C13_1D",
                            "DEPT135": "DEPT135",
                            "H1_pureshift": "H1_pureshift",
                            "DDEPT_CH3_ONLY": "DDEPT_CH3_ONLY",
                            "HSQC_CLIPCOSY": "HSQC_CLIPCOSY",
                            "NOESY": "NOESY",
                };

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
        MessageBox.warning("Error opening file for writing\n" +  filename);
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
        MessageBox.warning("Error opening file for reading\n" + filename);
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
    const title = spec.getParam("Title");
    return nucleus + " " + experiment + " " + pulsesequence + " " + title;
}

///////////////////////////////////////////////////////////////////////////////
// Shared code for simplePREDICT_eeh.qs, simpleASSIGN_eeh.qs and collectPeaks.qs
//
// Use:   var su = new simpleUtils();
//        var spectra = su.collectSpectra(doc, {details: true});
///////////////////////////////////////////////////////////////////////////////

// spectrum labels accepted by the export scripts (simplePREDICT_eeh, collectPeaks)
simpleUtils.prototype.exportLabels = ["HSQC",
                                      "HMBC",
                                      "COSY",
                                      "NOESY",
                                      "H1_1D",
                                      "C13_1D",
                                      "DEPT135",
                                      "PureShift",
                                      "DDEPTCH3ONLY",
                                      "SKIP",
                                      "HSQC_CLIPCOSY"];

// labels whose spectra are always 2D
simpleUtils.prototype.twoDLabels = ["HSQC", "HMBC", "COSY", "NOESY", "HSQC_CLIPCOSY",
                                    "DDEPT_CH3_ONLY", "DDEPTCH3ONLY"];

// keys in the exported JSON that are not spectra
simpleUtils.prototype.nonSpectrumKeys = ["smiles",
                                         "molfile",
                                         "nmrAssignments",
                                         "c13predictions",
                                         "carbonAtomsInfo",
                                         "allAtomsInfo",
                                         "hostname"];

// { "datatype": datatype, "count": 1, "data": { "0": value } }
simpleUtils.prototype.makeBlock = function (datatype, value) {
    var block = {};
    block["datatype"] = datatype;
    block["count"] = 1;
    block["data"] = {};
    block["data"]["0"] = value;
    return block;
}

// { "datatype": datatype, "count": n, "data": { "0": list[0], "1": list[1], ... } }
simpleUtils.prototype.makeListBlock = function (datatype, list) {
    var block = {};
    block["datatype"] = datatype;
    block["count"] = list.length;
    block["data"] = {};
    for (var i = 0; i < list.length; i++) {
        block["data"][i] = list[i];
    }
    return block;
}

// read a JSON file if it exists, otherwise return defaultValue (no warning)
simpleUtils.prototype.readJsonIfExists = function (filename, defaultValue) {
    if (!new Dir(filename).fileExists(filename)) {
        return defaultValue;
    }
    var fin = new File(filename);
    var rtn = defaultValue;
    if (fin.open(File.ReadOnly)) {
        var sin = new TextStream(fin, 'UTF-8');
        rtn = JSON.parse(sin.readAll());
    }
    fin.close();
    return rtn;
}

// write obj as JSON (indent = undefined for one line, 4 for the export files)
simpleUtils.prototype.writeJson = function (filename, obj, indent) {
    var fout = new File(filename);
    if (fout.open(File.WriteOnly)) {
        var sout = new TextStream(fout, 'UTF-8');
        sout.writeln(indent === undefined ? JSON.stringify(obj) : JSON.stringify(obj, null, indent));
    }
    else {
        print("Error opening file for writing\n", filename);
    }
    fout.close();
}

// smiles and molfile blocks (count 0 and empty data when there is no molecule)
simpleUtils.prototype.addMoleculeBlocks = function (spectra, mol) {
    if (mol === undefined) {
        spectra["smiles"] = {"datatype": "smiles", "count": 0, "data": {}};
        spectra["molfile"] = {"datatype": "molfile", "count": 0, "data": {}};
        return;
    }
    spectra["smiles"] = this.makeBlock("smiles", mol.generateSMILES());
    spectra["molfile"] = this.makeBlock("molfile", mol.getMolfile());
}

simpleUtils.prototype.hostnameBlock = function () {
    return this.makeBlock("hostname", SysInfo.hostID);
}

// allAtomsInfo and carbonAtomsInfo blocks
simpleUtils.prototype.addAtomInfoBlocks = function (spectra, mol) {
    var all = {"datatype": "allAtomsInfo", "data": {}};
    var carbons = {"datatype": "carbonAtomsInfo", "data": {}};
    var allCount = 0;
    var carbonCount = 0;
    for (var i = 1; i <= mol.atomCount; i++) {
        var atom = new Atom(mol.atom(i));
        var j = i - 1;
        var info = {};
        info["atom_idx"] = j;
        info["id"] = j;
        info["atomNumber"] = atom.number;
        info["symbol"] = atom.elementSymbol;
        info["numProtons"] = atom.nH;
        all["data"][j] = info;
        allCount++;
        if (atom.elementSymbol == "C") {
            carbons["data"][j] = JSON.parse(JSON.stringify(info));
            carbonCount++;
        }
    }
    all["count"] = allCount;
    carbons["count"] = carbonCount;
    spectra["allAtomsInfo"] = all;
    spectra["carbonAtomsInfo"] = carbons;
}

// carbon and attached-proton shifts assigned in Mnova (nmrAssignments block)
simpleUtils.prototype.nmrAssignmentsBlock = function (mol) {
    var block = {"datatype": "nmrAssignments", "count": 0, "data": {}};
    var assignments = NMRAssignments(mol);
    if (assignments.hasAssignments() != true) {
        return block;
    }
    var count = 0;
    for (var i = 1; i <= mol.atomCount; i++) {
        var atom = new Atom(mol.atom(i));
        if (atom.elementSymbol != "C") {
            continue;
        }
        block["data"][i] = {};
        var cshift = assignments.chemShift(i, 0);
        if (cshift == undefined) {
            continue;
        }
        var a = block["data"][i];
        a["atom_idx"] = i;
        a["atomNumber"] = atom.number;
        a["numProtons"] = atom.nH;
        a["f1_ppm"] = cshift.shift;
        a["label"] = atom.Label;
        a["iupacLabel"] = "";
        count++;
        for (var k = 1; k <= 3; k++) {
            var hshift = assignments.chemShift(i, k);
            a["f2_ppm" + k] = (hshift === undefined) ? "undefined" : hshift.shift;
        }
    }
    block["count"] = count;
    return block;
}

// One block per NMR spectrum in the document, keyed "<title>_<page index>".
// options.details = true  (collectPeaks): every peak, with kind, widths and flags,
//                          and the spectral width
// options.details = false (simplePREDICT): compound peaks without flags only
simpleUtils.prototype.collectSpectra = function (doc, options) {
    var details = options && options.details;
    var spectra = {};
    for (var i = 0, pageCount = doc.pageCount(); i < pageCount; i++) {
        var page = doc.page(i);
        for (var j = 0, itemCount = page.itemCount(); j < itemCount; j++) {
            var spec = new NMRSpectrum(page.item(j));
            if (!spec.isValid()) {
                continue;
            }
            var spectitle = (spec.title + "_" + i).replace(/\s/g, "_");
            spectra[spectitle] = this.spectrumBlock(spec, details);
        }
    }
    return spectra;
}

simpleUtils.prototype.spectrumBlock = function (spec, details) {
    var s = {};
    s["datatype"] = "nmrspectrum";
    s["origin"] = spec.originalFormat;
    s["type"] = spec.type;
    s["subtype"] = spec.subtype;
    s["experimenttype"] = spec.experimentType;
    s["experiment"] = spec.getParam("Experiment");
    s["class"] = spec.getParam("Class");
    s["spectype"] = spec.getParam("Spectrum Type");
    s["pulsesequence"] = spec.getParam("Pulse Sequence");
    s["intrument"] = spec.getParam("Instrument");
    s["probe"] = spec.getParam("Probe");
    s["datafilename"] = spec.getParam("Data File Name");
    s["nucleus"] = spec.getParam("Nucleus");
    s["specfrequency"] = spec.getParam("Spectrometer Frequency");
    if (details) {
        s["sw_Hz"] = spec.getParam("Spectral Width");
    }
    s["acquired_size"] = spec.getParam("Acquired Size");
    s["spectral_size"] = spec.getParam("Spectral Size");
    s["solvent"] = spec.getParam("Solvent");
    s["temperature"] = spec.getParam("Temperature");

    // a "1D" spectrum with two nuclei is really 2D
    if (s["type"] == "1D" && s["nucleus"].indexOf(",") != -1) {
        s["type"] = "2D";
    }

    // multiplets (compound multiplets only)
    var multiplets = spec.multiplets();
    s["multiplets"] = {"datatype": "multiplets", "normValue": multiplets.normValue, "data": {}};
    var goodMultiplets = 0;
    for (var m = 0; m < multiplets.count; m++) {
        var multiplet = multiplets.at(m);
        if (multiplet.type != 0) {
            continue;
        }
        var jlist = multiplet.jList();
        var jvals = [];
        for (var q = 0; q < jlist.count; q++) {
            jvals.push(jlist.at(q));
        }
        s["multiplets"]["data"][m] = {
            "delta1": multiplet.delta,
            "nH": multiplet.nH,
            "realH": multiplet.realH,
            "integralValue": multiplet.integralValue(),
            "category": multiplet.category,
            "type": multiplet.type,
            "jlistcount": jlist.count,
            "jvals": jvals
        };
        goodMultiplets++;
    }
    s["multiplets"]["count"] = goodMultiplets;

    // peaks
    var peaks = spec.peaks();
    s["peaks"] = {"datatype": "peaks", "data": {}};
    var goodPeaks = 0;
    for (var p = 0; p < peaks.count; p++) {
        var pk = peaks.at(p);
        if (!details && !(pk.type == 0 && pk.flagsToString() == "None")) {
            continue;
        }
        var d = {};
        d["delta1"] = pk.delta(1);
        d["delta2"] = pk.delta(2);
        d["intensity"] = pk.intensity;
        d["type"] = pk.type;
        if (details) {
            d["kind"] = pk.kind;
            d["width1"] = pk.width(1);
            d["width2"] = pk.width(2);
            d["flags"] = pk.flags;
        }
        d["annotation"] = pk.annotation;
        s["peaks"]["data"][p] = d;
        goodPeaks++;
    }
    s["peaks"]["count"] = goodPeaks;

    // integrals
    var integrals = spec.integrals();
    s["integrals"] = {"datatype": "integrals", "count": integrals.count,
                      "normValue": integrals.normValue, "type": integrals.type, "data": {}};
    for (var n = 0; n < integrals.count; n++) {
        var integral = integrals.at(n);
        s["integrals"]["data"][n] = {
            "intensity": integral.integralValue(),
            "rangeMin1": integral.rangeMin(1),
            "rangeMin2": integral.rangeMin(2),
            "rangeMax1": integral.rangeMax(1),
            "rangeMax2": integral.rangeMax(2),
            "delta1": (integral.rangeMin(1) + integral.rangeMax(1)) / 2,
            "delta2": (integral.rangeMin(2) + integral.rangeMax(2)) / 2,
            "type": 0
        };
    }
    return s;
}

// "<nucleus> <experiment> <pulse sequence> <key>" for every spectrum with peaks
simpleUtils.prototype.spectraWithPeaks = function (spectra) {
    var list = [];
    for (var key in spectra) {
        if (this.nonSpectrumKeys.indexOf(key) != -1) {
            continue;
        }
        if (spectra[key]["peaks"]["count"] == 0) {
            continue;
        }
        list.push(spectra[key]["nucleus"] + " " + spectra[key]["experiment"] + " " +
                  spectra[key]["pulsesequence"] + " " + key);
    }
    return list;
}

// Rename the chosen spectra to "<label>_<n>" (HSQC_0, HMBC_0, ...), store the
// original key in "filename", mark 2D spectra and add the chosenSpectra block.
// Returns {"<label>_<n>": "<original key>"}, or undefined after a warning if a
// label is not in this.exportLabels.
simpleUtils.prototype.applyChosenSpectra = function (spectra, chosen_spectra) {
    spectra["chosenSpectra"] = this.makeListBlock("chosenSpectra", chosen_spectra["checked_spectra"]);

    var ids = {};
    var counts = {};
    for (var i = 0; i < this.exportLabels.length; i++) {
        counts[this.exportLabels[i]] = -1;
    }
    var checked = chosen_spectra["checked_spectra"];
    for (var c = 0; c < checked.length; c++) {
        var parts = checked[c].split(" ");
        var label = parts[parts.length - 1];
        var original = parts[parts.length - 2];
        if (counts[label] === undefined) {
            MessageBox.warning("Error: Spectrum type " + label + " unknown. Please check the spectra in the document.");
            return undefined;
        }
        counts[label] += 1;
        ids[label + "_" + counts[label]] = original;
    }

    for (var key in ids) {
        spectra[key] = spectra[ids[key]];
        spectra[key]["filename"] = ids[key];
    }
    for (var key2 in ids) {
        if (this.twoDLabels.indexOf(key2.split("_")[0]) != -1) {
            spectra[key2]["type"] = "2D";
        }
    }
    for (var key3 in ids) {
        delete spectra[ids[key3]];
    }
    return ids;
}

// true if one of the chosen spectra has this label (e.g. "HSQC")
simpleUtils.prototype.hasChosenLabel = function (ids, label) {
    for (var key in ids) {
        if (key.split("_")[0] == label) {
            return true;
        }
    }
    return false;
}

// For every HSQC_<n>: if it has as many integrals as peaks, replace each peak's
// intensity with the integral whose region contains it.  Returns the last HSQC block.
simpleUtils.prototype.hsqcIntensitiesFromIntegrals = function (spectra) {
    var hsqc = undefined;
    for (var key in spectra) {
        var split_key = key.split("_");
        if (split_key[0] != "HSQC" || split_key.length != 2) {
            continue;
        }
        hsqc = spectra[key];
        var pks = hsqc["peaks"]["data"];
        var pks_count = hsqc["peaks"]["count"];
        var integrals = hsqc["integrals"]["data"];
        var integrals_count = hsqc["integrals"]["count"];

        if (pks_count == integrals_count) {
            for (var i = 0; i < pks_count; i++) {
                var pk = pks[i];
                for (var j = 0; j < integrals_count; j++) {
                    var integral = integrals[j];
                    if ((pk["delta1"] >= integral["rangeMax1"]) &&
                        (pk["delta1"] <= integral["rangeMin1"]) &&
                        (pk["delta2"] >= integral["rangeMin2"]) &&
                        (pk["delta2"] <= integral["rangeMax2"])) {
                        pk["intensity"] = integral["intensity"];
                        break;
                    }
                }
            }
        }
        else {
            MessageBox.warning("Number of peaks and integrals do not match in HSQC data. Therefore using peaks only. Program works best with integrals");
        }
    }
    return hsqc;
}

// Show the messages and open the registration page if the server says this
// computer is not (or no longer) registered.  Returns true if the script must stop.
simpleUtils.prototype.registrationProblem = function (responseObj) {
    if (responseObj.status !== "unregistered" && responseObj.status !== "registration_expired") {
        return false;
    }
    if (responseObj.status === "registration_expired") {
        MessageBox.warning("Host registration has expired. Please reregister to use this feature.");
    }
    else {
        MessageBox.warning("Host is not registered. Please register to use simpleNMR prediction.");
    }
    if (typeof Application.openUrl === 'function') {
        Application.openUrl(responseObj.registration_url);
    } else {
        print("Application.openUrl is not a function");
    }
    return true;
}

// Ask the server whether this computer is registered.
// Returns {"ml_consent": ...} if it is, undefined if the script must stop.
simpleUtils.prototype.checkRegistration = function (server, web_utils) {
    var hostname = SysInfo.hostID;
    print("hostname ", hostname);
    print("server ", server);

    var entry_point = server + "check_machine_learning";
    var json_obj_string = JSON.stringify({"hostname": hostname}, null, 4);
    print("Entry point", entry_point);
    var rtn = web_utils.jsonRequest(entry_point, json_obj_string, "", "", false);

    print("rtn.exitCode", rtn.exitCode);
    print("rtn.allErrorOutput", rtn.allErrorOutput);
    print("\nrtn.response\n");
    print(rtn.response);

    var responseObj;
    try {
        responseObj = JSON.parse(rtn.response);
    } catch (e) {
        print("Error parsing JSON response: ", e);
        MessageBox.warning("Error parsing JSON response: " + e + "\n" + rtn.response +
                           "\nrtn.exitCode" + rtn.exitCode + "\nrtn.allErrorOutput" + rtn.allErrorOutput);
        return undefined;
    }
    if (this.registrationProblem(responseObj)) {
        return undefined;
    }
    print("ml_consent ", responseObj.ml_consent);
    return {"ml_consent": responseObj.ml_consent};
}

// Local file path -> file:// URL that tolerates spaces etc.
simpleUtils.prototype.pathToFileUrl = function (path) {
    var p = String(path).replace(/\\/g, "/");   // backslashes -> forward slashes
    var parts = p.split("/");
    for (var i = 0; i < parts.length; i++) {
        // leave a Windows drive letter such as "C:" unencoded
        if (i === 0 && /^[A-Za-z]:$/.test(parts[i])) {
            continue;
        }
        parts[i] = encodeURIComponent(parts[i]);  // space -> %20, # -> %23, etc.
    }
    p = parts.join("/");
    if (p.charAt(0) === "/") {
        return "file://" + p;      // Mac/Linux: /Users/... -> file:///Users/...
    }
    return "file:///" + p;         // Windows: C:/... -> file:///C:/...
}

// Write the simpleNMR HTML page to <directory>/html/<name>_d3.html; returns the file name.
simpleUtils.prototype.writeHtml = function (htmlString, directoryPath, name) {
    var htmlfilename = directoryPath + "/html/" + name + "_d3.html";
    print("htmlfilename ", htmlfilename);

    var dir = new Dir(directoryPath + "/html");
    if (!dir.exists) {
        dir.mkdir(directoryPath + "/html");
    }

    var fout = new File(htmlfilename);
    if (fout.open(File.WriteOnly)) {
        var sout = new TextStream(fout, 'UTF-8');
        sout.write(htmlString);
    }
    fout.close();
    return htmlfilename;
}

// Open a local file in the default web browser.
simpleUtils.prototype.openInBrowser = function (filename) {
    if (typeof Application.openUrl === 'function') {
        var url = this.pathToFileUrl(filename);
        print("Opening: " + url);
        Application.openUrl(url);
    }
    else {
        print("Application.openUrl is not a function");
    }
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





