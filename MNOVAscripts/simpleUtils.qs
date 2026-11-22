
// Function to check if one string is a substring of another
function isSubstring(sub, str) {
    return str.indexOf(sub) !== -1;
}

///////////////////////////////////////////////////////////////////////////////
// THE list of NMR experiment types - defined here once, used by every script.
//
// The names are the same as simpleNMR's NMREXPERIMENTS (simpleNMRtools
// config/globals.py) and simpleNMRbuilder, so the JSON keys match the server.
// Spectra in the JSON are keyed <name>_<n>: build them with exptKey(EXPT.HSQC, 0)
// ("HSQC_0") and read them back with exptFromKey("HSQC_0") ("HSQC").  Refer to an
// experiment as EXPT.HSQC, EXPT.H1_pureshift, ... - never type the name in a script.
//
// Each entry: dims (1 or 2), the Mnova subtype a 1D spectrum must have (null = any),
// and the pulse-sequence names (substrings) used to identify the experiment.  The
// order is the order in which the identification tries them (first match wins).
///////////////////////////////////////////////////////////////////////////////
var SIMPLENMR_EXPERIMENTS = [
    {"name": "H1_1D", "dims": 1, "subtype": "1H",
     "pulseSequences": ["zg30", "s2pul", "zg", "zgcppr"]},
    {"name": "C13_1D", "dims": 1, "subtype": "13C",
     "pulseSequences": ["zgdc30", "s2pul", "zgpg30", "zgzrse", "zg0dc.fr", "js_carbon.jxp"]},
    {"name": "H1_pureshift", "dims": 1, "subtype": null,
     "pulseSequences": ["ja_PSYCHE_pr_03b", "reset_psychetse.ptg"]},
    {"name": "DEPT135", "dims": 1, "subtype": "13C",
     "pulseSequences": ["dept135.wu", "DEPT", "deptsp135"]},
    {"name": "HSQC", "dims": 2, "subtype": null,
     "pulseSequences": ["hsqcedetgpsisp2.3.ptg", "hsqcedetgpsisp2.3", "gHSQCAD", "hsqcedetgpsp.3", "gHSQC", "inv4gp.wu", "hsqcetgp", "gns_noah3-BSScc.eeh", "hsqcedetgpsisp2.4", "hsqcedetgpsisp.2", "js_hsqc.jxp"]},
    {"name": "HMBC", "dims": 2, "subtype": null,
     "pulseSequences": ["ghmbc.wu", "gHMBC", "hmbcetgpl3nd", "hmbcetgpl3nd.ptg", "gHMBCAD", "hmbcgpndqf", "gns_noah3-BSScc.eeh", "shmbcctetgpl2nd", "js_hmbcet.jxp", "hmbcgplpndqf"]},
    {"name": "COSY", "dims": 2, "subtype": null,
     "pulseSequences": ["cosygpqf", "cosygp", "gcosy", "cosygpmfppqf", "cosygpmfqf", "gCOSY", "js_cosy.jxp", "cosygpppqf_ptype", "cosyqf45", "cosygpmfphpp"]},
    {"name": "DDEPTCH3ONLY", "dims": 2, "subtype": null,
     "pulseSequences": ["hcdeptedetgpzf"]},
    {"name": "HSQC_CLIPCOSY", "dims": 2, "subtype": null,
     "pulseSequences": ["hsqc_clip_cosy_mc_notation.eeh", "gns_noah3-BSScc.eeh"]},
    {"name": "NOESY", "dims": 2, "subtype": null,
     "pulseSequences": ["noesygpphppzs"]}
];

// EXPT.H1_1D == "H1_1D", ... plus the special labels used by the dialogs and the
// spectrum identification (SKIP = do not use this spectrum)
var EXPT = {"SKIP": "SKIP", "UNIDENTIFIED": "unidentified", "PREDICTED": "Predicted"};
for (var iExpt = 0; iExpt < SIMPLENMR_EXPERIMENTS.length; iExpt++) {
    EXPT[SIMPLENMR_EXPERIMENTS[iExpt].name] = SIMPLENMR_EXPERIMENTS[iExpt].name;
}

// names of all experiment types, in table order
function exptNames() {
    var names = [];
    for (var i = 0; i < SIMPLENMR_EXPERIMENTS.length; i++) {
        names.push(SIMPLENMR_EXPERIMENTS[i].name);
    }
    return names;
}

// the table entry for one experiment name (undefined if the name is not known)
function exptInfo(name) {
    for (var i = 0; i < SIMPLENMR_EXPERIMENTS.length; i++) {
        if (SIMPLENMR_EXPERIMENTS[i].name === name) {
            return SIMPLENMR_EXPERIMENTS[i];
        }
    }
    return undefined;
}

// the choices offered in the spectrum dialogs: SKIP followed by every experiment
function exptDialogList() {
    return [EXPT.SKIP].concat(exptNames());
}

// "HSQC", 0 -> "HSQC_0".  Throws on an unknown name, so a typo fails loudly.
function exptKey(name, n) {
    if (exptInfo(name) === undefined && name !== EXPT.SKIP) {
        throw "exptKey: unknown experiment type '" + name + "'";
    }
    return name + "_" + (n === undefined ? 0 : n);
}

// "HSQC_0" -> "HSQC", "HSQC_CLIPCOSY_1" -> "HSQC_CLIPCOSY"; undefined if not a spectrum key
function exptFromKey(key) {
    var m = /^(.*)_([0-9]+)$/.exec(String(key));
    if (m === null) {
        return undefined;
    }
    return (exptInfo(m[1]) !== undefined || m[1] === EXPT.SKIP) ? m[1] : undefined;
}

// true for 2D experiment types
function exptIsTwoD(name) {
    var info = exptInfo(name);
    return info !== undefined && info.dims === 2;
}

// true if the spectrum's pulse sequence contains one of validPulseSequences
function comparePulseSequence(validPulseSequences, pulseSequence){
    for (var i = 0; i < validPulseSequences.length; i++) {
        if (isSubstring(validPulseSequences[i], pulseSequence)) {
            return true;
        }
    }
    return false;
}

// true for a 2D spectrum whose pulse sequence is one of the HSQC pulse sequences
function isHSQC(spectrum) {
    print("spectrum dimCount " + spectrum.dimCount);
    if (spectrum.dimCount != 2) {
        return false;
    }
    return comparePulseSequence(exptInfo(EXPT.HSQC).pulseSequences, spectrum.getParam("Pulse Sequence"));
}

// Lists kept for the older scripts (identify_spectrum.qs, idspectra_dialog.qs);
// both are derived from SIMPLENMR_EXPERIMENTS.
var spectra_keys = [EXPT.UNIDENTIFIED, EXPT.PREDICTED, EXPT.SKIP].concat(exptNames());
var experiments_dict = {};
for (var iKey = 0; iKey < spectra_keys.length; iKey++) {
    experiments_dict[spectra_keys[iKey]] = spectra_keys[iKey];
}

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

// ----------------------------------------------------------------------------
// Reading and writing files - use these in every script: readJson / writeJson for
// JSON, writeText for plain text.
// ----------------------------------------------------------------------------

// Read a JSON file.  Returns defaultValue if the file does not exist, cannot be
// opened or is not valid JSON; the reason is printed, and also shown in a message
// box when warn is true.
simpleUtils.prototype.readJson = function (filename, defaultValue, warn) {
    var problem = null;
    var rtn = defaultValue;
    if (!new Dir(filename).fileExists(filename)) {
        problem = "File not found";
    }
    else {
        var fin = new File(filename);
        if (fin.open(File.ReadOnly)) {
            var text = new TextStream(fin, 'UTF-8').readAll();
            fin.close();
            try {
                rtn = JSON.parse(text);
            } catch (e) {
                problem = "Not a valid JSON file (" + e + ")";
            }
        }
        else {
            problem = "Could not open the file";
        }
    }
    if (problem !== null) {
        print(problem + ": " + filename);
        if (warn) {
            MessageBox.warning(problem + ":\n" + filename);
        }
    }
    return rtn;
}

// Write plain text to a file (UTF-8), exactly as given.  Returns true if written.
simpleUtils.prototype.writeText = function (filename, text) {
    var fout = new File(filename);
    if (!fout.open(File.WriteOnly)) {
        print("Error opening file for writing\n", filename);
        return false;
    }
    var sout = new TextStream(fout, 'UTF-8');
    sout.write(text);
    fout.close();
    return true;
}

// Write obj to a JSON file: indent undefined = one line, a number (2, 4) = indented.
// Returns true if the file was written.
simpleUtils.prototype.writeJson = function (filename, obj, indent) {
    var fout = new File(filename);
    if (!fout.open(File.WriteOnly)) {
        print("Error opening file for writing\n", filename);
        return false;
    }
    var sout = new TextStream(fout, 'UTF-8');
    sout.writeln(indent === undefined ? JSON.stringify(obj) : JSON.stringify(obj, null, indent));
    fout.close();
    return true;
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
// label is not SKIP or one of SIMPLENMR_EXPERIMENTS.
simpleUtils.prototype.applyChosenSpectra = function (spectra, chosen_spectra) {
    spectra["chosenSpectra"] = this.makeListBlock("chosenSpectra", chosen_spectra["checked_spectra"]);

    var ids = {};
    var counts = {};
    var labels = exptDialogList();
    for (var i = 0; i < labels.length; i++) {
        counts[labels[i]] = -1;
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
        ids[exptKey(label, counts[label])] = original;
    }

    for (var key in ids) {
        spectra[key] = spectra[ids[key]];
        spectra[key]["filename"] = ids[key];
    }
    for (var key2 in ids) {
        if (exptIsTwoD(exptFromKey(key2))) {
            spectra[key2]["type"] = "2D";
        }
    }
    for (var key3 in ids) {
        delete spectra[ids[key3]];
    }
    return ids;
}

// true if one of the chosen spectra is of this type (e.g. EXPT.HSQC)
simpleUtils.prototype.hasChosenLabel = function (ids, label) {
    for (var key in ids) {
        if (exptFromKey(key) === label) {
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
        if (exptFromKey(key) !== EXPT.HSQC) {
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





