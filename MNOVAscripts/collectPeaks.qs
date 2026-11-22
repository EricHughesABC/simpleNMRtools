/*globals Application, simpleUtils, peakRefineDialog, MessageBox, print*/

// Export the picked peaks of the current Mnova document as JSON for automatic peak
// refinement (simpleNMRclaudePeakRefine).  The spectra are chosen in peakRefineDialog.qs;
// the choices are remembered in peakRefine_exptIdentifiers.json in the document's folder,
// separate from the simpleNMR verify step's exptIdentifiers.json.
// Output: <document>_claudePeakPicked_mresnova.json next to the Mnova document.
// No server is contacted.  Shared code is in simpleUtils.qs.

function collectPeaks() {

    var su = new simpleUtils();
    var doc = Application.mainWindow.activeDocument;
    var result = su.splitPathAndFilename(doc.name);

    // molecule (optional: only its formula is used by the refinement) and hostname
    var spectra = {};
    su.addMoleculeBlocks(spectra, su.getActiveMolecule(doc, Application.molecule));
    spectra["hostname"] = su.hostnameBlock();

    // every spectrum, with all peaks and their details
    var collected = su.collectSpectra(doc, {details: true});
    for (var key in collected) {
        spectra[key] = collected[key];
    }
    var spectra_with_peaks = su.spectraWithPeaks(spectra);

    // the labels chosen last time for this folder (falls back once to the verify step's choices)
    var exptIdentifiersJsonFilename = result.directoryPath + "/peakRefine_exptIdentifiers.json";
    var previousChoices = su.readJson(exptIdentifiersJsonFilename, null);
    if (previousChoices === null) {
        previousChoices = su.readJson(result.directoryPath + "/exptIdentifiers.json", null);
    }

    var chosen_spectra = peakRefineDialog(spectra_with_peaks, previousChoices, result["name"]);
    if (chosen_spectra === undefined) {
        MessageBox.warning("No spectra chosen");
        return;
    }
    print("chosen_spectra\n", chosen_spectra);

    // HSQC_0, HMBC_0, ... and the chosenSpectra block
    if (su.applyChosenSpectra(spectra, chosen_spectra) === undefined) {
        return;
    }
    su.hsqcIntensitiesFromIntegrals(spectra);

    // not used by the peak refinement; written with neutral values so the file has
    // every field of the shared simpleNMR JSON format (simpleNMRbuilder contract)
    spectra["MNOVAcalcMethod"] = su.makeBlock("MNOVAcalcMethod", "");
    spectra["exptIdentifiers"] = su.makeListBlock("exptIdentifiers", chosen_spectra["exptIdentifiers"]);

    // remember the choices for next time (by spectrum name) - never written to the
    // verify step's exptIdentifiers.json
    su.writeJson(exptIdentifiersJsonFilename, chosen_spectra["choices"], 2);

    spectra["workingDirectory"] = su.makeBlock("workingDirectory", result["directoryPath"]);
    spectra["workingFilename"] = su.makeBlock("workingFilename", result["name"]);
    spectra["simulatedAnnealing"] = su.makeBlock("simulatedAnnealing", false);   // placeholder
    spectra["ml_consent"] = su.makeBlock("ml_consent", false);                   // placeholder
    spectra["spectraWithPeaks"] = su.makeListBlock("spectraWithPeaks", spectra_with_peaks);

    var jsonfilename = result.directoryPath + "/" + result.name + "_claudePeakPicked_mresnova.json";
    su.writeJson(jsonfilename, spectra, 4);

    print("Peaks exported to " + jsonfilename);
    MessageBox.information("Peaks of " + chosen_spectra["checked_spectra"].length +
                           " spectra exported for refinement to\n" + jsonfilename);
}
