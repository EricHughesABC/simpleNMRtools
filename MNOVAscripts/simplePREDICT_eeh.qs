/*globals Application, simpleUtils, server_address, WebUtilsQT, xy3_dialog2, predictCarbon_eeh,
          ProgressDialog, MessageBox, File, Dir, TextStream, print, MnUi*/

// simpleNMR prediction: export the molecule, the Mnova assignments and the chosen spectra,
// send them to the simpleNMR server and open the HTML page it returns.
// Shared code is in simpleUtils.qs.

function simplePREDICT_eeh() {

    // Look for a previous export file (written by the "Export" button in
    // the HTML page -- see exportToMnova() in export.js) sitting next to
    // the mnova document, and pull out the "which experimental peak is
    // attached to which fixed atom slot" state (nodes_now) at export time.
    // atomNumber/x/y are fixed structural properties of the molecule
    // (from the molfile and its RDKit depiction, same every run); what
    // varies run to run is which experimental ppm gets attached to a
    // given atomNumber -- decided fresh each time by the server's
    // peak-to-atom matching, with no memory of a previous run unless we
    // feed it back in here. The path is constructed the same way the HTML
    // page names the file (workingDirectory + workingFilename +
    // "_assignments_from_simplemnova.json"), so no extra input from the
    // user is needed to find it.
    //
    // Returns an array of {atomNumber, ppm, x, y} objects (x/y are carried
    // along for completeness but not needed server-side, since they
    // follow deterministically from atomNumber) -- empty if no previous
    // export is found, or if it cannot be read/parsed.
    function readPreviousNodePositions(directoryPath, name) {

        var positions = [];
        var previousExportFilename = directoryPath + "/" + name + "_assignments_from_simplemnova.json";
        var previousExport = new simpleUtils().readJson(previousExportFilename, {});

        var nodesNow = previousExport["nodes_now"];
        if (nodesNow === undefined) {
            print("Previous export file has no nodes_now entry");
            return positions;
        }

        for (var i = 0; i < nodesNow.length; i++) {
            var node = nodesNow[i];
            if (node["atomNumber"] === undefined || node["ppm"] === undefined) {
                continue;
            }
            positions.push({
                "atomNumber": node["atomNumber"],
                "ppm": node["ppm"],
                "x": node["x"],
                "y": node["y"]
            });
        }

        return positions;
    }

    var su = new simpleUtils();
    var doc = Application.mainWindow.activeDocument;

    // server, and is this computer registered?
    var server = server_address();
    var web_utils = new WebUtilsQT();
    var registration = su.checkRegistration(server, web_utils);
    if (registration === undefined) {
        return;
    }
    var ml_consent = registration.ml_consent;

    // molecule (required here)
    var mol = su.getActiveMolecule(doc, Application.molecule);
    if (mol === undefined) {
        MessageBox.warning("No molecule found!");
        return;
    }

    var spectra = {};
    su.addMoleculeBlocks(spectra, mol);
    spectra["hostname"] = su.hostnameBlock();
    su.addAtomInfoBlocks(spectra, mol);
    spectra["nmrAssignments"] = su.nmrAssignmentsBlock(mol);

    // predict the C13 chemical shifts of the molecule from MNOVA
    spectra["c13predictions"] = predictCarbon_eeh();

    var result = su.splitPathAndFilename(doc.name);

    // every spectrum, compound peaks only
    var collected = su.collectSpectra(doc, {details: false});
    for (var key in collected) {
        spectra[key] = collected[key];
    }
    print("\nfinding spectra with peaks\n");
    var spectra_with_peaks = su.spectraWithPeaks(spectra);

    // dialog settings and spectrum labels chosen last time
    var dialogParametersJsonFilename = result.directoryPath + "/dialogPrameters.json";
    var dialogParameters = su.readJson(dialogParametersJsonFilename, {});
    dialogParameters["ml_consent"] = ml_consent;
    print("dialogParameters\n", dialogParameters);

    var exptIdentifiersJsonFilename = result.directoryPath + "/exptIdentifiers.json";
    var expt_identifiers = su.readJson(exptIdentifiersJsonFilename, []);

    print("server", server);
    var chosen_spectra = xy3_dialog2(spectra_with_peaks,
                                     expt_identifiers,
                                     dialogParameters,
                                     result["name"],
                                     server);
    if (chosen_spectra === undefined) {
        MessageBox.warning("No spectra chosen");
        return;
    }
    print("chosen_spectra\n", chosen_spectra);

    // HSQC_0, HMBC_0, ... and the chosenSpectra block
    var chosen_spectra_ids = su.applyChosenSpectra(spectra, chosen_spectra);
    if (chosen_spectra_ids === undefined) {
        return;
    }
    if (!su.hasChosenLabel(chosen_spectra_ids, EXPT.HSQC)) {
        MessageBox.warning("No HSQC found in the chosen spectra. Please choose a HSQC spectrum.");
        return;
    }

    // NOTE: the keys are "H1_1D_0", ... so this lookup never finds the 1H spectrum and the
    // HSQC-peaks-vs-1H-multiplets check further down never runs (unchanged behaviour,
    // under review).
    var h1_1d = spectra[exptKey(EXPT.H1_1D, 0)];

    // HSQC peak intensities from the HSQC integrals (when the counts match)
    var hsqc = su.hsqcIntensitiesFromIntegrals(spectra);

    spectra["MNOVAcalcMethod"] = su.makeBlock("MNOVAcalcMethod", chosen_spectra["calcSimpleMNOVA"]);
    spectra["exptIdentifiers"] = su.makeListBlock("exptIdentifiers", chosen_spectra["exptIdentifiers"]);

    // remember the labels and the dialog settings for next time
    su.writeJson(exptIdentifiersJsonFilename, chosen_spectra["exptIdentifiers"]);
    var dialogParams = {};
    dialogParams["simulatedAnnealing"] = chosen_spectra["simulatedAnnealing"];
    dialogParams["calcSimpleNMR"] = chosen_spectra["calcSimpleNMR"];
    dialogParams["calcSimpleMNOVA"] = chosen_spectra["calcSimpleMNOVA"];
    dialogParams["usePreviousPositions"] = chosen_spectra["usePreviousPositions"];
    su.writeJson(dialogParametersJsonFilename, dialogParams);

    spectra["workingDirectory"] = su.makeBlock("workingDirectory", result["directoryPath"]);
    spectra["workingFilename"] = su.makeBlock("workingFilename", result["name"]);
    spectra["simulatedAnnealing"] = su.makeBlock("simulatedAnnealing", chosen_spectra["simulatedAnnealing"]);
    spectra["ml_consent"] = su.makeBlock("ml_consent", chosen_spectra["ml_consent"]);
    spectra["spectraWithPeaks"] = su.makeListBlock("spectraWithPeaks", spectra_with_peaks);

    // optionally carry forward the previous iteration's peak-to-atom
    // assignment (which experimental ppm was attached to which fixed atom
    // slot) so the server can restore it, ahead of its normal matching --
    // see readPreviousNodePositions() above and _apply_previous_atom_assignments
    // in core/nmrsolution.py
    if (chosen_spectra["usePreviousPositions"]) {
        var previousPositions = readPreviousNodePositions(result.directoryPath, result.name);
        if (previousPositions.length > 0) {
            spectra["previousNodePositions"] = su.makeListBlock("previousNodePositions", previousPositions);
        }
        else {
            MessageBox.warning("No previous node positions found to reuse - continuing without them.");
        }
    }

    // save the JSON sent to the server next to the Mnova document
    var jsonfilename = result.directoryPath + "/" + result.name + "_assignments_mresnova.json";
    var json_spectra_string = JSON.stringify(spectra, null, 4);
    su.writeJson(jsonfilename, spectra, 4);

    ////////////////////
    // Do some checking
    ////////////////////

    // check if there are any assignments if the user has chosen to use MNOVA manually assigned
    if (spectra["nmrAssignments"]["count"] == 0 && chosen_spectra["calcSimpleMNOVA"] == "MNova Manually Assigned") {
        MessageBox.warning("No MNOVA assignment data found");
        return;
    }

    // check if there are any C13 predictions if the user has chosen to use MNOVA Predict
    if (spectra["c13predictions"]["count"] == 0 && chosen_spectra["calcSimpleMNOVA"] == "MNOVA Predict") {
        MessageBox.warning("No C13 predictions found");
        return;
    }

    // check if there are no integrals and no peaks in the HSQC spectrum
    if (hsqc["peaks"]["count"] == 0 && hsqc["integrals"]["count"] == 0) {
        MessageBox.warning("No peaks or integrals found in HSQC spectrum");
        return;
    }

    // warn the user that there are no integrals in the HSQC spectrum
    if (hsqc["integrals"]["count"] == 0) {
        MessageBox.warning("No integrals found in HSQC spectrum, using peaks only. Program works best with integrals");
    }

    // check if the number of peaks in hsqc is == number of multiplets in h1_1d
    if (h1_1d !== undefined) {
        if (hsqc["peaks"]["count"] != h1_1d["multiplets"]["count"]) {
            MessageBox.warning("Number of peaks in HSQC not equal to number of multiplets in H1_1D. Do not use H1_1D");
            return;
        }
    }

    ///////////////////////////////////////
    // Send to the server
    //////////////////////////////////////

    // MNOVA Predict, NMRSHIFTDB2 Predict and MNova Manually Assigned all use simpleMNOVA
    var entry_point = server + "simpleMNOVA";

    var progressDialog = new ProgressDialog();
    progressDialog.labelText = "Performing Prediction ...";
    progressDialog.setRange(0, 100);
    progressDialog.value = 50;
    progressDialog.visible = true;

    var rtn = web_utils.jsonRequest(entry_point, json_spectra_string, "", "", false);

    progressDialog.visible = false;

    var htmloutputstr;
    try {
        var responseObj = JSON.parse(rtn.response || rtn.allStrOutput);
        if (su.registrationProblem(responseObj)) {
            return;
        }
        // response is normal JSON
        htmloutputstr = responseObj.html || responseObj.content;
    } catch (e) {
        // not JSON: the response is the HTML page itself
        htmloutputstr = rtn.response.replace("dummy_title", result.name);
    }

    MessageBox.warning("Prediction complete. Please check the output file.");

    var htmlfilename = su.writeHtml(htmloutputstr, result.directoryPath, result.name);
    su.openInBrowser(htmlfilename);
}


if (this.MnUi && MnUi.simpleNMRtools) {
	MnUi.simpleNMRtools.simpleNMRtools_simplePREDICT_eeh = simplePREDICT_eeh;
}
