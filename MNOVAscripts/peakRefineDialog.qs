/*globals Dialog, GroupBox, CheckBox, ComboBox, print*/

// Spectrum-selection dialog for exporting the raw peaks for automatic peak refinement
// (simpleNMRclaudePeakRefine).  Called by collectPeaks.qs.
//
// Unlike the simpleNMR verify dialog (xy3_dialog2.qs) it only asks which spectra to use
// and what each one is.  The prediction method, simulated annealing, machine-learning
// consent, previous node positions and the server are not needed for the refinement.
//
// Arguments
//   spectra_available : ["<nucleus> <experiment> <pulse sequence> <spectrum key>", ...]
//   previousChoices   : the labels chosen last time, either
//                         {"<spectrum string>": "<label>", ...}  (saved by this dialog), or
//                         ["<label>", ...]  (old positional list, e.g. the verify step's
//                                            exptIdentifiers.json; used only when it has one
//                                            entry per spectrum)
//   exptName          : name of the Mnova document, shown in the title
//
// Returns undefined if cancelled, otherwise
//   { "checked_spectra": ["<spectrum string> <label>", ...],   // ticked spectra only
//     "exptIdentifiers": ["<label>", ...],                      // one per spectrum, in order
//     "choices": {"<spectrum string>": "<label>", ...} }         // saved for next time

function peakRefineDialog(spectra_available, previousChoices, exptName) {

    // SKIP and every experiment type, from SIMPLENMR_EXPERIMENTS in simpleUtils.qs
    var NMRexpts = exptDialogList();

    var dialog = new Dialog();
    dialog.title = "Export peaks for automatic refinement - " + exptName;

    var checkboxes = [];
    var comboboxes = [];

    for (var i = 0; i < spectra_available.length; i++) {
        var cb = new CheckBox();
        cb.text = spectra_available[i];
        cb.checked = true;

        var combo = new ComboBox();
        for (var j = 0; j < NMRexpts.length; j++) {
            combo.addItem(NMRexpts[j]);
        }

        // restore the previous choice: by spectrum name, or by position for an old list
        var previous = undefined;
        if (previousChoices && !(previousChoices instanceof Array)) {
            previous = previousChoices[spectra_available[i]];
        } else if (previousChoices instanceof Array && previousChoices.length == spectra_available.length) {
            previous = previousChoices[i];
        }
        if (previous !== undefined && NMRexpts.indexOf(previous) != -1) {
            combo.currentText = previous;
        }

        checkboxes.push(cb);
        comboboxes.push(combo);
    }

    var gbSpectra = new GroupBox();
    gbSpectra.title = "Choose what each spectrum is from the drop-downs (SKIP to leave it out)";
    for (var k = 0; k < checkboxes.length; k++) {
        gbSpectra.add(checkboxes[k], comboboxes[k]);
    }
    dialog.add(gbSpectra);

    if (!dialog.exec()) {
        return undefined;
    }

    var rtn = {"checked_spectra": [], "exptIdentifiers": [], "choices": {}};
    for (var m = 0; m < spectra_available.length; m++) {
        var label = comboboxes[m].currentText;
        if (checkboxes[m].checked) {
            rtn["checked_spectra"].push(spectra_available[m] + " " + label);
            rtn["exptIdentifiers"].push(label);
        } else {
            // an unticked spectrum is not exported, so it must not take a <label>_<n> number
            rtn["exptIdentifiers"].push("SKIP");
        }
        rtn["choices"][spectra_available[m]] = label;     // remembered even when unticked
    }
    print("peakRefineDialog: ", JSON.stringify(rtn));
    return rtn;
}
