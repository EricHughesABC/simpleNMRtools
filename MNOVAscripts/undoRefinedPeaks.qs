/*globals Application, NMRSpectrum, Peak, File, Dir, TextStream, MessageBox, FileDialog, print*/

// Undo the refined peaks: restore the peaks that were in the spectra before
// readRevisedPeaks.qs replaced them with the refined (simpleNMRclaudePeakRefine) peaks.
//
// readRevisedPeaks.qs saves the original peaks of every spectrum it changes to
//     <document>_peaksBeforeRefine.json
// in the folder of the Mnova document.  This script puts those peaks back - positions,
// intensities and, where Mnova allows it, the peak type (compound, impurity, solvent ...),
// kind, flags and annotation - into the spectra with the same data file name.  Spectra
// that were not changed by the read-back are left alone.
//
// Afterwards the snapshot is marked as restored (it is kept on disk), so the next
// read-back of refined peaks saves a fresh snapshot.  Running this script twice does no
// harm: it asks before restoring an already-restored snapshot.

function undoRefinedPeaks() {

    var su = new simpleUtils();




    // set a peak attribute if this Mnova version allows it
    function trySet(pk, attr, value) {
        if (value === undefined || value === null) return;
        try {
            pk[attr] = value;
        } catch (e) {
            print("could not set peak " + attr + ": " + e);
        }
    }

    var doc = Application.mainWindow.activeDocument;

    // find the snapshot next to the Mnova document; ask for it if it is not there
    var snapshotFilename = "";
    if (doc.name) {
        var details = su.splitPathAndFilename(doc.name);
        snapshotFilename = details.directoryPath + "/" + details.name + "_peaksBeforeRefine.json";
    }
    if (!snapshotFilename || !new Dir(snapshotFilename).fileExists(snapshotFilename)) {
        var startDir = doc.name ? su.splitPathAndFilename(doc.name).directoryPath : "";
        snapshotFilename = FileDialog.getOpenFileName("*_peaksBeforeRefine.json",
                                                      "Choose the saved peaks to restore", startDir);
        if (!snapshotFilename) {
            print("No snapshot chosen - nothing undone.");
            return;
        }
    }

    var snap = su.readJson(snapshotFilename, null);
    if (!snap || snap.datatype !== "peaksBeforeRefine" || !snap.spectra) {
        MessageBox.warning("Not a saved-peaks file:\n" + snapshotFilename);
        return;
    }
    if (snap.restored) {
        var answer = MessageBox.question("These peaks were already restored on " + snap.restoredOn +
                                         ".\nRestore them again?",
                                         MessageBox.Yes, MessageBox.No);
        if (answer != MessageBox.Yes) {
            return;
        }
    }

    var restored = [];
    for (var i = 0; i < doc.pageCount(); i++) {
        var page = doc.page(i);
        for (var j = 0; j < page.itemCount(); j++) {
            var spec = new NMRSpectrum(page.item(j));
            if (!spec.isValid()) continue;
            var df = spec.getParam("Data File Name");
            var saved = snap.spectra[df];
            if (!saved) continue;

            var newPeaks = spec.peaks();
            newPeaks.clear();
            for (var k = 0; k < saved.peaks.length; k++) {
                var p = saved.peaks[k];
                var pk = saved.dimCount == 2 ? new Peak(p.delta1, p.delta2, p.intensity)
                                             : new Peak(p.delta1, p.intensity);
                trySet(pk, "type", p.type);
                trySet(pk, "kind", p.kind);
                trySet(pk, "flags", p.flags);
                trySet(pk, "annotation", p.annotation);
                newPeaks.append(pk);
            }
            spec.setPeaks(newPeaks);
            spec.update();
            restored.push((saved.key || saved.title) + " (" + saved.peaks.length + " peaks)");
            print("Restored " + saved.peaks.length + " peaks in " + df);
        }
    }

    if (restored.length === 0) {
        MessageBox.warning("None of the saved spectra were found in this document.\n" +
                           "Is the right Mnova document open?\n\n" + snapshotFilename);
        return;
    }

    snap.restored = true;
    snap.restoredOn = new Date().toString();
    su.writeJson(snapshotFilename, snap, 2);

    var missing = [];
    for (var key in snap.spectra) {
        var found = false;
        for (var r = 0; r < restored.length; r++) {
            if (restored[r].indexOf(snap.spectra[key].key || snap.spectra[key].title) === 0) found = true;
        }
        if (!found) missing.push(snap.spectra[key].key || key);
    }
    MessageBox.information("Unrefined peaks restored:\n  " + restored.join("\n  ") +
                           (missing.length ? "\n\nNot found in this document:\n  " + missing.join("\n  ") : ""));
}
