// Read the refined peaks (simpleNMRclaudePeakRefine output) back into Mnova.
//
// Before anything is changed, the current peaks of every spectrum that will be
// replaced are saved to <document>_peaksBeforeRefine.json next to the Mnova document
// (see savePeaksBeforeRefine below).  undoRefinedPeaks.qs restores them.  The snapshot
// is only created when there is none, so reading refined peaks several times in a row
// still keeps the original, unrefined peaks; spectra not yet in the snapshot are added.

function readRevisedPeaks() {

    var su = new simpleUtils();

    // The spectra this script replaces, in this order (experiment names from
    // SIMPLENMR_EXPERIMENTS in simpleUtils.qs; the refined JSON keys are <name>_0).
    //   negativeOnly: keep only compound peaks with negative intensity (HSQC-CLIP-COSY)
    //   process:      re-process the spectrum after setting the peaks (COSY)
    var REPLACED = [{"name": EXPT.HSQC},
                    {"name": EXPT.HMBC},
                    {"name": EXPT.C13_1D},
                    {"name": EXPT.H1_1D},
                    {"name": EXPT.H1_pureshift},
                    {"name": EXPT.DDEPTCH3ONLY},
                    {"name": EXPT.HSQC_CLIPCOSY, "negativeOnly": true},
                    {"name": EXPT.COSY, "process": true}];
    var REPLACED_KEYS = [];
    for (var r = 0; r < REPLACED.length; r++) {
        REPLACED_KEYS.push(exptKey(REPLACED[r].name, 0));
    }

    // Replace the peaks of every spectrum in the document whose data file is datafilename
    // with the compound peaks of block (one refined spectrum from the JSON).
    function replacePeaks(doc, block, opts) {
        var twoD = exptIsTwoD(opts.name);
        for (var i = 0, pageCount = doc.pageCount(); i < pageCount; i++) {
            var page = doc.page(i);
            for (var j = 0, itemCount = page.itemCount(); j < itemCount; j++) {
                var spec = new NMRSpectrum(page.item(j));
                if (!spec.isValid() || spec.getParam("Data File Name") !== block.datafilename) {
                    continue;
                }
                print("Found matching " + opts.name + " spectrum: " + block.datafilename);

                var newPeaks = spec.peaks();
                newPeaks.clear();
                if (!opts.negativeOnly) {
                    spec.update();
                    newPeaks = spec.peaks();
                }
                var nKept = 0, nSkipped = 0;
                for (var k = 0; k < block.peaks.count; k++) {
                    var peak = block.peaks.data[k];
                    var keep = (Number(peak.type) === 0) && (!opts.negativeOnly || Number(peak.intensity) < 0);
                    if (!keep) {
                        nSkipped++;
                        continue;
                    }
                    newPeaks.append(twoD ? Peak(Number(peak.delta1), Number(peak.delta2), Number(peak.intensity))
                                         : Peak(Number(peak.delta1), Number(peak.intensity)));
                    nKept++;
                }
                if (opts.negativeOnly) {
                    print("Kept " + nKept + " negative peaks, skipped " + nSkipped);
                }
                spec.setPeaks(newPeaks);
                if (opts.process) {
                    spec.process();
                }
                spec.update();
            }
        }
    }



    // all peaks of one Mnova spectrum, with every attribute needed to restore them
    function peaksOf(spec) {
        var out = [];
        var pks = spec.peaks();
        var twoD = spec.dimCount == 2 || String(spec.getParam("Nucleus")).indexOf(",") != -1;
        for (var p = 0; p < pks.count; p++) {
            var pk = pks.at(p);
            out.push({
                "delta1": pk.delta(1),
                "delta2": twoD ? pk.delta(2) : 0,
                "intensity": pk.intensity,
                "type": pk.type,
                "kind": pk.kind,
                "flags": pk.flags,
                "annotation": pk.annotation
            });
        }
        return out;
    }

    // Save the current peaks of the spectra about to be replaced (matched by data file name).
    function savePeaksBeforeRefine(doc, refined, snapshotFilename) {
        var dfs = {};
        for (var r = 0; r < REPLACED_KEYS.length; r++) {
            var blk = refined[REPLACED_KEYS[r]];
            if (blk && blk.datafilename) {
                dfs[blk.datafilename] = REPLACED_KEYS[r];
            }
        }
        var snap = null;
        if (new Dir(snapshotFilename).fileExists(snapshotFilename)) {
            snap = su.readJson(snapshotFilename, null);
            if (snap && snap.restored) {
                snap = null;            // already undone: start a fresh snapshot
            }
        }
        if (!snap) {
            snap = {"datatype": "peaksBeforeRefine", "version": 1, "document": doc.name,
                    "created": new Date().toString(), "restored": false, "spectra": {}};
        }
        var added = 0;
        for (var i = 0; i < doc.pageCount(); i++) {
            var page = doc.page(i);
            for (var j = 0; j < page.itemCount(); j++) {
                var spec = new NMRSpectrum(page.item(j));
                if (!spec.isValid()) continue;
                var df = spec.getParam("Data File Name");
                if (!(df in dfs) || (df in snap.spectra)) continue;
                snap.spectra[df] = {"key": dfs[df], "title": spec.title,
                                    "dimCount": (spec.dimCount == 2 || String(spec.getParam("Nucleus")).indexOf(",") != -1) ? 2 : 1,
                                    "peaks": peaksOf(spec)};
                added++;
            }
        }
        if (added > 0) {
            if (su.writeJson(snapshotFilename, snap, 2)) {
                print("Saved the peaks of " + added + " spectra before refinement to " + snapshotFilename);
            } else {
                MessageBox.warning("Could not write " + snapshotFilename + "\nThe peaks can not be undone.");
            }
        } else {
            print("Snapshot " + snapshotFilename + " already holds the unrefined peaks");
        }
    }




    var doc = Application.mainWindow.activeDocument;
    var nmrspec = new NMRSpectrum(nmr.activeSpectrum());
    var specFileName = nmrspec.getParam("Data File Name");
    print(specFileName);

    var path = doc.name;
    print(path);
    var filePathDetails = su.splitPathAndFilename(path);
    print(filePathDetails);

    // // split specFilename to get the path upto and including craft
    // var specPath = specFileName.substring(0, specFileName.lastIndexOf("craft") + 5);
    // print(specPath);

    // // add the rest of the path to where revised peaks json file is located
    // specPath += "/craftdbox5b/CRAFTpublish";

    var revisedPeaksDialog = FileDialog.getOpenFileName("*.json", "Read revised peaks file", filePathDetails.directoryPath);

    print(revisedPeaksDialog);

    var revisedPeaks = null;

    if (revisedPeaksDialog) {
        print("Selected revised peaks file: " + revisedPeaksDialog);
        revisedPeaks = su.readJson(revisedPeaksDialog, null, true);
    } else {
        print("No file selected.");
        return;
    }

    // return if revised peaks is null
    if (!revisedPeaks) {
        return;
    }

    // the HSQC is required
    var hsqc_data = revisedPeaks[exptKey(EXPT.HSQC, 0)];
    if (!hsqc_data) {
        print("No HSQC data found.");
        return;
    }

    // keep the current (unrefined) peaks so that undoRefinedPeaks.qs can restore them
    var snapDetails = doc.name ? su.splitPathAndFilename(doc.name)
                               : su.splitPathAndFilename(revisedPeaksDialog.replace(/\.json$/i, "") + ".x");
    var snapshotFilename = snapDetails.directoryPath + "/" + snapDetails.name + "_peaksBeforeRefine.json";
    savePeaksBeforeRefine(doc, revisedPeaks, snapshotFilename);

    // print out the number of HSQC peaks
    print("Number of HSQC peaks: " + hsqc_data.peaks.count);

    // replace the peaks of each refined spectrum that is in the file
    for (var s = 0; s < REPLACED.length; s++) {
        var block = revisedPeaks[exptKey(REPLACED[s].name, 0)];
        if (block) {
            replacePeaks(doc, block, REPLACED[s]);
        }
        else {
            print("No " + REPLACED[s].name + " data found");
        }
    }
}