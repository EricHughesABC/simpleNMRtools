function readRevisedPeaks() {

    // read in the revised peaks json file
    function readRevisedPeaksFile(jsonfilename){
        var peaksTxt = "";
        var fin = new File(jsonfilename);
        if (fin.open(File.ReadOnly)) {
            var sin = new TextStream(fin, 'UTF-8');
            peaksTxt = sin.readAll();
        }
        else{
            MessageBox.warning("Could not open json file\n" + jsonfilename);
        }
        fin.close();
        return JSON.parse(peaksTxt);
    }


    function splitPathAndFilename(path) {
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

    var doc = Application.mainWindow.activeDocument;
    var nmrspec = new NMRSpectrum(nmr.activeSpectrum());
    var specFileName = nmrspec.getParam("Data File Name");
    print(specFileName);

    var path = doc.name;
    print(path);
    var filePathDetails = splitPathAndFilename(path);
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
        revisedPeaks = readRevisedPeaksFile(revisedPeaksDialog);
    } else {
        print("No file selected.");
        return;
    }

    // return if revised peaks is null
    if (!revisedPeaks) {
        return;
    }

    var hsqc_data = revisedPeaks["HSQC_0"];
    var hmbc_data = revisedPeaks["HMBC_0"];
    var c13_data = revisedPeaks["C13_1D_0"];
    var h1_data = revisedPeaks["H1_1D_0"];
    var pureshift_data = revisedPeaks["PureShift_0"];
    var ddeptCH3_data = revisedPeaks["DDEPTCH3ONLY_0"];
    var hsqcclipcosy_data = revisedPeaks["HSQC_CLIPCOSY_0"];
    var cosy_data = revisedPeaks["COSY_0"];


    // check  if HSQC exists
    if (!hsqc_data) {
        print("No HSQC data found.");
        return;
    }

    // print out the number of HSQC peaks
    print("Number of HSQC peaks: " + hsqc_data.peaks.count);

    var _hsqc_filename = hsqc_data.datafilename;


    // search for HSQC spectrum in mnova document
    for (i = 0, pageCount = doc.pageCount(); i < pageCount; i++) {
        page = doc.page(i);
        for (j = 0, itemCount = page.itemCount(); j < itemCount; j++) {
            spec = new NMRSpectrum(page.item(j));
      
            if( spec.isValid() ){
                var specFile = spec.getParam("Data File Name");
                if (specFile === _hsqc_filename) {
                    print("Found matching HSQC spectrum: " + specFile);
                    // You can add further processing for the matching spectrum here

                    // delete peaks in spec
                    var oldPeaks = spec.peaks();
                    // remove peaks from the spectrum using clear method
                    oldPeaks.clear();
                    // update the spectrum to reflect the changes
                    spec.update();


                    var hsqcPeaks = spec.peaks();
                    for (var k = 0; k < hsqc_data.peaks.count; k++) {
                        var peak = hsqc_data.peaks.data[k];
                        // only assign valid peaks
                        if (peak.type == 0) {
                            // print(peak.delta1 + ", " + peak.delta2 + ", " + peak.intensity);
                            hsqcPeaks.append(Peak(peak.delta1, peak.delta2, peak.intensity));
                        }
                    }
                    spec.setPeaks(hsqcPeaks);
                    // spec.process()
                    spec.update();

                }
            }
        }
    }

    // replace the HMBC peaks
        // search for HMBC spectrum in mnova document
    var _hmbc_filename = hmbc_data.datafilename;
    for (i = 0, pageCount = doc.pageCount(); i < pageCount; i++) {
        page = doc.page(i);
        for (j = 0, itemCount = page.itemCount(); j < itemCount; j++) {
            spec = new NMRSpectrum(page.item(j));
      
            if( spec.isValid() ){
                var specFile = spec.getParam("Data File Name");
                if (specFile === _hmbc_filename) {
                    print("Found matching HMBC spectrum: " + specFile);
                    // You can add further processing for the matching spectrum here

                    // delete peaks in spec
                    var oldPeaks = spec.peaks();
                    // remove peaks from the spectrum using clear method
                    oldPeaks.clear();
                    // update the spectrum to reflect the changes
                    spec.update();


                    var hmbcPeaks = spec.peaks();
                    for (var k = 0; k < hmbc_data.peaks.count; k++) {
                        var peak = hmbc_data.peaks.data[k];
                        // only assign valid peaks
                        if (peak.type == 0) {
                            // print(peak.delta1 + ", " + peak.delta2 + ", " + peak.intensity);
                            hmbcPeaks.append(Peak(peak.delta1, peak.delta2, peak.intensity));
                        }
                    }
                    spec.setPeaks(hmbcPeaks);
                    // spec.process()
                    spec.update();

                }
            }
        }
    }

    if (c13_data) {
        var _c13_filename = c13_data.datafilename;
        for (i = 0, pageCount = doc.pageCount(); i < pageCount; i++) {
            page = doc.page(i);
            for (j = 0, itemCount = page.itemCount(); j < itemCount; j++) {
                spec = new NMRSpectrum(page.item(j));
                if( spec.isValid() ){
                    var specFile = spec.getParam("Data File Name");
                    if (specFile === _c13_filename) {
                        print("Found matching C13 spectrum: " + specFile);
                        var oldPeaks = spec.peaks();
                        oldPeaks.clear();
                        spec.update();

                        var c13Peaks = spec.peaks();
                        for (var k = 0; k < c13_data.peaks.count; k++) {
                            var peak = c13_data.peaks.data[k];
                            if (peak.type == 0) {
                                // print(peak.delta1 + ", " + peak.intensity);
                                c13Peaks.append(Peak(peak.delta1, peak.intensity));
                            }
                        }
                        spec.setPeaks(c13Peaks);
                        // spec.process()
                        spec.update();
                    }
                }
            }
        }
    }

    if (h1_data) {
        var _h1_filename = h1_data.datafilename;
        for (i = 0, pageCount = doc.pageCount(); i < pageCount; i++) {
            page = doc.page(i);
            for (j = 0, itemCount = page.itemCount(); j < itemCount; j++) {
                spec = new NMRSpectrum(page.item(j));
                if( spec.isValid() ){
                    var specFile = spec.getParam("Data File Name");
                    if (specFile === _h1_filename) {
                        print("Found matching H1 spectrum: " + specFile);
                        var oldPeaks = spec.peaks();
                        oldPeaks.clear();
                        spec.update();

                        var h1Peaks = spec.peaks();
                        for (var k = 0; k < h1_data.peaks.count; k++) {
                            var peak = h1_data.peaks.data[k];
                            if (peak.type == 0) {
                                // print(peak.delta1 + ", " + peak.intensity);
                                h1Peaks.append(Peak(peak.delta1, peak.intensity));
                            }
                        }
                        spec.setPeaks(h1Peaks);
                        // spec.process()
                        spec.update();
                    }
                }
            }
        }
    }

    if (pureshift_data) {
        var _pureshift_filename = pureshift_data.datafilename;
        for (i = 0, pageCount = doc.pageCount(); i < pageCount; i++) {
            page = doc.page(i);
            for (j = 0, itemCount = page.itemCount(); j < itemCount; j++) {
                spec = new NMRSpectrum(page.item(j));
                if( spec.isValid() ){
                    var specFile = spec.getParam("Data File Name");
                    if (specFile === _pureshift_filename) {
                        print("Found matching Pure Shift spectrum: " + specFile);
                        var oldPeaks = spec.peaks();
                        oldPeaks.clear();
                        spec.update();

                        var pureshiftPeaks = spec.peaks();
                        for (var k = 0; k < pureshift_data.peaks.count; k++) {
                            var peak = pureshift_data.peaks.data[k];
                            if (peak.type == 0) {
                                // print(peak.delta1 + ", " + peak.intensity);
                                pureshiftPeaks.append(Peak(peak.delta1, peak.intensity));
                            }
                        }
                        spec.setPeaks(pureshiftPeaks);
                        // spec.process()
                        spec.update();
                    }
                }
            }
        }
    }

    if (ddeptCH3_data) {
        var _ddeptCH3_filename = ddeptCH3_data.datafilename;
        for (i = 0, pageCount = doc.pageCount(); i < pageCount; i++) {
            page = doc.page(i);
            for (j = 0, itemCount = page.itemCount(); j < itemCount; j++) {
                spec = new NMRSpectrum(page.item(j));
                if( spec.isValid() ){
                    var specFile = spec.getParam("Data File Name");
                    if (specFile === _ddeptCH3_filename) {
                        print("Found matching DDEPT CH3 spectrum: " + specFile);
                        var oldPeaks = spec.peaks();
                        oldPeaks.clear();
                        spec.update();

                        var ddeptCH3Peaks = spec.peaks();
                        for (var k = 0; k < ddeptCH3_data.peaks.count; k++) {
                            var peak = ddeptCH3_data.peaks.data[k];
                            if (peak.type == 0) {
                                // print(peak.delta1 + ", " + peak.intensity);
                                ddeptCH3Peaks.append(Peak(peak.delta1, peak.delta2,peak.intensity));
                            }
                        }
                        spec.setPeaks(ddeptCH3Peaks);
                        // spec.process()
                        spec.update();
                    }
                }
            }
        }
    }


    if (hsqcclipcosy_data) {
        var _hsqcclipcosy_filename = hsqcclipcosy_data.datafilename;
        print("Looking for HSQC CLIP COSY spectrum: " + _hsqcclipcosy_filename);
        for (var i = 0, pageCount = doc.pageCount(); i < pageCount; i++) {
            var page = doc.page(i);
            for (var j = 0, itemCount = page.itemCount(); j < itemCount; j++) {
                var spec = new NMRSpectrum(page.item(j));
                if (!spec.isValid()) continue;

                var specFile = spec.getParam("Data File Name");
                if (specFile !== _hsqcclipcosy_filename) continue;

                print("Found matching HSQC CLIP COSY spectrum: " + specFile);

                // Take one copy, empty it, and build the new list in it
                var newPeaks = spec.peaks();
                newPeaks.clear();

                var nKept = 0, nSkipped = 0;
                for (var k = 0; k < hsqcclipcosy_data.peaks.count; k++) {
                    var peak = hsqcclipcosy_data.peaks.data[k];
                    var inten = Number(peak.intensity);
                    if (Number(peak.type) === 0 && inten < 0) {
                        newPeaks.append(Peak(Number(peak.delta1), Number(peak.delta2), inten));
                        nKept++;
                    } else {
                        nSkipped++;
                    }
                }
                print("Kept " + nKept + " negative peaks, skipped " + nSkipped);

                spec.setPeaks(newPeaks);   // no spec.process() after this
                spec.update();
            }
        }
    }
    else {
        print("No HSQC CLIP COSY data found");
    }

    if (cosy_data) {
        var _cosy_filename = cosy_data.datafilename;
        for (i = 0, pageCount = doc.pageCount(); i < pageCount; i++) {
            page = doc.page(i);
            for (j = 0, itemCount = page.itemCount(); j < itemCount; j++) {
                spec = new NMRSpectrum(page.item(j));
                if( spec.isValid() ){
                    var specFile = spec.getParam("Data File Name");
                    if (specFile === _cosy_filename) {
                        print("Found matching COSY spectrum: " + specFile);
                        var oldPeaks = spec.peaks();
                        oldPeaks.clear();
                        spec.update();

                        var cosyPeaks = spec.peaks();
                        for (var k = 0; k < cosy_data.peaks.count; k++) {
                            var peak = cosy_data.peaks.data[k];
                            if (peak.type == 0) {
                                // print(peak.delta1 + ", " + peak.delta2 + ", " + peak.intensity);
                                cosyPeaks.append(Peak(peak.delta1, peak.delta2, peak.intensity));
                            }
                        }
                        spec.setPeaks(cosyPeaks);
                        spec.process()
                        spec.update();
                    }
                }
            }
        }
    }
}