function snapToCarbonHSQC_eeh()
{
    var doc = Application.mainWindow.activeDocument;
    var simpleutils = new simpleUtils();

    // identify the spectra in the document
    var spectra_lst = identify_spectrum();

    print("spectra_lst", spectra_lst);

    // find the carbon spectrum
    // var simpleutils = new simpleUtils();
    // var spectra = simpleutils.identifySpectra(doc);
    var c13spec = simpleutils.findC13spectrum(spectra_lst);

    if( c13spec == undefined){
        MessageBox.warning("No C13 1D spectrum found");
        return;
    }
    if(c13spec.peaks().count == 0){
        MessageBox.warning("No peaks found in C13 1D spectrum");
        return;
    }

    // identify the HSQC spectrum
    var hsqc_idx = -1;
    var hsqcspec = undefined;
    for (var i = 0; i < spectra_lst.length; i++) {
        if (spectra_lst[i].experimentEEH == "HSQC") {
            hsqc_idx = i;
            break;
        }
    }
    if (hsqc_idx != -1) {
        hsqcspec = spectra_lst[hsqc_idx];
    }
    else {
        MessageBox.warning("No HSQC spectrum found");
        return;
    }

    if (hsqcspec.peaks().count == 0) {
        MessageBox.warning("No peaks found in HSQC spectrum");
        return;
    }

    // collect the c13 peaks

    var c13peaksList = [];
    
    var c13peaks = c13spec.peaks();
    for (var i = 0; i < c13peaks.count; i++) {
        var pk = c13peaks.at(i);
        c13peaksList.push(pk.delta(1));
    }

    // keep original HSQC peaks both proton and carbon
    var hsqcpeaks = hsqcspec.peaks();
    var hsqcC13peaksList = [];
    var hsqcH1peaksList = [];
    for (var i = 0; i < hsqcpeaks.count; i++) {
        var pk = hsqcpeaks.at(i);
        hsqcC13peaksList.push(pk.delta(1));
        hsqcH1peaksList.push(pk.delta(2));
    }

    // update the HSQC C13 values based on the the C13 values
    var newHSQCpeaks = new Peaks();

    for (var i = 0; i < hsqcpeaks.count; i++) {
        var pk = hsqcpeaks.at(i);
        var delta1 = pk.delta(1);
        var delta2 = pk.delta(2);

        // find the closest peak in the C13 list
        var closestPeakDelta1 = c13peaksList[0];
        
        var smallestC13Delta = Math.abs(delta1 - closestPeakDelta1);
        var idx = 0

        // check if the peak is in the carbon peaks list
        for (var j=1; j<c13peaksList.length; j++){
            var c13delta = Math.abs(delta1 - c13peaksList[j]);
            if (c13delta < smallestC13Delta){
                closestPeakDelta1 = c13peaksList[j];
                smallestC13Delta = c13delta;
                idx = j;
            }
        }
        newHSQCpeaks.append(new Peak(closestPeakDelta1, delta2, 999));
    }

    // var hsqcpeaks = hsqcspec.peaks();
    for(var i = 0; i < hsqcpeaks.count; i++) {
        hsqcpeaks.clear(i)
    }
    // hsqcspec.process();
    hsqcspec.setPeaks(newHSQCpeaks);

    // clear hsqcspec integrals
    for(var i = 0; i < hsqcspec.integrals().count; i++) {
        var pk = hsqcspec.integrals().at(i);
        hsqcspec.integrals().clear(i);
    }

    // add the integrals using simpleintegrate so that we can get the correct 
    // intensities on the HSQC peak table
    simpleIntegrate(hsqcspec);

    hsqcspec.process();
    doc.update();

}

if (this.MnUi && MnUi.simpleNMRtools) {
	MnUi.simpleNMRtools.simpleNMRtools_snapToCarbonHSQC_eeh = snapToCarbonHSQC_eeh;
}