function snapToCarbonHSQC()
{
    var doc = Application.mainWindow.activeDocument;
    var simpleutils = new simpleUtils();

    // identify the spectra in the document
    var spectra_lst = identify_spectrum();

    // find the carbon spectrum
    var simpleutils = new simpleUtils();
    // var spectra = simpleutils.identifySpectra(doc);
    var c13spec = simpleutils.findC13spectrum(spectra_lst);

    if( c13spec == undefined){
        MessageBox.warning("No C13 1D spectrum found");
        return;
    }
    else if(c13spec.peaks().count == 0){
        MessageBox.warning("No peaks found in C13 1D spectrum");
        return;
    }
    else {
        var c13peaksList = [];
        
        var c13peaks = c13spec.peaks();
        for (var i = 0; i < c13peaks.count; i++) {
            var pk = c13peaks.at(i);
            c13peaksList.push(pk.delta(1));
        }
    }
    
    // //  print out the carbon peaks
    // for (var i = 0; i < c13peaksList.length; i++) {
    //     print("C13 Peak: " + c13peaksList[i].toFixed(2));
    // }


    // identify the HSQC spectrum
    var hsqc_idx = -1;
    var hsqcspec = undefined;
    for (var i = 0; i < spectra_lst.length; i++) {
        print(spectra_lst[i].experimentEEH);
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


    if (hsqcspec.integrals().count == 0) {
        MessageBox.warning("No integrals found in HSQC spectrum");
        return;
    }

    if (hsqcspec.peaks().count == 0) {
        MessageBox.warning("No peaks found in HSQC spectrum");
        return;
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



    for(var i = 0; i < hsqcspec.integrals().count; i++) {
        var pk = hsqcspec.integrals().at(i);
        if(pk.integralValue() > 0)
            continue;
        else

            var pk_centre1 = pk.rangeMin(1) + (pk.rangeMax(1) - pk.rangeMin(1)) / 2;
            var pk_centre2 = pk.rangeMin(2) + (pk.rangeMax(2) - pk.rangeMin(2)) / 2;
            print("HSQC Integral: " + pk.id + ", " +  pk_centre1.toFixed(2) + ", " +  pk_centre2.toFixed(2));
    }

    var newPeakList = new Peaks();

    for(var i = 0; i < hsqcspec.integrals().count; i++) {
        var pk = hsqcspec.integrals().at(i);

        var newPeakIntensity = pk.integralValue();


        var newPeakDelta1 = pk.rangeMin(1) + (pk.rangeMax(1) - pk.rangeMin(1)) / 2;
        var newPeakDelta2 = pk.rangeMin(2) + (pk.rangeMax(2) - pk.rangeMin(2)) / 2;

        var newPeak = new Peak(newPeakDelta1, newPeakDelta2, newPeakIntensity);

        // add the new peak to the list of peaks
        newPeakList.append(newPeak);
    }

    print("New Peak List: " + newPeakList.count);

    // remove the old peaks from the HSQC spectrum

    var peaks = hsqcspec.peaks();

    for(var i = 0; i < peaks.count; i++) {
        peaks.clear(i)
    }
    // hsqcspec.process();
    hsqcspec.setPeaks(newPeakList);
    doc.update();


    // snap the HSQC peaks to the C13 peaks
    var hsqcPeakList = new Peaks();
    var hsqcpeaks = hsqcspec.peaks();
    for (var i=0; i<hsqcpeaks.count; i++){
        var pk = hsqcpeaks.at(i);
        var delta1 = pk.delta(1);
        var delta2 = pk.delta(2);
        print(idx + "," + delta1.toFixed(2) + "," + delta2.toFixed(2));


        var closestPeakDelta1 = c13peaks.at(0).delta(1);
        var closestPeakDelta2 = c13peaks.at(0).delta(2);
        var smallestDelta = Math.abs(delta1 - closestPeakDelta1);
        var idx = 0
        

        // check if the peak is in the carbon peaks list
        for (var j=1; j<c13peaks.count; j++){
            var delta = Math.abs(delta1 - c13peaks.at(j).delta(1));
            if (delta < closestPeakDelta1){
                closestPeakDelta1 = delta;
                idx = j;
            }
        }
        print(idx + "," + c13peaks.at(idx).delta(1).toFixed(2) + "," + pk.delta(2).toFixed(2));
        hsqcPeakList.append(new Peak(c13peaks.at(idx).delta(1), pk.delta(2), pk.intensity));
    }

    var hsqcpeaks = hsqcspec.peaks();

    for(var i = 0; i < hsqcpeaks.count; i++) {
        hsqcpeaks.clear(i)
    }
    // hsqcspec.process();
    hsqcspec.setPeaks(hsqcPeakList);
    hsqcspec.process();
    doc.update();
}

if (this.MnUi && MnUi.simpleNMRtools) {
	MnUi.simpleNMRtools.simpleNMRtools_snapToCarbonHSQC = snapToCarbonHSQC;
}