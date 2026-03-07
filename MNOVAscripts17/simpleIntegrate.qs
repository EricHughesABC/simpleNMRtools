/*globals  MnUi*/


function simpleIntegrate(nmrspecpassed) {

    // check if a spectrum has been passed across or is undefined
    if (nmrspecpassed === undefined) {
        print("smpleIntegrate: nmrspec is undefined");
        // nmrspec = nmr.activeSpectrum();
    }
    else{
        print("smpleIntegrate: nmrspec is defined");
    }

    // iterate through pages and print out what is on the page
    var doc = Application.mainWindow.activeDocument;
    var simpleutils = new simpleUtils();

    if (nmrspecpassed === undefined) {

        var nmrspec = new NMRSpectrum(nmr.activeSpectrum());
    }
    else{
        var nmrspec = new NMRSpectrum(nmrspecpassed);
    }
    // check if the spectrum is a 2D spectrum
    if(!simpleutils.is2DSpectrum(nmrspec)){
        return;
    }

    // get the peak list
    var peaks = nmrspec.peaks();
    // save peak coordinates
    var delta1 = [];
    var delta2 = [];
    for( var i = 0; i < peaks.length; i++ ) {
        var pk = peaks.at(i);
        delta1.push(pk.delta(1));
        delta2.push(pk.delta(2));
    }
    var peaks_bckup = JSON.parse(JSON.stringify(peaks));
    var peakCount = peaks.count;
    if(peakCount == 0){
        MessageBox.warning("No peaks found in 2D spectrum");
        return;
    }



    var nucleus = nmrspec.getParam("Nucleus");
    // split string on ,
    var nucleusArray = nucleus.split(",");
    // remove square brackets from the strings in the array
    // and save the string back in the array
    for (var i=0; i<nucleusArray.length; i++){
        nucleusArray[i] = nucleusArray[i].replace("[", "");
        nucleusArray[i] = nucleusArray[i].replace("]", "");
    }

    // remove any whitespace from the strings in the array
    for (var i=0; i<nucleusArray.length; i++){
        nucleusArray[i] = nucleusArray[i].replace(/\s/g, "");
    }

    // set the default peak width for the integrals
    var pk_pm1 = 0.1;
    var pk_pm2 = 0.1;

    if( nucleusArray[0] == "1H"){
        pk_pm1 = 0.025;
    }
    else {
        pk_pm1 = 0.25;
    }

    if( nucleusArray[1] == "1H"){
        pk_pm2 = 0.025;
    }
    else {
        pk_pm2 = 0.25;
    }   

    // create a new integral list
    var intList = new Integrals();
    for( var i = 0; i < peaks.length; i++ ) {
        var pk = peaks.at(i);
        
        var sReg = new SpectrumRegion(pk.delta(2) - pk_pm2,
                                      pk.delta(2) + pk_pm2,
                                      pk.delta(1) - pk_pm1,
                                      pk.delta(1) + pk_pm1);
        var newInt = new Integral(nmrspec, sReg, false);

        intList.append(newInt);
    }

    //Replace spectrum integrals
	nmrspec.setIntegrals(intList);
    // nmrspec.process();
    // get the peaks from the spectrum again
    var peaks = nmrspec.peaks();
    // clear the peaks from the spectrum
    for( var i = 0; i < peaks.count; i++ ) {
        peaks.clear(i);
    }

    // replace them with original peaks as the integrals may add new peaks
    // and definitely change the position of the peaks slightly
    for( var i = 0; i < peaks_bckup.count; i++ ) {
        var p = new Peak(delta1[i], delta2[i], nmrspec);//single peak
        peaks.append(p);
    }
    nmrspec.setPeaks(peaks);
	nmrspec.process();
}

if (this.MnUi && MnUi.simpleNMRtools) {
    MnUi.simpleNMRtools.simpleNMRtools_simpleIntegrate = simpleIntegrate;
}