function pickPeaks2Dregion()
{
    // get simpleutils functions
    var simpleutils = new simpleUtils();
    var doc = Application.mainWindow.activeDocument;
    // get active spectrum
    var activeSpectrum = new NMRSpectrum(nmr.activeSpectrum());

    // if not 2D spectrum, exit with message
    if(!simpleutils.is2DSpectrum(activeSpectrum)){
        MessageBox.warning("Current spectrum is not a 2D spectrum. Please select a 2D spectrum and try again.");
        return;
    }

    // get zoomed region

    var zoomLimits = activeSpectrum.getScaleLimits();
    


    // get peaks in zoomed region
    var peaksRegion = SpectrumRegion( zoomLimits.fromX, zoomLimits.toX, zoomLimits.fromY, zoomLimits.toY);
    var allPeaksFound = Peaks(activeSpectrum, peaksRegion);

    print(allPeaksFound);

        // extract delta(1) and delta(2) values from the peaks
    for( var i = 0; i < allPeaksFound.count; i++ ) {
        var pk = allPeaksFound.at(i);
        print(pk.delta(2), pk.delta(1));
    }


    // add peaks to the 2D spectrum
    activeSpectrum.setPeaks(allPeaksFound);
    activeSpectrum.process();

    print(zoomLimits);
    // print number of peaks found
    MessageBox.information("Number of peaks found in the zoomed region: " + allPeaksFound.length);



}