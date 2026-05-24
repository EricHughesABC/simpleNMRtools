function peakPickHSQC_fromC131D(){

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
    if(c13spec.peaks().count == 0){
        MessageBox.warning("No peaks found in C13 1D spectrum");
        return;
    }
    else if(c13spec.peaks().count > 0){
        print("C13 1D spectrum found with " + c13spec.peaks().count + " peaks");
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


    // print("HSQC spectrum found with " + hsqcspec.peaks().count + " peaks");
    // print("hsqcspec.isValid() = " + hsqcspec.isValid());
    // print("hsqcspec.cout(1) = " + hsqcspec.count(1));
    // print("hsqcspec.cout(2) = " + hsqcspec.count(2));
    // print("hsqcspec.scaleWidth(1) = " + hsqcspec.scaleWidth(1));
    // print("hsqcspec.scaleWidth(2) = " + hsqcspec.scaleWidth(2));
    // print("hsqcspec.hz(1) = " + hsqcspec.hz(1));
    // print("hsqcspec.hz(2) = " + hsqcspec.hz(2));
    var f1npts = hsqcspec.count(1);
    var f2npts = hsqcspec.count(2);

    var scalelimits = hsqcspec.getFullScaleLimits();
    // print(scalelimits)
    var f1ppm_min = parseFloat(scalelimits.fromY);
    var f1ppm_max = parseFloat(scalelimits.toY);

    // print("f1ppm_min = " + f1ppm_min);
    // print("f1ppm_max = " + f1ppm_max);

    var ppmTopoints  = f1npts / (f1ppm_max - f1ppm_min);

    // print("ppmTopoints = " + ppmTopoints);

    function c13ppmToF1Point(c13ppm, f1ppm_min, f1ppm_max, ppmTopoints, f1npts) {
        // c13ppm: the C13 chemical shift in ppm
        // Returns the corresponding point index in the F1 dimension of the HSQC spectrum

        // Ensure c13ppm is within the F1 ppm range
        if (c13ppm < f1ppm_min || c13ppm > f1ppm_max) {
            print("Warning: c13ppm value " + c13ppm + " is outside the F1 ppm range.");
            return -1; // Return -1 or handle as needed
        }

        // Calculate the point index (0-based)
        var point = f1npts - Math.round((c13ppm - f1ppm_min) * ppmTopoints);

        // Clamp to valid range
        if (point < 0) point = 0;
        if (point >= f1npts) point = f1npts - 1;

        return point;
    }


    var idx = 7; // Example index, change as needed
    var baseline_pnts = 20; // Number of points to consider for baseline
    var pk = c13spec.peaks().at(idx);
    var c13ppm = pk.delta(1);
    var f1point = c13ppmToF1Point(c13ppm, f1ppm_min, f1ppm_max, ppmTopoints, f1npts);
  

    print("C13 peak " + idx + " at " + c13ppm + " ppm corresponds to F1 point " + f1point);


    var pk_offset = -4 //; Offset to include surrounding points
	var rrr = hsqcspec.real({from:f1point-pk_offset, to:f1point-pk_offset}, {from:0,to:f2npts});
    var rrr_2D = hsqcspec.real({from:f1point-pk_offset, to:f1point+pk_offset}, {from:0,to:f2npts});

    print("rrr_2D = " + JSON.stringify(rrr_2D));

    for (var i=f1point-pk_offset+1; i <= f1point + pk_offset; i++) {
        var rrr_temp = hsqcspec.real({from:i, to:i}, {from:0,to:f2npts});
        print("i = " + i);
        // add the values in rrr_temp to rrr
        for (var j=0; j < f2npts; j++) {
            rrr[j] += rrr_temp[j];
        }

    }
    var maxY = Math.max.apply(Math, rrr);
    var maxIndex = rrr.indexOf(maxY);

    print("Max value in the real spectrum at F1 point " + f1point + " is " + maxY + " at index " + maxIndex);

    var yyy = rrr.slice(maxIndex - baseline_pnts, maxIndex + baseline_pnts);

    var xxx = [];
    for(var i=maxIndex - baseline_pnts; i < maxIndex + baseline_pnts; i++){
        xxx.push(i);
    }

    // print the length of yyy and xxx
    print("Length of yyy: " + yyy.length);
    print("Length of xxx: " + xxx.length);

    var fitter = GaussianFitter();
    var initialGuess = fitter.generateInitialGuess(xxx, yyy);

    var fitResult = fitter.fit(xxx, yyy);

    print('Initial Guess:');
    print('--------------------');
    print(initialGuess)
    // print('Amplitude: ' + initialGuess.amplitude.toFixed(2));
    // print('Center: ' + initialGuess.center.toFixed(2));
    // print('Sigma: ' + initialGuess.sigma.toFixed(2));
    // print('Baseline: ' + initialGuess.baseline.toFixed(2));
    print('--------------------');

    var init_yyy = [];
    for (var i = 0; i < xxx.length; i++) {
        init_yyy.push(fitter.gaussian( xxx[i], initialGuess));
    }


    print('Fitted Gaussian Peak:');
    print('Amplitude: ' + fitResult.amplitude.toFixed(2));
    print('Center: ' + fitResult.center.toFixed(2));
    print('Sigma: ' + fitResult.sigma.toFixed(2));
    print('Baseline: ' + fitResult.baseline.toFixed(2));
    print('R-squared: ' + fitResult.rSquared.toFixed(4));
    print('Iterations: ' + fitResult.iterations);
    print('Converged: ' + fitResult.converged);
    print('xxx = ' + JSON.stringify(xxx));
    print('yyy = ' + JSON.stringify(yyy));
    print('yyy_res = ' + JSON.stringify(fitResult.fittedY));
    print('init_plot_x = ' + JSON.stringify(xxx));
    print('init_plot_y = ' + JSON.stringify(init_yyy));
    



	

}