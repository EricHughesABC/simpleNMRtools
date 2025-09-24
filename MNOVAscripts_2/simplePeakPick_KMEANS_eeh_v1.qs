/*globals WebUtilsQT NMRAssignments xy3_dialog Dir brukerExperiments String iupac_dialog server_address MnUi*/


function simplePeakPick_KMEANS_eeh_v1() {

    function find2DPeaks(f2ppm_list, f1ppm_list, spec2D) {

        var simpleutils = new simpleUtils();

        print("f2ppm_list ", f2ppm_list);
        print("f1ppm_list ", f1ppm_list);
     
        // var f2peaks = f2spec.peaks();
        // var f1peaks = f1spec.peaks();
    
        // // get the peak ppm values from the f2 peaks
        // var f2ppm_list = [];
        // for( var i = 0; i < f2peaks.length; i++ ) {
        //     var pk = f2peaks.at(i);
        //     f2ppm_list.push(pk.delta(1));
        // }
    
        // var f1ppm_list = [];
        // for( var i = 0; i < f1peaks.length; i++ ) {
        //     var pk = f1peaks.at(i);
        //     f1ppm_list.push(pk.delta(1));
        // }
     
        // create an array of f1ppm and f2ppm coordinates
        var possiblePeaksData = [];
        for( var i = 0; i < f1ppm_list.length; i++ ) {
            for( var j = 0; j < f2ppm_list.length; j++ ) {
                possiblePeaksData.push([f2ppm_list[j], f1ppm_list[i]]);
            }
        }

        print("Number of possible peaks ", possiblePeaksData.length);
        print();
    
        var nmr2Dspectrum = spec2D;
        var sLimits = spec2D.getScaleLimits();
    
        var twoDppeaksfound = nmr2Dspectrum.peaks();
        var twoDintegalsfound = nmr2Dspectrum.integrals();
        print("the number of peaks found in the 2D spectrum ", twoDppeaksfound.count);
        print("The number of integrals found in the 2D spectrum ", twoDintegalsfound.count);
        print();


    
         // findout the range of the peak intensiities in order to determine the cutoff height
         // create a peaks region
        var peaksRegion = SpectrumRegion( sLimits.fromX, sLimits.toX, sLimits.fromY, sLimits.toY);
        var allPeaksFound = Peaks(nmr2Dspectrum, peaksRegion);
        var allIntegralsFound = Integral(nmr2Dspectrum, peaksRegion);
    
        print("The numbr of integrals found ", allIntegralsFound.count);
        print("The numbr of peaks found ", allPeaksFound.count);
        print();

        // print out the peaks in json format using stringify
        print("All peaks found in the region in json format \n", JSON.stringify(allPeaksFound));
        print();
        print(Integral.calculationParams);



                        // add peaks to the 2D spectrum
        var newPeaks = new Peaks();
        for( var i = 0; i < allPeaksFound.count; i++ ) {
            var pk = allPeaksFound.at(i);
            newPeaks.append(pk);
        }
        nmr2Dspectrum.setPeaks(newPeaks);
        nmr2Dspectrum.process();

        return;
    
         // extract delta(1) and delta(2) values from the peaks
        var allPeaksData = [];
        for( var i = 0; i < allPeaksFound.count; i++ ) {
            var pk = allPeaksFound.at(i);
            allPeaksData.push([pk.delta(2), pk.delta(1)]);
        }

            print("All peaks data \n", allPeaksData);
            print();
    
         // find the possible peaks enclosed in the SLimits
        var possiblePeaksInRegion = [];
        for( var i = 0; i < possiblePeaksData.length; i++ ) {
            var pk = possiblePeaksData[i];
            if( (pk[0] >= sLimits.fromX) && (pk[0] <= sLimits.toX) && (pk[1] >= sLimits.fromY )&& (pk[1] <= sLimits.toY) ) {
                possiblePeaksInRegion.push(pk);
            }
        }

        print("Possible peaks in region \n", possiblePeaksInRegion);
        print();
        return;
         var result = kMeans2D(allPeaksData, possiblePeaksInRegion);
    
         // create a list of peaks found from the cluster results
        for (var i = 0; i < result.clusters.length; i++) {
            if( result.clusters[i].length > 0 ) {
                // peak found
                var pk = new Peak( possiblePeaksInRegion[i][1], possiblePeaksInRegion[i][0], nmr2Dspectrum);
                twoDppeaksfound.append(pk);
            }
        }
    
        // if the spectrum is a clipcosyhsqc only keep the peaks that have a negative intensity
        if( simpleutils.isHSQCCLIPCOSY(nmr2Dspectrum) ) {
            var newPeaks = new Peaks();
            for( var i = 0; i < twoDppeaksfound.count; i++ ) {
                var pk = twoDppeaksfound.at(i);
                if( pk.intensity < 0.0 ) {
                    newPeaks.append(pk);
                }
            }
            twoDppeaksfound = newPeaks;
        }
        else if( simpleutils.isCOSY(nmr2Dspectrum) ) {
            var newPeaks = new Peaks();
            for( var i = 0; i < twoDppeaksfound.count; i++ ) {
                var pk = twoDppeaksfound.at(i);
                if( pk.delta(1) != pk.delta(2) ) {
                    newPeaks.append(pk);
                }
            }
            twoDppeaksfound = newPeaks;
        }
    
        print("Number of peaks found: ", twoDppeaksfound.count);
    
        // add integrals to the peaks
    
    
        // add peaks to the 2D spectrum
        nmr2Dspectrum.setPeaks(twoDppeaksfound);
        nmr2Dspectrum.process();
    
    }

    // Function to calculate the Euclidean distance between two points in 2D
    function euclideanDistance2D(point1, point2) {
        var sum = 0;
        for (var i = 0; i < point1.length; i++) {
            sum += Math.pow(point1[i] - point2[i], 2);
        }
        return Math.sqrt(sum);
    }

    // Function to assign each data point to the nearest centroid
    function assignClusters2D(data, centroids) {
        var clusters = [];
        for (var i = 0; i < centroids.length; i++) {
            clusters.push([]);
        }

        for (var i = 0; i < data.length; i++) {
            var point = data[i];
            var distances = [];
            for (var j = 0; j < centroids.length; j++) {
                distances.push(euclideanDistance2D(point, centroids[j]));
            }
            var closestCentroidIndex = distances.indexOf(Math.min.apply(Math, distances));
            clusters[closestCentroidIndex].push(point);
        }
        return clusters;
    }

    // Function to update centroids by calculating the mean of each cluster
    function updateCentroids2D(clusters) {
        var newCentroids = [];
        for (var i = 0; i < clusters.length; i++) {
            var cluster = clusters[i];
            if (cluster.length === 0) {
                newCentroids.push([0, 0]); // Placeholder for empty clusters
                continue;
            }
            var centroid = [0, 0];
            for (var j = 0; j < cluster.length; j++) {
                var point = cluster[j];
                centroid[0] += point[0];
                centroid[1] += point[1];
            }
            centroid[0] /= cluster.length;
            centroid[1] /= cluster.length;
            newCentroids.push(centroid);
        }
        return newCentroids;
    }

    // K-Means clustering algorithm with predefined initial centroids
    function kMeans2D(data, initialCentroids, maxIterations) {
        maxIterations = maxIterations || 100;
        var centroids = initialCentroids.slice(); // Copy the initial centroids to avoid modifying the input array
        var clusters;

        for (var i = 0; i < maxIterations; i++) {
            clusters = assignClusters2D(data, centroids);
            var newCentroids = updateCentroids2D(clusters);

            // Check for convergence (if centroids do not change)
            if (JSON.stringify(newCentroids) === JSON.stringify(centroids)) {
                break;
            }
            centroids = newCentroids;
        }
        return { clusters: clusters, centroids: centroids };
    }

    // iterate through pages and print out what is on the page
    var doc = Application.mainWindow.activeDocument;
 
    // var pageCount = doc.pageCount();

    // MessageBox.information("Number of pages in the document: " + pageCount);

    // identify the spectra in the document
    var spectra_lst = identify_spectrum();

    print("Number of spectra in the document: " + spectra_lst.length);

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
    

    // MessageBox.information("Carbon Spectrum Found");

    var pshift = simpleutils.findPureShiftSpectrum(spectra_lst);
    if(pshift == undefined && hsqcspec == undefined){
        MessageBox.warning("No Proton peaks defined");
        return
    }

    var protonPeaksList = [];
    if(pshift != undefined && pshift.peaks().count > 0){
        if (pshift.peaks().count > 0 ) {
            // creat protons peaks list from the pshift spectrum
            var protonsPeaks = pshift.peaks();

            for (var i = 0; i < protonsPeaks.count; i++) {
                var pk = protonsPeaks.at(i);
                protonPeaksList.push(pk.delta(1));
            }
        }
    }
    else if (hsqcspec != undefined && hsqcspec.peaks().count > 0) {
        // creat protons peaks list from the hsqc spectrum
        var protonsPeaks = hsqcspec.peaks();

        for (var i = 0; i < protonsPeaks.count; i++) {
            var pk = protonsPeaks.at(i);
            protonPeaksList.push(pk.delta(2));
        }
    }
    else {
        MessageBox.warning("No Proton peaks defined");
        return
    }

    // MessageBox.information("Pureshift Spectrum Found");

    var nmrspec = new NMRSpectrum(nmr.activeSpectrum());
    // check if the spectrum is a 2D spectrum
    if(!simpleutils.is2DSpectrum(nmrspec)){
        MessageBox.warning("Please select a 2D spectrum");
        return;
    }

    if(isSubstring("13C", nmrspec.getParam("Nucleus"))){
        find2DPeaks( protonPeaksList, c13peaksList, nmrspec);
    }
    else {
        find2DPeaks( protonPeaksList, protonPeaksList, nmrspec);
    }
    // MessageBox.information("Find 2D Peaks Function finished");
}
 
 



 if (this.MnUi && MnUi.simpleNMRtools) {
	MnUi.simpleNMRtools.simpleNMRtools_simplePeakPick_eeh = simplePeakPick_eeh;
}





 
