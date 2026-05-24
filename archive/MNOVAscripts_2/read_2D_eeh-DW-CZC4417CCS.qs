function read_2D_eeh() {

    /**
     * Gaussian filter implementation for 2D arrays
     */
    function gaussianFilter(matrix, sigma) {
        if (sigma === 0) return matrix;
        
        var rows = matrix.length;
        var cols = matrix[0].length;
        
        var kernelSize = Math.max(3, Math.ceil(6 * sigma) | 1);
        var kernelRadius = Math.floor(kernelSize / 2);
        
        var kernel = [];
        var sum = 0;
        for (var i = 0; i < kernelSize; i++) {
            var x = i - kernelRadius;
            var value = Math.exp(-(x * x) / (2 * sigma * sigma));
            kernel[i] = value;
            sum += value;
        }
        
        for (var i = 0; i < kernelSize; i++) {
            kernel[i] /= sum;
        }
        
        var tempResult = [];
        for (var r = 0; r < rows; r++) {
            tempResult[r] = [];
            for (var c = 0; c < cols; c++) {
                var value = 0;
                for (var k = 0; k < kernelSize; k++) {
                    var sourceCol = c + k - kernelRadius;
                    sourceCol = Math.max(0, Math.min(cols - 1, sourceCol));
                    value += matrix[r][sourceCol] * kernel[k];
                }
                tempResult[r][c] = value;
            }
        }
        
        var result = [];
        for (var r = 0; r < rows; r++) {
            result[r] = [];
            for (var c = 0; c < cols; c++) {
                var value = 0;
                for (var k = 0; k < kernelSize; k++) {
                    var sourceRow = r + k - kernelRadius;
                    sourceRow = Math.max(0, Math.min(rows - 1, sourceRow));
                    value += tempResult[sourceRow][c] * kernel[k];
                }
                result[r][c] = value;
            }
        }
        
        return result;
    }

    /**
     * Maximum filter implementation
     */
    function maximumFilter(matrix, size) {
        var rows = matrix.length;
        var cols = matrix[0].length;
        var result = [];
        var radius = Math.floor(size / 2);
        
        for (var r = 0; r < rows; r++) {
            result[r] = [];
            for (var c = 0; c < cols; c++) {
                var maxValue = -Infinity;
                
                for (var dr = -radius; dr <= radius; dr++) {
                    for (var dc = -radius; dc <= radius; dc++) {
                        var nr = r + dr;
                        var nc = c + dc;
                        
                        if (nr >= 0 && nr < rows && nc >= 0 && nc < cols) {
                            maxValue = Math.max(maxValue, matrix[nr][nc]);
                        }
                    }
                }
                
                result[r][c] = maxValue;
            }
        }
        
        return result;
    }

    function abs2D(matrix) {
        var result = [];
        for (var r = 0; r < matrix.length; r++) {
            result[r] = [];
            for (var c = 0; c < matrix[r].length; c++) {
                result[r][c] = Math.abs(matrix[r][c]);
            }
        }
        return result;
    }

    function findCoordinates(mask) {
        var coords = [];
        for (var r = 0; r < mask.length; r++) {
            for (var c = 0; c < mask[r].length; c++) {
                if (mask[r][c]) {
                    coords.push([r, c]);
                }
            }
        }
        return coords;
    }

    function extractValues(matrix, coords) {
        var values = [];
        for (var i = 0; i < coords.length; i++) {
            var r = coords[i][0];
            var c = coords[i][1];
            values.push(matrix[r][c]);
        }
        return values;
    }

    function logicalAnd(mask1, mask2) {
        var result = [];
        for (var r = 0; r < mask1.length; r++) {
            result[r] = [];
            for (var c = 0; c < mask1[r].length; c++) {
                result[r][c] = mask1[r][c] && mask2[r][c];
            }
        }
        return result;
    }

    function compareEqual(matrix1, matrix2) {
        var result = [];
        for (var r = 0; r < matrix1.length; r++) {
            result[r] = [];
            for (var c = 0; c < matrix1[r].length; c++) {
                result[r][c] = matrix1[r][c] === matrix2[r][c];
            }
        }
        return result;
    }

    function greaterThan(matrix, threshold) {
        var result = [];
        for (var r = 0; r < matrix.length; r++) {
            result[r] = [];
            for (var c = 0; c < matrix[r].length; c++) {
                result[r][c] = matrix[r][c] > threshold;
            }
        }
        return result;
    }

    // Standard scaler function
    function standardScale(data) {
        if (data.length === 0) return data;
        
        var scaled = [];
        var means = [];
        var stds = [];
        
        // Calculate means for each dimension
        for (var dim = 0; dim < data[0].length; dim++) {
            var sum = 0;
            for (var i = 0; i < data.length; i++) {
                sum += data[i][dim];
            }
            means[dim] = sum / data.length;
        }
        
        // Calculate standard deviations
        for (var dim = 0; dim < data[0].length; dim++) {
            var sumSquaredDiffs = 0;
            for (var i = 0; i < data.length; i++) {
                sumSquaredDiffs += Math.pow(data[i][dim] - means[dim], 2);
            }
            stds[dim] = Math.sqrt(sumSquaredDiffs / data.length);
            // Avoid division by zero
            if (stds[dim] === 0) stds[dim] = 1;
        }
        
        // Scale the data
        for (var i = 0; i < data.length; i++) {
            var scaledPoint = [];
            for (var dim = 0; dim < data[i].length; dim++) {
                scaledPoint[dim] = (data[i][dim] - means[dim]) / stds[dim];
            }
            scaled.push(scaledPoint);
        }
        
        return scaled;
    }

    // Calculate cluster centers (centroids)
    function calculateClusterCenters(data, result) {
        var centers = {};
        var clusterKeys = Object.keys(result.clusters);
        
        for (var i = 0; i < clusterKeys.length; i++) {
            var clusterId = clusterKeys[i];
            var cluster = result.clusters[clusterId];
            var center = [0, 0];
            
            // Sum all points in the cluster
            for (var j = 0; j < cluster.length; j++) {
                var point = data[cluster[j].index];
                center[0] += point[0];
                center[1] += point[1];
            }
            
            // Calculate mean (centroid)
            center[0] /= cluster.length;
            center[1] /= cluster.length;
            
            centers[clusterId] = center;
        }
        
        return centers;
    }

    /**
     * Peak detection function with enhanced debugging
     */
    function findPeakCandidatesReal(spectrum,  useLocalMaxima, minDistance) {
        if (typeof useLocalMaxima === 'undefined') useLocalMaxima = true;
        if (typeof minDistance === 'undefined') minDistance = 3;
        
        print('Running peak detection with minDistance =', minDistance);
        
        var coords = [];
        var intensities = [];
        
        if (useLocalMaxima) {
            var absSpectrum = abs2D(spectrum);
            var smoothed = gaussianFilter(absSpectrum, 1);
            
            // Find all local maxima first (without distance filtering)
            var maxFiltered = maximumFilter(smoothed, 3); // Use fixed small size first
            var localMaxima = compareEqual(maxFiltered, smoothed);
            var nonZeroMask = greaterThan(smoothed, 0);
            var peakMask = logicalAnd(localMaxima, nonZeroMask);
            
            var allCoords = findCoordinates(peakMask);
            var allIntensities = extractValues(smoothed, allCoords);
            
            print('Found', allCoords.length, 'initial local maxima');
            
            // Sort peaks by intensity (descending) to prioritize stronger peaks
            var peakData = [];
            for (var i = 0; i < allCoords.length; i++) {
                peakData.push({
                    coord: allCoords[i],
                    intensity: allIntensities[i]
                });
            }
            peakData.sort(function(a, b) { return b.intensity - a.intensity; });
            
            // Apply distance filtering
            var finalCoords = [];
            var finalIntensities = [];
            
            for (var i = 0; i < peakData.length; i++) {
                var currentPeak = peakData[i];
                var tooClose = false;
                
                // Check distance to all already accepted peaks
                for (var j = 0; j < finalCoords.length; j++) {
                    var distance = Math.sqrt(
                        Math.pow(currentPeak.coord[0] - finalCoords[j][0], 2) + 
                        Math.pow(currentPeak.coord[1] - finalCoords[j][1], 2)
                    );
                    
                    if (distance < minDistance) {
                        tooClose = true;
                        break;
                    }
                }
                
                if (!tooClose) {
                    finalCoords.push(currentPeak.coord);
                    finalIntensities.push(currentPeak.intensity);
                }
            }
            
            coords = finalCoords;
            intensities = finalIntensities;

            // threshold the peaks using 50% of maximum peak found
            var maxIntensity = -Infinity;
            for (var i = 0; i < intensities.length; i++) {
                if (intensities[i] > maxIntensity) {
                    maxIntensity = intensities[i];
                }
            }
            var intensityThreshold = maxIntensity * 0.5;
            var thresholdedCoords = [];
            var thresholdedIntensities = [];

            for (var i = 0; i < coords.length; i++) {
                if (intensities[i] >= intensityThreshold) {
                    thresholdedCoords.push(coords[i]);
                    thresholdedIntensities.push(intensities[i]);
                }
            }

            coords = thresholdedCoords;
            intensities = thresholdedIntensities;

            print('After distance filtering with minDistance=' + minDistance + ':', coords.length, 'peaks remain');

        } else {
            var absSpectrum = abs2D(spectrum);
            var mask = greaterThan(absSpectrum, 0);
            
            coords = findCoordinates(mask);
            intensities = extractValues(absSpectrum, coords);
            
            print('Found', coords.length, 'non-zero points in thresholded spectrum');
        }
        
        return {
            coords: coords,
            intensities: intensities,
            smoothedSpec: smoothed
        };
    }

    var spec = nmr.activeSpectrum();
  
  	if (!spec.isValid()) {
		return;
	}
	
	
	var zoomLimits = spec.getScaleLimits();

    print("\nzoomLimits\n");
    print(zoomLimits);
    print();
    	
	
	print("dimCount " + spec.dimCount );	
	print("spec.count " + spec.count(1) + " " + spec.count(2) );
	print("spec Threshold " + spec.threshold );
	print(" spec realMax " + spec.realMax );
	print("spec realMin " + spec.realMin );
    print();

	// print out the the sLimits
	//  sLimits.fromX, sLimits.toX, sLimits.fromY, sLimits.toY
	
    var fullLimits = spec.getFullScaleLimits();

    var f2_full_ppm_from = fullLimits.fromX;
    var f2_full_ppm_to = fullLimits.toX;
    var f1_full_ppm_from = fullLimits.fromY;
    var f1_full_ppm_to = fullLimits.toY;


    // convert the limits from ppm to pts
    // unitsToUnits (String aFromUnits, String aToUnits, Number aValue, Number aDim) : Number
    var ptsfromX = spec.unitsToUnits("ppm", "pt", zoomLimits.fromX, 2);
    var ptsfromY = spec.unitsToUnits("ppm", "pt", zoomLimits.fromY, 1);
    var ptstoX = spec.unitsToUnits("ppm", "pt", zoomLimits.toX, 2);
    var ptstoY = spec.unitsToUnits("ppm", "pt", zoomLimits.toY, 1);

    var f2_zoom_pts_from = ptsfromX;
    var f2_zoom_pts_to = ptstoX;
    var f1_zoom_pts_from = ptsfromY;
    var f1_zoom_pts_to = ptstoY;

    var f2_zoom_ppm_from = zoomLimits.fromX;
    var f2_zoom_ppm_to = zoomLimits.toX;
    var f1_zoom_ppm_from = zoomLimits.fromY;
    var f1_zoom_ppm_to = zoomLimits.toY;


	// var rrr = spec.real({from:ptstoY, to:ptsfromY}, {from:ptstoX,to:ptsfromX});
    var rrr = spec.real({from:f1_zoom_pts_to, to:f1_zoom_pts_from}, 
                        {from:f2_zoom_pts_to, to:f2_zoom_pts_from});

    print("rrr.length " + rrr.length);
    print("rrr[0].length " + rrr[0].length);
    print();

    // convert rrr to json
    var rrr_json = JSON.stringify(rrr);
    // print("rrr_json.length " + rrr_json.length);
    // print("rrr_json " + rrr_json);



    results = findPeakCandidatesReal(rrr);

    // print results
    print("Found", results.coords.length, "peaks");
    print("Peak intensities:\n", results.intensities);
    print("Peak coordinates:\n", results.coords);
    print();

    for (var i = 0; i < results.coords.length; i++) {
        print("Peak " + (i + 1) + " coordinates: " + results.coords[i] + " intensity: " + results.intensities[i]);
    }

    var smoothed_json = JSON.stringify(results.smoothedSpec);
    var coords_json = JSON.stringify(results.coords);

    // save "exptIdentifiers" as a json file
    var fout = new File("/Users/vsmw51/OneDrive - Durham University/projects/programming/2025/python/awh/simpleNMRtools2/MNOVAscripts_2/rrr.json");
    if (fout.open(File.WriteOnly)) {
        sout = new TextStream(fout, 'UTF-8');
        sout.writeln(smoothed_json);
    }
    fout.close();

    // saved coords_json to coords.json
    var fout_coords = new File("/Users/vsmw51/OneDrive - Durham University/projects/programming/2025/python/awh/simpleNMRtools2/MNOVAscripts_2/coords.json");
    if (fout_coords.open(File.WriteOnly)) {
        sout = new TextStream(fout_coords, 'UTF-8');
        sout.writeln(coords_json);
    }
    fout_coords.close();

    // Example usage:
    var dbscan = new DBSCAN(0.5, 1); // eps=1.0, minPts=3



    // Scale the data for DBSCAN
    scaledData = standardScale(results.coords);
    var result = dbscan.fitPredict(scaledData);

    var clusterCenters = calculateClusterCenters(results.coords, result);

    print("Num clusters: " + Object.keys(result.clusters).length);
    print("Num noise points: " + result.noise.length);
    print("Num labels: " + result.labels.length);

    print('Clusters:\n', result.clusters);
    print('Noise points:\n', result.noise);
    print('Labels:\n', result.labels);

    print('Cluster centers:\n', clusterCenters);

    var clusterCenters_json = JSON.stringify(clusterCenters);



    // saved clusterCenters_json to clusterCenters.json
    var fout_centers = new File("/Users/vsmw51/OneDrive - Durham University/projects/programming/2025/python/awh/simpleNMRtools2/MNOVAscripts_2/clusterCenters.json");
    if (fout_centers.open(File.WriteOnly)) {
        sout = new TextStream(fout_centers, 'UTF-8');
        sout.writeln(clusterCenters_json);
    }
    fout_centers.close();

    // convert clusterCenters to ppm
    var clusterCenters_ppm = {};
    for (var k in clusterCenters) {
        clusterCenters_ppm[k] = [
            spec.unitsToUnits("pt", "ppm", clusterCenters[k][0], 1),
            spec.unitsToUnits("pt", "ppm", clusterCenters[k][1], 2)
        ];
    }

    print('Cluster centers (ppm):\n', clusterCenters_ppm);

    // get the peaks list and append the cluster centers
    var peaks = spec.peaks();
    if (!peaks) {
        peaks = new Peaks();
    }
    print("peaks.count " + peaks.count);
    print("clusterCenters_ppm.count " + clusterCenters_ppm.length);
    for( var k in clusterCenters_ppm ) {
        print("k: " + k);
        var p = new Peak(clusterCenters_ppm[k][1], clusterCenters_ppm[k][0], spec);//single peak
        peaks.append(p);
    }
    spec.setPeaks(peaks);
    spec.process();
}