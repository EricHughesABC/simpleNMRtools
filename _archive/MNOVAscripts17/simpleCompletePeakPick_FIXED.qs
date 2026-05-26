


/**
 * Finds and returns the first 1D 13C spectrum object from a list of spectra.
 *
 * @param {Array<Object>} spectra_lst - List of spectrum objects to search.
 * @param {string} experiment - The experiment type to search for.
 * @returns {Object|undefined} The first spectrum object with "experimentEEH" equal to the specified experiment, or undefined if not found.
 */
function find_NMR_spectrum( spectra_lst, experiment ) {
    var nmr_spectrum = undefined;

    // search field "experimentEEH" for the specified experiment
    for (var i=0; i<spectra_lst.length; i++){
        if (spectra_lst[i]["experimentEEH"] == experiment) {
            nmr_spectrum = spectra_lst[i];
            break;
        }
    }
    return nmr_spectrum;
}


function getZoomed2Dspec_rrr(spec) {

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

    // print("spec");
    // print(spec);

    return {"rrr": rrr,
            "num_F2_pts": rrr.length,
            "num_F1_pts": rrr[0].length,
            "F2": [f2_zoom_ppm_to, f2_zoom_ppm_from],
            "F1": [f1_zoom_ppm_to, f1_zoom_ppm_from],
            "F2_pts": [f2_zoom_pts_to, f2_zoom_pts_from],
            "F1_pts": [f1_zoom_pts_to, f1_zoom_pts_from]
    };
}

function find_max_min_in_matrix(matrix) {

    // find the max and min in the zoomed region zoomed_region.rrr
    var max_val_rrr = -Infinity;
    var min_val_rrr = Infinity;

    for (var i=0; i<matrix.length; i++){
        for (var j=0; j<matrix[i].length; j++){
            if (matrix[i][j] > max_val_rrr){
                max_val_rrr  = matrix[i][j];
            }
            if (matrix[i][j] < min_val_rrr){
                min_val_rrr = matrix[i][j];
            }
        }
    }
    return { "max": max_val_rrr, "min": min_val_rrr };  
}

function determine_peaks_pos_and_neg( max_val, min_val ) {

    // check if we have both positive and negative peaks or just positive or just negative
    var pos_only = false;
    var neg_only = false;
    var both_pos_neg = false;

    var pos_max = max_val;
    var neg_max = Math.abs(min_val);
    if (pos_max > neg_max)
        var ratio = pos_max / neg_max;
    else
        var ratio = neg_max / pos_max;

    if( ratio > 10 && pos_max > neg_max){
        pos_only = true;
    }
    else if( ratio > 10 && neg_max > pos_max){
        neg_only = true;
    }
    else{
        both_pos_neg = true;
    }

    return { "pos_only": pos_only, "neg_only": neg_only, "both_pos_neg": both_pos_neg };
}


function find_pks( matrix, peak_types, max_val, min_val ) {
    var pk_coords = [];

    var threshold_neg = min_val - 0.9 * min_val;
    var threshold_pos = max_val - 0.9 * max_val;

    var pos_only = peak_types.pos_only;
    var neg_only = peak_types.neg_only;
    var both_pos_neg = peak_types.both_pos_neg;

    var nrows = matrix.length;
    var ncols = matrix[0].length;

    for (var i=0; i<matrix.length; i++){
        for (var j=0; j<matrix[i].length; j++){
            var coord_index = i * ncols + j;

            if (pos_only){
                if (matrix[i][j] > threshold_pos){
                    pk_coords.push( [i,j, coord_index, matrix[i][j]] );       
                }
            }
            else if (neg_only){
                if (matrix[i][j] < threshold_neg){
                    pk_coords.push( [i,j, coord_index, matrix[i][j]] );       
                }
            }
            else if (both_pos_neg){
                if (matrix[i][j] > threshold_pos || matrix[i][j] < threshold_neg){
                    pk_coords.push( [i,j, coord_index, matrix[i][j]] );       
                }
            }
        }
    }
    return pk_coords;
} 



function create_mask_index(index_val, nrows, ncols) {
    var mask = [];

    for (var i=-1; i<=1; i++){
        for (var j=-ncols; j<=ncols; j+=ncols){
            var new_index = index_val + j + i;
            mask.push(new_index);
        }
    }
    return mask;
}

function group_peaks(peak_info, nrows, ncols) {
    var peaks_found = [];
    if (peak_info.length == 0 || peak_info.length == 1) {
        return peaks_found;
    }

    peaks_found.push( [ peak_info[0][2] ] );

    for (var p=1; p<peak_info.length; p++){
        var signal = peak_info[p];
        var new_pk_found = true;
        var index = signal[2];

        print("Processing index: " + index);

        for (var g=0; g<peaks_found.length; g++){
            var peak_group = peaks_found[g];
            // CRITICAL FIX: Capture the initial length to avoid infinite loop
            // when we modify the array during iteration
            var initial_group_length = peak_group.length;
            print(" Checking against peak group " + g + ": " + peak_group);
            
            // check if the index is within the mask of any existing peak group
            for (var k=0; k<initial_group_length; k++){
                var pk = peak_group[k];
                var mask_indices = create_mask_index(pk, nrows, ncols);
                print("  Checking peak " + pk);
                
                // Check if index is in this peak's mask
                var found_in_mask = false;
                for (var m=0; m<mask_indices.length; m++){
                    if ( mask_indices[m] === index ){
                        found_in_mask = true;
                        break;
                    }
                }
                
                // If found, add to group and stop checking this group
                if (found_in_mask){
                    print("   -> Match! Adding index " + index + " to group " + g);
                    peaks_found[g].push(index);
                    new_pk_found = false;
                    break;  // Break from peak_group loop - don't check other peaks in this group
                }
            }
            
            // If we found a match in this group, don't check other groups
            if (!new_pk_found){
                break;
            }
        }
        
        if (new_pk_found){
            print("   -> New peak group created for index " + index);
            peaks_found.push( [ index ] );
        }
        
    }

    // remove any peaks groups that only have one point
    var filtered_peaks_found = [];
    for (var m=0; m<peaks_found.length; m++){
        if (peaks_found[m].length > 1){
            filtered_peaks_found.push(peaks_found[m]);
        }
    }

    return filtered_peaks_found;
}   



function simpleCompletePeakPick() {

    // get simpleutils functions
    var simpleutils = new simpleUtils();



    // //  get the active 2D spectrum
    // var nmr2Dspec = new NMRSpectrum(nmr.activeSpectrum());

    var doc = Application.mainWindow.activeDocument;
    var simpleutils = new simpleUtils();
    var spectra_lst = simpleutils.identifySpectra(doc);
    spectra_lst = identify_spectraType(spectra_lst);

    var nmr2Dspec = new NMRSpectrum(nmr.activeSpectrum());

    // if not 2D spectrum, exit with message
    if(!simpleutils.is2DSpectrum(nmr2Dspec)){
        return;
    }
    
    // find all the spectra in the document
    // identify the spectra in the document
    var currentSpectrum = identify_spectraType(nmr2Dspec);

    // print out the active spectrum type
    print("\nActive spectrum type: " + currentSpectrum["experimentEEH"]);
    print();

    const expt_types = ["H1_1D", "C13_1D", "PURESHIFT", "DEPT135", "HSQC", "HMBC", "COSY", "NOESY", "DOUBLEDEPTCH3", "HSQCCLIPCOSY"];

    var spectra_dict = {};

    for( var i=0; i<expt_types.length; i++){
        var spectrum = find_NMR_spectrum(spectra_lst, expt_types[i]);
        if (spectrum) {
          spectra_dict[expt_types[i]] = spectrum;
        }
    }

    // print(spectra_dict);

    print("\nSpectra summary:\n");
    for ( var key in spectra_dict){
      print("\t" + key + ": has peaks? " + spectra_dict[key]["hasPeaks"] + ", has multiplets? " + spectra_dict[key]["hasMultiplets"] + ", has integrals? " + spectra_dict[key]["hasIntegrals"]);
    }

    // dictionary of spectra with peaks
    var spectra_with_peaks = {};
    for (var key in spectra_dict) {
        if (spectra_dict[key]["hasPeaks"]) {
            spectra_with_peaks[key] = spectra_dict[key];
        }
    }

    // print(spectra_with_peaks);
    print("\nSpectra with peaks:\n");
    for( var key in spectra_with_peaks){
      print("\t" + key);
    }

    print( spectra_with_peaks.C13_1D.hasPeaks );

    // var peakPickingPath = 0;

    // if (spectra_with_peaks.C13_1D.hasPeaks) peakPickingPath |= 1;
    // if (spectra_with_peaks.PURESHIFT.hasPeaks) peakPickingPath |= 2;

    // if (spectra_with_peaks.HSQC && spectra_with_peaks.HSQC.hasPeaks) peakPickingPath |= 4;

    // print( "peakPickingPath: " + peakPickingPath );

    // switch(peakPickingPath){

    //     case 1:
    //         print("Using 1D 13C spectrum only for peak picking.");
    //         break;
    //     case 2:
    //         print("Using 1D PureShift spectrum only for peak picking.");
    //         break;
    //     case 3:
    //         print("Using both 1D 13C and PureShift spectra for peak picking.");
    //         break;

    //     case 4:
    //         print("Using 2D HSQC spectrum only for peak picking.");
    //         break;

    //     case 5:
    //         print("Using 1D 13C and 2D HSQC spectra for peak picking.");
    //         break;
    // // default:
    //     default:
    //         print("No case matched for peakPickingPath: " + peakPickingPath);
    //         break;
    // }


    // get active spectrum
    var activeSpectrum = new NMRSpectrum(nmr.activeSpectrum());

    // check if active spectrum is 2D
    if(!simpleutils.is2DSpectrum(activeSpectrum)){
        print("Active spectrum is not 2D. Please select a 2D spectrum and try again.");
        return;
    }
    else{
        print("Active spectrum is 2D. Continuing....")
    }

    // get the coordinates of the 2D spectrum and print

    var sLimits = activeSpectrum.getScaleLimits();

    print(sLimits);

    var zoomed_region = getZoomed2Dspec_rrr(activeSpectrum);

    print(zoomed_region.num_F2_pts + ", " + zoomed_region.num_F1_pts);

    print("F2: " + zoomed_region.F2);
    print("F1: " + zoomed_region.F1);

    print();

    print("F2 pts: " + zoomed_region.F2_pts);
    print("F1 pts: " + zoomed_region.F1_pts);

    var max_min_values = find_max_min_in_matrix(zoomed_region.rrr);
    var max_val_rrr = max_min_values.max;
    var min_val_rrr = max_min_values.min;

    // check if we have both positive and negative peaks or just positive or just negative
    var peak_types = determine_peaks_pos_and_neg( max_val_rrr, min_val_rrr );

    // find the peak coordinates
    var pk_coords = find_pks( zoomed_region.rrr, peak_types, max_val_rrr, min_val_rrr );


    print("Number of peak coordinates found: " + pk_coords.length );
    print();
    for (var i=0; i < pk_coords.length; i++){
        print( "\t" + pk_coords[i] );
    }

    var peaks_found = group_peaks( pk_coords, zoomed_region.num_F2_pts, zoomed_region.num_F1_pts );

    print("\nNumber of peak groups found: " + peaks_found.length );
    print();
    for (var i=0; i < peaks_found.length; i++){
        print( "\tPeak group " + (i+1) + ": " + peaks_found[i] );
    }
    // create info to send to flask server

    // # Prepare the data payload
    // payload = {
    //     "source": "MNOVA",
    //     "destination": "MNOVA",
    //     "type": "HSQC",
    //     "matrix": dummy_spectrum,
    //     "axis_info": {
    //         "F2": {
    //             "start": 10.0,   # 1H dimension starting at 10 ppm
    //             "end": 6.0,      # 1H dimension ending at 6 ppm
    //             "label": "1H"
    //         },
    //         "F1": {
    //             "start": 130.0,  # 15N dimension starting at 130 ppm
    //             "end": 100.0,    # 15N dimension ending at 100 ppm
    //             "label": "15N"
    //         }
    //     }
    // }

    var payload = {

        "source": "MNOVA",
        "destination": "MNOVA",
        "type": "HSQC",
        "matrix": zoomed_region.rrr,
        "axis_info": {
            "F2": {
                "start": zoomed_region.F2[0],
                "end": zoomed_region.F2[1],
                "label": "F2"
            },
            "F1": {
                "start": zoomed_region.F1[0],
                "end": zoomed_region.F1[1],
                "label": "F1"
            }
        }
    }

    // var url = "http://localhost:5000/process_spectrum";
    // var web_utils = new WebUtilsQT();
    // var json_obj_string = JSON.stringify(payload,null,4);
    // rtn = web_utils.jsonRequest(url, json_obj_string, "", "", false);

    // print(rtn);
    
}


