


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

    var iii = spec.imag({from:f1_zoom_pts_to, to:f1_zoom_pts_from}, 
                        {from:f2_zoom_pts_to, to:f2_zoom_pts_from});

    // calc magnitude of spectrum using real and imaginary values ie rrr and iii

    var mmm = [];

    //  to be completed


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
            "iii": iii,
            "mmm": mmm,
            "num_F2_pts": rrr.length,
            "num_F1_pts": rrr[0].length,
            "F2": [f2_zoom_ppm_to, f2_zoom_ppm_from],
            "F1": [f1_zoom_ppm_to, f1_zoom_ppm_from],
            "F2_pts": [f2_zoom_pts_to, f2_zoom_pts_from],
            "F1_pts": [f1_zoom_pts_to, f1_zoom_pts_from],
            'F2_ppm_delta': Math.abs(f2_zoom_ppm_to - f2_zoom_ppm_from) / rrr[0].length,
            'F1_ppm_delta': Math.abs(f1_zoom_ppm_to - f1_zoom_ppm_from) / rrr.length
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

    var threshold_neg = min_val - 0.85 * min_val;
    var threshold_pos = max_val - 0.85 * max_val;

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
                    pk_coords.push( {"F1": i, 
                                     "F2": j, 
                                     "index": coord_index, 
                                     "intensity": matrix[i][j],
                                     "pk_id": -1} );       
                }
            }
            else if (neg_only){
                if (matrix[i][j] < threshold_neg){
                    pk_coords.push( {"F1": i, 
                                     "F2": j, 
                                     "index": coord_index, 
                                     "intensity": matrix[i][j],
                                     "pk_id": -1} );       
                }
            }
            else if (both_pos_neg){
                if (matrix[i][j] > threshold_pos || matrix[i][j] < threshold_neg){
                    pk_coords.push( {"F1": i, 
                                     "F2": j, 
                                     "index": coord_index, 
                                     "intensity": matrix[i][j],
                                     "pk_id": -1} );       
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




// 

function group_pts(peak_info, nrows, ncols) {
    var peaks_found = [];
    if (peak_info.length == 0 || peak_info.length == 1) {
        return peaks_found;
    }

    // Initialize first group as an object with nested properties
    var first_group = {};
    var first_index = peak_info[0].index;
    first_group[first_index] = {
        'peak_height': peak_info[0]["intensity"],
        'grp_id': 1,
        'F1_ppm': peak_info[0].F1_ppm,
        'F2_ppm': peak_info[0].F2_ppm
    };
    peaks_found.push(first_group);

    var current_grp_id = 1;

    for (var p=1; p<peak_info.length; p++){
        var signal = peak_info[p];
        var new_pk_found = true;
        var index = signal.index;
        var peak_height = signal["intensity"];

        for (var g=0; g<peaks_found.length; g++){
            var peak_group = peaks_found[g];
            
            // Get all current indices in this group as an array
            var group_indices = Object.keys(peak_group);
            
            // Get the grp_id for this group (from first entry)
            var grp_id = peak_group[group_indices[0]].grp_id;
            
            // check if the index is within the mask of any existing peak in the group
            for (var k=0; k<group_indices.length; k++){
                var pk = parseInt(group_indices[k]);  // Convert string key back to integer
                var mask_indices = create_mask_index(pk, nrows, ncols);
                
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
                    peaks_found[g][index] = {
                        'peak_height': peak_height,
                        'grp_id': grp_id,
                        'F1_ppm': signal.F1_ppm,
                        'F2_ppm': signal.F2_ppm
                    };
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
            current_grp_id++;
            var new_group = {};
            new_group[index] = {
                'peak_height': peak_height,
                'grp_id': current_grp_id,
                'F1_ppm': signal.F1_ppm,
                'F2_ppm': signal.F2_ppm
            };
            peaks_found.push(new_group);
        }
        
    }

    // remove any peak groups that only have one point
    var filtered_peaks_found = [];
    for (var m=0; m<peaks_found.length; m++){
        var group_size = Object.keys(peaks_found[m]).length;
        if (group_size > 1){
            filtered_peaks_found.push(peaks_found[m]);
        }
    }

    return filtered_peaks_found;
}

function F1_from_index( index, ncols) {
    var F1 = Math.floor( index / ncols );
    return F1;
}

function F2_from_index( index, ncols) {
    var F2 = index % ncols;
    return F2;
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

    const expt_types = ["H1_1D", "C13_1D", "H1_pureshift", "DEPT135", "HSQC", "HMBC", "COSY", "NOESY", "DDEPT_CH3_ONLY", "HSQC_CLIPCOSY"];

    var spectra_dict = {};

    for( var i=0; i<expt_types.length; i++){
        var spectrum = find_NMR_spectrum(spectra_lst, expt_types[i]);
        if (spectrum) {
          spectra_dict[expt_types[i]] = spectrum;
        }
    }


    // dictionary of spectra with peaks
    var spectra_with_peaks = {};
    for (var key in spectra_dict) {
        if (spectra_dict[key]["hasPeaks"]) {
            spectra_with_peaks[key] = spectra_dict[key];
        }
    }


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

    

    var max_min_values = find_max_min_in_matrix(zoomed_region.rrr);
    var max_val_rrr = max_min_values.max;
    var min_val_rrr = max_min_values.min;

    // check if we have both positive and negative peaks or just positive or just negative
    var peak_types = determine_peaks_pos_and_neg( max_val_rrr, min_val_rrr );

    // find the peak coordinates
    var pk_coords = find_pks( zoomed_region.rrr, peak_types, max_val_rrr, min_val_rrr );


    print("Number of pt coordinates found: " + pk_coords.length );



    for (var i=0; i < pk_coords.length; i++){
        var F1_ppm = zoomed_region.F1[0] - pk_coords[i].F1 * zoomed_region.F1_ppm_delta;
        var F2_ppm = zoomed_region.F2[0] - pk_coords[i].F2 * zoomed_region.F2_ppm_delta;
        pk_coords[i].F1_ppm = F1_ppm;
        pk_coords[i].F2_ppm = F2_ppm;
    }

    print();

    // print out the pk_coords F1_ppm and F2_ppm values

    // for (var i=0; i < pk_coords.length; i++){
    //     print( "\tPeak coord " + (i+1) + ": F1_ppm: " + pk_coords[i].F1_ppm + ", F2_ppm: " + pk_coords[i].F2_ppm + ", intensity: " + pk_coords[i].intensity );
    // }
    
    // print();

    // // display all peak coordinates found on the spectrum
    // if (pk_coords.length == 0)
    //     print("No peak coordinates found.");

    // var mnova_peaks = new Peaks();
    // for (var i=0; i < pk_coords.length; i++){
    //     var peak = new Peak(pk_coords[i].F1_ppm, pk_coords[i].F2_ppm, pk_coords[i].intensity);
    //     mnova_peaks.append(peak);
    // }
    // activeSpectrum.setPeaks(mnova_peaks);
    // activeSpectrum.process();
    

    var peaks_found = group_pts( pk_coords, zoomed_region.num_F2_pts, zoomed_region.num_F1_pts );

    // // print out the peaks_found details
    // print("From Groups:\t")
    // print();
    // for (var i=0; i < peaks_found.length; i++){
    //     var peak_group = peaks_found[i];
    //     var group_indices = Object.keys(peak_group);
        
    //     // Print the group with its details
    //     print( "\tPeak group " + (i+1) + " (" + group_indices.length + " peaks):");
    //     for (var idx=0; idx < group_indices.length; idx++){
    //         var key = group_indices[idx];
    //         var peak_data = peak_group[key];
    //         print("\t  Index " + key + ": height=" + peak_data.peak_height + ", grp_id=" + peak_data.grp_id + ", F1_ppm=" + peak_data.F1_ppm + ", F2_ppm=" + peak_data.F2_ppm   );
    //     }
    // }

    

    print("\nNumber of peak groups found: " + peaks_found.length );
    print();
    for (var i=0; i < peaks_found.length; i++){
        var peak_group = peaks_found[i];
        var group_indices = Object.keys(peak_group);
        
        // Print the group with its details
        print( "\tPeak group " + (i+1) + " (" + group_indices.length + " peaks):");
        for (var idx=0; idx < group_indices.length; idx++){
            var key = group_indices[idx];
            var peak_data = peak_group[key];
            print("\t  Index " + key + ": height=" + peak_data.peak_height + ", grp_id=" + peak_data.grp_id + ", F1_ppm=" + peak_data.F1_ppm + ", F2_ppm=" + peak_data.F2_ppm   );
        }
        
        // // match the index value in peaks_found to the index in pk_coords and append the group index value
        // for (var j=0; j < group_indices.length; j++){
        //     var index_val = parseInt(group_indices[j]);  // Convert string key to integer
        //     for (var k=0; k < pk_coords.length; k++){
        //         if (pk_coords[k].index == index_val){
        //             pk_coords[k].pk_id = i+1;  // append group index (1-based)
        //             break;
        //         }
        //     }
        // }       
    }




    var peaks_found_centroid_in_pts = [];

    const num_peaks = peaks_found.length;
    for( var i=0; i<num_peaks; i++){
        var peak_group = peaks_found[i];
        var group_indices = Object.keys(peak_group);
        var num_pts_in_pk = group_indices.length;

        var sum_F1 = 0;
        var sum_F2 = 0;
        var peak_height = 0;

        for( var j=0; j<num_pts_in_pk; j++){
            var index = parseInt(group_indices[j]);  // Convert string key to integer
            var peak_data = peak_group[group_indices[j]];

            // sum_F1 += F1_from_index( index, zoomed_region.num_F1_pts);
            // sum_F2 += F2_from_index( index, zoomed_region.num_F1_pts);

            sum_F1 += peak_data.F1_ppm;
            sum_F2 += peak_data.F2_ppm;
            
            peak_height += peak_data.peak_height;
            // // Track maximum peak height in the group
            // if (Math.abs(peak_data.peak_height) > Math.abs(max_peak_height)){
            //     max_peak_height = peak_data.peak_height;
            // }
        }
        var centroid_F1 = sum_F1 / num_pts_in_pk;
        var centroid_F2 = sum_F2 / num_pts_in_pk;

        // var sum_F1 = 0;
        // var sum_F2 = 0;
        // var sum_weights = 0;
        // var peak_height = 0;

        // for( var j=0; j<num_pts_in_pk; j++){
        //     var index = parseInt(group_indices[j]);
        //     var peak_data = peak_group[group_indices[j]];
        //     var weight = Math.abs(peak_data.peak_height);
            
        //     sum_F1 += F1_from_index( index, zoomed_region.num_F1_pts) * weight;
        //     sum_F2 += F2_from_index( index, zoomed_region.num_F1_pts) * weight;
        //     sum_weights += weight;
        //     peak_height += peak_data.peak_height;
        // }

        // var centroid_F1 = sum_F1 / sum_weights;
        // var centroid_F2 = sum_F2 / sum_weights;


        print("From Centroids:\t")
        print( "peak " + i + ", F1: " + centroid_F1 + " F2: " + centroid_F2 + " peak_height: " + peak_height/num_pts_in_pk); 
        peaks_found_centroid_in_pts.push( { "F1": centroid_F1,
                                           "F2": centroid_F2,
                                           "max_height": peak_height/num_pts_in_pk } );
        // print( "peak " + i + ", F1: " + sum_F1/num_pts_in_pk + " F2: " + sum_F2/num_pts_in_pk + " peak_height: " + peak_height/num_pts_in_pk); 
        // peaks_found_centroid_in_pts.push( { "F1": sum_F1/num_pts_in_pk,
        //                                    "F2": sum_F2/num_pts_in_pk,
        //                                    "max_height": peak_height/num_pts_in_pk } );
    }

    // convert centroid positions from pts to ppm
    // first add the offset of the zoomed region
    // for( var i=0; i<peaks_found_centroid_in_pts.length; i++){
    //     peaks_found_centroid_in_pts[i].F1 += zoomed_region.F1_pts[0]-0.0;
    //     peaks_found_centroid_in_pts[i].F2 += zoomed_region.F2_pts[0]-0.0;
    // }   
    
    var peaks_found_centroid_in_ppm = peaks_found_centroid_in_pts;
    // for( var i=0; i<peaks_found_centroid_in_pts.length; i++){
    //     // var F1_ppm = activeSpectrum.unitsToUnits("pt", "ppm", peaks_found_centroid_in_pts[i].F1, 1);
    //     // var F2_ppm = activeSpectrum.unitsToUnits("pt", "ppm", peaks_found_centroid_in_pts[i].F2, 2);

    //     peaks_found_centroid_in_ppm.push( { "F1": F1_ppm,
    //                                         "F2": F2_ppm,
    //                                         "max_height": peaks_found_centroid_in_pts[i].max_height } );
    // }

    print();
    var mnova_peaks = new Peaks();
    for( var i=0; i<peaks_found_centroid_in_ppm.length; i++){
        print( "peak " + i + ", F1 ppm: " + peaks_found_centroid_in_ppm[i].F1 + " F2 ppm: " + peaks_found_centroid_in_ppm[i].F2 + " max_height: " + peaks_found_centroid_in_ppm[i].max_height); 
    
        // create MNOVA peak object and add to document
        var new_peak = new Peak(peaks_found_centroid_in_ppm[i].F1, peaks_found_centroid_in_ppm[i].F2, peaks_found_centroid_in_ppm[i].max_height);
        mnova_peaks.append(new_peak);
    }

    // add peaks to the active spectrum
    // add peaks to the 2D spectrum
    activeSpectrum.setPeaks(mnova_peaks);
    activeSpectrum.process();


    // print("Spectra keys:");
    // print(spectra_keys);

    // print( "Experiments dict");
    // print(experiments_dict);

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



