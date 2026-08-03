/* Globals Dialog, GroupBox, CheckBox */

// 01/05/2025 Updated to add wgich server we are using

function idspectra_dialog(spec_lst ) {

    // spectra_available is a list of spectra available in the document
    // create a dialog with checkboxes for each spectrum

    print("inside idspectra_dialog");

    // return if spec_lst is undefined
    if (spec_lst == undefined) {
        print("spec_lst is undefined");
        return undefined;
    }


    var dialog = new Dialog();
    dialog.title = "Identify Spectra";

    // create checkboxes and store them in a dictionary

    var checkboxes1 = {};
    var comboboxes1 = {};

    var simpleutils = new simpleUtils();

    // var NMRexpts = ["unidentified",
    //                 "H1_1D", 
    //                 "C13_1D", 
    //                 "DEPT135", 
    //                 "H1_pureshift", 
    //                 "COSY", "HSQC", 
    //                 "HMBC", 
    //                 "HSQC_CLIPCOSY", 
    //                 "DDEPT_CH3_ONLY", 
    //                 "NOESY", 
    //                 "Predicted"];

    var NMRexpts = spectra_keys;
        // var NMRexpts = ["undefined", "H1_1D", "C13_1D", "DEPT135", "PURESHIFT", "COSY", "HSQC", "HMBC", "HSQCCLIPCOSY", "DOUBLEDEPTCH3", "NOESY", "Predicted"];

    // create checkboxes for each spectrum in the list spectra_available
    print("\nCreating checkboxes for spectra:\n");
    for (var i = 0; i < spec_lst.length; i++) {
        var spectrum = spec_lst[i];

        var spec_id = simpleutils.exptUniqueIdentifier(spectrum);
        checkboxes1[spec_id] = new CheckBox();
        checkboxes1[spec_id].text = spec_id;

        comboboxes1[spec_id] = new ComboBox();
        for( var j = 0; j < NMRexpts.length; j++ ){        
            comboboxes1[spec_id].addItem(NMRexpts[j]);                           
        }
        // set current text to spectrum.experimentEEH
        print(i + " spectrum.experimentEEH: " + spectrum.experimentEEH);
        comboboxes1[spec_id].currentText = spectrum.experimentEEH;
    }
    print();




    // create a groupbox and add the checkboxes to it
    var gbSpectra = new GroupBox();
    gbSpectra.title = "Choose Spectra to use in Analysis";

    for (var key in checkboxes1) {
        gbSpectra.add(checkboxes1[key], comboboxes1[key]);
    }

    dialog.add(gbSpectra);



    // display the dialog and return the list of checked checkboxes
    if (!dialog.exec()) {
        return undefined;
    }
    else {
        var i_expt = 0;
        for( var ky in comboboxes1){
            print("Selected experiment for " + ky + " is " + comboboxes1[ky].currentText);
            spec_lst[i_expt]["experimentEEH"] = comboboxes1[ky].currentText;
            // set experimentIdentified to true if not undefined
            if( comboboxes1[ky].currentText != "undefined"){
                spec_lst[i_expt]['experimentIdentified'] = true;
            }
            else{
                spec_lst[i_expt]['experimentIdentified'] = false;
            }   
            // set Notes field to currentText
            spec_lst[i_expt]['page']['notes'] = comboboxes1[ky].currentText;
            i_expt = i_expt + 1;
        }
    }
}

// var spec_lst = identify_spectrum();

// print("spec_list.length", spec_lst.length);
// idspectra_dialog( spec_lst );