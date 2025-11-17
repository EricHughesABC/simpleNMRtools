/* Globals Dialog, GroupBox, CheckBox */

// 01/05/2025 Updated to add wgich server we are using

function xy3_dialog2(spectra_available, exptIdentifierIndex, dialogParameters, exptName, servername) {

    // spectra_available is a list of spectra available in the document
    // create a dialog with checkboxes for each spectrum

    print("inside xy3_dialog2");

    print("dialogParameters\n", dialogParameters);
    print("exptName ", exptName);
    print("servername ", servername);

    var dialog = new Dialog();
    print("servername ", servername);
    dialog.title = exptName + " " + servername;

    // create checkboxes and store them in a dictionary

    var checkboxes1 = {};
    var comboboxes1 = {};

    // var NMRexpts = ["SKIP", "H1_1D", "C13_1D", "DEPT135", "PureShift", "COSY", "HSQC", "HMBC", "HSQC_CLIPCOSY", "DDEPT_CH3_ONLY", "NOESY"];
    var NMRexpts = ["SKIP", "H1_1D", "C13_1D", "DEPT135", "PureShift", "COSY", "HSQC", "HMBC", "HSQC_CLIPCOSY", "DDEPT_CH3_ONLY", "NOESY"];

    // create checkboxes for each spectrum in the list spectra_available
    for (var i = 0; i < spectra_available.length; i++) {
        checkboxes1[spectra_available[i]] = new CheckBox();
        checkboxes1[spectra_available[i]].text = spectra_available[i];

        comboboxes1[spectra_available[i]] = new ComboBox();
        for( var j = 0; j < NMRexpts.length; j++ ){        
            comboboxes1[spectra_available[i]].addItem(NMRexpts[j]);                           
        }            
    }

    // set the expt identifier to the combo box

    if ( exptIdentifierIndex.length == spectra_available.length ){
        for( var i =0; i < spectra_available.length; i++){
            comboboxes1[spectra_available[i]].currentText = exptIdentifierIndex[i];
        }
    }

    // create a groupbox and add the checkboxes to it
    var gbSpectra = new GroupBox();
    gbSpectra.title = "Choose Spectra to use in Analysis";

    for (var key in checkboxes1) {
        gbSpectra.add(checkboxes1[key], comboboxes1[key]);
    }

    dialog.add(gbSpectra);

    // set the default values for the checkboxes based on  the list of spectra  in the document


    for (var key in checkboxes1) {
        if (spectra_available.indexOf(key) != -1) {
            checkboxes1[key].checked = true;
        } else {
            checkboxes1[key].checked = false;
            // set enabled to false
            checkboxes1[key].enabled = false;
        }
    }

    // add a radio button group to predict the assignment of the peaks
    var rButtonGroupBoxPredictAssign = new GroupBox();
    rButtonGroupBoxPredictAssign.objectName="rButtonsGroupBoxPredictAssign";
    rButtonGroupBoxPredictAssign.title = "Predict C13 Chemical Shifts Using... ";

    var rButton1_MNOVApredictAssign = new RadioButton();
    rButton1_MNOVApredictAssign.objectName = "radioButton1_MNOVApredictAssign";
    rButton1_MNOVApredictAssign.text = "MNOVA Predict and Assign data";

    if( dialogParameters["calcSimpleMNOVA"] == "MNOVA Predict" ){
        rButton1_MNOVApredictAssign.checked = true;
    }
    else{
        rButton1_MNOVApredictAssign.checked = false;
    }

    var rButton2_MNOVAmanuallyAssign = new RadioButton();
    rButton2_MNOVAmanuallyAssign.objectName = "radioButton2_MNOVAmanuallyAssign";
    rButton2_MNOVAmanuallyAssign.text = "MNova Manually Assigned data";
    if( dialogParameters["calcSimpleMNOVA"] == "MNova Manually Assigned" ){
        rButton2_MNOVAmanuallyAssign.checked = true;
    }
    else{
        rButton2_MNOVAmanuallyAssign.checked = false;
    }

    var rButton3_NMRSHIFTDB2predictAssign = new RadioButton();
    rButton3_NMRSHIFTDB2predictAssign.objectName = "radioButton3_NMRSHIFTDB2predictAssign";
    rButton3_NMRSHIFTDB2predictAssign.text = "NMRSHIFTDB2 Predict and Assign data";
    if( dialogParameters["calcSimpleMNOVA"] == "NMRSHIFTDB2 Predict" ){
        rButton3_NMRSHIFTDB2predictAssign.checked = true;
    }
    else{
        rButton3_NMRSHIFTDB2predictAssign.checked = false;
    }

    // add the radio buttons to the group box
    rButtonGroupBoxPredictAssign.add(rButton1_MNOVApredictAssign);
    // rButtonGroupBoxPredictAssign.add(rButton2_MNOVAmanuallyAssign);
    rButtonGroupBoxPredictAssign.add(rButton3_NMRSHIFTDB2predictAssign);
    rButtonGroupBoxPredictAssign.add(tickButton_SimulatedAnnealing);

    // add the group box to the dialog
    dialog.add(rButtonGroupBoxPredictAssign);

    // add a simulated Annealing group box
    var GroupBoxSimulatedAnnealing = new GroupBox();

    var tickButton_SimulatedAnnealing = new CheckBox();
    // set the text for the tick button
    tickButton_SimulatedAnnealing.text = "Optimize Assignments using HMBC and COSY Information";
    // set the checked property to True
    if (dialogParameters["simulatedAnnealing"]) {

        tickButton_SimulatedAnnealing.checked = true;
    } else {
        tickButton_SimulatedAnnealing.checked = false;
    }

    // add simulated annealing widgets to the group box
    GroupBoxSimulatedAnnealing.add(tickButton_SimulatedAnnealing);

    // add the group box to the dialog
    dialog.add(GroupBoxSimulatedAnnealing);

    var GroupBox_MachineLearningKeep = new GroupBox();
    var tickButton_MachineLearningKeep = new CheckBox();
    // set the text for the tick button
    tickButton_MachineLearningKeep.text = "Allow Results to be kept in Machine Learning Database";
    // set the checked property to True
    print("dialogParameters[\"ml_consent\"]", dialogParameters["ml_consent"]);
    print("dialogParameters[\"simulatedAnnealing\"]", dialogParameters["simulatedAnnealing"]);
    if (dialogParameters["ml_consent"]) {
        print("Machine Learning consent is true");
        tickButton_MachineLearningKeep.checked = true;
    } else {
        print("Machine Learning consent is false");
        tickButton_MachineLearningKeep.checked = false;
    }

    // add simulated annealing widgets to the group box
    GroupBox_MachineLearningKeep.add(tickButton_MachineLearningKeep);
    // add the group box to the dialog
    dialog.add(GroupBox_MachineLearningKeep);

    // display the dialog and return the list of checked checkboxes
    if (!dialog.exec()) {
        return undefined;
    }
    else {
        var rtn_dict = {};
        rtn_dict["checked_spectra"] = [];
        rtn_dict["exptIdentifiers"] = []
        // var checked = ;
        var i_expt = 0;
        for (var key in checkboxes1) {
            if (checkboxes1[key].checked) {
                rtn_dict["checked_spectra"].push(key + " " + comboboxes1[spectra_available[i_expt]].currentText);
            }
            rtn_dict["exptIdentifiers"].push(comboboxes1[spectra_available[i_expt]].currentText);
            i_expt = i_expt + 1;
        }

        print( "find out predict or use assigned data");
        if (rButton1_MNOVApredictAssign.checked) {
            print("MNOVA Predict and Assign data");
            rtn_dict["calcSimpleMNOVA"] = "MNOVA Predict";
        }
        else if (rButton3_NMRSHIFTDB2predictAssign.checked) {
            print("NMRSHIFTDB2 Predict and Assign data");
            rtn_dict["calcSimpleMNOVA"] = "NMRSHIFTDB2 Predict";
        }
        else{
            print("No radio button selected");
        }

        // return the simulated annealing boolean value
        rtn_dict["simulatedAnnealing"] = tickButton_SimulatedAnnealing.checked;
        rtn_dict["ml_consent"] = tickButton_MachineLearningKeep.checked;

        print(rtn_dict);
        return rtn_dict;
    }
}