/*globals WebUtilsQT NMRAssignments xy3_dialog Dir brukerExperiments String iupac_dialog server_address MnUi*/

function simpleASSIGN_eeh() {

    function splitPathAndFilename(path) {
        var separatorIndex = path.lastIndexOf('/');
        var directoryPath = path.substring(0, separatorIndex);
        var filename = path.substring(separatorIndex + 1);

        var dotIndex = filename.lastIndexOf('.');
        var name = filename.substring(0, dotIndex);
        var extension = filename.substring(dotIndex + 1);

        return {
            directoryPath: directoryPath,
            filename: filename,
            name: name,
            extension: extension
        };
    }

    var doc = Application.mainWindow.activeDocument;

    //////////////////////////////////////
    // choose which server we are using
    /////////////////////////////////////

    var server = server_address();
    var web_utils = new WebUtilsQT();

    // get the hostname
    var hostname = SysInfo.hostID;

    print("hostname ", hostname);

    // see if the hostname is registered
    // if it is then the server will return the machine_learning status
    // if is is not registered then the server will return a redirect to the registration page
    var entry_point = server + "check_machine_learning";
    var json_obj = {};
    json_obj["hostname"] = hostname;
    var json_obj_string = JSON.stringify(json_obj,null,4);
    // print("json_obj_string ", json_obj_string);
    rtn = web_utils.jsonRequest(entry_point, json_obj_string, "", "", false);

    // print("entry_point ", entry_point);

    // print(rtn.url);
    // print(rtn.exitCode);
    // print(rtn.allErrorOutput);
    // print(rtn.allStrOutput);
    print("\nrtn.response\n");
    print(rtn.response);

    // Check if response is registration redirect
    try {
        // var responseObj = JSON.parse(rtn.response || rtn.allStrOutput);
        var responseObj = JSON.parse(rtn.response);
        
        if (responseObj.status === "unregistered" || responseObj.status === "registration_expired") {
            // Handle unregistered host
            // display a message box
            if( responseObj.status === "registration_expired") {
                MessageBox.warning("Host registration has expired. Please reregister to use this feature.");
            }
            
            if (responseObj.status === "unregistered") {
                MessageBox.warning("Host is not registered. Please register to use simpleNMR prediction.");
            }
            
            if (typeof Application.openUrl === 'function') {
                Application.openUrl(responseObj.registration_url);
            } else {
                print("Application.openUrl is not a function");
            }
            
            return; // Exit early
        }
        
        var ml_consent = responseObj.ml_consent;
        print("ml_consent ", ml_consent);
        // If we get here, response is normal JSON
        // htmloutputstr = responseObj.html || responseObj.content;
    } catch (e) {
        // Not JSON or parsing error - continue with original response handling
        // replace the dummy title in the html output with result.name
        // MessageBox.warning("No JSON response found. Please check the server response.");
        // htmloutputstr = rtn.response.replace("dummy_title", result.name);
        print("Error parsing JSON response: ", e);
        MessageBox.warning("Error parsing JSON response: " + e);
        return;
    }



    var path = doc.name;
    var result = splitPathAndFilename(path);

    // read in experiment identifiers from file expt_identifiers
    var jsonfilename = result.directoryPath + "/" + result.name + "_assignments_from_simplemnova.json";
    print(jsonfilename);
        
    var fin = new File(jsonfilename);
    if (fin.open(File.ReadOnly)) {
        print("Reading json file");
        var sin = new TextStream(fin, 'UTF-8');
        var jsonstr = sin.readAll();
        var jsondata = JSON.parse(jsonstr);
        // var jsondata_orig = JSON.parse(jsonstr);
        
    }
    else{
        print("Could not open json file");
        print(jsonfilename);
        MessageBox.warning("Could not open json file");
        return;
    }
    fin.close();
        
    // make a copy of the jsondata

    var carbonatoms_orig = JSON.parse(JSON.stringify(jsondata["nodes_orig"]));
    var carbonatoms_moved = JSON.parse(JSON.stringify(jsondata["nodes_now"]));
    var oldjsondata = JSON.parse(JSON.stringify(jsondata["oldjsondata"]));
    // var carbonatoms_orig = jsondata["nodes_orig"];
    // var carbonatoms_moved = jsondata["nodes_now"];
    // convert string to json object
    // var oldjsondata = jsondata["oldjsondata"];
    // convert the string oldjsondata to a json object
    // var oldjsondata = jsondata["oldjsondata"];

    // replace the atom_idx in carbons_moved with the atom_idx in carbons_orig using the atomNumber to match
    var carbonatoms = [];
    carbonatoms_moved.forEach(function(atom){
        var atomNumber = atom["atomNumber"];

        for (var i = 0; i < carbonatoms_orig.length; i++) {
            if (carbonatoms_orig[i]["atomNumber"] == atomNumber) {
                atom["id"] = carbonatoms_orig[i]["id"];
                carbonatoms.push(atom);
                break;
            }
        }
    });
    

	"use strict";
	var moleculePlugin = new MoleculePlugin(Application.molecule);
    var assign = new NMRAssignments(moleculePlugin.activeMolecule());

	if( moleculePlugin === undefined) {

        MessageBox.warning("Could not find a molecule in the current document");
		return;
	}

	var mol = new Molecule(moleculePlugin.activeMolecule());
	if( mol.isValid() === false ) {

        MessageBox.warning("Could not create a mol object from the molecule in the current document");
		throw "Invalid Molecule";
	}


	// atom = new Atom(mol.atom(atomNumber));
	// assign = new NMRAssignments(mol);	
	// shift = assign.chemShift(atomNumber, h);
	// if( shift !== undefined ) {
	// 	print(JSON.stringify(shift, null, 2));
	// }

    carbonatoms.forEach(function(atom){
        //convert atomNumber to integer
        var atomNumber = parseInt(atom["atomNumber"]);
        var c13_ppm = parseFloat(atom["ppm"]);
        var numProtons = parseInt(atom["numProtons"]);

        var atom_idx = parseInt(atom["id"]);

        // add 1 to atom_idx if dataFrom == "nmrshiftdb2"
        print("jsondata[\"dataFrom\"] = ", jsondata["dataFrom"]);
        if( jsondata["dataFrom"] == "nmrshiftdb2"){
            print("Adding 1 to atom_idx");
            atom_idx += 1;
        }
        // var atom_idx = parseInt(atom["id"]); // atom index starts from 1
        print(atom_idx, atomNumber, c13_ppm);

        assign.setChemShift(atom_idx, 0, c13_ppm);



        // assign the protons
        if( numProtons == 1){
            if( atom["H1_ppm"] === undefined ){
                print("No protons found for atom ", atomNumber);
                return;
            }
            var h1_ppm = parseFloat(atom["H1_ppm"][0]);
            assign.setChemShift(atom_idx, 1, h1_ppm);
        }
        else if( numProtons == 3){
            if( atom["H1_ppm"] === undefined ){
                print("No protons found for atom ", atomNumber);
                return;
            }
            var h1_ppm = parseFloat(atom["H1_ppm"][0]);
            assign.setChemShift(atom_idx, 1, h1_ppm);
            assign.setChemShift(atom_idx, 2, h1_ppm);
            assign.setChemShift(atom_idx, 3, h1_ppm);
        }
        else if( numProtons == 2){
            if( atom["H1_ppm"] === undefined ){
                print("No protons found for atom ", atomNumber);
                return;
            }
            if( atom["H1_ppm"].length == 1){
                var h1_ppm = parseFloat(atom["H1_ppm"][0]);
                assign.setChemShift(atom_idx, 1, h1_ppm);
                assign.setChemShift(atom_idx, 2, h1_ppm);
            }
            else{
                var h1_ppm = parseFloat(atom["H1_ppm"][0]);
                var h2_ppm = parseFloat(atom["H1_ppm"][1]);
                assign.setChemShift(atom_idx, 1, h1_ppm);
                assign.setChemShift(atom_idx, 2, h2_ppm);
            }
        }
    });

    // rebuild html file after assigning molecule

    var spectra = {};

    // copy smiles from oldjsondata

    spectra["smiles"] = oldjsondata["smiles"];

    // copy molfile from oldjsondata
    spectra["molfile"] = oldjsondata["molfile"];

    // copy carbonAtomsInfo from oldjsondata
    spectra["carbonAtomsInfo"] = oldjsondata["carbonAtomsInfo"];


    // // add hostname to spectra object

    spectra["hostname"] = {};
    spectra["hostname"]["datatype"] = "hostname";
    spectra["hostname"]["count"] = 1;
    spectra["hostname"]["data"] = {};
    spectra["hostname"]["data"]["0"] = SysInfo.hostID;

    // copy nmrAssignments from mnova file

    // reget molecule from mnova file as it has been updated with assignments

    var moleculePlugin = new MoleculePlugin(Application.molecule);
    var assign = new NMRAssignments(moleculePlugin.activeMolecule());

	if( moleculePlugin === undefined) {
		return;
	}

	var mol = new Molecule(moleculePlugin.activeMolecule());
	if( mol.isValid() === false ) {
		throw "Invalid Molecule";
	}

    var assignments = NMRAssignments(mol);
    print("hasAssignments ", assignments.hasAssignments());
    spectra["nmrAssignments"] = {};
    spectra["nmrAssignments"]["datatype"] = "nmrAssignments";
    spectra["nmrAssignments"]["count"] = 0;
    spectra["nmrAssignments"]["data"] = {};

    // loop through atoms in molecule and extract proton and chemical shifts and atom numbering
    
    if ( assignments.hasAssignments() == true){

        var count = 0;
        for ( var i=1; i <= mol.atomCount; i++) {    
            var atom = new Atom(mol.atom(i));
    
            if( atom.elementSymbol == "C" ){
                spectra["nmrAssignments"]["data"][i] = {};
                cshift = assignments.chemShift(i,0); 
                if( cshift != undefined ){
                    print("cshift ", cshift ); 

                	spectra["nmrAssignments"]["data"][i]["atom_idx"] = i;      
                	spectra["nmrAssignments"]["data"][i]["atomNumber"] = atom.number;
                	spectra["nmrAssignments"]["data"][i]["numProtons"] = atom.nH;
                	spectra["nmrAssignments"]["data"][i]["f1_ppm"] = cshift.shift;
                	spectra["nmrAssignments"]["data"][i]["label"] = atom.Label;
                	spectra["nmrAssignments"]["data"][i]["iupacLabel"] = "";
                	count++;
    
                	for( j=1; j<=3; j++ ){
                    	hshift = assignments.chemShift(i,j);	    
                    	if( hshift === undefined ){	      
                        	spectra["nmrAssignments"]["data"][i]["f2_ppm"+j] = "undefined";
                    	}
                    	else {	      
                        	spectra["nmrAssignments"]["data"][i]["f2_ppm"+j] = hshift.shift
                    	}
                	}
                }	      
            }
            else{
                print(i, atom.elementSymbol);
                print(atom);
            }

        }
        spectra["nmrAssignments"]["count"] = count;
    }

    // copy c13predictions from oldjsondata
    spectra["c13predictions"] = oldjsondata["c13predictions"];

    // copy over spectra info from oldjsondata
    // loop through a list of the spectra "HSQC", "HMBC", "COSY", "H1_1D", "C13_1D", ""

    // print out the keys in oldjsondata
    // print("oldjsondata keys\n", Object.keys(oldjsondata));

    // print(oldjsondata);


    for( var idx=0; idx<oldjsondata["exptIdentifiers"]["count"]; idx++){
        var ky = oldjsondata["exptIdentifiers"]["data"][idx];
        print(ky);
        spectra[ky] = oldjsondata[ky];
    }

    // copy over exptIdentifiers from oldjsondata
    spectra["exptIdentifiers"] = oldjsondata["exptIdentifiers"];

    // copy over MNOVAcalcMethod from oldjsondata
    spectra["MNOVAcalcMethod"] = oldjsondata["MNOVAcalcMethod"];

    // change the data to use assign

    spectra["MNOVAcalcMethod"]["data"][0] = "MNova Manually Assigned data";

    // copy over carbonCalcPositionsMethod from oldjsondata
    spectra["carbonCalcPositionsMethod"] = oldjsondata["carbonCalcPositionsMethod"];

    // copy over workingDirectory from oldjsondata
    spectra["workingDirectory"] = oldjsondata["workingDirectory"];

    // copy over workingFilename from oldjsondata
    spectra["workingFilename"] = oldjsondata["workingFilename"];

    // copy over chosenSpectra from oldjsondata
    spectra["chosenSpectra"] = oldjsondata["chosenSpectra"];

    // save spectra JSON string to file in data directory
    var jsonfilename = result.directoryPath + "/" + result.name + "_assignments_mresnova.json";

    // save spectra json string to file
    var json_spectra_string = JSON.stringify(spectra,null,4);
    var json_data_string = JSON.stringify(jsondata,null,4);
    // var json_data_string = jsonstr;
    var fout1 = new File(jsonfilename);
    if (fout1.open(File.WriteOnly)) {
        var sout1 = new TextStream(fout1, 'UTF-8');
        sout1.writeln(json_spectra_string);
    }
    fout1.close(); 

    // post json string to server
    print("Post to the server");
    // var web_utils = new WebUtilsQT();

    ////////////////////
    // Do some checking
    ////////////////////

    var chosen_spectra = spectra["chosenSpectra"]["data"][0];
    var hsqc = spectra["HSQC"];
    var h1_1d = spectra["H1_1D"];

    // check if there are any assignments if the user has chosen to use MNOVA manually assigned
    if( spectra["nmrAssignments"]["count"] == 0 && chosen_spectra == "MNova Manually Assigned"){
        print("No assignments found");
        MessageBox.warning("No MNOVA assignment data found");  
        return;
    }

    // check if there are any C13 predictions if the user has chosen to use MNOVA Predict
    if( spectra["c13predictions"]["count"] == 0 && chosen_spectra == "MNOVA Predict"){
        print("No C13 predictions found");
        MessageBox.warning("No C13 predictions found");  
        return;
    }

    // check if there are no integrals and no peaks in the HSQC spectrum
    if( hsqc["peaks"]["count"] == 0 && hsqc["integrals"]["count"] == 0){
        print("No peaks or integrals found in HSQC spectrum");
        MessageBox.warning("No peaks or integrals found in HSQC spectrum");  
        return;
    }

    // warn the user that there are no integrals in the HSQC spectrum
    if(  hsqc["integrals"]["count"] == 0){
        print("No peaks or integrals found in HSQC spectrum");
        MessageBox.warning("No integrals found in HSQC spectrum, using peaks only. Program works best with integrals");  
    }

    // check if the number of peaks in hsqc is == number of multiplets in h1_1d
    if( h1_1d !== undefined ){
        if( hsqc["peaks"]["count"] != h1_1d["multiplets"]["count"]){
            print("Number of peaks in HSQC not equal to number of multiplets in H1_1D");
            MessageBox.warning("Number of peaks in HSQC not equal to number of multiplets in H1_1D. Do not use H1_1D");  
            return;
        }
    }

    //////////////////////////////////////
    // choose which server we are using
    /////////////////////////////////////

    var server = server_address();
    // var server = "http://simplenmr.pythonanywhere.com/";
    // var server = "http://simplenmr.awh.durham.ac.uk/";

    ///////////////////////////////////////
    // Decide which entry point to call
    //////////////////////////////////////


    // build entry point
    // var entry_point = server + "simpleMNOVAmanuallyassigned"
    if( spectra["nmrAssignments"]["count"] == 0 ){
        print("No assignments found");
        MessageBox.warning("No MNOVA assignment data found");  
        return;
    }
    print(spectra["nmrAssignments"]["count"])
    var entry_point = server + "/simpleMNOVAfinalHTML"

    

    print("entry_point ", entry_point);
    // json_spectra_string = JSON.encode(json_spectra_string.encode, 'utf-8');
    // var rtn = web_utils.jsonRequest(entry_point, json_spectra_string, "", "", false);
    var rtn = web_utils.jsonRequest(entry_point, json_data_string,  "", "", false);

    print(rtn.url);


    
    // write rtn.response string to a html file
    var htmlfilename = result.directoryPath + "/html/" + result.name + "_d3.html";

    print("htmlfilename ", htmlfilename);

    // create directory if it does not exist
    var dir = new Dir(result.directoryPath + "/html");
    if( !dir.exists ){
        dir.mkdir(result.directoryPath + "/html");
    }

    // replace the dummy title in the html ouput with result.name
    htmloutputstr = rtn.response.replace( "dummy_title", result.name );

    var fout = new File(htmlfilename);
    if (fout.open(File.WriteOnly)) {
        sout = new TextStream(fout, 'UTF-8');
        sout.write(htmloutputstr);
    }
    fout.close();

    // display message dialog to tell the user the html file has been updated
    MessageBox.information("Molecule has been assigned and  the html file has been updated");

    // // write rtn.response string to a html file
    // var htmlfilename = result.directoryPath + "/html/" + result.name + "_d3.html";

    // // create directory if it does not exist
    // var dir = new Dir(result.directoryPath + "/html");
    // if( !dir.exists ){
    //     dir.mkdir(result.directoryPath + "/html");
    // }

    // // replace the dummy title in the html ouput with result.name
    // htmloutputstr = rtn.response.replace( "dummy_title", result.name );

    // var fout = new File(htmlfilename);
    // if (fout.open(File.WriteOnly)) {
    //     sout = new TextStream(fout, 'UTF-8');
    //     sout.write(htmloutputstr);
    // }
    // fout.close();

    // open html file in browser
    if (typeof Application.openUrl === 'function') {
        if(SysInfo.isWin){
                print("htmlfilename\n\t", htmlfilename );
            Application.openUrl(htmlfilename);
        }
        else if (SysInfo.isMac){
            print("htmlfilename\n\tfile://", htmlfilename );
            Application.openUrl("file://" + htmlfilename);
        }       	 
    }
    else{
        print("Applcation.openUrl is not a function");
    }

}

if (this.MnUi && MnUi.simpleNMRtools) {
	MnUi.simpleNMRtools.simpleNMRtools_simpleASSIGN_eeh = simpleASSIGN_eeh;
}