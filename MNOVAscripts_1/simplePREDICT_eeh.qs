function simplePREDICT_eeh(){


    // var spectra_keys = ["HSQC", 
    //                     "HMBC", 
    //                     "COSY", 
    //                     "NOESY", 
    //                     "H1_1D", 
    //                     "C13_1D", 
    //                     "DEPT135", 
    //                     "PureShift", 
    //                     "DDEPT_CH3_ONLY", 
    //                     "SKIP", 
    //                     "HSQC_CLIPCOSY"];

    var spectra_keys = ["HSQC", 
                        "HMBC", 
                        "COSY", 
                        "NOESY", 
                        "H1_1D", 
                        "C13_1D", 
                        "DEPT135", 
                        "PureShift", 
                        "DDEPTCH3ONLY", 
                        "SKIP", 
                        "HSQCCLIPCOSY"];

    function stringInArray(needle, haystack){
        for( var i=0; i<haystack.length; i++){
            if( haystack[i] === needle){
                return true;
          }
        }
        return false;
    }

    function isObjectEmpty(obj) {
        return Object.keys(obj).length === 0;
    }  

	function getActiveMolecule(aDocWin, aMolPlugin) {
	    var molec = aMolPlugin.activeMolecule();
		
		if (molec.isValid()) {
		    return molec;
		}
		if (aDocWin.itemCount("Molecule") == 1) {
			molec = new Molecule(aDocWin.item(0, "Molecule"));
			return molec;
		}
		return undefined;
	}

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
    
    // read in the json file and return the json object
    function readJsonFile(jsonfilename){

        var jsonobj = {};
        var fin = new File(jsonfilename);
        if (fin.open(File.ReadOnly)) {
            var sin = new TextStream(fin, 'UTF-8');
            var jsonstr = sin.readAll();
            jsonobj = JSON.parse(jsonstr);
        }
        else{
            MessageBox.warning("Could not open json file\n" + jsonfilename);
        }
        fin.close();
        return jsonobj;
    }



    // iterate through pages and print out what is on the page
    var doc = Application.mainWindow.activeDocument;

    //////////////////////////////////////
    // choose which server we are using
    /////////////////////////////////////

    var server = server_address();
    var web_utils = new WebUtilsQT();

    // get the hostname
    var hostname = SysInfo.hostID;

    print("hostname ", hostname);
    print("server ", server);

    // see if the hostname is registered
    // if it is then the server will return the machine_learning status
    // if is is not registered then the server will return a redirect to the registration page
    var entry_point = server + "check_machine_learning";
    var json_obj = {};
    json_obj["hostname"] = hostname;
    var json_obj_string = JSON.stringify(json_obj,null,4);
    // print("json_obj_string ", json_obj_string);
    print( "Entry point", entry_point);
    rtn = web_utils.jsonRequest(entry_point, json_obj_string, "", "", false);

    // print("entry_point ", entry_point);

    // print(rtn.url);
    print("rtn.exitCode", rtn.exitCode);
    print("rtn.allErrorOutput", rtn.allErrorOutput);
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
        print("Error parsing JSON response: ", e);
        MessageBox.warning("Error parsing JSON response: " + e + "\n" + rtn.response + "\nrtn.exitCode" + rtn.exitCode + "\nrtn.allErrorOutput" + rtn.allErrorOutput);
        return;
    }
  
    // return molecule from document
    var mol = getActiveMolecule(doc, Application.molecule);
    
    var spectra = {};
    if( mol === undefined ){
        MessageBox.warning("No molecule found!")
        return;  
    }

    var smiles = mol.generateSMILES()
    spectra["smiles"] = {};
    spectra["smiles"]["datatype"] = "smiles";
    spectra["smiles"]["count"] = 1;
    spectra["smiles"]["data"] = {};
    spectra["smiles"]["data"]["0"] = smiles;

    spectra["molfile"] = {};
    spectra["molfile"]["datatype"] = "molfile";
    spectra["molfile"]["count"] = 1;
    spectra["molfile"]["data"] = {};
    spectra["molfile"]["data"]["0"] = mol.getMolfile(); 

    // add hostname to spectra object

    spectra["hostname"] = {};
    spectra["hostname"]["datatype"] = "hostname";
    spectra["hostname"]["count"] = 1;
    spectra["hostname"]["data"] = {};
    spectra["hostname"]["data"]["0"] = SysInfo.hostID;    

    // get atom numbers for all atoms in molecule
    spectra["allAtomsInfo"] = {};
    spectra["allAtomsInfo"]["datatype"] = "allAtomsInfo";

    var allAcomCount = 0;
    spectra["allAtomsInfo"]["data"] = {};
    for ( var i=1; i <= mol.atomCount; i++) {    
        var atom = new Atom(mol.atom(i));
        var j = i-1;

        spectra["allAtomsInfo"]["data"][j] = {};
        spectra["allAtomsInfo"]["data"][j]["atom_idx"] = j;  
        spectra["allAtomsInfo"]["data"][j]["id"] = j;  
        spectra["allAtomsInfo"]["data"][j]["atomNumber"] = atom.number;
        spectra["allAtomsInfo"]["data"][j][["symbol"]] = atom.elementSymbol;
        spectra["allAtomsInfo"]["data"][j]["numProtons"] = atom.nH; 
        allAcomCount++;                      
    }
    spectra["allAtomsInfo"]["count"] = allAcomCount;

    // get atom numbers for carbon atoms in molecule
    spectra["carbonAtomsInfo"] = {};
    spectra["carbonAtomsInfo"]["datatype"] = "carbonAtomsInfo";

    var carbonAcomCount = 0;
    spectra["carbonAtomsInfo"]["data"] = {};
    for ( var i=1; i <= mol.atomCount; i++) {    
        var atom = new Atom(mol.atom(i));
        var j = i-1;

        if( atom.elementSymbol == "C" ){
            spectra["carbonAtomsInfo"]["data"][j] = {};
            spectra["carbonAtomsInfo"]["data"][j]["atom_idx"] = j;      
            spectra["carbonAtomsInfo"]["data"][j]["id"] = j;      
            spectra["carbonAtomsInfo"]["data"][j]["atomNumber"] = atom.number;
            spectra["carbonAtomsInfo"]["data"][j][["symbol"]] = atom.elementSymbol;
            spectra["carbonAtomsInfo"]["data"][j]["numProtons"] = atom.nH; 
            carbonAcomCount++;                      
        }
    }
    spectra["carbonAtomsInfo"]["count"] = carbonAcomCount;

    var assignments = NMRAssignments(mol);
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

                	spectra["nmrAssignments"]["data"][i]["atom_idx"] = i;      
                	spectra["nmrAssignments"]["data"][i]["atomNumber"] = atom.number;
                	spectra["nmrAssignments"]["data"][i]["numProtons"] = atom.nH;
                	spectra["nmrAssignments"]["data"][i]["f1_ppm"] = cshift.shift;
                	spectra["nmrAssignments"]["data"][i]["label"] = atom.Label;
                	spectra["nmrAssignments"]["data"][i]["iupacLabel"] = "";
                	count++;
                
    
                	for( var k=1; k<=3; k++ ){
                    	hshift = assignments.chemShift(i,k);	    
                    	if( hshift === undefined ){	      
                        	spectra["nmrAssignments"]["data"][i]["f2_ppm"+k] = "undefined";
                    	}
                    	else {	 
                        	spectra["nmrAssignments"]["data"][i]["f2_ppm"+k] = hshift.shift;
                    	}
                	}
                }	      
            }
            else{
                print(i, atom.elementSymbol );
                print(atom);
            }

        }
        spectra["nmrAssignments"]["count"] =count;
    }
    // predict the C13 chemical shifts of the molecule from MNOVA
    var c13predictions = predictCarbon_eeh();
    spectra["c13predictions"] = c13predictions;

    // Example usage:
    var path = doc.name;
    var result = splitPathAndFilename(path);


    for (i = 0, pageCount = doc.pageCount(); i < pageCount; i++) {
        page = doc.page(i);
        for (j = 0, itemCount = page.itemCount(); j < itemCount; j++) {
            spec = new NMRSpectrum(page.item(j));
      
            if( spec.isValid() ){
      	        var spectitle = spec.title + "_" + i;
                // remove any spaces from the title
                spectitle = spectitle.replace(/\s/g, "_");
                spectra[spectitle] = {};
                spectra[spectitle]["datatype"] = "nmrspectrum";
                spectra[spectitle]["origin"] = spec.originalFormat;
                spectra[spectitle]["type"] = spec.type;
                spectra[spectitle]["subtype"] = spec.subtype;
                spectra[spectitle]["experimenttype"] = spec.experimentType;
                spectra[spectitle]["experiment"] = spec.getParam("Experiment");
                spectra[spectitle]["class"] = spec.getParam("Class");
                spectra[spectitle]["spectype"] = spec.getParam("Spectrum Type");
                spectra[spectitle]["pulsesequence"] = spec.getParam("Pulse Sequence");
                spectra[spectitle]["intrument"] = spec.getParam("Instrument");
                spectra[spectitle]["probe"] = spec.getParam("Probe");
                spectra[spectitle]["datafilename"] = spec.getParam("Data File Name");
                spectra[spectitle]["nucleus"] = spec.getParam("Nucleus");
                spectra[spectitle]["specfrequency"] = spec.getParam("Spectrometer Frequency");
                // add solvent
                spectra[spectitle]["solvent"] = spec.getParam("Solvent");
                // add temperature
                spectra[spectitle]["temperature"] = spec.getParam("Temperature");

                // double check the type, if type is 1D and the nucleus has a comma in it then it is a 2D spectrum then set type to "2D"
                if( spectra[spectitle]["type"] == "1D" && spectra[spectitle]["nucleus"].indexOf(",") != -1 ){
                    spectra[spectitle]["type"] = "2D";
                }
        
                // process multiplets
        
                var multiplets = spec.multiplets();
                
                spectra[spectitle]["multiplets"] = {};
                spectra[spectitle]["multiplets"]["datatype"] = "multiplets";
                spectra[spectitle]["multiplets"]["normValue"] = multiplets.normValue;
        
                // loop through multiplets and add information
                var good_multiplets = 0;
                spectra[spectitle]["multiplets"]["data"] = {};
                for( var m=0; m<multiplets.count; m++){
                    var multiplet = multiplets.at(m);
                    if (multiplet.type == 0){
                        spectra[spectitle]["multiplets"]["data"][m] = {};
                        spectra[spectitle]["multiplets"]["data"][m]["delta1"] = multiplet.delta;
                        spectra[spectitle]["multiplets"]["data"][m]["nH"] = multiplet.nH;
                        spectra[spectitle]["multiplets"]["data"][m]["realH"] = multiplet.realH;
                        spectra[spectitle]["multiplets"]["data"][m]["integralValue"] = multiplet.integralValue();
                        spectra[spectitle]["multiplets"]["data"][m]["category"] = multiplet.category;
                        spectra[spectitle]["multiplets"]["data"][m]["type"] = multiplet.type;
            
                        var jlist = multiplet.jList();
                        spectra[spectitle]["multiplets"]["data"][m]["jlistcount"] = jlist.count;
            
                        var jvalslist = [];
            
                        for( q=0; q<jlist.count; q++){
                            jvalslist.push(jlist.at(q));
                        }
                        spectra[spectitle]["multiplets"]["data"][m]["jvals"] = jvalslist;
                        good_multiplets++;
                    }
                              
                }
                spectra[spectitle]["multiplets"]["count"] = good_multiplets;

                // loop over peaks in spectrum and add information
                var peaks = spec.peaks()
                spectra[spectitle]["peaks"] = {};
                spectra[spectitle]["peaks"]["datatype"] = "peaks";
                spectra[spectitle]["peaks"]["data"] = {};
                var good_pks = 0;
                for( var p=0; p<peaks.count; p++ ){
                    var pk = peaks.at(p);
                    if((pk.type == 0) & (pk.flagsToString() == "None" )){
                        spectra[spectitle]["peaks"]["data"][p] = {};
                        spectra[spectitle]["peaks"]["data"][p]["delta1"] = pk.delta(1);
                        spectra[spectitle]["peaks"]["data"][p]["delta2"] = pk.delta(2);
                        spectra[spectitle]["peaks"]["data"][p]["intensity"] = pk.intensity;
                        spectra[spectitle]["peaks"]["data"][p]["type"] = pk.type;
                        //   add annotation to peaks
                        spectra[spectitle]["peaks"]["data"][p]["annotation"] = pk.annotation;

                        good_pks++; 
                    }
                    
                }
                spectra[spectitle]["peaks"]["count"] = good_pks;

                // loop over integrals in spectrum and add information
                var integrals = spec.integrals();
                spectra[spectitle]["integrals"] = {};
                spectra[spectitle]["integrals"]["datatype"] = "integrals"
                spectra[spectitle]["integrals"]["count"] = integrals.count;  
                spectra[spectitle]["integrals"]["normValue"] = integrals.normValue;
                spectra[spectitle]["integrals"]["type"] = integrals.type;    
            
                spectra[spectitle]["integrals"]["data"] = {};
                for( p=0; p<integrals.count; p++){
                    var integral = integrals.at(p);
                    spectra[spectitle]["integrals"]["data"][p] = {};
                    spectra[spectitle]["integrals"]["data"][p]["intensity"] = integral.integralValue();
                    spectra[spectitle]["integrals"]["data"][p]["rangeMin1"] = integral.rangeMin(1);
                    spectra[spectitle]["integrals"]["data"][p]["rangeMin2"] = integral.rangeMin(2);
                    spectra[spectitle]["integrals"]["data"][p]["rangeMax1"] = integral.rangeMax(1);
                    spectra[spectitle]["integrals"]["data"][p]["rangeMax2"] = integral.rangeMax(2);
                    spectra[spectitle]["integrals"]["data"][p]["delta1"] = (integral.rangeMin(1) + integral.rangeMax(1))/2;
                    spectra[spectitle]["integrals"]["data"][p]["delta2"] = (integral.rangeMin(2) + integral.rangeMax(2))/2;
                    spectra[spectitle]["integrals"]["data"][p]["type"] = 0;

                }
            }        
        }
    } 

//   // create a list of spectra to send to the server if there is data present
    const keys_to_ignore = ["smiles", 
                            "molfile", 
                            "nmrAssignments", 
                            "c13predictions", 
                            "carbonAtomsInfo", 
                            "allAtomsInfo",
                             "hostname"];


    // var brukerExpts = new brukerExperiments();

    // search through the spectra using type and subtype and pulsesequence to see which 
    // experiments are present in the document

    var spectra_with_peaks = [];

    print("\nfinding spectra with peaks\n")
    for ( var spectrum in spectra){
        
        if( keys_to_ignore.indexOf(spectrum) == -1 ){
        // if( spectrum !== "smiles" && spectrum !== "molfile" && spectrum !== "nmrAssignments" ){
            // if there are no peaks picked then skip the spectrum
            if( spectra[spectrum]["peaks"]["count"] == 0 ){
                continue;
            }
            const type = spectra[spectrum]["type"];
            const subtype = spectra[spectrum]["subtype"].split(",")[0];
            const pulsesequence = spectra[spectrum]["pulsesequence"];
            const nucleus = spectra[spectrum]["nucleus"];
            const experiment = spectra[spectrum]["experiment"];
            spectra_with_peaks.push(nucleus + " " + experiment + " " + pulsesequence + " " + spectrum);
        }
    }

    // read in saved dialog values from file
    var dialogParametersJsonFilename = result.directoryPath + "/dialogPrameters.json";
    var dialogParmetersJsonFilenamePath = new Dir(dialogParametersJsonFilename);
    var dialogParameters = null;

    if (dialogParmetersJsonFilenamePath.fileExists(dialogParametersJsonFilename)){
        var fin = new File(dialogParametersJsonFilename);
        if (fin.open(File.ReadOnly)) {
            var sin = new TextStream(fin, 'UTF-8');
            
            var jsonStr = sin.readAll();
            dialogParameters = JSON.parse(jsonStr);
            
           
        }
        fin.close();
    }
    else{
        dialogParameters = {}
    }

    dialogParameters["ml_consent"] = ml_consent;

    print("dialogParameters\n", dialogParameters);

    // read in experiment identifiers from file expt_identifiers
    var exptIdentifiersJsonFilename = result.directoryPath + "/exptIdentifiers.json";
    var exptIdentifiersJsonFilenamePath = new Dir(exptIdentifiersJsonFilename);
    var expt_identifiers = [];

    if (exptIdentifiersJsonFilenamePath.fileExists(exptIdentifiersJsonFilename)){
        var fin = new File(exptIdentifiersJsonFilename);
        if (fin.open(File.ReadOnly)) {
            var sin = new TextStream(fin, 'UTF-8');
            
            var jsonStr = sin.readAll();
            expt_identifiers = JSON.parse(jsonStr);
            print("server", server);
            var chosen_spectra = xy3_dialog2( spectra_with_peaks, 
                                              expt_identifiers, 
                                              dialogParameters, 
                                              result["name"], 
                                              server );
           
        }
        fin.close();
    }
    else{
        print("server", server);
        var chosen_spectra = xy3_dialog2( spectra_with_peaks, 
                                          expt_identifiers, 
                                          dialogParameters,  
                                          result["name"], 
                                          server );
    }

    if( chosen_spectra === undefined ){
        MessageBox.warning("No spectra chosen");
        return;
    }



    print("chosen_spectra\n", chosen_spectra);

    spectra["chosenSpectra"] = {};
    spectra["chosenSpectra"]["datatype"] = "chosenSpectra";
    spectra["chosenSpectra"]["count"] = chosen_spectra["checked_spectra"].length;
    spectra["chosenSpectra"]["data"] = {};
    for(var i=0; i<chosen_spectra["checked_spectra"].length; i++){
        spectra["chosenSpectra"]["data"][i] = chosen_spectra["checked_spectra"][i];
    }




    // loop through the chosen spectra split the string and createa a dictionary of the last two items 
    // with the last item as the key
    var chosen_spectra_ids = {};
    var spectra_inrements = {};
    // set spectra_Increments keys to experiment types and values to -1
    for( var i=0; i<spectra_keys.length; i++){
        spectra_inrements[spectra_keys[i]] = -1;
    }
    var skip_increment = 0
    for( var i=0; i<chosen_spectra["checked_spectra"].length; i++){
        var split_str = chosen_spectra["checked_spectra"][i].split(" ");
        var key = split_str[split_str.length-1];
        var value = split_str[split_str.length-2];
        print("key ", key, " value ", value);
        // if the key is in the spectra_increments then increment the value by 1
        if( spectra_inrements[key] !== undefined ){
            spectra_inrements[key] += 1;
            key = key + "_" + spectra_inrements[key];
            print("\tkey ", key, " value ", value);
            chosen_spectra_ids[key] = value;
        }
        else{
            // Big error if key is not in spectra_increments
            // warn user and end script
            MessageBox.warning("Error: Spectrum type " + key + " unknown. Please check the spectra in the document.");
            return;
        }
    }

    // check if HSQC is present in the chosen_spectra_ids2
    // loop through the keys and in chosen_spectr_ids2, split the key by "_" and check if the first part is "HSQC"
    var foundHSQC = false;
    for( var key in chosen_spectra_ids){
        var split_key = key.split("_");
        if( split_key[0] == "HSQC" ){
            foundHSQC = true;
            break;
        }
    }
    //  warn the user and end the script if HSQC is not found
    if( !foundHSQC ){
        MessageBox.warning("No HSQC found in the chosen spectra. Please choose a HSQC spectrum.");
        return;
    }


    // loop through the chosen spectra and replace the key with the value in chosen_spectra_ids
    for( var key in chosen_spectra_ids){
        spectra[key] = spectra[chosen_spectra_ids[key]];
        spectra[key]["filename"] = chosen_spectra_ids[key];
    }



    for( var key in chosen_spectra_ids){
        // split the key by "_" and check
        var split_key = key.split("_");
        if( split_key[0] == "HSQC" || split_key[0] == "HMBC" || split_key[0] == "COSY" || split_key[0] == "NOESY" || split_key[0] == "HSQC_CLIPCOSY" || split_key[0] == "DDEPT_CH3_ONLY"){
            spectra[key]["type"] = "2D";
        }
    }


    // loop through the chosen spectra and remove the key from the spectra object
    for( var key in chosen_spectra_ids){
        delete spectra[chosen_spectra_ids[key]];
    }

    // print the chosen spectra ids
    for( var key in chosen_spectra_ids){
        print("chosen_spectra_ids ", key);
    }

    // with the hsqc data replace the intensity from the peaks with the integrals intensity values.
    // only do it if the number of integrals equals the number of peaks
    // match the peaks to the integrals by the delta values being in the range of the integrals in the two dimensions
    // if there is a match then replace the intensity value of the peak with the integral value
    
    // assign H1_1D to a variable if it exists
    if( spectra["H1_1D"] !== undefined ){
        var h1_1d = spectra["H1_1D"];
    }
    else{
        var h1_1d = undefined;
    }

    for( var key in spectra){
        var split_key = key.split("_");

        // check if the first par of the key is HSQC and the second part is numeric
        // if it is not then continue to the next key
        if( split_key[0] != "HSQC" || split_key.length != 2 ){
            continue;
        }
        
        var hsqc = spectra[key];
        var pks = hsqc["peaks"]["data"];
        var pks_count = hsqc["peaks"]["count"];
        var integrals = hsqc["integrals"]["data"];
        var integrals_count = hsqc["integrals"]["count"];

        if (pks_count == integrals_count){
            for ( var i=0; i<pks_count; i++){
                var pk = pks[i];
                for( var j=0; j<integrals_count; j++){
                    var integral = integrals[j];

                    if( (pk["delta1"] >= integral["rangeMax1"]) && 
                        (pk["delta1"] <= integral["rangeMin1"]) && 
                        (pk["delta2"] >= integral["rangeMin2"]) && 
                        (pk["delta2"] <= integral["rangeMax2"])){

                        pk["intensity"] = integral["intensity"];
                        break;
                    }
                }
            }
        }
        
        else if( pks_count != integrals_count){
            MessageBox.warning("Number of peaks and integrals do not match in HSQC data. Therefore using peaks only. Program works best with integrals");
        }
    }

    spectra["MNOVAcalcMethod"] = {};
    spectra["MNOVAcalcMethod"]["datatype"] = "MNOVAcalcMethod";
    spectra["MNOVAcalcMethod"]["count"] = 1;
    spectra["MNOVAcalcMethod"]["data"] = {};
    spectra["MNOVAcalcMethod"]["data"]["0"] = chosen_spectra["calcSimpleMNOVA"];

    spectra["exptIdentifiers"] = {};
    spectra["exptIdentifiers"]["count"] = chosen_spectra["exptIdentifiers"].length;
    spectra["exptIdentifiers"]["datatype"] = "exptIdentifiers";
    spectra["exptIdentifiers"]["data"] = {};
    for( var i=0; i<chosen_spectra["exptIdentifiers"].length; i++ ){
        spectra["exptIdentifiers"]["data"][i] = chosen_spectra["exptIdentifiers"][i];
    }

    // save "exptIdentifiers" as a json file
    var fout = new File(exptIdentifiersJsonFilename);
    if (fout.open(File.WriteOnly)) {
        sout = new TextStream(fout, 'UTF-8');
        sout.writeln(JSON.stringify(chosen_spectra["exptIdentifiers"]));
    }
    fout.close();

    // save dialog parameters as a json file
    // parameters include the simulated annealing parameters,
    // the exptIdentifiers and the calcSimpleNMR and calcSimpleMNOVA

    var dialogParams = {};
    dialogParams["simulatedAnnealing"] = chosen_spectra["simulatedAnnealing"];

    dialogParams["calcSimpleNMR"] = chosen_spectra["calcSimpleNMR"];
    dialogParams["calcSimpleMNOVA"] = chosen_spectra["calcSimpleMNOVA"];

    var fout = new File(dialogParametersJsonFilename);
    if (fout.open(File.WriteOnly)) {
        sout = new TextStream(fout, 'UTF-8');
        sout.writeln(JSON.stringify(dialogParams));
    }

    spectra["workingDirectory"] = {};
    spectra["workingDirectory"]["datatype"] = "workingDirectory";
    spectra["workingDirectory"]["count"] = 1;
    spectra["workingDirectory"]["data"] = {};

    spectra["workingDirectory"]["data"]["0"] = result["directoryPath"];

    spectra["workingFilename"] = {};
    spectra["workingFilename"]["datatype"] = "workingFilename";
    spectra["workingFilename"]["count"] = 1;
    spectra["workingFilename"]["data"] = {};
    spectra["workingFilename"]["data"]["0"] = result["name"];

    // add simulated annealing to the spectra object
    spectra["simulatedAnnealing"] = {};
    spectra["simulatedAnnealing"]["datatype"] = "simulatedAnnealing";
    spectra["simulatedAnnealing"]["count"] = 1;
    spectra["simulatedAnnealing"]["data"] = {};
    spectra["simulatedAnnealing"]["data"]["0"] = chosen_spectra["simulatedAnnealing"];

    //add the ml_consent to the spectra object
    spectra["ml_consent"] = {};
    spectra["ml_consent"]["datatype"] = "ml_consent";
    spectra["ml_consent"]["count"] = 1;
    spectra["ml_consent"]["data"] = {};
    spectra["ml_consent"]["data"]["0"] = chosen_spectra["ml_consent"];


    spectra["spectraWithPeaks"] = {};
    spectra["spectraWithPeaks"]["datatype"] = "spectraWithPeaks";
    spectra["spectraWithPeaks"]["count"] = spectra_with_peaks.length;
    spectra["spectraWithPeaks"]["data"] = {};
    for( var i=0; i<spectra_with_peaks.length; i++){
        spectra["spectraWithPeaks"]["data"][i] = spectra_with_peaks[i];

    }

    if( isObjectEmpty(spectra) ){
        MessageBox.warning("No spectra found");
    }
    else{
        // save spectra JSON string to file in data directory
        var jsonfilename = result.directoryPath + "/" + result.name + "_assignments_mresnova.json";
    
        // save spectra json string to file
        var json_spectra_string = JSON.stringify(spectra,null,4);
        var fout1 = new File(jsonfilename);
        if (fout1.open(File.WriteOnly)) {
            sout1 = new TextStream(fout1, 'UTF-8');
            sout1.writeln(json_spectra_string);
        }
        fout1.close();        

        // post json string to server

        ////////////////////
        // Do some checking
        ////////////////////

        // check if there are any assignments if the user has chosen to use MNOVA manually assigned
        if( spectra["nmrAssignments"]["count"] == 0 && chosen_spectra["calcSimpleMNOVA"] == "MNova Manually Assigned"){
            MessageBox.warning("No MNOVA assignment data found");  
            return;
        }

        // check if there are any C13 predictions if the user has chosen to use MNOVA Predict
        if( spectra["c13predictions"]["count"] == 0 && chosen_spectra["calcSimpleMNOVA"] == "MNOVA Predict"){
            MessageBox.warning("No C13 predictions found");  
            return;
        }

        // check if there are no integrals and no peaks in the HSQC spectrum
        if( hsqc["peaks"]["count"] == 0 && hsqc["integrals"]["count"] == 0){
            MessageBox.warning("No peaks or integrals found in HSQC spectrum");  
            return;
        }

        // warn the user that there are no integrals in the HSQC spectrum
        if(  hsqc["integrals"]["count"] == 0){
            MessageBox.warning("No integrals found in HSQC spectrum, using peaks only. Program works best with integrals");  
        }

        // check if the number of peaks in hsqc is == number of multiplets in h1_1d
        if( h1_1d !== undefined ){
            if( hsqc["peaks"]["count"] != h1_1d["multiplets"]["count"]){
                MessageBox.warning("Number of peaks in HSQC not equal to number of multiplets in H1_1D. Do not use H1_1D");  
                return;
            }
        }

        // //////////////////////////////////////
        // // choose which server we are using
        // /////////////////////////////////////

        // var server = server_address();

        ///////////////////////////////////////
        // Decide which entry point to call
        //////////////////////////////////////

        if(chosen_spectra["calcSimpleMNOVA"] == "MNOVA Predict"){
            // build entry point
            entry_point = server + "simpleMNOVA"
        }
        else if (chosen_spectra["calcSimpleMNOVA"] == "NMRSHIFTDB2 Predict"){
            // build entry point
            entry_point = server + "simpleMNOVA"
        }
        else if (chosen_spectra["calcSimpleMNOVA"] == "MNova Manually Assigned"){
            // build entry point
            if( spectra["nmrAssignments"]["count"] == 0 ){
                MessageBox.warning("No MNOVA assignment data found");  
                return;
            }
            entry_point = server + "simpleMNOVA"
        }

        progressDialog = new ProgressDialog();
        progressDialog.labelText = "Performing Prediction ...";
        progressDialog.setRange(0, 100);
        progressDialog.value = 50;
        progressDialog.visible = true;


        var rtn = web_utils.jsonRequest(entry_point, json_spectra_string, "", "", false);

        progressDialog.visible = false;


        // Check if response is registration redirect
        try {
            var responseObj = JSON.parse(rtn.response || rtn.allStrOutput);
            
            if (responseObj.status === "unregistered" || responseObj.status === "registration_expired") {
                // Handle unregistered host
                // display a message box
                if( responseObj.status === "registration_expired") {
                    MessageBox.warning("Host registration has expired. Please reregister to use this feature.");
                }
                else if (responseObj.status === "unregistered") {
                    MessageBox.warning("Host is not registered. Please register to use simpleNMR prediction.");
                }
                
                if (typeof Application.openUrl === 'function') {
                    Application.openUrl(responseObj.registration_url);
                } else {
                    print("Application.openUrl is not a function");
                }
                
                return; // Exit early
            }
            
            // If we get here, response is normal JSON
            htmloutputstr = responseObj.html || responseObj.content;
        } catch (e) {
            // Not JSON or parsing error - continue with original response handling
            // replace the dummy title in the html output with result.name
            // MessageBox.warning("No JSON response found. Please check the server response.");
            htmloutputstr = rtn.response.replace("dummy_title", result.name);
        }

        MessageBox.warning("Prediction complete. Please check the output file.");

        // Continue with file saving logic
        var htmlfilename = result.directoryPath + "/html/" + result.name + "_d3.html";

        // create directory if it does not exist
        var dir = new Dir(result.directoryPath + "/html");
        if (!dir.exists) {
            dir.mkdir(result.directoryPath + "/html");
        }

        var fout = new File(htmlfilename);
        if (fout.open(File.WriteOnly)) {
            sout = new TextStream(fout, 'UTF-8');
            sout.write(htmloutputstr);
        }
        fout.close();

        // open html file in browser
        if (typeof Application.openUrl === 'function') {
            if (SysInfo.isWin) {
                Application.openUrl(htmlfilename);
            }
            else if (SysInfo.isMac) {
                Application.openUrl("file://" + htmlfilename);
            }
        }
        else {
            print("Application.openUrl is not a function");
        }

      
    }
}


if (this.MnUi && MnUi.simpleNMRtools) {
	MnUi.simpleNMRtools.simpleNMRtools_simplePREDICT_eeh = simplePREDICT_eeh;
}