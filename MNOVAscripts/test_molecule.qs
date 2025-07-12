function test_molecule()
{

    var spectra_keys = ["HSQC", 
                    "HMBC", 
                    "COSY", 
                    "NOESY", 
                    "H1_1D", 
                    "C13_1D", 
                    "DEPT135", 
                    "PureShift", 
                    "DDEPT_CH3_ONLY", 
                    "SKIP", 
                    "HSQC_CLIPCOSY"];

	function getActiveMolecule(aDocWin, aMolPlugin) {
	    var molec = aMolPlugin.activeMolecule();
		
		if (molec.isValid()) {
            print("Active molecule is valid.");
		    return molec;
		}
		if (aDocWin.itemCount("Molecule") == 1) {
            print("Only one molecule found in the document, returning it.");
			molec = new Molecule(aDocWin.item(0, "Molecule"));
			return molec;
		}
        else {
            print("No valid molecule found in the active document.");
            return molec;
        }
	}

    // iterate through pages and print out what is on the page
    var doc = Application.mainWindow.activeDocument;
    
    print("Application.molecule.activeMolecule().isValid()")
    print(Application.molecule.activeMolecule().isValid())

    // return molecule from document
    var mol = getActiveMolecule(doc, Application.molecule);

    // // check if the molecule is undefined
    // if (mol === undefined) {
    //     MessageBox.warning("No molecule found in the active document.");
    //     return;
    // }

    // if no molecule is found, return
    if (!mol.isValid()) {
        MessageBox.warning("<span style='color:red'>No valid molecule found in the active document.</span>");
        return;
    }
    

    // get graphic properties of molecule
    var props = mol.graphicProperties();

    // print("props: " + props);

    print(mol.paramNames());
    print("Molecule name: " + mol.name);
    print("Molecule ID: " + mol.id);
    print("Molecule atom count: " + mol.atomCount);

    var molHasSymmetry = false;
    // create a dictionary to hold symmetric atoms
    var symmetricAtoms = {};
    var symmetricCarbonAtoms = {};


    for ( var i=1; i <= mol.atomCount; i++) {    
        var atom = new Atom(mol.atom(i));
        var j = i-1;

        var symAtoms = mol.symmetricalAtoms(i);


        if (symAtoms.length != 0) {
            molHasSymmetry = true;

            var symAtom = new Atom(mol.atom(symAtoms[0]));
            print("atom idx " + i + ", atomNumber " + atom.number + ", sym atomNumber  " + symAtom.number);

            // add atom to symmetric atoms dictionary
            symmetricAtoms[i] = symAtoms;
            // if the atom is carbon, add it to the symmetric carbon atoms dictionary
            if (atom.elementSymbol == "C") {
                symmetricCarbonAtoms[i] = symAtoms;
            }

        

            // spectra["allAtomsInfo"]["data"][j] = {};
            // spectra["allAtomsInfo"]["data"][j]["atom_idx"] = j;  
            // spectra["allAtomsInfo"]["data"][j]["id"] = j;  
            // spectra["allAtomsInfo"]["data"][j]["atomNumber"] = atom.number;
            // spectra["allAtomsInfo"]["data"][j][["symbol"]] = atom.elementSymbol;
            // spectra["allAtomsInfo"]["data"][j]["numProtons"] = atom.nH; 
            // allAcomCount++;                      
        }
    }
    print("Molecule has symmetry: " + molHasSymmetry);

    if (molHasSymmetry) {
        print("Symmetric atoms: ");
        for (var key in symmetricAtoms) {
            if (symmetricAtoms.hasOwnProperty(key)) {
                print("Atom " + key + " has symmetric atoms: " + symmetricAtoms[key].join(", "));
            }
        }
    } else {
        print("No symmetric atoms found.");
    }

    // print out symmetric carbon atoms
    if (Object.keys(symmetricCarbonAtoms).length > 0) {
        print("Symmetric carbon atoms: ");
        for (var key in symmetricCarbonAtoms) {
            if (symmetricCarbonAtoms.hasOwnProperty(key)) {
                print("Carbon Atom " + key + " has symmetric atoms: " + symmetricCarbonAtoms[key].join(", "));
            }
        }
    } else {
        print("No symmetric carbon atoms found.");
    }

    // count the number of different carbon types: CH3, CH2, CH, C
    var carbonInformation = {};

    for (var i = 1; i <= mol.atomCount; i++) {
        var atom = new Atom(mol.atom(i));
        if (atom.elementSymbol == "C") {   
            var carbonType = "C";
            if (atom.nH == 3) {
                carbonType = "CH3";
            } else if (atom.nH == 2) {
                carbonType = "CH2";
            } else if (atom.nH == 1) {
                carbonType = "CH";
            }
            carbonInformation[i] = carbonType;
        }
    }
    print("Carbon types in the molecule: ");
    for (var key in carbonInformation) {
        if (carbonInformation.hasOwnProperty(key)) {
            print("Atom " + key + " is of type: " + carbonInformation[key]);
        }
    }
    print("Total number of carbon atoms: " + Object.keys(carbonInformation).length);
    print("Total number of atoms in the molecule: " + mol.atomCount);

    // count the number of carbonTypes in the molecule
    var carbonTypesCount = {};
    for (var key in carbonInformation) {
        if (carbonInformation.hasOwnProperty(key)) {
            var carbonType = carbonInformation[key];
            if (carbonTypesCount[carbonType]) {
                carbonTypesCount[carbonType]++;
            } else {
                carbonTypesCount[carbonType] = 1;
            }
        }
    }
    print("Count of carbon types in the molecule: ");
    for (var carbonType in carbonTypesCount) {
        if (carbonTypesCount.hasOwnProperty(carbonType)) {
            print(carbonType + ": " + carbonTypesCount[carbonType]);
        }
    }

    // get spectra_keys from simpleUtils file and print them
    var simpleutils = new simpleUtils();
    // var spectra_keys = simpleUtils.spectra_keys;
    for (var i = 0; i < spectra_keys.length; i++) {
        var key = spectra_keys[i];
        print("Spectrum key: " + key);
    }

    // identify the spectra in the document
    var spectra_lst = identify_spectrum();

    print("spectra_lst\n\t" + spectra_lst.length + " spectra found in the document.");

    for( var i=0; i<spectra_lst.length; i++){
        var spectrum = spectra_lst[i];
        print("spectrum.experimentEEH " + spectrum.experimentEEH + " " + spectrum.getParam("Pulse Sequence"));

    }

    // create a dictionary of NMR experiments to list of spectra indices
    var nmrExperiments = {};
    for (var i = 0; i < spectra_lst.length; i++) {
        var spectrum = spectra_lst[i];
        var experimentName = spectrum.experimentEEH;
        
        if (!nmrExperiments[experimentName]) {
            nmrExperiments[experimentName] = [];
        }
        nmrExperiments[experimentName].push(i);
    }

    print("NMR Experiments found in the document: ");
    for (var experiment in nmrExperiments) {
        if (nmrExperiments.hasOwnProperty(experiment)) {
            print("Experiment: " + experiment + ", Spectra indices: " + nmrExperiments[experiment].join(", "));
        }
    }

    // search for all the C13 spectra and count the number of peaks in each and report the number of peaks corresponding to each spectrum found
    var c13Spectra = [];
    if (nmrExperiments["C13_1D"]) {
        for (var i = 0; i < nmrExperiments["C13_1D"].length; i++) {
            var spectrumIndex = nmrExperiments["C13_1D"][i];
            var spectrum = spectra_lst[spectrumIndex];
            if (spectrum.subtype == "13C") {
                c13Spectra.push(spectrum);
                print("C13 Spectrum found: " + spectrum.name + ", Peaks count: " + spectrum.peaks().count);
            }
        }
    }

    // print out the number of carbon atoms in the molecule minus half the number of symmetric carbon atoms
    var carbonAtomCount = Object.keys(carbonInformation).length;
    var symmetricCarbonCount = Object.keys(symmetricCarbonAtoms).length;
    var effectiveCarbonCount = carbonAtomCount - Math.floor(symmetricCarbonCount / 2);
    print("Total number of carbon atoms in the molecule: " + carbonAtomCount);
    print("Number of symmetric carbon atoms: " + symmetricCarbonCount);
    print("Effective number of carbon atoms in the molecule: " + effectiveCarbonCount); 

    // what type are the symmetic carbon atoms?
    if (symmetricCarbonCount > 0) {
        print("Types of symmetric carbon atoms: ");
        for (var key in symmetricCarbonAtoms) {
            if (symmetricCarbonAtoms.hasOwnProperty(key)) {
                var atom = new Atom(mol.atom(key));
                var carbonType = carbonInformation[key];
                print("Atom " + key + " is of type: " + carbonType + ", Element Symbol: " + atom.elementSymbol);
            }
        }
    }

    // calculate the number of symmetric carbon atoms of each type
    var symmetricCarbonTypesCount = {};
    // set them to zero
    symmetricCarbonTypesCount["CH3"] = 0;
    symmetricCarbonTypesCount["CH2"] = 0;
    symmetricCarbonTypesCount["CH"] = 0;
    symmetricCarbonTypesCount["C"] = 0;
    for (var key in symmetricCarbonAtoms) {
        if (symmetricCarbonAtoms.hasOwnProperty(key)) {
            var atom = new Atom(mol.atom(key));
            var carbonType = carbonInformation[key];

            symmetricCarbonTypesCount[carbonType]++;

        }
    }
    // divide the count by 2 for each symmetric carbon atom
    print("Count of symmetric carbon types in the molecule: ");
    for (var carbonType in symmetricCarbonTypesCount) {
        if (symmetricCarbonTypesCount.hasOwnProperty(carbonType)) {
            print(carbonType + ": " + symmetricCarbonTypesCount[carbonType]/2);
        }
    }
    
    // search for all the HSQC spectra and count the number of peaks in each and report the number of peaks corresponding to each spectrum found
    var hsqcSpectra = [];
    if (nmrExperiments["HSQC"]) {
        for (var i = 0; i < nmrExperiments["HSQC"].length; i++) {
            var spectrumIndex = nmrExperiments["HSQC"][i];
            var spectrum = spectra_lst[spectrumIndex];
            hsqcSpectra.push(spectrum);
            print("HSQC Spectrum found: " + spectrum.name + ", Peaks count: " + spectrum.peaks().count);
            print("Integrals: " + spectrum.integrals().count);
        }
    }

    // count the number of positive integrals in the HSQC spectra and print out
    var totalPositiveIntegrals = 0;
    var totalNegativeIntegrals = 0;
    if (nmrExperiments["HSQC"]) {
        for (var i = 0; i < nmrExperiments["HSQC"].length; i++) {
            var spectrumIndex = nmrExperiments["HSQC"][i];
            var spectrum = spectra_lst[spectrumIndex];
            var integrals = spectrum.integrals();
            for (var j = 0; j < integrals.count; j++) {
                var integral = integrals.at(j);
                if (integral.integralValue() > 0) {
                    totalPositiveIntegrals++;
                }
                else if (integral.integralValue() < 0) {
                    totalNegativeIntegrals++;
                }
            }
        }
    }
    print("Total number of positive integrals in HSQC spectra: " + totalPositiveIntegrals);
    print("Total number of negative integrals in HSQC spectra: " + totalNegativeIntegrals);



    if (hsqcSpectra.length == 0) {
        MessageBox.warning("No HSQC spectra found in the document.");
        return;
    }

    //  output message as html in a message box
    var message = "<h2>Test Molecule Results</h2>";
    message += "<p>Molecule name: <strong>" + mol.name + "</strong></p>";
    message += "<p>Molecule ID: <strong>" + mol.id + "</strong></p>";
    message += "<p>Molecule atom count: <strong>" + mol.atomCount + "</strong></p>";
    message += "<p>Total number of carbon atoms: <strong>" + carbonAtomCount + "</strong></p>";

    message += "<p>Molecule has symmetry: <strong>" + molHasSymmetry + "</strong></p>";
    if (molHasSymmetry) {
        message += "<p>Symmetric atoms: <strong>" + Object.keys(symmetricAtoms).length + "</strong></p>";
        for (var key in symmetricAtoms) {
            if (symmetricAtoms.hasOwnProperty(key)) {
                message += "<p>Atom " + key + " has symmetric atoms: " + symmetricAtoms[key].join(", ") + "</p>";
            }
        }
    }

    if (Object.keys(symmetricCarbonAtoms).length > 0) {
        message += "<p>Symmetric carbon atoms: <strong>" + Object.keys(symmetricCarbonAtoms).length + "</strong></p>";
        for (var key in symmetricCarbonAtoms) {
            if (symmetricCarbonAtoms.hasOwnProperty(key)) {
                message += "<p>Carbon Atom " + key + " has symmetric atoms: " + symmetricCarbonAtoms[key].join(", ") + "</p>";
            }
        }
    }

// add carbon types information, ie the number of each type of carbon in the molecule

// add CH3, CH2, CH, C counts
    message += "<p>Count of carbon types in the molecule: </p>";
    message += "<ul>";
    for (var carbonType in carbonTypesCount) {
        if (carbonTypesCount.hasOwnProperty(carbonType)) {
            message += "<li>" + carbonType + ": " + carbonTypesCount[carbonType] + "</li>";
        }
    }
    message += "</ul>";

    // add experiment information in terms of the number of carbon peaks in the C13 spectra
    // the number of CH2 (-ve), combined total of CH3 and CH1 (+ve) peaks in HSQC spectra
    message += "<p>Number of C13 spectra found: <strong>" + c13Spectra.length + "</strong></p>";
    if (c13Spectra.length > 0) {
        message += "<p>Number of peaks in C13 spectra: </p>";
        message += "<ul>";
        for (var i = 0; i < c13Spectra.length; i++) {
            var spectrum = c13Spectra[i];
            message += "<li>" + spectrum.name + ": " + spectrum.peaks().count + " peaks</li>";
        }
        message += "</ul>";
    }
    // add HSQC spectra information
    message += "<p>Number of HSQC spectra found: <strong>" + hsqcSpectra.length + "</strong></p>";
    if (hsqcSpectra.length > 0) {
        message += "<p>Number of peaks in HSQC spectra: </p>";
        message += "<ul>";
        for (var i = 0; i < hsqcSpectra.length; i++) {
            var spectrum = hsqcSpectra[i];
            message += "<li>" + spectrum.name + ": " + spectrum.peaks().count + " peaks</li>";
        }
        message += "</ul>";
    }
    // add total number of positive and negative integrals in HSQC spectra
    message += "<p>Total number of positive integrals in HSQC spectra: <strong>" + totalPositiveIntegrals + "</strong></p>";
    message += "<p>Total number of negative integrals in HSQC spectra: <strong>" + totalNegativeIntegrals + "</strong></p>";    
    
    MessageBox.information(message, "Test Molecule Results");
    
}

