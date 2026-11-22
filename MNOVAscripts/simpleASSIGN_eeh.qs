/*globals Application, simpleUtils, server_address, WebUtilsQT, MoleculePlugin, Molecule,
          NMRAssignments, SysInfo, MessageBox, print, MnUi*/

// simpleNMR assignment: read the assignments made in the simpleNMR HTML page
// (<document>_assignments_from_simplemnova.json), put the carbon and proton shifts on the
// molecule in Mnova, and ask the server for the final HTML page.
// Shared code is in simpleUtils.qs.

function simpleASSIGN_eeh() {

    var su = new simpleUtils();
    var doc = Application.mainWindow.activeDocument;

    // server, and is this computer registered?
    var server = server_address();
    var web_utils = new WebUtilsQT();
    if (su.checkRegistration(server, web_utils) === undefined) {
        return;
    }

    var result = su.splitPathAndFilename(doc.name);

    // the assignments exported from the HTML page
    var jsonfilename = result.directoryPath + "/" + result.name + "_assignments_from_simplemnova.json";
    print(jsonfilename);
    var jsondata = su.readJson(jsonfilename, undefined);
    if (jsondata === undefined) {
        print("Could not open json file");
        print(jsonfilename);
        MessageBox.warning("Could not open json file");
        return;
    }
    print("Reading json file");

    var carbonatoms_orig = JSON.parse(JSON.stringify(jsondata["nodes_orig"]));
    var carbonatoms_moved = JSON.parse(JSON.stringify(jsondata["nodes_now"]));
    var oldjsondata = JSON.parse(JSON.stringify(jsondata["oldjsondata"]));

    // replace the atom_idx in carbons_moved with the atom_idx in carbons_orig using the atomNumber to match
    var carbonatoms = [];
    carbonatoms_moved.forEach(function (atom) {
        var atomNumber = atom["atomNumber"];
        for (var i = 0; i < carbonatoms_orig.length; i++) {
            if (carbonatoms_orig[i]["atomNumber"] == atomNumber) {
                atom["id"] = carbonatoms_orig[i]["id"];
                carbonatoms.push(atom);
                break;
            }
        }
    });

    // the molecule in the current document
    var moleculePlugin = new MoleculePlugin(Application.molecule);
    if (moleculePlugin === undefined) {
        MessageBox.warning("Could not find a molecule in the current document");
        return;
    }
    var assign = new NMRAssignments(moleculePlugin.activeMolecule());
    var mol = new Molecule(moleculePlugin.activeMolecule());
    if (mol.isValid() === false) {
        MessageBox.warning("Could not create a mol object from the molecule in the current document");
        throw "Invalid Molecule";
    }

    // put the shifts on the molecule
    carbonatoms.forEach(function (atom) {
        var atomNumber = parseInt(atom["atomNumber"]);
        var c13_ppm = parseFloat(atom["ppm"]);
        var numProtons = parseInt(atom["numProtons"]);

        // atom index starts from 1 in Mnova (for all data sources)
        print("jsondata[\"dataFrom\"] = ", jsondata["dataFrom"]);
        var atom_idx = parseInt(atom["id"]) + 1;
        print(atom_idx, atomNumber, c13_ppm);

        assign.setChemShift(atom_idx, 0, c13_ppm);

        if (!(numProtons == 1 || numProtons == 2 || numProtons == 3)) {
            return;                                         // no attached protons to assign
        }
        if (atom["H1_ppm"] === undefined) {
            print("No protons found for atom ", atomNumber);
            return;
        }
        var h1_ppm = parseFloat(atom["H1_ppm"][0]);
        if (numProtons == 1) {
            assign.setChemShift(atom_idx, 1, h1_ppm);
        }
        else if (numProtons == 3) {
            assign.setChemShift(atom_idx, 1, h1_ppm);
            assign.setChemShift(atom_idx, 2, h1_ppm);
            assign.setChemShift(atom_idx, 3, h1_ppm);
        }
        else {                                              // CH2
            var h2_ppm = (atom["H1_ppm"].length == 1) ? h1_ppm : parseFloat(atom["H1_ppm"][1]);
            assign.setChemShift(atom_idx, 1, h1_ppm);
            assign.setChemShift(atom_idx, 2, h2_ppm);
        }
    });

    // rebuild the export with the new assignments (saved next to the document)
    var spectra = {};
    spectra["smiles"] = oldjsondata["smiles"];
    spectra["molfile"] = oldjsondata["molfile"];
    spectra["carbonAtomsInfo"] = oldjsondata["carbonAtomsInfo"];
    spectra["hostname"] = su.hostnameBlock();

    // re-read the molecule from the document, as it now carries the assignments
    mol = new Molecule(moleculePlugin.activeMolecule());
    if (mol.isValid() === false) {
        throw "Invalid Molecule";
    }
    spectra["nmrAssignments"] = su.nmrAssignmentsBlock(mol);
    print("nmrAssignments count ", spectra["nmrAssignments"]["count"]);

    spectra["c13predictions"] = oldjsondata["c13predictions"];

    // NOTE: exptIdentifiers holds labels such as "HSQC", while the spectra in oldjsondata are
    // keyed "HSQC_0", ... so nothing is copied here (unchanged behaviour; this file is only a
    // local record - the server is sent jsondata below).
    for (var idx = 0; idx < oldjsondata["exptIdentifiers"]["count"]; idx++) {
        var ky = oldjsondata["exptIdentifiers"]["data"][idx];
        print(ky);
        spectra[ky] = oldjsondata[ky];
    }

    spectra["exptIdentifiers"] = oldjsondata["exptIdentifiers"];
    spectra["MNOVAcalcMethod"] = oldjsondata["MNOVAcalcMethod"];
    spectra["MNOVAcalcMethod"]["data"][0] = "MNova Manually Assigned data";
    spectra["carbonCalcPositionsMethod"] = oldjsondata["carbonCalcPositionsMethod"];
    spectra["workingDirectory"] = oldjsondata["workingDirectory"];
    spectra["workingFilename"] = oldjsondata["workingFilename"];
    spectra["chosenSpectra"] = oldjsondata["chosenSpectra"];

    su.writeJson(result.directoryPath + "/" + result.name + "_assignments_mresnova.json", spectra, 4);

    // ask the server for the final HTML page
    if (spectra["nmrAssignments"]["count"] == 0) {
        print("No assignments found");
        MessageBox.warning("No MNOVA assignment data found");
        return;
    }
    var entry_point = server + "/simpleMNOVAfinalHTML";
    print("Post to the server");
    print("entry_point ", entry_point);
    var rtn = web_utils.jsonRequest(entry_point, JSON.stringify(jsondata, null, 4), "", "", false);
    print(rtn.url);

    var htmlfilename = su.writeHtml(rtn.response.replace("dummy_title", result.name),
                                    result.directoryPath, result.name);

    MessageBox.information("Molecule has been assigned and  the html file has been updated");

    su.openInBrowser(htmlfilename);
}

if (this.MnUi && MnUi.simpleNMRtools) {
	MnUi.simpleNMRtools.simpleNMRtools_simpleASSIGN_eeh = simpleASSIGN_eeh;
}
