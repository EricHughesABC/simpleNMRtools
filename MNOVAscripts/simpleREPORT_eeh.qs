/* globals MnUi, print*/

function simpleREPORT_eeh() {

    var simpleutils = new simpleUtils();

    // Get the current document
    var doc = Application.mainWindow.activeDocument;

    // get the assignments

    // return molecule from document
    var mol = simpleutils.getActiveMolecule(doc, Application.molecule);

    var spectra = {};
    if( mol === undefined ){
        MessageBox.warning("No molecule found!")
        print("molecule is undefined");
        return;  
    }

    var assignments = NMRAssignments(mol);
    // print("hasAssignments ", assignments.hasAssignments());
    spectra["nmrAssignments"] = {};
    spectra["nmrAssignments"]["datatype"] = "nmrAssignments";
    spectra["nmrAssignments"]["count"] = 0;
    spectra["nmrAssignments"]["data"] = {};

    if ( assignments.hasAssignments() == false){
        MessageBox.warning("No assignments found!")
        print("No assignments found");
        return;
    }

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
                // spectra["nmrAssignments"]["data"][i]["iupacLabel"] = mnova_to_iupac[i][1];
                spectra["nmrAssignments"]["data"][i]["iupacLabel"] = "";
                count++;
                print("atom.label ", atom.Label);
                print("atom properties\n", atom);
            

                for( var k=1; k<=3; k++ ){
                    hshift = assignments.chemShift(i,k);	    
                    if( hshift === undefined ){	      
                        spectra["nmrAssignments"]["data"][i]["f2_ppm"+k] = "undefined";
                    }
                    else {	 
                        print("hshift ", hshift );     
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
    spectra["nmrAssignments"]["count"] = count;

    print("assignments\n");

    // copy the data to a new object
    var adata = JSON.parse(JSON.stringify(spectra["nmrAssignments"]["data"]));

    // remove the data abject if f1_ppm is undefined
    for( var i in adata ){
        if( adata[i]["f1_ppm"] === undefined ){
            delete adata[i];
        }
    }
    for( var i in adata ){
        print(i, adata[i]["f1_ppm"], adata[i]["atomNumber"]);
    }

    // // numbers
    // var numbers = [1,2,5,3,4,8,7];
    // // sort numbers
    // numbers.sort(function(a, b){return a-b});
    // print("numbers ", numbers);

    //sort adata based on f1_ppm
    var sorted = Object.keys(adata).sort(function(a, b) {
        return adata[b]["f1_ppm"] - adata[a]["f1_ppm"];
    });

    // sor the adat based on f1_ppm in descending order
    // var sorted = Object.keys(adata).sort(function(a, b) {


    // print("\nsorted.length\n", sorted.length);

    // print("\nsorted\n");
    // for( var i=0; i<sorted.length; i++ ){
    //     print(i + "\t" + adata[sorted[i]]["atom_idx"]+ "\t" + adata[sorted[i]]["atomNumber"] + "\t" + adata[sorted[i]]["f1_ppm"].toFixed(2));
    // }

    // print the f2_ppm values from the sorted array
    // only print values if they are defined and different 

    // for( var i=0; i<sorted.length; i++ ){
    //     // check if f1_ppm is undefined continue
    //     if( adata[sorted[i]]["f2_ppm1"] == "undefined" ){
    //         continue;
    //     }
    //     else{
    //         print(i, adata[sorted[i]]["f2_ppm1"], adata[sorted[i]]["numProtons"]);
    //     }
    //     print("\n");
    //     print(i + "\t" + adata[sorted[i]]["atom_idx"]+ "\t" + adata[sorted[i]]["atomNumber"] + "\t" + adata[sorted[i]]["f1_ppm"].toFixed(2));
        
    //     var numProtons = adata[sorted[i]]["numProtons"];
    //     if( numProtons == 1 || numProtons == 3 ){
    //         print("f2_ppm1 ", adata[sorted[i]]["f2_ppm1"].toFixed(2));
    //     }
    //     else if( numProtons == 2 ){
    //         var f2_ppm1 = adata[sorted[i]]["f2_ppm1"];
    //         var f2_ppm2 = adata[sorted[i]]["f2_ppm2"];

    //         if( f2_ppm1 != f2_ppm2 ){
    //             print("f2_ppm1 ", f2_ppm1.toFixed(2));
    //             print("f2_ppm2 ", f2_ppm2.toFixed(2));
    //         }
    //         else{
    //             print("f2_ppm1 ", f2_ppm1.toFixed(2));
    //         }
    //     }
    // }

    print("\n");
    print(adata[sorted[0]]);

    carbonReportString = "13C NMR (ppm): ";
    for( var i=0; i<sorted.length; i++ ){
        var carbon = adata[sorted[i]];
        carbonReportString +=  carbon["f1_ppm"].toFixed(2) + " (C" + carbon["atomNumber"] + "), ";
        if( i == sorted.length-1 ){
            // remove the last comma and space and replace with a period
            carbonReportString = carbonReportString.slice(0, -2) + ".";
        }
        
    }
    print("carbonReportString\n", carbonReportString);

    // create the protonReportString
    protonReportString = "1H NMR (ppm): ";
    
    // create a list  from adata keeping only items with numProtons > 0
    var protondata = JSON.parse(JSON.stringify(spectra["nmrAssignments"]["data"]));
    for( var i in protondata ){
        if( protondata[i]["numProtons"] == 0 ){
            delete protondata[i];
        }
    }
    // keep only the items where f1_ppm is defined
    for( var i in protondata ){
        if( protondata[i]["f1_ppm"] === undefined ){
            delete protondata[i];
        }
    }

    // keep only the items where f2_ppm1 is defined
    for( var i in protondata ){
        if( protondata[i]["f2_ppm1"] === "undefined" ){
            delete protondata[i];
        }
    }

    // sort protondata based on f1_ppm in descending order
    var sortedProtons = Object.keys(protondata).sort(function(a, b) {
        return protondata[b]["f2_ppm1"] - protondata[a]["f2_ppm1"];
    });

    // continue with the protonReportString
    for( var i=0; i<sortedProtons.length; i++ ){
        var proton = protondata[sortedProtons[i]];
        // add first proton f2_ppm1
        if (proton["f2_ppm1"] == "undefined"){
            continue;
        }
        protonReportString +=  proton["f2_ppm1"].toFixed(2);
        if( proton["numProtons"] == 2){
            if(proton["f2_ppm1"] != proton["f2_ppm2"]){

                protonReportString += ", " + proton["f2_ppm2"].toFixed(2);
            }
        }
        protonReportString += " (C" + proton["atomNumber"] + "), ";
        if( i == sortedProtons.length-1 ){
            // remove the last comma and space and replace with a period
            protonReportString = protonReportString.slice(0, -2) + ".";
        }

    }
    print("protonReportString\n", protonReportString);

    var reportText = carbonReportString + "\n\n" + protonReportString;

    Application.clipboard.text = reportText;

    MessageBox.information(reportText, "NMR Report");

    // save the report as plain text next to the Mnova document: <document>_NMR_Report.txt
    if (!doc.name) {
        MessageBox.warning("The document has not been saved yet, so the report was only copied to the clipboard.");
        return;
    }
    var result = simpleutils.splitPathAndFilename(doc.name);
    var reportFilename = result.directoryPath + "/" + result.name + "_NMR_Report.txt";
    if (simpleutils.writeText(reportFilename, reportText + "\n")) {
        print("Report saved to " + reportFilename);
    }
    else {
        MessageBox.warning("Could not save the report to\n" + reportFilename);
    }
}

if (this.MnUi && MnUi.simpleNMRtools) {
	MnUi.simpleNMRtools.simpleNMRtools_simpleREPORT_eeh = simpleREPORT_eeh;
}