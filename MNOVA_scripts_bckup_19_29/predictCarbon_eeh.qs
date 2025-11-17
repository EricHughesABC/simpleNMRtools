function predictCarbon_eeh() {
	'use strict';
	var mol = new Molecule(Application.molecule.activeMolecule());
	var predictionResult = undefined;

    var c13predictions = {};
    c13predictions["datatype"] = "c13predictions";
    c13predictions["count"] = 0;
    c13predictions["data"] = {}
	
	print( Application.NMRPredictor === undefined);

	if (Application.NMRPredictor === undefined) {
		print( "Predictor Plugin not loaded");
        return c13predictions
	}
	else {
		print("predictor plugin found");
	}					
	if (!mol.isValid()){
		print( "Invalid Molecule" );
        return c13predictions
	}
	else {
		print("Molecule is valid" );
	}

	predictionResult = mol.nmrPrediction("C13", true);
    print("Prediction Result:\n", predictionResult);

    // if the prediiction result is null or undefined, return

    if (predictionResult === null || predictionResult === undefined) {
    	print("Prediction Result is null or undefined");
    	return c13predictions;
    }

    //for (var i = 0; i < predictionResult.length; i++) {
    for (var i = 0; i < predictionResult.length; i++) {
        for (var j = 0; j < predictionResult[i].atom.length; j++) {
            print(i, ' ', predictionResult[i].atom[j].index, ' ', predictionResult[i].shift.value);
        }
    	// print(i, ' ', predictionResult[i].atom[0].index, ' ', predictionResult[i].shift.value);
    }


    // for ( var i=1; i <= mol.atomCount; i++) {    
    //     var atom = new Atom(mol.atom(i));
    //     print(i, " ", atom.elementSymbol, ' ', atom.number);            
    // }  
    
    //create a dictionary of the atom number to atom index
    var atomNumberToIndex = {};
    var atomIndexToNumber = {};
    var atomNumberToNumProtons = {};
    for ( i = 1; i <= mol.atomCount; i++) {    
        var atom = new Atom(mol.atom(i));
        atomNumberToIndex[atom.number] = i-1; // zero based index changed EEH 2025-sep-06
        atomIndexToNumber[i-1] = atom.number; // zero based index changed EEH 2025-sep-06
        atomNumberToNumProtons[atom.number] = atom.nH;
    }

    var k = 0;
    for (var i = 0; i < predictionResult.length; i++) {

        for (var j = 0; j < predictionResult[i].atom.length; j++) {
            print(k, ' ', predictionResult[i].atom[j].index, ' ', predictionResult[i].shift.value);
            c13predictions["data"][k] = {
                                    atomNumber:  predictionResult[i].atom[j].index, 
                                    atom_idx: atomNumberToIndex[predictionResult[i].atom[j].index],
                                    numProtons: atomNumberToNumProtons[predictionResult[i].atom[j].index], 
                                    ppm: parseFloat(predictionResult[i].shift.value.toFixed(3))};
            k++;
        }
    }
    c13predictions["count"] = k;


    // dump the c13predictions to a json string
    var c13predictionsStr = JSON.stringify(c13predictions, null, 4);

    return c13predictions


    //

}