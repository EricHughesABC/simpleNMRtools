function test_mnovaC13predictions(){
	'use strict';
	var mol = new Molecule(Application.molecule.activeMolecule());
	var predictionResult = undefined;

    var c13predictions = {};
    c13predictions["datatype"] = "c13predictions";
    c13predictions["count"] = 0;
    c13predictions["data"] = {}

    predictionResult = mol.nmrPrediction("C13", false);
    print("Prediction Result:\n", predictionResult);
    
    print("predictionResult.length " + predictionResult.length );

    // for (var i = 0; i < predictionResult.length; i++) {
        
    //     for (var j = 0; j < predictionResult[i].atom.length; j++) {
    //         var atom_idx = predictionResult[i].atom[j].index;
    //         print("atom_idx " + atom_idx);
    //         var atom = new Atom(mol.atom(atom_idx));
    //         print(i + ' ' + predictionResult[i].atom[j].index + '\t' + atom.number + '\t' + predictionResult[i].shift.value.toFixed(2));
    //     }
    // }

    for (var i = 0; i < predictionResult.length; i++) {
        var atom_idx = predictionResult[i].atom.index;
        var atom = new Atom(mol.atom(atom_idx));
        print(predictionResult[i].atom.index + "\t" + atom.number + "\t" + atom.nH + "\t" + predictionResult[i].shift.value.toFixed(2));
        
        c13predictions["data"][i] = {
            "atom_idx": atom_idx-1, // zero based index changed EEH 2025-sep-06
            "atomNumber": atom.number,
            "numProtons": atom.nH,
            "shift": predictionResult[i].shift.value
        }
        c13predictions["count"] += 1;
    }
    print("c13predictions:\n", c13predictions);
}