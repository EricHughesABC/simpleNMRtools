function test_molecule()
{

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

    // iterate through pages and print out what is on the page
    var doc = Application.mainWindow.activeDocument;
    
    

    // return molecule from document
    var mol = getActiveMolecule(doc, Application.molecule);
    

    // get graphic properties of molecule
    var props = mol.graphicProperties();

    // print("props: " + props);

    print(mol.paramNames());

    for ( var i=1; i <= mol.atomCount; i++) {    
        var atom = new Atom(mol.atom(i));
        var j = i-1;

        var symAtoms = mol.symmetricalAtoms(i);


        if (symAtoms.length != 0) {

            var symAtom = new Atom(mol.atom(symAtoms[0]));
            print("atom idx " + i + ", atomNumber " + atom.number + ", sym atomNumber  " + symAtom.number);
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

