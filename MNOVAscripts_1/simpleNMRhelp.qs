
/*globals WebUtilsQT NMRAssignments xy3_dialog Dir brukerExperiments String iupac_dialog server_address MnUi*/

function simpleNMRhelp() {

    var server = server_address();
    var entry_point = server + "documentation"

    Application.openUrl(entry_point);

}

 if (this.MnUi && MnUi.simpleNMRtools) {
	MnUi.simpleNMRtools.simpleNMRtools_simpleNMRhelp = simpleNMRhelp;
}