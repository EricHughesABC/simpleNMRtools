/*
 * export.js
 * ─────────────────────────────────────────────────────────────────────────
 * Getting assignments out of the tool: reportAssignments (clipboard/alert
 * summary) and exportToMnova (the JSON file MNova reads back in). Both are
 * pure consumers of state.js's data — no dependency on view/nodes/links/
 * interactions/moves.
 */

    function reportAssignments() {
        // using the current nodes print out the assignments seperately for carbons and protons
        // create a string to print out in a single text block do not use carriege retruns to separate lines
        var carbonAssignments = "Carbon Assignments (ppm): ";
        var protonAssignments = "Proton Assignments (ppm): ";

        // sort nodes by atomNumber
        nodes.sort((a, b) => a.atomNumber - b.atomNumber);

        nodes.forEach(node => {
            carbonAssignments += `C${node.atomNumber} ${node.ppm.toFixed(2)}, `;
        });
        // replace the last comma with a full stop
        carbonAssignments = carbonAssignments.slice(0, -2) + ".\n";

        // create the proton assignments
        nodes.forEach(node => {
            if (node.numProtons > 0) {
                protonAssignments += `C${node.atomNumber} `;
                for (var i = 0; i < node.H1_ppm.length; i++) {
                    protonAssignments += `${node.H1_ppm[i].toFixed(2)}`;
                    if (i < node.H1_ppm.length - 1) {
                        protonAssignments += " & ";
                    }
                }
                protonAssignments += ", ";
            }
        });
        // replace the last comma with a full stop
        protonAssignments = protonAssignments.slice(0, -2) + ".\n";


        // display the assignments in an alert box
        alert(carbonAssignments + "\n" + protonAssignments + "\n\n(The assignments have also been copied to the clipboard)");

        // is it possible to copy the assignments to the clipboard
        var assignmentText = carbonAssignments + "\n" + protonAssignments;
        navigator.clipboard.writeText(assignmentText).then(function() {
        }, function(err) {
            console.error('Could not copy text: ', err);
        });

    }

function exportToMnova() {

    const jsonFilename = workingFilename + '_assignments_from_simplemnova.json';

    var nodes_export = JSON.parse(JSON.stringify(nodes));

    unscaleNodes(nodes_export);
    unscaleNodes(catoms_orig);
    unscaleNodes(catoms);
    unscaleNodes(nodes_orig);

    var nodes_data = {};
    nodes_data['nodes_orig']       = catoms_orig;
    nodes_data['nodes_now']        = nodes_export;
    nodes_data['links']            = links;
    nodes_data['smilesString']     = smilesString;
    nodes_data['molfile']          = molfile;
    nodes_data['dataFrom']         = dataFrom;
    nodes_data['oldjsondata']      = oldjsondata;
    nodes_data['molgraph']         = molgraph;
    nodes_data['shortest_paths']   = shortest_paths;
    nodes_data['svg']              = svg_bckgrnd_image_str;
    nodes_data['catoms']           = catoms;
    nodes_data['catoms_orig']      = catoms_orig;
    nodes_data['best_results']     = best_results;
    nodes_data['workingDirectory'] = workingDirectory;
    nodes_data['workingFilename']  = workingFilename;
    nodes_data['title']            = title;
    nodes_data['number_of_hmbc_cosy_subgraphs'] = number_of_hmbc_cosy_subgraphs;   

    // Copy working directory to clipboard. Fire-and-forget: do NOT await
    // this. Some embedding environments (e.g. a QWebEngineView with no
    // clipboard permission handler wired up) never resolve or reject this
    // promise, which would otherwise block everything below indefinitely
    // and silently prevent the export from ever happening.
    navigator.clipboard.writeText(workingDirectory).then(function () {
    }).catch(function (err) {
        console.warn('Clipboard write failed (non-fatal):', err);
    });

    // Trigger download immediately — no longer gated on the clipboard call above
    const data = "text/json;charset=utf-8," + encodeURIComponent(JSON.stringify(nodes_data, null, 2));
    const a = document.createElement('a');
    a.href = 'data:' + data;
    a.download = jsonFilename;
    a.click();

    scaleNodes(catoms_orig);
    scaleNodes(catoms);
    scaleNodes(nodes_orig);
}

