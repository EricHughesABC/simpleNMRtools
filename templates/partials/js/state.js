/*
 * state.js
 * ─────────────────────────────────────────────────────────────────────────
 * The application's data model: server-derived arrays (catoms, nodes,
 * links, mapping, ...) and the pure helpers that build/transform them.
 * The raw Jinja-rendered values these depend on (catoms, rawNodes, links,
 * best_results, etc.) are declared in the small inline <script> block in
 * the HTML template, which loads immediately before this file.
 *
 * findIsolatedNodes and the isolatedCarbons/isolatedCarbonsCount it derives
 * were moved here as a unit from where they originally sat, textually,
 * inside the Jinja data block — the call needs the function already
 * defined, and keeping the inline bootstrap script pure data (no logic)
 * keeps that block a straight list of server-value substitutions.
 */

    function findIsolatedNodes(nodes, links) {
        const connectedIds = new Set();
        
        links.forEach(link => {
            connectedIds.add(link.source);
            connectedIds.add(link.target);
        });
        
        return nodes.filter(node => 
            node.symbol === 'C' && !connectedIds.has(node.id)
        );
    }

    const isolatedCarbons = findIsolatedNodes(rawNodes, links);
    const isolatedCarbonsCount = isolatedCarbons.length;

    /*
     * Coordinate scales — normalised server space → SVG viewBox space
     * ─────────────────────────────────────────────────────────────────────────
     * The server stores all atom positions as normalised values in [0, 1].
     * These scales map them into the SVG viewBox coordinate space (0–1000,
     * 0–600) so D3 can place nodes and links at the correct positions
     * relative to the background molecule image.
     *
     * These scales are *constant* — they never change on resize.  The browser
     * handles all further scaling (viewBox → screen pixels) automatically via
     * preserveAspectRatio, so no JS recalculation is needed when the window
     * is resized or the device is rotated.  Attempting to recompute these
     * scales in screen pixels on resize would cause a mismatch between the
     * background SVG (which the browser scales in viewBox space) and the D3
     * overlay (which would then be in a different coordinate space).
     */
    const xScale = d3.scaleLinear().domain([0, 1]).range([0, svg_width_eeh]);
    const yScale = d3.scaleLinear().domain([0, 1]).range([0, svg_height_eeh]);

    // Scale an array of nodes from normalised [0,1] space into viewBox space (in-place)
    function scaleNodes(arr) {
        arr.forEach(node => {
            node.x = xScale(node.x);
            node.y = yScale(node.y);
        });
    }

    // Inverse: convert viewBox coordinates back to normalised [0,1] space (in-place)
    function unscaleNodes(arr) {
        arr.forEach(node => {
            node.x = xScale.invert(node.x);
            node.y = yScale.invert(node.y);
        });
    }

function create_tooltip_string(weight, mae, lae, lae_atomNumber, isolatedCarbonsCount, totalNodes, totalEdges, number_of_hmbc_cosy_subgraphs) {
    // if weight is between 0 and 5, use a green color for the weight value, if it is between 5 and 15,use an  orange colour, if it is above 15, use a red colour and in bold font  
    var weight_color = "green";
    var weight_fontWeight = "bold";
    var fontSize = "20px";
    if (weight > 0 && weight <= 5) {
        weight_color = "orange";
    } else if (weight > 5) {
        weight_color = "red";
    }
    // if the isolatedCarbonsCount is 0, use a green color for the isolatedCarbonsCount value, if it is between 1 and 3, use an orange colour, if it is above 3, use a red colour and in bold font
    var isolatedCarbonsColor = "green";
    var isolatedCarbonsFontWeight = "bold";
    if (isolatedCarbonsCount > 0 && isolatedCarbonsCount <= 3) {
        isolatedCarbonsColor = "orange";
    } else if (isolatedCarbonsCount > 3) {
        isolatedCarbonsColor = "red";
    }

    var number_of_hmbc_cosy_subgraphs_color = "green";
    var number_of_hmbc_cosy_subgraphs_fontWeight = "bold";
    if (number_of_hmbc_cosy_subgraphs > 2 && number_of_hmbc_cosy_subgraphs <= 4) {
        number_of_hmbc_cosy_subgraphs_color = "orange";
    } else if (number_of_hmbc_cosy_subgraphs == 0 || number_of_hmbc_cosy_subgraphs > 4) {
        number_of_hmbc_cosy_subgraphs_color = "red";
    }

    // highlight the MAE and LAE values depending on their values, if the MAE is between 0 and 1, use a green color for the MAE value, if it is between 1 and 2, use an orange colour, if it is above 2, use a red colour and in bold font. If the LAE is between 0 and 2, use a green color for the LAE value, if it is between 2 and 4, use an orange colour, if it is above 4, use a red colour and in bold font
    var mae_color = "green";
    var mae_fontWeight = "bold";
    if (mae > 1 && mae <= 2) {
        mae_color = "orange";
    } else if (mae > 2) {
        mae_color = "red";
    }
    // highlight the LAE values depending on their values
    var lae_color = "green";
    var lae_fontWeight = "bold";
    if (lae > 2 && lae <= 10) {
        lae_color = "orange";
    } else if (lae > 10) {
        lae_color = "red";
    }
    if ( totalEdges === 0 ||  totalEdges === undefined) {
        return `Correlation Penalty:&nbsp;<span style="color: red; font-weight:${weight_fontWeight}; font-size:${fontSize};">Unused</span><br>Uncorrelated&nbsp;Nodes:&nbsp;<span style="color:${isolatedCarbonsColor}; font-weight:${isolatedCarbonsFontWeight}; font-size:${fontSize};">${isolatedCarbonsCount}</span>&nbsp;out&nbsp;of&nbsp;${totalNodes}<br>MAE:&nbsp;<span style="color:${mae_color}; font-weight:${mae_fontWeight}; font-size:${fontSize};">${mae.toFixed(2)}</span>,&nbsp;LAE:&nbsp;<span style="color:${lae_color}; font-weight:${lae_fontWeight}; font-size:${fontSize};">${lae.toFixed(2)}</span>&nbsp;(NODE&nbsp;${lae_atomNumber})`;
    }
    else{
        return `Correlation Penalty:&nbsp;<span style="color:${weight_color}; font-weight:${weight_fontWeight}; font-size:${fontSize};">${weight}</span><br>Uncorrelated&nbsp;Nodes:&nbsp;<span style="color:${isolatedCarbonsColor}; font-weight:${isolatedCarbonsFontWeight}; font-size:${fontSize};">${isolatedCarbonsCount}</span>&nbsp;out&nbsp;of&nbsp;${totalNodes}<br>Sub-networks:&nbsp;<span style="color:${number_of_hmbc_cosy_subgraphs_color}; font-weight:${number_of_hmbc_cosy_subgraphs_fontWeight}; font-size:${fontSize};">${number_of_hmbc_cosy_subgraphs}</span><br>MAE:&nbsp;<span style="color:${mae_color}; font-weight:${mae_fontWeight}; font-size:${fontSize};">${mae.toFixed(2)}</span>,&nbsp;LAE:&nbsp;<span style="color:${lae_color}; font-weight:${lae_fontWeight}; font-size:${fontSize};">${lae.toFixed(2)}</span>&nbsp;(NODE&nbsp;${lae_atomNumber})`;
    }
}

    // Modify INFO tooltip to display the best results from the optimization with only 2 decimal places
    // d3.select("#info_tooltip").html(`Correlation Penalty:&nbsp;${best_results.best_weight}<br>Uncorrelated&nbsp;Nodes:&nbsp;${isolatedCarbonsCount}&nbsp;out&nbsp;of&nbsp;${nodes.length}<br>MAE:&nbsp;${best_results.best_mae.toFixed(2)},&nbsp;LAE:&nbsp;${best_results.best_lae.toFixed(2)}&nbsp;(NODE&nbsp;${best_results.best_lae_atomNumber})`);
    var info_tooltip_str = create_tooltip_string(best_results.best_weight, best_results.best_mae, best_results.best_lae, best_results.best_lae_atomNumber, isolatedCarbonsCount, rawNodes.length, links.length, number_of_hmbc_cosy_subgraphs);
    d3.select("#info_tooltip").html(info_tooltip_str);

    // Set the visible property of the catoms array to false
    // loop through the nodes
    catoms.forEach(node => {
        node.visible = false;
    });

    // Populate catoms with NMR data from the server-rendered node set
    rawNodes.forEach(node => {
        const catom = catoms.find(catom => catom.id === node.id);
        if (catom) {
            catom.visible      = true;
            catom.ppm          = node.ppm;
            catom.iupacLabel   = node.iupacLabel;
            catom.jCouplingVals = node.jCouplingVals;
            catom.jCouplingClass = node.jCouplingClass;
            catom.H1_ppm       = node.H1_ppm;
        }
    })

    var catoms_orig = JSON.parse(JSON.stringify(catoms));

    // Map all node arrays from normalised [0,1] space into SVG viewBox space
    scaleNodes(rawNodes);
    scaleNodes(nodes_orig);
    scaleNodes(catoms);
    scaleNodes(catoms_orig);

    // define mapping
    var mapping = {};
    for(var i=0; i<catoms_orig.length; i++) {
        mapping[catoms_orig[i]["id"]] = catoms_orig[i]["id"];
    }

    // Display nodes: the filtered, working set used for all D3 rendering.
    // Distinct from rawNodes (the raw server data used only during init above).
    var nodes = JSON.parse(JSON.stringify(catoms.filter(node => node.visible)));

