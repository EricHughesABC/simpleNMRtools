/*
 * nodes.js
 * ─────────────────────────────────────────────────────────────────────────
 * Carbon-atom node creation and appearance: the node <g>/<circle>/<text>
 * elements, the experimental/calculated shift toggle (nodeFillColor,
 * nodePpmText, toggleCalculatedShifts, C keyboard shortcut target), and
 * the grey "not yet added" placeholder-node layer for un-assigned carbons.
 *
 * Depends on constants.js, state.js (nodes, catoms), moves.js
 * (dragstarted/dragged/dragended — must already be defined; see moves.js's
 * own comment) and interactions.js (handleMouseOver/handleMouseOut — these
 * are only referenced inside deferred arrow-function callbacks, so
 * interactions.js does NOT need to load before this file; nodes.js loads
 * before interactions.js instead, since interactions.js's own setup
 * queries the .node DOM elements this file creates).
 */

    /*
     * Calculated (predicted) chemical shift toggle
     * ─────────────────────────────────────────────────────────────────────────
     * Server-rendered nodes carry both an experimental shift (d.ppm) and a
     * calculated/predicted shift (d.ppm_calculated). When no calculated value
     * has been added yet for a node, ppm_calculated comes through as an empty
     * array ([]) rather than a number — see the same check used in the
     * tooltip builder below.
     *
     * "C" toggles the node labels between the two, and — while calculated
     * shifts are showing — nodes without a calculated value yet are drawn as
     * plain grey disks instead of their usual numProtons colour, flagging
     * them as "not yet added".
     */
    let showCalculatedPpm = false;

    function nodeHasCalculatedPpm(d) {
        return !(Array.isArray(d.ppm_calculated) && d.ppm_calculated.length === 0);
    }

    function nodeFillColor(d) {
        if (showCalculatedPpm && !nodeHasCalculatedPpm(d)) {
            return grey;
        }
        return colorScale(d.numProtons);
    }

    function nodePpmText(d) {
        if (showCalculatedPpm) {
            return nodeHasCalculatedPpm(d) ? d.ppm_calculated.toFixed(2) : "–";
        }
        return d.ppm.toFixed(2);
    }

    // Redraw every node's fill and ppm label to reflect the current showCalculatedPpm state
    function toggleCalculatedShifts() {
        showCalculatedPpm = !showCalculatedPpm;

        allnodes.each(function(d) {
            d3.select(this).select("circle").style("fill", nodeFillColor(d));
            d3.select(this).select(".node-text-ppm").text(nodePpmText(d));
        });

        // Show/hide grey placeholder nodes for carbons with no HSQC-assigned node
        renderCalcOnlyNodes();

        // Warn that the displayed shifts are predicted, not experimental
        d3.select("#calculated-shift-banner")
            .style("display", showCalculatedPpm ? "block" : "none");
    }


    const svg_display = d3.select('svg')
    // Create/select nodesGroup AFTER
    const nodesGroup = d3.select('svg').select('.nodes-group').empty()
        ? d3.select('svg').append('g').attr('class', 'nodes-group') 
        : d3.select('svg').select('.nodes-group');

    // Create the nodes
    const allnodes = nodesGroup.selectAll(".nodes-group")
      .data(nodes)
      .enter().append("g")
      .attr("class", "node")
      .attr("transform", d => `translate(${d.x},${d.y})`)
      .attr("atomNumber", d => d.atomNumber)
      .attr("id", d => d.id)
      .attr("ppm", d => d.ppm)
      .attr("ppm_calculated", d => d.ppm_calculated)
      .attr("numProtons", d => d.numProtons)
      .attr("x", d => d.x)
      .attr("y", d => d.y)
      .on("mouseenter", (event, d) => handleMouseOver(event, d))
      .on("mouseleave", (event, d) => handleMouseOut(event, d))
      .call(d3.drag()
        .on("start", dragstarted)
        .on("drag", dragged)
        .on("end", dragended));


    // Add circles to represent nodes
    allnodes.append("circle")
      .attr("r", nodeRadius)
      .attr("fill", d => nodeFillColor(d)) // Set fill color based on numProtons (grey if showing an as-yet-uncalculated shift)

    // Append atomNumber text
    // font-size comes from CSS (.node-text-atomNumber, .node-text-ppm in
    // molplot.css), driven by the --inverse-scale custom property — see
    // render() in view.js.
    allnodes.append("text")
        .attr("class", "node-text-atomNumber")
        .text(d => `${d.atomNumber}`) // Display the id
        .attr("dy", -textOffsetInt) // Center the text vertically
        .attr("text-anchor", "middle"); // Center the text horizontally

    // Append ppm text
    allnodes.append("text")
        .attr("class", "node-text-ppm")
        .text(d => nodePpmText(d)) // Display the experimental or calculated shift, per showCalculatedPpm
        .attr("dy", +textOffsetInt) // Center the text vertically
        .attr("text-anchor", "middle"); // Center the text horizontally

    /*
     * Placeholder nodes for carbons with no HSQC-assigned node
     * ─────────────────────────────────────────────────────────────────────────
     * A separate group, kept empty until showCalculatedPpm is toggled on, so
     * it never interferes with the index-based .data(nodes) rebinds that
     * updateMovedAtoms() does against nodesGroup/allnodes elsewhere.
     * renderCalcOnlyNodes() keys the join on atom id, so repeated toggling
     * just adds/removes the same elements rather than rebuilding them.
     */
    const calcOnlyNodesGroup = d3.select('svg').select('.calc-only-nodes-group').empty()
        ? d3.select('svg').append('g').attr('class', 'calc-only-nodes-group')
        : d3.select('svg').select('.calc-only-nodes-group');

    function renderCalcOnlyNodes() {
        // Recomputed from catoms (not cached) so a reassignment via the atom-move
        // feature — which flips a catom's visible flag — is picked up immediately.
        const calcOnlyNodes = calcOnlyNodesGroup.selectAll(".calc-only-node")
            .data(showCalculatedPpm ? catoms.filter(node => !node.visible) : [], d => d.id);

        calcOnlyNodes.exit().remove();

        const calcOnlyNodesEnter = calcOnlyNodes.enter().append("g")
            .attr("class", "node calc-only-node")
            .attr("transform", d => `translate(${d.x},${d.y})`)
            .attr("atomNumber", d => d.atomNumber)
            .attr("id", d => `calc-${d.id}`);

        calcOnlyNodesEnter.append("circle")
            .attr("r", nodeRadius / viewState.scaleLevel)
            .attr("fill", grey);

        // font-size comes from CSS here too — same rule, inherited from
        // #svg-image, so a newly-created placeholder node picks up the
        // current zoom level automatically with no JS calculation needed.
        calcOnlyNodesEnter.append("text")
            .attr("class", "node-text-atomNumber")
            .text(d => `${d.atomNumber}`)
            .attr("dy", `${-textOffsetInt / viewState.scaleLevel}px`)
            .attr("text-anchor", "middle");

        calcOnlyNodesEnter.append("text")
            .attr("class", "node-text-ppm")
            .text(d => d.ppm_calculated.toFixed(2))
            .attr("dy", `${textOffsetInt / viewState.scaleLevel}px`)
            .attr("text-anchor", "middle");
    }

