/*
 * moves.js
 * ─────────────────────────────────────────────────────────────────────────
 * Reassigning a peak by dragging its node onto another: the drag handlers,
 * the pure findValidMoves/commitMoves pair (detection separated from
 * mutation — see their own comments), the updateMovedAtoms3 orchestrator,
 * updateMovedAtoms (redraws + recomputes MAE/LAE/correlation penalty after
 * a committed move), and the 2D correlation-penalty weight function.
 *
 * Must load before nodes.js: node creation wires dragstarted/dragged/
 * dragended directly (`.call(d3.drag().on("start", dragstarted)...)`),
 * and that reference is resolved immediately, not deferred to a later
 * event — so these functions need to already be defined as globals by
 * the time nodes.js runs.
 */

// Drag handlers
function dragstarted(event, d) {
    // Only hide COSY links if they are currently visible (COSY button not toggled off)
    if (!d3.select("#sunkenButton").classed("sunken")) {
        link.style("opacity", 0);
    }
    d3.select(this).raise().classed("active", true);
    viewState.dragStarted = true;
}

function dragged(event, d) {
    d.x = event.x;
    d.y = event.y;
    link.filter(l => l.source === d.id).attr("x1", d.x).attr("y1", d.y);
    link.filter(l => l.target === d.id).attr("x2", d.x).attr("y2", d.y);
    hmbclinks_svg.filter(l => l.source === d.id).attr("x1", d.x).attr("y1", d.y);
    hmbclinks_svg.filter(l => l.target === d.id).attr("x2", d.x).attr("y2", d.y);

    d3.select(this).attr("transform", `translate(${d.x},${d.y})`);
}

function dragended(event, d) {
    d3.select(this).classed("active", false);
    viewState.dragStarted = false;
    handleMouseOver(event, d);
}

    function updateMovedAtoms() {
        // Update the moved atoms

        if( !updateMovedAtoms3()){
            return
        };

        nodes = JSON.parse(JSON.stringify(catoms.filter(node => node.visible)));

        // Keep the "not yet added" grey placeholder layer in sync — a move can
        // flip which catoms are visible, so refresh it if it's currently shown.
        if (showCalculatedPpm) {
            renderCalcOnlyNodes();
        }

        // Select all circle elements with the class 'node' in the nodesGroup
        var circlenodes = nodesGroup.selectAll('.node circle');

        // Add circles to represent nodes
        // svg.selectAll(".node circle")

        // update node colour
        nodesGroup.selectAll(".node circle")
            .data(nodes)
            .style("fill", d => nodeFillColor(d))
            .raise(); // Display the id

        // update atomNumber
        nodesGroup.selectAll(".node text.node-text-atomNumber")
            .data(nodes)
            .text(d => `${d.atomNumber}`)
            .raise(); // Display the id

        // update the ppm text
        nodesGroup.selectAll(".node text.node-text-ppm")
            .data(nodes)
            .text(d => nodePpmText(d))
            .raise(); // Display the ppm

        // move the nodes
        nodesGroup.selectAll(".node")
            .data(nodes)
            .attr("class", "node")
            .transition()
            .duration(1000)
            .attr("transform", d => `translate(${d.x},${d.y})`);

        // Update COSY link endpoints using the typed named selection
        link
            .attr("x1", d => nodes.find(node => node.id === d.source).x)
            .attr("y1", d => nodes.find(node => node.id === d.source).y)
            .attr("x2", d => nodes.find(node => node.id === d.target).x)
            .attr("y2", d => nodes.find(node => node.id === d.target).y);

        // Update HMBC link endpoints using the typed named selection
        hmbclinks_svg
            .attr("x1", d => nodes.find(node => node.id === d.source).x)
            .attr("y1", d => nodes.find(node => node.id === d.source).y)
            .attr("x2", d => nodes.find(node => node.id === d.target).x)
            .attr("y2", d => nodes.find(node => node.id === d.target).y);

        // recalculate the MAE and LAE

        var mae = 0.0;
        var lae_best = 0.0;
        var lae_best_atomNumber = -1;
        
        for(var i=0; i<nodes.length; i++) {
            var mae_latest = Math.abs(nodes[i].ppm - nodes[i].ppm_calculated);

            if (mae_latest > lae_best) {
                lae_best = mae_latest;
                lae_best_atomNumber = nodes[i].atomNumber;
            }
            mae += mae_latest           
        }
        mae = mae/nodes.length;


        // Rebuild mapping: NMR-peak id → the original catom id at the same structural position.
        // Used by computeTotalWeight2D for HMBC correlation scoring.
        // Step 1: reset all entries to identity (invisible atoms map to themselves)
        catoms_orig.forEach(c => { mapping[c.id] = c.id; });
        // Step 2: build atomNumber → original-id lookup once
        const atomNumberToOrigId = Object.fromEntries(
            catoms_orig.map(c => [c.atomNumber, c.id])
        );
        // Step 3: update mapping for each currently displayed node
        nodes.forEach(node => {
            mapping[node.id] = atomNumberToOrigId[node.atomNumber];
        });


        // recaclulate correlation penalty
        var best_weight = computeTotalWeight2D(links, mapping, shortest_paths, 0);

        // update the best results
        best_results.best_weight = best_weight;
        best_results.best_mae = mae;
        best_results.best_lae = lae_best;
        best_results.best_lae_atomNumber = lae_best_atomNumber;

        // update the INFO tooltip
        // d3.select("#infoButton")
        // .attr("title", `Correlation Penalty: ${best_results.best_weight}, MAE: ${best_results.best_mae.toFixed(2)}, LAE: ${best_results.best_lae.toFixed(2)} (NODE ${best_results.best_lae_atomNumber})`);
        var info_tooltip_str = create_tooltip_string(best_results.best_weight, best_results.best_mae, best_results.best_lae, best_results.best_lae_atomNumber, isolatedCarbonsCount, nodes.length, links.length, number_of_hmbc_cosy_subgraphs);
        // d3.select("#info_tooltip").html(`Correlation Penalty:&nbsp;${best_results.best_weight}<br>Uncorrelated&nbsp;Nodes:&nbsp;${isolatedCarbonsCount}&nbsp;out&nbsp;of&nbsp;${nodes.length}<br>MAE:&nbsp;${best_results.best_mae.toFixed(2)},&nbsp;LAE:&nbsp;${best_results.best_lae.toFixed(2)}&nbsp;(NODE&nbsp;${best_results.best_lae_atomNumber})`);
        d3.select("#info_tooltip").html(info_tooltip_str);
        colorHmbcLinksByDistance(mapping, 0);
        reorderElements();
    }


    /*
     * findValidMoves — pure detection + validation, no side effects
     * ─────────────────────────────────────────────────────────────────────────
     * Determines which dragged nodes can legitimately be committed to a new
     * structural position.  Returns an array of move descriptors; an empty
     * array means nothing should be committed.
     *
     * Each returned object has the shape:
     *   { movedAtomNumber, targetAtomNumber, snapshot }
     *
     * where snapshot is a deep copy of the node taken before any mutations,
     * used later by commitMoves to read original NMR values regardless of order.
     *
     * Validation rules (a move is rejected if ANY of these are true):
     *   1. Source and destination have different numProtons (type mismatch)
     *   2. Source and destination are the same atom (no-op)
     *   3. Destination is NOT also moving AND is already visible (occupied)
     *   4. Destination IS also moving but that move is itself invalid (broken chain)
     */
    function findValidMoves(nodes, catoms, catoms_orig) {

        // --- Step 1: snapshot nodes that have moved from their original positions ---
        // Use atomNumber (structural identity) not array index for the catoms_orig lookup
        const movedSnapshots = [];
        nodes.forEach(node => {
            const orig = catoms_orig.find(c => c.atomNumber === node.atomNumber);
            if (orig && (node.x !== orig.x || node.y !== orig.y)) {
                movedSnapshots.push(JSON.parse(JSON.stringify(node)));
            }
        });

        if (movedSnapshots.length === 0) return [];

        // --- Step 2: find the closest catoms_orig position for each moved node ---
        const candidates = movedSnapshots.map(snapshot => {
            const distancesSquared = catoms_orig.map(c =>
                (snapshot.x - c.x) ** 2 + (snapshot.y - c.y) ** 2
            );
            const minDist  = Math.min(...distancesSquared);
            const target   = catoms_orig[distancesSquared.indexOf(minDist)];
            return {
                movedAtomNumber:  snapshot.atomNumber,
                targetAtomNumber: target.atomNumber,
                distanceSquared:  minDist,
                snapshot,
                invalid: false
            };
        });

        // --- Step 3: validation pass 1 — type mismatch and no-op ---
        candidates.forEach(c => {
            const src = catoms.find(a => a.atomNumber === c.movedAtomNumber);
            const dst = catoms.find(a => a.atomNumber === c.targetAtomNumber);
            if (src.numProtons !== dst.numProtons) c.invalid = true;
            if (src.atomNumber === dst.atomNumber)  c.invalid = true;
        });

        // --- Step 4: validation pass 2 — destination occupancy / chain moves ---
        // A destination that is itself being moved is valid if that move is valid.
        // A destination that is not moving must be unoccupied (visible === false).
        candidates.forEach(c => {
            if (c.invalid) return;
            const dst          = catoms.find(a => a.atomNumber === c.targetAtomNumber);
            const chainIdx     = candidates.findIndex(other => other.movedAtomNumber === c.targetAtomNumber);
            const destIsMoving = chainIdx !== -1;

            if (!destIsMoving) {
                c.invalid = dst.visible;           // occupied position → reject
            } else {
                c.invalid = candidates[chainIdx].invalid;  // inherit chain validity
            }
        });

        return candidates.filter(c => !c.invalid);
    }

    /*
     * commitMoves — apply validated moves to catoms (side effects only)
     * ─────────────────────────────────────────────────────────────────────────
     * Reads NMR data exclusively from each move's snapshot (the pre-commit state)
     * so swap and chain moves are order-independent in the data-copy phase.
     *
     * Visibility updates are split into two sub-passes:
     *   pass A: mark all source positions invisible  (clears old occupancy)
     *   pass B: mark all destination positions visible (sets new occupancy)
     * This order is required for swaps: without it, a source cleared in pass A
     * would be immediately re-set visible in pass B before its pair is cleared.
     */
    function commitMoves(validMoves) {

        // Pass 1 — copy NMR data to destination (reads from snapshot, no order issues)
        validMoves.forEach(move => {
            const dst = catoms.find(c => c.atomNumber === move.targetAtomNumber);
            dst.ppm             = move.snapshot.ppm;
            dst.H1_ppm          = move.snapshot.H1_ppm;
            dst.jCouplingVals   = move.snapshot.jCouplingVals;
            dst.jCouplingClass  = move.snapshot.jCouplingClass;
            dst.id              = move.snapshot.id;
        });

        // Pass 2a — clear all source positions
        validMoves.forEach(move => {
            catoms.find(c => c.atomNumber === move.movedAtomNumber).visible = false;
        });

        // Pass 2b — activate all destination positions
        validMoves.forEach(move => {
            catoms.find(c => c.atomNumber === move.targetAtomNumber).visible = true;
        });
    }

    /*
     * updateMovedAtoms3 — thin orchestrator
     * Calls findValidMoves then commitMoves; returns false when there is
     * nothing to commit so that updateMovedAtoms() can exit early.
     */
    function updateMovedAtoms3() {
        const validMoves = findValidMoves(nodes, catoms, catoms_orig);
        if (validMoves.length === 0) return false;
        commitMoves(validMoves);
        return true;
    }

    function computeTotalWeight2D(links, mapping, shortestPaths, node_offset) {
        let totalWeight = 0;

        links.forEach(link => {
            const u = link.source;
            const v = link.target;
            const d = link;

            if (d.cosy || d.hmbc) {
                const uu = mapping[u] - node_offset;
                const vv = mapping[v] - node_offset;


                if(uu == vv) {
                    return;
                }

                if (d.cosy) {
                    let weight = shortestPaths[uu][vv];
                    totalWeight += Math.pow(weight - 1, 3);
                }
                if (d.hmbc) {
                    let weight = shortestPaths[uu][vv];
                    if (weight < 3) {
                        weight = 2;
                    }
                    totalWeight += Math.pow(weight - 2, 3);
                }
            }
        });
        return totalWeight;
    }

