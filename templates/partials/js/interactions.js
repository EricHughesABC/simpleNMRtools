/*
 * interactions.js
 * ─────────────────────────────────────────────────────────────────────────
 * Hover/tooltip behaviour and the three toggle buttons that share state
 * with it (COSY, HMBC, INFO) — grouped together because handleMouseOut
 * forcibly resets laeVisible/hmbcVisible when the pointer leaves a node
 * while INFO/HMBC mode is active, so those flags and their toggle
 * functions need to live in one file rather than being split apart.
 *
 * Must load after nodes.js — reorderElements() (called below, defined in
 * links.js) requires nodesGroup to already exist in the DOM — and after
 * links.js itself, for that same reorderElements() call plus the HMBC
 * toggle's references to hmbclinks_svg/link/colorHmbcLinksByDistance.
 */

    // Tooltip for node hover
    const tooltip = d3.select("#tooltip");

    tooltip.style("opacity", 0);

    let hmbcLinks = []; // Track hmbc links separately

    reorderElements();

    // Function to find nearest neighbors and non-nearest neighbours of a node with a specific attribute
    function findNearestNeighbors(nodeId, attribute) {
        const neighbors = [];
        const nonNeighbors = [];

        // Iterate through the links
        links.forEach(link => {
            if (link[attribute] && (link.source === nodeId || link.target === nodeId)) {
                neighbors.push(link.source === nodeId ? link.target : link.source);
            }
        });

        // Iterate through all nodes to find non-neighbors
        nodes.forEach(node => {
            if (node.id !== nodeId && !neighbors.includes(node.id)) {
                nonNeighbors.push(node.id);
            }
        });

        return { neighbors, nonNeighbors };
    }

    // Function to display tooltip
    function displayTooltip(event, d){
        // generate the html sub string for d.H1 ppm which is an array

        // When numProtons is 0
        if (d.numProtons == 0) {
            if (d.ppm_calculated.length == 0) {
                var htmlstr =  `Node ${d.atomNumber} ${d.iupacLabel} <b>${generateHTMLnumProtonsString(d.numProtons)}</b>: <sup>13</sup>C&nbsp;=&nbsp;${d.ppm.toFixed(2)}&nbsp;ppm`;
            }
            else {
                var htmlstr =  `Node ${d.atomNumber} ${d.iupacLabel} <b>${generateHTMLnumProtonsString(d.numProtons)}</b>: <sup>13</sup>C&nbsp;=&nbsp;${d.ppm.toFixed(2)}&nbsp;(${d.ppm_calculated.toFixed(2)})&nbsp;ppm`;
            }
        } 
        else if (d.H1_ppm.length == 1) {
            // Expect d.H1 to be an array of length 1
            if (d.ppm_calculated.length == 0) {
                var htmlstr =  `Node ${d.atomNumber} ${d.iupacLabel} <b>${generateHTMLnumProtonsString(d.numProtons)}</b>: <sup>13</sup>C&nbsp;=&nbsp;${d.ppm.toFixed(2)}&nbsp;ppm, <sup>1</sup>H&nbsp;=&nbsp;${d.H1_ppm[0].toFixed(2)}&nbsp;ppm`;
            }
            else {
                var htmlstr =  `Node ${d.atomNumber} ${d.iupacLabel} <b>${generateHTMLnumProtonsString(d.numProtons)}</b>: <sup>13</sup>C&nbsp;=&nbsp;${d.ppm.toFixed(2)}&nbsp;(${d.ppm_calculated.toFixed(2)})&nbsp;ppm, <sup>1</sup>H&nbsp;=&nbsp;${d.H1_ppm[0].toFixed(2)}&nbsp;ppm`;
            }
                // var htmlstr =  `C${d.atomNumber} ${d.iupacLabel} <b>${generateHTMLnumProtonsString(d.numProtons)}</b>: <sup>13</sup>C = ${d.ppm.toFixed(2)} ppm, <sup>1</sup>H = ${d.H1_ppm[0].toFixed(2)} ppm`;
            // if jCouplingVals is not empty, add the J coupling values
            if (d.jCouplingClass.length > 0) {
                if (d.jCouplingClass === "m" || d.jCouplingClass === "s") {
                    htmlstr = htmlstr + `<br>J's : ${d.jCouplingClass}`;
                } 
                else {
                // htmlstr = htmlstr + `<br>J's : ${d.jCouplingClass[0]} : ${d.jCouplingVals[0]} Hz`;
                    htmlstr = htmlstr + `<br>J's : ${d.jCouplingClass} : ${d.jCouplingVals} Hz`;
                }
            }
        } 
        else {
            // Expect d.H1 to be an array of length 2
            if (d.ppm_calculated.length == 0) {
                var htmlstr =  `Node ${d.atomNumber} ${d.iupacLabel} <b>${generateHTMLnumProtonsString(d.numProtons)}</b>: <sup>13</sup>C&nbsp;=&nbsp;${d.ppm.toFixed(2)}&nbsp;ppm, <sup>1</sup>H&nbsp;=&nbsp;${d.H1_ppm[0].toFixed(2)} & ${d.H1_ppm[1].toFixed(2)}&nbsp;ppm`;
            }
            else {
                var htmlstr =  `Node ${d.atomNumber} ${d.iupacLabel} <b>${generateHTMLnumProtonsString(d.numProtons)}</b>: <sup>13</sup>C&nbsp;=&nbsp;${d.ppm.toFixed(2)}&nbsp;(${d.ppm_calculated.toFixed(2)})&nbsp;ppm, <sup>1</sup>H&nbsp;=&nbsp;${d.H1_ppm[0].toFixed(2)} & ${d.H1_ppm[1].toFixed(2)}&nbsp;ppm`;
            }
            // var htmlstr =  `C${d.atomNumber} ${d.iupacLabel} <b>${generateHTMLnumProtonsString(d.numProtons)}</b>: <sup>13</sup>C&nbsp;=&nbsp;${d.ppm.toFixed(2)}&nbsp;ppm, <sup>1</sup>H&nbsp;=&nbsp;${d.H1_ppm[0].toFixed(2)} & ${d.H1_ppm[1].toFixed(2)}&nbsp;ppm`;
            // add the J coupling values if not empty
            // print("d.jCouplingClass.length, d.jCouplingVals.length", d.jCouplingClass.length, d.jCouplingVals.length);
            if(d.jCouplingClass.length > 0) {
                // htmlstr = htmlstr + `<br>J's : ${d.jCouplingClass[0]} : ${d.jCouplingVals[0]} Hz, ${d.jCouplingClass[1]} : ${d.jCouplingVals[1]} Hz`;
                if (d.jCouplingClass === "m" || d.jCouplingClass === "s") {
                    htmlstr = htmlstr + `<br>J's : ${d.jCouplingClass}`;
                }
                else {
                    htmlstr = htmlstr + `<br>J's : ${d.jCouplingClass} : ${d.jCouplingVals} Hz`;
                }
            }
        }

        // Render off-screen first to measure actual dimensions, then reposition.
        tooltip.html(htmlstr)
            .style("left", "-9999px")
            .style("top",  "-9999px")
            .style("display", "inline-block")
            .style("opacity", 9.75);

        const tooltipNode = tooltip.node();
        const pageWidth   = document.documentElement.clientWidth;
        const tooltipX    = (pageWidth - tooltipNode.clientWidth) / 5;

        tooltip
            .style("left", `${tooltipX}px`)
            .style("top",  "25px");

        // Increase the radius of the hovered over node using d.id
        allnodes.filter(nodeDatum => nodeDatum.id === d.id)
            .select("circle")
            .attr("r", nodeRadiusLarge / viewState.scaleLevel);
    }

    // Function to highlight HMBC neighbors (updated)
    function highlightHMBCneighbors(event, d) {

        // check if the HMBC button is sunken and if it is release it and  reset the HMBC links to hidden and default colour
        if (d3.select("#hmbcButton").classed("sunken")){
            hmbcVisible = false;
            d3.select("#hmbcButton").classed("sunken", false);
            hmbclinks_svg.attr("opacity", 0);
            // hmbclinks_svg.attr("stroke", hmbcEdgeColor);
            colorHmbcLinksByDistance(mapping, 0);

            // check if the COSY button is not sunken, if it is not, set the opacity of the cosy links to 1
            if (!d3.select("#sunkenButton").classed("sunken")){
                // change the opacity of the cosy links to 1
                link.style("opacity", 1);
            }
        }
        
        // Define the attribute to find nearest neighbors
        const attribute = "hmbc";
        const attribute2 = "noesy";
        
        // Call the function to find nearest neighbors
        const hmbcNearestNodes = findNearestNeighbors(d.id, attribute);
        const noesyNearestNodes = findNearestNeighbors(d.id, attribute2);

        const nearestNeighborIDs = hmbcNearestNodes.neighbors;
        const nearestNeighborIDs2 = noesyNearestNodes.neighbors;

        // change the opacity of the nodes to 0.5
        const nonNeighborIds = new Set(hmbcNearestNodes.nonNeighbors);
        const dimmedNodes = allnodes.filter(nodeDatum => nonNeighborIds.has(nodeDatum.id));
        dimmedNodes.select("circle")
            .style("opacity", nodeHoverOpacity)
            .attr("r", nodeRadiusSmall / viewState.scaleLevel);

        // set the font colour to grey for all the text nodes
        dimmedNodes.selectAll(".node-text-atomNumber, .node-text-ppm")
            .style("fill", nodeHoverColor)
            .style("opacity", textHoverOpacity);

        // increase the size of the circles for the nearest neighbors
        const neighborIds = new Set(nearestNeighborIDs);
        allnodes.filter(nodeDatum => neighborIds.has(nodeDatum.id))
            .select("circle")
            .attr("r", nodeRadiusLarge / viewState.scaleLevel);

        // change the opacity of the cosy links to 0
        link.style("opacity", 0);

        // First reset all HMBC links to invisible if the global toggle is off
        if (!hmbcVisible) {
            hmbclinks_svg.attr("opacity", 0);
        }

        // Now show only connections to the current node, regardless of global toggle
        hmbclinks_svg.filter(l => l.source === d.id).attr("opacity", 1);
        hmbclinks_svg.filter(l => l.target === d.id).attr("opacity", 1);
    }

    // Function to handle mouseover event
    function handleMouseOver(event, d) {

        if (viewState.dragStarted) {
            return;
        }   
        // Display tooltip
        if (event) {
            displayTooltip(event, d);
            highlightHMBCneighbors(event,d);
        }
    }

    // Function to handle mouseout event
    function handleMouseOut(event, d) {

        if (viewState.dragStarted) {
            return;
        }   
        // change the opacity of the text nodes to 1
        // reset the font colour of the text in the nodes to black
        d3.selectAll(".node text")
            .style("opacity", 1)
            .style("fill", "black");

        // set the opacity of all the nodes to 1
        d3.selectAll(".node circle")
            .style("opacity", 1);

        // reset the fill color and radius of every real node's circle
        allnodes.each(function(nodeDatum) {
            d3.select(this).select("circle")
                .style("fill", nodeFillColor(nodeDatum))
                .attr("r", nodeRadius / viewState.scaleLevel);
        });

        // Hide the tooltip
        tooltip.transition()
                .duration(100)
                .style("opacity", 0);
        // Revert the radius of this node's circle
        // (event.currentTarget is the <g> node group; event.target would be whichever
        //  child the pointer was over, which is unreliable with mouseleave)
        d3.select(event.currentTarget).select("circle")
            .transition()
            .duration(100)
            .attr("r", nodeRadius/viewState.scaleLevel);

        // Remove hmbc links
        hmbcLinks.forEach(link => link.remove());
        hmbcLinks = []; // Clear the hmbc links array

        // change the opacity of the cosy links to 1
        link.style("opacity", 1);

        // hmbclinks_svg.style("opacity", 0);
        hmbclinks_svg.filter(l => l.source === d.id).attr("opacity", 0);
        hmbclinks_svg.filter(l => l.target === d.id).attr("opacity", 0);

        if (d3.select("#infoButton").classed("sunken")){
            // Forcibly deactivate INFO mode on mouseout — keep laeVisible in sync
            laeVisible = false;
            d3.select("#infoButton").classed("sunken", false);
            // Reset node fill colours and isolated-node strokes
            allnodes.each(function(d) {
                d3.select(this).select("circle").style("fill", nodeFillColor(d));
            });
            applyIsolatedNodeStrokes(false);
        }
    }

    function generateHTMLnumProtonsString(index) {
        switch (index) {
            case 0:
            return "-C-";
            case 1:
            return "-CH";
            case 2:
            return "-CH<sub>2</sub>";
            case 3:
            return "-CH<sub>3</sub>";
            default:
            return "<p>Invalid index</p>";
        }
    }

    const cosybutton = document.getElementById('sunkenButton');

    // Toggle between 'sunken' and 'normal' styles when clicked
    let linksVisible = true;
    cosybutton.addEventListener('click', function() {
        cosybutton.classList.toggle('sunken');

        if (linksVisible) {
            d3.selectAll(".cosy-link").style("display", "none");
        } else {
            d3.selectAll(".cosy-link").style("display", "inline");
        }
        linksVisible = !linksVisible;
    });


    // Toggle code
    const hmbcbutton = document.getElementById('hmbcButton');
    let hmbcVisible = false;

    function toggleHmbcVisibility() {
        hmbcbutton.classList.toggle('sunken');
        hmbcVisible = !hmbcVisible;
        
        if (hmbcVisible) {
            // Show all HMBC edges with opacity
            d3.selectAll(".hmbc-link").attr("opacity", 1);
            // If the cosy button is not sunken, set the opacity of the cosy links to 0
            if (!d3.select("#sunkenButton").classed("sunken")) {
                // change the opacity of the cosy links to 0
                link.style("opacity", 0);
            }
        } else {
            // Hide all HMBC edges with opacity
            d3.selectAll(".hmbc-link").attr("opacity", 0);
            // If the cosy button is not sunken, set the opacity of the cosy links to 1
            if (!d3.select("#sunkenButton").classed("sunken")) {
                // change the opacity of the cosy links to 1
                link.style("opacity", 1);       
            }
        }
    }


    hmbcbutton.addEventListener('click', toggleHmbcVisibility);

    const infobutton = document.getElementById('infoButton');

    // Toggle between 'sunken' and 'normal' styles when clicked
    // Starts false to match the button's initial unsunken visual state.
    let laeVisible = false;
    infobutton.addEventListener('click', function() {
        infobutton.classList.toggle('sunken');
        laeVisible = !laeVisible;

        if (laeVisible) {
            // Single pass: colour each node by its assignment quality.
            // Priority: worst-LAE node → red, large ppm error (>10) → pink, else unchanged.
            allnodes.each(function(d) {
                const circle = d3.select(this).select("circle");
                if (d.atomNumber == best_results.best_lae_atomNumber) {
                    circle.style("fill", "red");
                } else if (Math.abs(d.ppm - d.ppm_calculated) > 10) {
                    circle.style("fill", "pink");
                }
            });

            // Highlight isolated nodes (no links) with a black circumference
            applyIsolatedNodeStrokes(true);

        } else {

            // Assuming allnodes is already defined as in your code
            allnodes.each(function(d) {
                // Select the circle within the current node group and reset its fill color
                d3.select(this).select("circle").style("fill", nodeFillColor(d));
            });

            // Remove isolated node stroke highlighting
            applyIsolatedNodeStrokes(false);

        }
    });


    // Build a Set of isolated node IDs for O(1) lookup
    const isolatedNodeIds = new Set(isolatedCarbons.map(n => n.id));

    /*
     * applyIsolatedNodeStrokes — add or remove the black circumference ring
     * on carbon nodes that have no COSY or HMBC links.
     *
     * show = true  → thick black stroke, scaled so it stays visually
     *                constant regardless of the current zoom level.
     * show = false → stroke removed (null resets the inline style so the
     *                CSS default of "none" takes over again).
     */
    function applyIsolatedNodeStrokes(show) {
        allnodes.each(function(d) {
            if (isolatedNodeIds.has(d.id)) {
                const circle = d3.select(this).select("circle");
                if (show) {
                    circle
                        .style("stroke", "black")
                        .style("stroke-width", `${3 / viewState.scaleLevel}px`);
                } else {
                    circle
                        .style("stroke", null)
                        .style("stroke-width", null);
                }
            }
        });
    }

    // Hover preview: show isolated-node rings while hovering the INFO button
    // (only when the button is not already in its sunken/active state)
    infobutton.addEventListener('mouseenter', function() {
        if (!d3.select("#infoButton").classed("sunken")) {
            applyIsolatedNodeStrokes(true);
        }
    });

    infobutton.addEventListener('mouseleave', function() {
        if (!d3.select("#infoButton").classed("sunken")) {
            applyIsolatedNodeStrokes(false);
        }
    });

