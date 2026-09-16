/*
 * links.js
 * ─────────────────────────────────────────────────────────────────────────
 * COSY/HMBC edge data and rendering: the <line> elements between nodes,
 * HMBC distance-based colouring, and SVG paint-order (links must sit
 * behind nodes). Depends on constants.js (edge colours/widths) and
 * state.js (nodes, mapping, links, shortest_paths) — both load first.
 */

    // Set initial positions of nodes and links

    const linksGroup = d3.select('svg').select('.links-group').empty() 
        ? d3.select('svg').insert('g', ':first-child').attr('class', 'links-group')
        : d3.select('svg').select('.links-group');

    // Filter links that have the attribute cosy as true
    const cosyLinks = links.filter(eachLink => eachLink.cosy);
    const hmbcLinksData = links.filter(eachLink => eachLink.hmbc);


    function colorHmbcLinksByDistance(mapping, offset) {
        d3.selectAll(".hmbc-link").attr("stroke", d => {
            const sourceId = mapping[d.source] - offset;
            const targetId = mapping[d.target] - offset;
            
            // Check if the path length is greater than 2 (long range)
            if (shortest_paths[sourceId] && 
                shortest_paths[sourceId][targetId] && 
                shortest_paths[sourceId][targetId] > 2) {
                return "#FF5733"; // Red color for long range HMBC links
            } else {
                return "#808080"; // Grey color for short range HMBC links
            }
        });
    }

function reorderElements() {
    const svg = d3.select('svg');
    
    // Get both groups
    const linksGroupNode = svg.select('.links-group').node();
    const nodesGroupNode = svg.select('.nodes-group').node();
    
    if (!linksGroupNode || !nodesGroupNode) {
        console.error("Could not find links or nodes group");
        return;
    }
    
    // CRITICAL: Ensure linksGroup comes before nodesGroup in SVG
    // Remove and re-insert linksGroup as the first child
    const svgNode = svg.node();
    svgNode.removeChild(linksGroupNode);
    svgNode.insertBefore(linksGroupNode, nodesGroupNode);
    
    // Now reorder links within the linksGroup (COSY → grey HMBC → red HMBC)
    d3.selectAll(".cosy-link").each(function() {
        linksGroupNode.appendChild(this);
    });
    
    d3.selectAll(".hmbc-link").each(function() {
        const color = d3.select(this).attr("stroke");
        if (color === "#808080" || color === "rgb(128, 128, 128)") {
            linksGroupNode.appendChild(this);
        }
    });
    
    d3.selectAll(".hmbc-link").each(function() {
        const color = d3.select(this).attr("stroke");
        if (color === "#FF5733" || color === "rgb(255, 87, 51)") {
            linksGroupNode.appendChild(this);
        }
    });
}

    // Create the COSY links to start
    // Create a group for links
    const link = linksGroup.selectAll(".links-group")
        .data(cosyLinks)
        .enter().append("line")
        .attr("class", "link cosy-link")
        .attr("stroke", cosyEdgeColor)
        .attr("source", d => d.source)
        .attr("target", d => d.target)
        .attr("x1", d => nodes.find(node => node.id === d.source).x)
        .attr("y1", d => nodes.find(node => node.id === d.source).y)
        .attr("x2", d => nodes.find(node => node.id === d.target).x)
        .attr("y2", d => nodes.find(node => node.id === d.target).y);
        // stroke-width comes from CSS (.cosy-link, .hmbc-link in molplot.css)


    // Edge initialization
    const hmbclinks_svg = linksGroup.selectAll(".links-group")
        .data(hmbcLinksData)
        .enter().append("line")
        .attr("class", "link hmbc-link")
        .attr("stroke", hmbcEdgeColor)
        .attr("x1", d => nodes.find(node => node.id === d.source).x)
        .attr("y1", d => nodes.find(node => node.id === d.source).y)
        .attr("x2", d => nodes.find(node => node.id === d.target).x)
        .attr("y2", d => nodes.find(node => node.id === d.target).y)
        .attr("source", d => d.source)
        .attr("target", d => d.target)
        .attr("opacity", 0); // Initially hidden with opacity

    // Color HMBC links by distance
    colorHmbcLinksByDistance(mapping, 0);

