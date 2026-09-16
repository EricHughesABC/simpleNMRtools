/*
 * view.js
 * ─────────────────────────────────────────────────────────────────────────
 * The SVG viewport: element setup (id/viewBox/preserveAspectRatio), the
 * zoom/pan/rotate transform pipeline (viewState, render, handleZoom,
 * rotateSVG), resize handling, the double-click view reset, and trackpad
 * rotation gestures.
 *
 * Depends on constants.js (svg_width_eeh/svg_height_eeh, nodeRadius,
 * textOffsetInt) and, in render()'s transform block only, nodesGroup/
 * calcOnlyNodesGroup (nodes.js) and linksGroup (links.js). render()'s
 * radius/dy updates reach node elements from both nodes.js groups via a
 * plain `.node ...` selector rather than referencing either group's
 * variable directly — see the comment at that code for why. None of this
 * runs until the first user interaction (zoom/pan/rotate/resize), well
 * after every script on the page has loaded, so there's no load-order
 * requirement between view.js and nodes.js/links.js.
 */

    // Select the SVG element and add an id to it
    d3.select("svg").attr("id", "svg-image");
        
    // Select the first <g> element and add an id to it
    d3.select("#svg-image").select("g").attr("id", "molimage")

    const svgElement = document.querySelector('svg.center');

    /*
     * SVG element responsive setup
     * ─────────────────────────────────────────────────────────────────────────
     * The server may have written fixed pixel values (e.g. width="1000"
     * height="600") directly into the <svg> element.  These are overridden here
     * so that CSS — not inline attributes — controls the rendered dimensions.
     *
     * width / height → "100%"
     *   Together with the position:absolute CSS rule, the SVG stretches to
     *   fill its .plot-container parent exactly, at any screen size.
     *
     * viewBox → "0 0 1000 600"
     *   Declares the internal coordinate system.  All atom positions, node
     *   x/y values, link endpoints, and transforms use coordinates in this
     *   0–1000 / 0–600 space.  The browser maps it to screen pixels
     *   automatically — no JavaScript recalculation is needed on resize.
     *
     * preserveAspectRatio → "xMidYMid meet"
     *   "meet": scale the viewBox uniformly (maintain aspect ratio) until it
     *   fits inside the element rectangle — equivalent to CSS object-fit:contain.
     *   "xMidYMid": centre the viewBox horizontally and vertically within the
     *   element when the element's aspect ratio differs from 5:3 (e.g. when
     *   the viewport height cap kicks in and the element is wider than it is
     *   5:3 tall).  This keeps the molecule centred rather than pinned to a
     *   corner.
     */
    svgElement.setAttribute('width',  '100%');
    svgElement.setAttribute('height', '100%');
    svgElement.setAttribute('viewBox', `0 0 ${svg_width_eeh} ${svg_height_eeh}`);
    svgElement.setAttribute('preserveAspectRatio', 'xMidYMid meet');

    // Initial D3 zoom/pan transform
    const initialTransform = d3.zoomIdentity.translate(0, 0);

    // Track the current transform
    let currentTransform = initialTransform;

    // Append a group element to hold the content
    const svg = d3.select("#svg-image").attr("transform", initialTransform);

    // Add zoom and pan functionality.
    // The filter excludes:
    //   dblclick   — replaced by the view-reset handler below
    //   shift+wheel — claimed by the Shift+scroll rotation handler below
    const myZoom = d3.zoom()
        .scaleExtent([0.5, 5])
        .filter(event => event.type !== 'dblclick' &&
                         !(event.type === 'wheel' && event.shiftKey))
        .on("zoom", handleZoom);

    svg.call(myZoom);

    const molimage = d3.select("#molimage");

    /*
     * Resize and orientation-change handling
     * ─────────────────────────────────────────────────────────────────────────
     * Because all D3 node coordinates are stored in SVG viewBox space
     * (0–1000 / 0–600) and the browser scales the entire SVG uniformly via
     * viewBox + preserveAspectRatio, there is *nothing to reposition* when the
     * window is resized or the device is rotated.  The background molecule
     * image and all D3 nodes/links scale together automatically.
     *
     * The only actions needed on resize are calling render(), which redraws
     * the current pan/zoom/rotation state and re-applies counter-scaling so
     * that node circles and text stay visually the same size at any zoom level.
     *
     * Debounce (150 ms):
     *   Both 'resize' and 'orientationchange' fire many times per second during
     *   a window drag or device rotation animation.  Debouncing ensures
     *   render() is called only once, 150 ms after the last event, avoiding
     *   unnecessary layout thrashing.
     *
     * Three event listeners for full browser coverage:
     *   'resize'            — window size change on desktop; also fires on
     *                         mobile when the soft keyboard appears/disappears.
     *   'orientationchange' — legacy event; fires on older iOS and Android
     *                         browsers that do not support screen.orientation.
     *   screen.orientation  — modern W3C API; used on current Chrome/Firefox/
     *                         Safari where 'orientationchange' may be deprecated.
     */
    let resizeTimer;
    function onResize() {
        clearTimeout(resizeTimer);
        resizeTimer = setTimeout(() => {
            render();
        }, 150);
    }

    window.addEventListener('resize', onResize);
    window.addEventListener('orientationchange', onResize);  // legacy Safari/Android
    if (screen.orientation) {
        screen.orientation.addEventListener('change', onResize); // modern API
    }

/*
 * viewState — single source of truth for all pan / zoom / rotation / drag state.
 * Grouping these prevents accidental omission when resetting or inspecting the
 * current view, and makes the relationship between the variables obvious.
 */
const viewState = {
    panX:          0,
    panY:          0,
    scaleLevel:    1,
    rotation:      0,
    lastTransform: d3.zoomIdentity,
    dragStarted:   false
};

/*
 * render — apply the current pan/zoom/rotation state and counter-scale elements
 * ───────────────────────────────────────────────────────────────────────────────
 * Called after every state change (zoom, pan, rotate, resize). Does two things
 * in a single pass so call sites can never forget one half:
 *
 * 1. Builds a transform string and applies it identically to all three SVG
 *    groups (molecule image, link lines, atom nodes) so they stay aligned.
 *    The transform sequence (applied right-to-left per SVG convention):
 *      translate(-cx,-cy)   move viewBox centre to origin
 *      rotate(viewState.rotation)     rotate around that origin
 *      scale(viewState.scaleLevel)    uniform scale around origin
 *      translate(cx,cy)     restore centre
 *      translate(viewState.panX,viewState.panY) apply accumulated pan
 *    Rotation centre is always the viewBox centre (500,300) — using viewBox
 *    coordinates here is critical so the centre doesn't drift on resize.
 *
 * 2. Counter-scales nodes, text, and edge strokes so they appear the same
 *    visual size regardless of the current zoom level.
 */
function render() {
    const cx = svg_width_eeh  / 2;
    const cy = svg_height_eeh / 2;
    const transformString =
        `translate(${viewState.panX}, ${viewState.panY}) ` +
        `translate(${cx}, ${cy}) ` +
        `scale(${viewState.scaleLevel}) ` +
        `rotate(${viewState.rotation}) ` +
        `translate(${-cx}, ${-cy})`;

    molimage.attr("transform", transformString);
    linksGroup.attr("transform", transformString);
    nodesGroup.attr("transform", transformString);
    calcOnlyNodesGroup.attr("transform", transformString);

    // font-size (node text) and stroke-width (cosy/hmbc links) are driven by
    // CSS via --inverse-scale (see molplot.css) — one property update here
    // replaces what used to be four separate selectAll().attr() passes.
    svg.style("--inverse-scale", 1 / viewState.scaleLevel);

    // Radius and text dy stay JS-driven — see molplot.css's comment on why
    // (radius has competing hover states a CSS rule would override; dy's
    // CSS support for SVG text is less consistent across browsers).
    //
    // One fresh selector reaches both the real (nodesGroup/allnodes) and
    // grey placeholder (calcOnlyNodesGroup) node elements at once — they
    // all share the .node class — instead of repeating each call for both
    // groups. This is a plain, unbound d3.selectAll() (no .data() join),
    // so it neither affects nor is affected by the index-keyed data rebind
    // updateMovedAtoms() does against nodesGroup specifically (moves.js).
    d3.selectAll(".node circle")
        .attr("r", nodeRadius / viewState.scaleLevel);

    d3.selectAll(".node-text-atomNumber")
        .attr("dy", `${-textOffsetInt / viewState.scaleLevel}px`);

    d3.selectAll(".node-text-ppm")
        .attr("dy", `${textOffsetInt / viewState.scaleLevel}px`);
}

function handleZoom(event) {
    const { x, y, k } = event.transform;

    // Only update pan position for drag events (not wheel/zoom)
    if (event.sourceEvent) {
        if (event.sourceEvent.type === 'mousemove' || event.sourceEvent.type === 'touchmove') {
            // This is a pan/drag - accumulate the delta
            const dx = x - viewState.lastTransform.x;
            const dy = y - viewState.lastTransform.y;
            viewState.panX += dx;
            viewState.panY += dy;
        }
    }

    // Update scale level
    viewState.scaleLevel = k;
    viewState.lastTransform = event.transform;

    render();
}

// Function to handle rotation
function rotateSVG(angle) {
    viewState.rotation = (viewState.rotation + angle) % 360;

    render();

    // Counter-rotate text labels so they stay horizontal after the group
    // rotates — one CSS custom property, inherited by every node's text
    // (both allnodes and calcOnlyNodesGroup), replaces four selectAll().attr()
    // passes.
    svg.style("--counter-rotate", `${-viewState.rotation}deg`);
}

/*
 * Double-click on the SVG — reset pan, zoom and rotation to the initial state.
 *
 * D3's built-in dblclick zoom is suppressed via the .filter() on myZoom above.
 * This handler replaces it with a full view reset:
 *   - viewState.panX/viewState.panY/viewState.scaleLevel/viewState.rotation are set back to their starting values
 *   - D3's internal zoom transform is synced via myZoom.transform so that the
 *     next wheel or pinch gesture starts from the correct baseline
 *   - Text counter-rotation is removed (viewState.rotation is now 0)
 *   - render() redraws all groups with the reset transform
 */
d3.select("#svg-image").on("dblclick", function(event) {
    event.preventDefault();

    viewState.panX      = 0;
    viewState.panY      = 0;
    viewState.scaleLevel = 1;
    viewState.rotation  = 0;
    viewState.lastTransform = d3.zoomIdentity;

    // Keep D3's internal zoom state in sync with the manual reset
    svg.call(myZoom.transform, d3.zoomIdentity);

    // Remove text counter-rotation (viewState.rotation is now 0)
    svg.style("--counter-rotate", "0deg");

    render();
});

/*
 * Trackpad rotation
 * ─────────────────────────────────────────────────────────────────────────────
 * Option A — Two-finger twist gesture (macOS, Chrome + Safari)
 *   Chrome on macOS fires WebKit GestureEvents when the user performs a
 *   two-finger twist on the trackpad.  event.rotation is the cumulative angle
 *   in degrees since gesturestart, so we track the previous value and pass
 *   the per-frame delta to rotateSVG().
 *
 * Option B — Shift + two-finger scroll (universal fallback)
 *   When Shift is held during a scroll, wheel deltaY drives rotation instead
 *   of zoom.  The myZoom filter above already excludes these events so D3
 *   never sees them.  passive:false is required so preventDefault() can stop
 *   the browser from scrolling the page while Shift is held.
 *   Multiplier 0.2 converts typical trackpad deltaY values (≈5–50 px per
 *   event) into a comfortable 1–10° of rotation per scroll tick.
 */
const svgDomNode = document.getElementById('svg-image');
let lastGestureRotation = 0;

// Option A — gesture events
svgDomNode.addEventListener('gesturestart', (event) => {
    event.preventDefault();
    lastGestureRotation = event.rotation;
});

svgDomNode.addEventListener('gesturechange', (event) => {
    event.preventDefault();
    const delta = event.rotation - lastGestureRotation;
    lastGestureRotation = event.rotation;
    rotateSVG(delta);
});

svgDomNode.addEventListener('gestureend', (event) => {
    event.preventDefault();
});

// Option B — Shift + scroll
svgDomNode.addEventListener('wheel', (event) => {
    if (event.shiftKey) {
        event.preventDefault();
        rotateSVG(event.deltaY * 0.2);
    }
}, { passive: false });

