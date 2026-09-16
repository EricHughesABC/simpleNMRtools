/*
 * constants.js
 * ─────────────────────────────────────────────────────────────────────────
 * Design-space dimensions, touch detection, colour palette, node/edge
 * sizing, and the numProtons -> colour scale. Pure constants only — nothing
 * here depends on server-rendered data or the DOM (beyond isTouch's
 * matchMedia check), so this loads first, before state.js and view.js.
 */

    /*
     * Design-space dimensions
     * ─────────────────────────────────────────────────────────────────────────
     * These constants define the SVG *viewBox* coordinate space — the fixed
     * 1000 × 600 grid in which the server generates all atom positions and in
     * which D3 places all nodes, links, and transforms.
     *
     * They are NOT the pixel dimensions of the rendered element on screen.
     * The browser maps this coordinate space to whatever pixel size the
     * .plot-container happens to be (via viewBox + preserveAspectRatio set
     * below), scaling everything uniformly without any JS involvement.
     *
     * Keeping these as named constants (rather than inline literals) means
     * the viewBox, xScale/yScale, and applyTransform() rotation centre all
     * stay in sync if the design dimensions ever change.
     */
    const svg_width_eeh  = 1000;
    const svg_height_eeh = 600;

    /*
     * Touch / coarse-pointer detection
     * ─────────────────────────────────────────────────────────────────────────
     * Evaluated once at page load.  matchMedia('(pointer: coarse)') returns
     * true for touchscreen devices (tablets, phones) and false for mouse-driven
     * devices.  The result is used below to increase node radii so that carbon
     * atoms are easier to tap and drag with a finger.
     */
    const isTouch = window.matchMedia('(pointer: coarse)').matches;

    // Colour-blind-friendly palette
    const orange = 'rgba(255, 165, 0, 1)';
    const blue   = 'rgba(173, 216, 230, 1)';
    const green  = 'rgba(152, 251, 152, 1)';
    const yellow = 'rgba(255, 255, 0, 1)';
    const cyan   = 'rgba(0, 255, 255, 1)';
    const grey   = 'rgba(128, 128, 128, 1)';
    const black  = 'rgba(0, 0, 0, 1)';
    const white  = 'rgba(255, 255, 255, 1)';

    // Assign colours to carbon multiplicity types
    const CH3 = cyan;
    const CH2 = yellow;
    const CH  = green;
    const C   = orange;

    // Edge colours
    const cosyEdgeColor = blue;
    const hmbcEdgeColor = grey;

    // Node hover appearance
    const nodeHoverColor   = "lightgrey";
    const nodeHoverOpacity = 0.4;
    const textHoverOpacity = 0;

    /*
     * Node radii — larger on touch screens for easier tap/drag interaction
     * ─────────────────────────────────────────────────────────────────────────
     * Values are in SVG viewBox units (not screen pixels), so they scale
     * correctly with the rest of the diagram at any screen size.
     *   nodeRadius      — default resting size of each node circle
     *   nodeRadiusSmall — size applied to non-hovered neighbours on mouseover
     *   nodeRadiusLarge — size applied to the hovered node on mouseover
     */
    const nodeRadius      = isTouch ? 28 : 22;
    const nodeRadiusSmall = isTouch ? 20 : 15;
    const nodeRadiusLarge = isTouch ? 32 : 26;

    // Text sizing
    // Note: font-size (12px) and cosy/hmbc edge stroke-width (8px) are no
    // longer JS constants — they're set directly in molplot.css, driven by
    // the --inverse-scale custom property (see render() in view.js). Kept
    // here only as a value textOffsetInt still needs.
    const textOffsetInt = 8;

    // Colour scale: maps numProtons (0–3) → C / CH / CH2 / CH3 colours
    const colorScale = d3.scaleOrdinal()
        .domain([0, 1, 2, 3])
        .range([C, CH, CH2, CH3]);

