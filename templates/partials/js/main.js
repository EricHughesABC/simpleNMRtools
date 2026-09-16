/*
 * main.js
 * ─────────────────────────────────────────────────────────────────────────
 * Loads last. The one piece of genuinely cross-cutting wiring — the global
 * keyboard shortcuts (R/L rotate, C toggle calculated shifts) — lives here
 * rather than in view.js/nodes.js individually, since it calls into both.
 * Escape-to-close for both modal dialogs is handled natively by the
 * browser now (see modals.js) and no longer needs a branch here.
 * Every reference below is inside the keydown callback itself, so it only
 * needs to resolve by the time a key is actually pressed — well after
 * every other file on this page has loaded — which is why this file has
 * no ordering requirements of its own beyond loading after the rest.
 */

// Add keydown event listener for rotation
d3.select("body").on("keydown", (event) => {
    if (event.key === "R" && event.shiftKey) {
        rotateSVG(90);
    }
    else if (event.key === "L" && event.shiftKey) {
        rotateSVG(-90);
    }
    else if (event.key === "r" || event.key === "R") {
        rotateSVG(10);
    }
    else if (event.key === "l" || event.key === "L") {
        rotateSVG(-10);
    }
    else if (event.key === "c" || event.key === "C") {
        // Toggle node labels between experimental and calculated shifts;
        // nodes without a calculated shift yet turn grey while it's showing.
        toggleCalculatedShifts();
    }
});

