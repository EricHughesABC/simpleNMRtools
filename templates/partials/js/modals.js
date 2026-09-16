/*
 * modals.js
 * ─────────────────────────────────────────────────────────────────────────
 * The two popups, both native <dialog> elements: the Help-button controls
 * reference (Shift/Ctrl+click) and the Solution Statistics popup.
 * showModal()/close() hand centering, the dimmed backdrop, and
 * Escape-to-close to the browser natively — the one thing dialogs don't
 * do on their own is close on a backdrop click, which wireModal() below
 * adds for both. No dependency on any other custom file — safe to load
 * anywhere after the DOM exists.
 */

// Wires a dialog's close button and backdrop-click-to-close, and returns
// the dialog element so callers can open it with .showModal().
function wireModal(dialogId, closeButtonId) {
    const dialog = document.getElementById(dialogId);
    document.getElementById(closeButtonId).addEventListener("click", () => dialog.close());
    // A click on the dimmed backdrop lands on the <dialog> element itself —
    // there's no separate backdrop element to target — while a click on
    // anything inside the dialog hits that inner element instead.
    dialog.addEventListener("click", (event) => {
        if (event.target === dialog) {
            dialog.close();
        }
    });
    return dialog;
}

const controlsModal = wireModal("controls-modal", "controls-modal-close");
const statsModal = wireModal("statsModal", "statsModalClose");

/*
 * Help button — a plain click opens the documentation page as before;
 * Shift+click or Ctrl+click instead shows the controls reference dialog.
 */
d3.select("#helpButton").on("click", (event) => {
    if (event.shiftKey || event.ctrlKey) {
        controlsModal.showModal();
    } else {
        window.open('http://simplenmr.pythonanywhere.com/documentation/index.html', '_blank');
    }
});

// Called from the "Solution Statistics" button's onclick in the template
function openStatsModal() {
    statsModal.showModal();
}
