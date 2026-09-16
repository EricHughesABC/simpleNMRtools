r"""
strip_comments.py
─────────────────────────────────────────────────────────────────────────
Strips comments (and, as a side effect, the incidental whitespace they
leave behind) from a fully-rendered HTML page before it's sent to a
client — without touching the source template/CSS/JS files, which keep
their comments untouched on disk.

Used in routes.py, at the end of the existing rtn_html post-processing
chain in both simpleMNOVAfinalHTML() and _render_and_save() — the two
places render_template("d3molplotmnova_template.html", ...) is called.
Runs after the True/False and np.float64/np.int64 cleanup already there,
so it only ever sees valid, final JS/JSON, never an intermediate Python
repr fragment:

    rtn_html = rtn_html.replace("True", "true").replace("False", "false")
    rtn_html = re.sub(r"np\.float64\(([\d\.]+)\)", r"\1", rtn_html)
    rtn_html = re.sub(r"np\.int64\(([\d]+)\)", r"\1", rtn_html)
    rtn_html = strip_comments.maybe_strip_for_delivery(rtn_html)

Toggling it off for debugging
─────────────────────────────────────────────────────────────────────────
maybe_strip_for_delivery() checks the SIMPLENMR_STRIP_HTML_COMMENTS
environment variable rather than a CLI flag directly, because the app
starts two different ways that don't share an argv:
  - locally, `python run.py` goes through normal argument parsing
  - on PythonAnywhere, the WSGI server imports wsgi.py -> application
    directly; nothing on that path ever sees a command line
An environment variable is the one thing both paths can read. run.py's
--debug-html flag sets it before the dev server starts; on PythonAnywhere
set it directly via the Web tab's Environment variables section (same
place DB_USERNAME etc. already live there) and Reload the web app.

Why rjsmin/rcssmin and not a hand-rolled regex
─────────────────────────────────────────────────────────────────────────
This page's JS legitimately contains "//" inside string literals — e.g.
window.open('http://simplenmr.pythonanywhere.com/...') and the embedded
background SVG's xmlns='http://www.w3.org/2000/svg' attributes. A naive
line-based "strip from // to end of line" approach would truncate those
lines and corrupt the page. rjsmin/rcssmin are mature, widely-used pure
-Python packages (no Node.js or other runtime needed) built specifically
to track string/regex context while stripping, so they don't fall into
that trap — verified against this exact file's risky lines before this
was written.

Install: pip install rjsmin rcssmin
"""

import os
import re
import rjsmin
import rcssmin

_STYLE_BLOCK = re.compile(r"<style>(.*?)</style>", re.S)
_SCRIPT_BLOCK = re.compile(r"<script>(.*?)</script>", re.S)
_HTML_COMMENT = re.compile(r"<!--.*?-->", re.S)

_ENV_VAR = "SIMPLENMR_STRIP_HTML_COMMENTS"
_ENV_VAR = "SIMPLENMR_STRIP_HTML_COMMENTS"


def comments_enabled() -> bool:
    """
    True unless SIMPLENMR_STRIP_HTML_COMMENTS is set to a falsy value
    ("0", "false", or "no", case-insensitive). Stripping is the default
    (production) behaviour — the variable only needs setting at all when
    you want to turn it off for local debugging.

    Reads the environment fresh on every call rather than caching at
    import time, so it doesn't matter whether run.py sets the variable
    before or after this module gets imported — only that it's set
    before the first request that should see the effect.
    """
    value = os.environ.get(_ENV_VAR, "")
    print(f"SIMPLENMR_STRIP_HTML_COMMENTS value: {value}")
    return value.strip().lower() not in ("0", "false", "no")


def strip_for_delivery(html: str) -> str:
    """
    Returns a copy of `html` with CSS comments, JS comments, and HTML
    comments removed, and the whitespace rjsmin/rcssmin naturally trim
    along with them. Only touches inline <style>...</style> and
    <script>...</script> blocks with no src attribute (i.e. this page's
    own bundled CSS/JS) — a <script src="..."> tag like the d3.js CDN
    include has no body and is left untouched.

    Unconditional — strips every time, regardless of
    SIMPLENMR_STRIP_HTML_COMMENTS. Call maybe_strip_for_delivery()
    instead unless you specifically want that.
    """
    html = _STYLE_BLOCK.sub(lambda m: f"<style>{rcssmin.cssmin(m.group(1))}</style>", html)
    html = _SCRIPT_BLOCK.sub(lambda m: f"<script>{rjsmin.jsmin(m.group(1))}</script>", html)
    html = _HTML_COMMENT.sub("", html)
    return html


def maybe_strip_for_delivery(html: str) -> str:
    """
    strip_for_delivery(html), unless SIMPLENMR_STRIP_HTML_COMMENTS has
    disabled it (see comments_enabled()) — in which case `html` is
    returned unchanged, comments and all, for debugging in an
    editor/DevTools. This is what routes.py calls.
    """
    print(f"comments_enabled: {comments_enabled()}")
    if comments_enabled():
        return strip_for_delivery(html)
    return html
