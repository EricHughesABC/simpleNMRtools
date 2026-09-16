import argparse
import os

from app import create_app

app = create_app()

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Run the simpleNMR Flask dev server.")
    parser.add_argument(
        "--debug-html",
        action="store_true",
        help=(
            "Serve rendered HTML (d3molplotmnova_template.html) uncompressed, "
            "with comments intact, for easier debugging in an editor/DevTools. "
            "Sets SIMPLENMR_STRIP_HTML_COMMENTS=false for this process."
        ),
    )
    args = parser.parse_args()
    print(f"[run.py] args.debug_html = {args.debug_html}")

    if args.debug_html:
        os.environ["SIMPLENMR_STRIP_HTML_COMMENTS"] = "false"

    print(f"[run.py] env var right before app.run() = {os.environ.get('SIMPLENMR_STRIP_HTML_COMMENTS')!r}")

    app.run(debug=True, use_reloader=False)
        
