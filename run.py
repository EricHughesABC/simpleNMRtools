"""
Local development entry point.

    python run.py

Do NOT use this for PythonAnywhere — use wsgi.py there.
"""

from logging_config import setup_logging

setup_logging(log_level="DEBUG")

# During Phase 1/2: import directly from the existing flat app file.
# Phase 3: replace with the app factory once routes are in app/routes.py.
#
#   from app import create_app
#   app = create_app()
#
from simpleNMRtest_app import app  # noqa: E402

if __name__ == "__main__":
    app.run(debug=True, host="0.0.0.0", port=5000)
