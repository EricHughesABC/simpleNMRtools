"""
Local development entry point.

    python run.py

Do NOT use this for PythonAnywhere — use wsgi.py there.
"""

from log_setup import setup_logging

setup_logging(log_level="DEBUG")

from app import create_app  # noqa: E402
from app.config import DevelopmentConfig

app = create_app(DevelopmentConfig)

if __name__ == "__main__":
    app.run(debug=True, host="0.0.0.0", port=5000)
