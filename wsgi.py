"""
WSGI entry point for PythonAnywhere.

Update your PythonAnywhere WSGI configuration file to point at this file::

    import sys
    sys.path.insert(0, '/home/<username>/simpleNMRtools')
    from wsgi import application

This replaces the previous direct import of simpleNMRtest_app.
Only switch PythonAnywhere to use this file after Phase 3 is complete
(routes extracted into app/routes.py and registered in app/__init__.py).
"""

from logging_config import setup_logging
from app import create_app

setup_logging(log_level="INFO")
application = create_app()
