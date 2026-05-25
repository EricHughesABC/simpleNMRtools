"""
WSGI entry point for PythonAnywhere.

Update your PythonAnywhere WSGI configuration file to replace:

    from simpleNMRtest_app import app as application

with:

    import sys
    sys.path.insert(0, '/home/<username>/simpleNMRtools')
    from wsgi import application
"""

from log_setup import setup_logging
from app import create_app

setup_logging(log_level="INFO")
application = create_app()
