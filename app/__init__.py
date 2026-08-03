"""
Flask application factory for simpleNMR.

Usage
-----
Development::

    python run.py

PythonAnywhere (wsgi.py)::

    from app import create_app
    application = create_app()
"""

from flask import Flask
from loguru import logger

from app.extensions import db, migrate
from app.config import Config


def create_app(config_class: type = Config) -> Flask:
    """Create and configure the Flask application.

    Parameters
    ----------
    config_class:
        Configuration class to use.  Defaults to ``Config`` (reads from
        environment variables).  Pass a subclass for testing.

    Returns
    -------
    Flask
        Fully configured application instance.
    """
    app = Flask(__name__, template_folder="../templates")
    app.config.from_object(config_class)

    # Initialise extensions
    db.init_app(app)
    migrate.init_app(app, db)

    # Import models so Flask-Migrate picks them up for schema tracking
    from app.models import User, Device, Result  # noqa: F401

    with app.app_context():
        db.create_all()
        logger.info("Database tables ready")

    # Register routes blueprint
    from app.routes import bp as main_bp
    app.register_blueprint(main_bp)

    # Shell context — makes db and models available in `flask shell`
    @app.shell_context_processor
    def make_shell_context() -> dict:
        return {"db": db, "User": User, "Device": Device, "Result": Result}

    # Backward-compatibility shim: templates use url_for('endpoint') but
    # routes now live in the 'main' blueprint, so endpoints are
    # 'main.endpoint'.  Override url_for in the Jinja2 globals to try
    # the blueprint-prefixed name automatically so no template changes
    # are needed.
    _original_url_for = app.jinja_env.globals["url_for"]

    def _compat_url_for(endpoint: str, **values):
        if "." not in endpoint and endpoint != "static":
            try:
                return _original_url_for(f"main.{endpoint}", **values)
            except Exception:
                pass
        return _original_url_for(endpoint, **values)

    app.jinja_env.globals["url_for"] = _compat_url_for

    return app
