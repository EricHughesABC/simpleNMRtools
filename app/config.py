"""
Flask configuration for simpleNMR.

Environment variables
---------------------
SECRET_KEY          Flask secret key (required in production)
DATABASE_URL        Full SQLAlchemy URI (takes priority if set)
DB_HOST             PostgreSQL host  ]
DB_NAME             PostgreSQL db    ]  used together when DATABASE_URL is absent
DB_USER             PostgreSQL user  ]
DB_PASSWORD         PostgreSQL pass  ]

If none of the above are set, falls back to a local SQLite file
``registrations.db`` in the project root.
"""

import os
from pathlib import Path

from dotenv import load_dotenv, find_dotenv
from loguru import logger

# Load .env from the project root (two levels up from this file)
_env_path = Path(__file__).resolve().parent.parent / ".env"
_loaded = load_dotenv(_env_path)
logger.debug(f"Config: .env path={_env_path}  loaded={_loaded}")


def _build_database_uri() -> str:
    """Resolve the database URI from environment variables."""
    # 1. Explicit full URI
    url = os.environ.get("DATABASE_URL", "")
    if url:
        logger.info(f"Database: using DATABASE_URL (host: {url.split('@')[-1]})")
        return url

    # 2. Individual PostgreSQL credentials
    host = os.environ.get("DB_HOST", "")
    name = os.environ.get("DB_NAME", "")
    user = os.environ.get("DB_USER", "")
    password = os.environ.get("DB_PASSWORD", "")
    if host and name:
        logger.info(f"Database: using DB_* variables (host: {host}/{name})")
        return f"postgresql://{user}:{password}@{host}/{name}"

    # 3. SQLite fallback
    sqlite_path = Path(__file__).resolve().parent.parent / "registrations.db"
    logger.info(f"Database: using SQLite fallback ({sqlite_path})")
    return f"sqlite:///{sqlite_path}"


class Config:
    """Base configuration — reads from environment variables."""

    SECRET_KEY: str = os.environ.get("SECRET_KEY", "dev-secret-key-change-in-production")
    SQLALCHEMY_DATABASE_URI: str = _build_database_uri()
    SQLALCHEMY_TRACK_MODIFICATIONS: bool = False
    JSON_AS_ASCII: bool = False


class DevelopmentConfig(Config):
    """Local development overrides."""
    DEBUG: bool = True


class ProductionConfig(Config):
    """PythonAnywhere / production overrides."""
    DEBUG: bool = False


class TestingConfig(Config):
    """In-memory SQLite for tests."""
    TESTING: bool = True
    SQLALCHEMY_DATABASE_URI: str = "sqlite:///:memory:"
