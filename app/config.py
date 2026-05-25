"""
Flask configuration for simpleNMR.

Environment variables
---------------------
SECRET_KEY                  Flask secret key (required in production)
DATABASE_URL                Full SQLAlchemy URI (takes priority if set)
DB_USERNAME                 MySQL/Postgres username  ]
DB_PASSWORD                 MySQL/Postgres password  ]  used together when
DB_HOST                     MySQL/Postgres host       ]  DATABASE_URL absent
DB_NAME                     MySQL/Postgres db name    ]
SQLALCHEMY_POOL_SIZE        Connection pool size        (default: 5)
SQLALCHEMY_POOL_TIMEOUT     Pool timeout in seconds     (default: 30)
SQLALCHEMY_POOL_RECYCLE     Recycle connections after N seconds (default: 1800)
SQLALCHEMY_MAX_OVERFLOW     Max overflow connections    (default: 2)

If none of the above database variables are set, falls back to a local
SQLite file ``registrations.db`` in the project root.
"""

import os
from pathlib import Path

from dotenv import load_dotenv
from loguru import logger

# Load .env from the project root (two levels up from app/config.py)
_env_path = Path(__file__).resolve().parent.parent / ".env"
_loaded = load_dotenv(_env_path)
logger.debug(f"[dotenv] .env path={_env_path}  loaded={_loaded}")


def _build_database_uri() -> str:
    """Resolve the database URI from environment variables.

    Priority order:
      1. DATABASE_URL  — full URI, works on any server
      2. DB_* variables — PythonAnywhere dashboard format
      3. SQLite fallback — used automatically for local development
    """
    # 1. Explicit full URI
    url = os.environ.get("DATABASE_URL", "")
    if url:
        logger.info(f"[db] Using DATABASE_URL (host: {url.split('@')[-1]})")
        return url

    # 2. Individual credentials (PythonAnywhere dashboard format)
    db_username = os.environ.get("DB_USERNAME", "")
    db_password = os.environ.get("DB_PASSWORD", "")
    db_host     = os.environ.get("DB_HOST", "")
    db_name     = os.environ.get("DB_NAME", "")

    if all([db_username, db_password, db_host, db_name]):
        logger.info(f"[db] Using DB_* variables (host: {db_host}/{db_name})")
        return f"mysql+pymysql://{db_username}:{db_password}@{db_host}/{db_name}"

    # 3. SQLite fallback
    sqlite_path = Path(__file__).resolve().parent.parent / "registrations.db"
    logger.info(f"[db] No DATABASE_URL/DB_* found — using SQLite: {sqlite_path}")
    return f"sqlite:///{sqlite_path}"


class Config:
    """Base configuration — reads from environment variables."""

    SECRET_KEY: str = os.environ.get("SECRET_KEY", "dev-secret-key-change-in-production")
    SQLALCHEMY_DATABASE_URI: str = _build_database_uri()
    SQLALCHEMY_TRACK_MODIFICATIONS: bool = False
    JSON_AS_ASCII: bool = False

    # Connection pool — only meaningful for MySQL/Postgres; ignored by SQLite
    SQLALCHEMY_POOL_SIZE: int    = int(os.environ.get("SQLALCHEMY_POOL_SIZE", 5))
    SQLALCHEMY_POOL_TIMEOUT: int = int(os.environ.get("SQLALCHEMY_POOL_TIMEOUT", 30))
    SQLALCHEMY_POOL_RECYCLE: int = int(os.environ.get("SQLALCHEMY_POOL_RECYCLE", 1800))
    SQLALCHEMY_MAX_OVERFLOW: int = int(os.environ.get("SQLALCHEMY_MAX_OVERFLOW", 2))


class DevelopmentConfig(Config):
    """Local development overrides."""
    DEBUG: bool = True


class ProductionConfig(Config):
    """PythonAnywhere / production overrides."""
    DEBUG: bool = False


class TestingConfig(Config):
    """In-memory SQLite for unit tests."""
    TESTING: bool = True
    SQLALCHEMY_DATABASE_URI: str = "sqlite:///:memory:"
