"""
Logging configuration for simpleNMR using loguru.

Usage
-----
In each module that needs logging::

    from loguru import logger

At application startup (e.g. in simpleNMRtest_app.py)::

    from logging_config import setup_logging
    setup_logging()

Log levels used across the project
------------------------------------
logger.debug()    — diagnostic detail useful during development
                    (ppm values, index loops, intermediate DataFrame states)
logger.info()     — normal operational messages
                    (DB connected, experiment loaded, assignment complete)
logger.warning()  — unexpected but recoverable situations
                    (no solution found, missing optional data, fallback used)
logger.error()    — errors that prevent normal processing
                    (NMR solution failed, required data missing, mismatch)
"""

import sys
from pathlib import Path
from loguru import logger


def setup_logging(
    log_level: str = "DEBUG",
    log_file: str | None = None,
) -> None:
    """Configure loguru handlers for the application.

    Parameters
    ----------
    log_level:
        Minimum level to emit. Use "INFO" in production,
        "DEBUG" during development.
    log_file:
        Optional path for a rotating log file.  If None,
        only stderr output is configured.
    """
    logger.remove()  # Remove the default handler

    logger.add(
        sys.stderr,
        level=log_level,
        format=(
            "<green>{time:HH:mm:ss}</green> | "
            "<level>{level: <8}</level> | "
            "<cyan>{name}</cyan>:<cyan>{function}</cyan>:<cyan>{line}</cyan> — "
            "<level>{message}</level>"
        ),
        colorize=True,
    )

    if log_file:
        logger.add(
            log_file,
            level=log_level,
            rotation="10 MB",
            retention="30 days",
            encoding="utf-8",
            format=(
                "{time:YYYY-MM-DD HH:mm:ss} | "
                "{level: <8} | "
                "{name}:{function}:{line} — "
                "{message}"
            ),
        )
