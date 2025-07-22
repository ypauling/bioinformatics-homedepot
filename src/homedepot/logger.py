"""Logging configuration module"""

import logging


def setup_logger(
    name: str,
    level: int = logging.INFO,
    fmt: str = "[%(asctime)s] %(levelname)s:%(name)s: %(message)s",
    datefmt: str = "%Y-%m-%d %H:%M:%S",
    stream: bool = True,
) -> logging.Logger:
    """
    Create and configure a logger.

    Args:
        name (str): The name of the logger, usually __name__ of the caller.
        level (int): Logging level.
        fmt (str): Log message format string.
        datefmt (str): Format of timestamps.
        stream (bool): Whether to add a StreamHandler.

    Returns:
        logging.Logger: Configured Logger instance
    """

    logger = logging.getLogger(name)
    logger.setLevel(level)
    if not logger.hasHandlers():
        formatter = logging.Formatter(fmt=fmt, datefmt=datefmt)
        if stream:
            console_handler = logging.StreamHandler()
            console_handler.setFormatter(formatter)
            logger.addHandler(console_handler)
    return logger
