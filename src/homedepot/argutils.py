"""Shared argument parsing utilities"""

import argparse
import logging
from typing import Dict

LOG_LEVELS: Dict[str, int] = {
    "critical": logging.CRITICAL,
    "error": logging.ERROR,
    "warning": logging.WARNING,
    "info": logging.INFO,
    "debug": logging.DEBUG,
    "notset": logging.NOTSET,
}


def set_log_level_argument(parser: argparse.ArgumentParser) -> None:
    """
    Add a --log-level argument to the parser

    Args:
        parser (argparser.ArgumentParser): the parser to be extended
    """

    parser.add_argument(
        "--log-level",
        choices=LOG_LEVELS.keys(),
        default="info",
        help="Set logging level (default: info)",
    )


def get_log_level_argument(args: argparse.Namespace) -> int:
    """
    Extract the logging level from the parser

    Args:
        parser (argparse.Namespace): the given parser
    """

    return LOG_LEVELS[args.log_level]
