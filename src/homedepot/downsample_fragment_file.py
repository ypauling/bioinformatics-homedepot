"""CLI for downsampling fragment file"""

import csv
import gzip
import random
import collections
import argparse
import logging
import pandas as pd
from concurrent.futures import ThreadPoolExecutor
from typing import List, Tuple

from homedepot.logger import setup_logger
from homedepot.argutils import set_log_level_argument, get_log_level_argument

LOGGER_NAME = "ds_frag"


def parse_args() -> argparse.Namespace:
    """Parser for downsample-frag"""

    parser = argparse.ArgumentParser(
        description="Downsample multiple fragment files to the same coverage."
    )

    parser.add_argument(
        "-i",
        "--input_file",
        type=str,
        required=True,
        help="A tab delimited file with two columns: <input_frag>\t<output_frag>.",
    )
    parser.add_argument(
        "-t", "--threads", type=int, default=8, help="Number of threads."
    )
    set_log_level_argument(parser)

    return parser.parse_args()


def load_input_table(input_file: str) -> List[Tuple[int, int]]:
    """Read in the input file"""

    pairs = []
    logger = logging.getLogger(LOGGER_NAME)

    with open(input_file, "r") as f:
        for line in f:
            infile, outfile = line.strip().split()
            logger.debug("Get input: {}; output: {}.".format(infile, outfile))
            pairs.append((infile, outfile))

    return pairs


def count_fragments(file_path: str) -> int:
    """Count the total number of fragments in a fragment file"""

    total = 0
    with gzip.open(file_path, "rt") as f:
        for line in f:
            if line.startswith("#"):
                continue
            total += int(line.strip().split("\t")[4])

    logger = logging.getLogger(LOGGER_NAME)
    logger.debug("File {} has {} fragments.".format(file_path, total))

    return total


def downsample_file(input_path: str, output_path: str, N: int) -> None:
    """
    Downsample the fragment file

    Args:
        input_path (str): the input file path
        output_path (str): the output file path
        N (int): the downsample target number
    """

    logger = logging.getLogger(LOGGER_NAME)

    fragments = []
    logger.info("Reading {}.".format(input_path))
    with gzip.open(input_path, "rt") as f:
        for line in f:
            if line.startswith("#"):
                continue
            parts = line.strip().split("\t")
            count = int(parts[4])
            fragments.extend(tuple(parts[:4]) * count)
    logger.info("Finished reading {}.".format(input_path))

    logger.info("Downsampling {}.".format(input_path))
    sampled = random.sample(fragments, N)
    collapsed = collections.Counter(sampled)
    logger.info("Finished downsampling {}.".format(input_path))

    logger.info("Writing {} to {}.".format(input_path, output_path))
    with gzip.open(output_path, "wt") as out:
        writer = csv.writer(out, delimiter="\t", lineterminator="\n")
        for frag, count in collapsed.items():
            writer.writerow(list(frag) + [count])
    logger.info("Finished writing to {}.".format(output_path))


def main():

    args = parse_args()
    log_level = get_log_level_argument(args)
    logger = setup_logger(LOGGER_NAME, level=log_level)

    logger.debug(
        """Parameters:
        Input file: {},
        Number of threads: {}
        """.format(
            args.input_file, args.threads
        )
    )

    try:
        file_pairs = load_input_table(args.input_file)
    except Exception as e:
        logger.info("Unexpected error when reading {}: {}.".format(args.input_file, e))

    with ThreadPoolExecutor(max_workers=args.threads) as executor:
        counts = list(executor.map(lambda p: count_fragments(p[0]), file_pairs))

    N = min(counts)
    logger.info("Downsample target N: {}.".format(N))

    with ThreadPoolExecutor(max_workers=args.threads) as executor:
        args = [(infile, outfile, N) for infile, outfile in file_pairs]
        executor.map(lambda p: downsample_file(*p), args)
