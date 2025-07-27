"""CLI extracting information from GTF/GFF to bed"""

import argparse
import logging
import sys
import time
import re
import pandas as pd
from typing import Dict, Tuple, List
from homedepot.logger import setup_logger
from homedepot.argutils import set_log_level_argument, get_log_level_argument

LOGGER_NAME = "gtf-to-bed"


def parse_args() -> argparse.Namespace:
    """Parser for gtf-to-bed"""

    parser = argparse.ArgumentParser(
        description="Extract and write sequence information from GTF/GFF to bed format"
    )

    parser.add_argument(
        "-i", "--input_file", type=str, required=True, help="Input GTF/GFF file name."
    )
    parser.add_argument(
        "--input_type",
        type=str,
        choices=["GTF", "GFF"],
        default="GTF",
        help="Inut format (default GTF)",
    )
    parser.add_argument(
        "-o", "--output_file", type=str, required=True, help="Output bed file name."
    )
    parser.add_argument(
        "-t",
        "--type",
        type=str,
        default="gene",
        help="What feature to extract (default: gene).",
    )
    parser.add_argument(
        "-u",
        "--upstream",
        type=int,
        default=2000,
        help="Distance to extract from upstream of TSS, can be negative (default: 2000).",
    )
    parser.add_argument(
        "--upstream_anchor",
        type=str,
        choices=["start", "end"],
        default="start",
        help="Which end to apply the upstream distance? (default: start, meaning 5' end)",
    )
    parser.add_argument(
        "-d",
        "--downstream",
        type=int,
        default=2000,
        help="Distance to extract from downstream of TSS (defalt: 2000)",
    )
    parser.add_argument(
        "--downstream_anchor",
        type=str,
        choices=["start", "end"],
        default="start",
        help="Which end to apply the downstream distance? (default: start, meaning 5' end)",
    )
    parser.add_argument(
        "-c",
        "--chrom_size_file",
        type=str,
        required=True,
        help="The size of each chromsome.",
    )
    parser.add_argument(
        "-r",
        "--attr_id",
        type=str,
        default="gene_id",
        help="The attr tag to extract record name (default: gene_id)",
    )
    set_log_level_argument(parser)

    return parser.parse_args()


def write_bed(df: pd.DataFrame, output_file: str, columns: List[str]) -> None:
    """
    Write the resulting dataframe to the output file

    Args:
        df (pd.DataFrame): the input dataframe
        output_file (str): the output file path
        columns (List[str]): the columns to write
    """

    logger = logging.getLogger(LOGGER_NAME)
    logger.debug("Writing to file {}.".format(output_file))
    df[columns].to_csv(output_file, header=False, sep="\t", index=False)
    logger.debug("Writing to bed finished.")


def read_gtf_input(input_file: str) -> pd.DataFrame:
    """
    Read the input GTF/GFF file

    Args
        input_file (str): the input file path

    Returns:
        pd.DataFrame: the input GFF/GTF dataframe
    """

    gtf_columns = [
        "chrom",
        "source",
        "feature",
        "start",
        "end",
        "score",
        "strand",
        "fname",
        "attribute",
    ]
    dtypes = {
        "chrom": str,
        "source": str,
        "feature": str,
        "start": int,
        "end": int,
        "score": str,
        "strand": str,
        "fname": str,
        "attribute": str,
    }
    try:
        df = pd.read_csv(
            input_file,
            sep="\t",
            header=None,
            comment="#",
            names=gtf_columns,
            dtype=dtypes,
            low_memory=False,
        )
    except Exception as e:
        raise e

    return df


def read_chrom_size_file(input_file: str) -> Dict[str, int]:
    """
    Args:
        input_file (str): the input file path

    Returns:
        Dict[str, str]: a dict containing mapping from chrom name to chrom size
    """

    try:
        df = pd.read_csv(
            input_file,
            sep="\t",
            header=None,
            comment="#",
            names=["name", "size"],
            dtype={"name": str, "size": int},
        )
    except Exception as e:
        raise e

    try:
        result_dict = dict(zip(df["name"], df["size"]))
    except Exception as e:
        raise e

    return result_dict


def process_records(
    row: pd.Series,
    updist: int,
    upanchor: str,
    downdist: int,
    downanchor: str,
    chrom_df: Dict[str, int],
) -> Tuple[int, int]:
    """
    Process each row to get the desired coordinates

    Args:
        row (pd.Seires): the record row
        updist (int): the upstream distance, can be negative
        upanchor (str): the upstream anchor, start means 5', end means 3'
        downdist (int): the downstream distance, can be negative
        downanchor (str): the downstream anchor, start means 5', end means 3'
        chrom_df (Dict[str, int]): the mapping from chrom name to size

    Returns:
        Tuple[int, int]: the new start ane end coordinates
    """

    min_end = 0
    max_end = chrom_df[row["chrom"]]

    start = 0
    end = 0
    if row["strand"] == "." or row["strand"] == "+":
        if upanchor == "start":
            if updist > 0:
                start = max(min_end, row["start"] - updist)
            else:
                start = min(row["start"] - updist, max_end)
        else:
            if updist > 0:
                start = max(min_end, row["end"] - updist)
            else:
                start = min(row["end"] - updist, max_end)

        if downanchor == "start":
            if downdist > 0:
                end = min(row["start"] + downdist, max_end)
            else:
                end = max(min_end, row["start"] + downdist)
        else:
            if downdist > 0:
                end = min(row["end"] + downdist, max_end)
            else:
                end = max(max_end, row["end"] + downdist)

    else:
        if upanchor == "start":
            if updist > 0:
                end = min(row["end"] + updist, max_end)
            else:
                end = max(min_end, row["end"] + updist)
        else:
            if updist > 0:
                end = min(row["start"] + updist, max_end)
            else:
                end = max(min_end, row["start"] + updist)

        if downanchor == "start":
            if downdist > 0:
                start = max(min_end, row["end"] - downdist)
            else:
                start = min(row["end"] - downdist, max_end)
        else:
            if downdist > 0:
                start = max(min_end, row["start"] - downdist)
            else:
                start = min(row["start"] - downdist, max_end)

    return start, end


def parse_attr(attrstr: str, format: str) -> Dict[str, str]:
    """
    Parse the attr string based on file format

    Args:
        attrstr (str): the attribute string
        format (str): the file format, can be either 'GTF' or 'GFF'

    Returns:
        Dict[str, str]: the resulting tag-value pair
    """
    attrs = {}

    if format == "GTF":
        pattern = re.compile(r'\s*(\S+)\s+"([^"]+)"\s*;?')
        for match in pattern.finditer(attrstr):
            key, val = match.groups()
            attrs[key] = val
    else:
        pattern = re.compile(r"(\S+?)=([^;]+);?")
        for match in pattern.finditer(attrstr):
            key, val = match.groups()
            attrs[key] = val

    return attrs


def main():
    args = parse_args()
    log_level = get_log_level_argument(args)
    logger = setup_logger(LOGGER_NAME, level=log_level)

    logger.debug(
        """Parameters:
        Input file: {},
        Input format: {}, 
        Output file: {},
        Type of information: {},
        Upstream distance: {},
        Upstream anchor: {},
        Downstream distance: {},
        Downstream anchor: {},
        Chrom size file: {},
        Attr tag: {}
        """.format(
            args.input_file,
            args.input_type,
            args.output_file,
            args.type,
            args.upstream,
            args.upstream_anchor,
            args.downstream,
            args.downstream_anchor,
            args.chrom_size_file,
            args.attr_id,
        )
    )

    stime = time.perf_counter()

    try:
        input_df = read_gtf_input(args.input_file)
    except Exception as e:
        logger.info(
            "Unexpected error when reading input file %s: %s.", args.input_file, e
        )
        sys.exit(1)
    logger.debug("Finish reading input: %s.", args.input_file)
    logger.debug(f"Reading input takes: {time.perf_counter() - stime:.3f}s.")

    stime = time.perf_counter()
    try:
        chrom_dict = read_chrom_size_file(args.chrom_size_file)
    except Exception as e:
        logger.info(
            "Unexpected error when reading chrom size file %s: %s.",
            args.chrom_size_file,
            e,
        )
        sys.exit(1)
    logger.debug("Finish reading chrom size file: %s.", args.chrom_size_file)
    logger.debug(f"Reading chrom size file takes: {time.perf_counter() - stime:.3f}s.")

    input_df = input_df[input_df["chrom"].isin(chrom_dict.keys())].copy()
    input_df = input_df[input_df["feature"] == args.type].copy()

    stime = time.perf_counter()
    try:
        input_df[["ns", "ne"]] = input_df.apply(
            lambda row: process_records(
                row,
                args.upstream,
                args.upstream_anchor,
                args.downstream,
                args.downstream_anchor,
                chrom_dict,
            ),
            axis=1,
            result_type="expand",
        )
    except Exception as e:
        logger.info("Error when processing coordinates: %s.", e)
        sys.exit(1)
    logger.debug("Finish transform coordinates.")
    logger.debug(f"Transform coordinates takes: {time.perf_counter() - stime:.3f}s.")

    stime = time.perf_counter()
    try:
        input_df["attr_dict"] = input_df["attribute"].apply(
            lambda x: parse_attr(x, args.input_type)
        )
        print(input_df["attribute"])
        input_df["name"] = input_df["attr_dict"].apply(lambda d: d.get(args.attr_id))
    except Exception as e:
        logger.info("Error when extract %s from attributes: %s.", args.attr_id, e)
        sys.exit(1)
    logger.debug("Finish extract %s information.", args.attr_id)
    logger.debug(f"Extract id information takes: {time.perf_counter() - stime:.3f}s.")

    input_df = input_df[input_df["ns"] < input_df["ne"]].copy()

    stime = time.perf_counter()
    try:
        write_bed(
            input_df, args.output_file, ["chrom", "ns", "ne", "name", "score", "strand"]
        )
    except Exception as e:
        logger.info("Error when writing result to %s.", args.output_file)
        sys.exit(1)
    logger.debug("Finished...")
    logger.debug(f"Writing takes: {time.perf_counter() - stime:.3f}s.")
