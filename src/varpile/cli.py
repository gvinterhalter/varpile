"""Main entry point into the cli application."""

import argparse
from pathlib import Path

from varpile.run import count, finalize
import re


class Region(argparse.Action):
    """Converts a string representation of a region into a list of strings."""

    def __call__(self, parser, namespace, values, option_string=None):
        regions: list[str] = values.split(",")

        # Validate regions are correct
        region_pattern = re.compile(r"^(\w+)(:\d+|:\d+-\d+)?$")
        for region in regions:
            m = region_pattern.match(region)
            if not m:
                raise argparse.ArgumentError(
                    self, f"Invalid region format: '{region}'. Expected format: contig[:start[-stop]]"
                )

            if m.group(2):
                numbers: list[str] = m.group(2)[1:].split("-")
                start = int(numbers[0])
                if len(numbers) > 1:
                    stop = int(numbers[1])
                    if stop < start:
                        raise argparse.ArgumentError(self, f"Invalid region format: '{region}', stop < start")

        setattr(namespace, self.dest, regions)


def make_parser() -> argparse.ArgumentParser:
    """Parses command-line arguments for the CLI."""
    parser = argparse.ArgumentParser(description="Varpile CLI Application")
    subparsers = parser.add_subparsers(dest="action", required=True, help="Subcommands for varpile")

    # Count subcommand
    count_parser = subparsers.add_parser("count", help="Count allele frequencies from VCF files")
    count_parser.add_argument("paths", nargs="+", type=Path, help="Provide one or several paths to file or directory")
    count_parser.add_argument("-o", "--output", type=Path, required=True, help="Specify the output directory path")
    count_parser.add_argument(
        "-r",
        "--regions",
        type=str,
        action=Region,
        help="comma separated regions of form contig[:begin[-end]] (1-based) ",
    )
    count_parser.add_argument("--min-DP", type=int, default=10, help="Variants with lower DP are discarded")
    count_parser.add_argument("-@", "--threads", type=int, default=1, help="Number of threads to use (default 1)")
    count_parser.add_argument("--debug", action="store_true", help="Enable debug mode that preserves per sample output")

    # Merge subcommand
    merge_parser = subparsers.add_parser("merge", help="Merge multiple count datasets into one")

    # Finalize subcommand
    finalize_parser = subparsers.add_parser("finalize", help="Finalize the output (calculate AN, DP_mean and DP_std")
    finalize_parser.add_argument("path", type=Path, help="Input counts")
    finalize_parser.add_argument("-o", "--output", type=Path, required=True, help="Output directory")
    finalize_parser.add_argument("-@", "--threads", type=int, default=1, help="Number of threads to use (default 1)")

    # we can use this to conform to a type
    # actions = {a.dest: a.type for a in parser._actions}
    # print(actions)
    return parser


def main():
    args = make_parser().parse_args()
    opt = args.__dict__

    action = opt.pop("action")

    if action == "count":
        count(opt)
    elif action == "finalize":
        finalize(opt["path"], opt["output"], opt["threads"])
    elif action == "merge":
        pass  # TOOD


if __name__ == "__main__":
    main()
