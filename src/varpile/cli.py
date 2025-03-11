"""Main entry point into the cli application."""

import argparse
from pathlib import Path

from varpile.errors import RegionError
from varpile.run import count, finalize

from varpile.utils import Region1


class ParseRegion(argparse.Action):
    """Converts a string representation of a region into a list of Region1 objects."""

    def __call__(self, parser, namespace, values, option_string=None):
        try:
            regions = [Region1.from_string(x) for x in values.split(",")]
        except RegionError as e:
            raise argparse.ArgumentError(self, str(e))
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
        action=ParseRegion,
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
