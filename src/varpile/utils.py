import os
import sys
import re
import shutil
from dataclasses import dataclass
from pathlib import Path
from typing import ClassVar

import duckdb

from varpile.errors import RegionError


class OutFile:
    """Class that abstracts a parquet file

    Output is the parquet file, but we don't write parquet directly,
    we first write a temporary tsv file which we then convert to parquet file.
    This was just easier.
    We don't have to convert to parquet, but it might help out during the merge step.
    """

    def __init__(self, file_path: Path, columns: dict) -> None:
        """

        Args:
            file_path: resulting parquet file
            columns: dict of the form name: type (type is duckdb SQL type)
        """
        self.output_path = file_path  # parquet file
        self.columns = columns
        self.tmp_path = Path(file_path.parent / "tmp.tsv")  # temporary file that we will later convert

        block_size = os.statvfs(file_path.parent).f_bsize
        self.file_handle = open(self.tmp_path, "w", buffering=block_size)

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        try:
            self.file_handle.close()
            self._convert_to_parquet()
        except Exception:
            raise
        finally:
            if not exc_type:  # Only delete temp file if no exception occurred
                self.tmp_path.unlink()

    def write_line(self, line) -> None:
        self.file_handle.write(line)

    def _convert_to_parquet(self) -> None:
        con = duckdb.connect()
        con.query("set threads=1")  # No need to multithread

        rel = con.query(
            f"""FROM read_csv('{str(self.tmp_path)}', columns={self.columns},
                               HEADER=FALSE, DELIM='\t', HIVE_PARTITIONING=FALSE, AUTO_DETECT=FALSE );
            """
        )

        # # # TODO: we might need to do some binning (this is how we could do it)
        # rel = rel.select("*, pos // 100_000_000 as bin")
        # rel.to_parquet(str(self.output_path), compression="ZSTD", partition_by=["bin"])
        # flatten_dir(self.output_path)

        rel.to_parquet(str(self.output_path), compression="ZSTD")


def flatten_dir(dir_path: Path) -> None:
    """Move the contents of dir_path to dir_path.parent and remove dir_path."""
    for path in dir_path.iterdir():
        new_path = dir_path.parent / path.name
        if new_path.exists():
            if new_path.is_dir():
                shutil.rmtree(new_path)
            else:
                path.unlink()
        shutil.move(path, dir_path.parent)
    dir_path.rmdir()


@dataclass(frozen=True)
class Region1:
    """Genomic region (1-based)."""

    contig: str
    begin: int | None
    end: int | None

    _REGION_PATTERN: ClassVar = re.compile(r"^(\w+)(:\d+|:\d+-\d+)?$")

    def __post_init__(self):
        self._validate_coordinates()

    def _validate_coordinates(self):
        if self.begin is None and self.end is not None:
            raise RegionError(f"Invalid region: {repr(self)}")
        if self.begin is not None and self.end is not None:
            if self.begin > self.end:
                raise RegionError(f"Invalid region: '{self}', begin < end")

    def __str__(self) -> str:
        s = self.contig
        if self.begin is not None:
            s += f":{self.begin}"
        if self.end is not None:
            s += f"-{self.end}"
        return s

    @classmethod
    def from_string(cls, region: str) -> "Region1":
        """Create a Region1 instance from a string.

        Args:
            s (str): Region string in the format 'contig[:start[-stop]]'.

        Returns:
            Region1: An instance of Region1 parsed from the string.

        Examples:
            >>> Region1.from_string("chr3")
            Region1(contig='chr3', begin=None, end=None)
            >>> Region1.from_string("chr2:150")
            Region1(contig='chr2', begin=150, end=None)
            >>> Region1.from_string("chr1:100-200")
            Region1(contig='chr1', begin=100, end=200)
        """

        # Validate region pattern is correct
        if not cls._REGION_PATTERN.match(region):
            raise RegionError(f"Invalid region: '{region}'. Expected format: contig[:start[-stop]]")

        try:
            begin, end = None, None
            if ":" in region:
                contig, coords = region.split(":", 1)
                if "-" in coords:
                    begin, end = map(int, coords.split("-"))
                else:
                    begin = int(coords)
            else:
                contig = region
        except Exception:
            raise RegionError(f"Invalid region: '{region}'. Expected format: contig[:start[-stop]]")

        return cls(contig.strip(), begin, end)

    def to_pysam_tuple(self) -> tuple:
        """Return a tuple and convert to 0-based coordinates."""
        begin0 = None if self.begin is None else self.begin - 1
        return (self.contig, begin0, self.end)


def is_test() -> bool:
    """Return True if the unit tests are executing."""
    return "pytest" in sys.modules
