import shutil
from dataclasses import dataclass
from io import StringIO
from pathlib import Path

import duckdb


class OutFile:
    """Class that abstracts a parquet file

    Output is the parquet file, but we don't write parquet directly,
    we first write a temporary tsv file which we then convert to parquet file.
    This was just easier.
    We don't have to convert to parquet, but it might help out during the merge step.
    """

    BUFFER_LIMIT = 1_000_000  # bytes

    def __init__(self, file_path: Path, columns: dict) -> None:
        """

        Args:
            file_path: resulting parquet file
            columns: dict of the form name: type (type is duckdb SQL type)
        """
        self.output_path = file_path  # parquet file
        self.columns = columns
        self.tmp_path = Path(file_path.parent / "tmp.tsv")  # temporary file that we will later convert

        self.buffer = StringIO()
        self.file_handle = open(self.tmp_path, "w")  # Open file in truncate mode

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        try:
            self.flush()
            self._convert_to_parquet()
        except Exception:
            raise
        finally:
            self.file_handle.close()
            if not exc_type:  # Only delete temp file if no exception occurred
                self.tmp_path.unlink()

    def write_line(self, line) -> None:
        self.buffer.write(line)
        if self.buffer.tell() > OutFile.BUFFER_LIMIT:
            self.flush()

    def flush(self):
        """Flush the internal buffer to the file."""
        if self.buffer.tell() > 0:
            self.buffer.seek(0)
            self.file_handle.write(self.buffer.read())
            self.buffer = StringIO()

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


def decode_region_string(region: str) -> tuple[str, int, int] | tuple[str, int, None] | tuple[str, None, None]:
    # Split the region. The pattern is contig[:begin[-end]]
    try:
        contig, begin, end = None, None, None
        if ":" in region:
            contig, coords = region.split(":", 1)
            if "-" in coords:
                begin, end = map(int, coords.split("-"))
            else:
                begin = int(coords)
        else:
            contig = region
    except Exception:
        raise ValueError(f"Invalid {region=}")

    return (contig, begin, end)
