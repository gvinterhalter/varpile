from io import StringIO
from pathlib import Path


class OutFile:
    BUFFER_LIMIT = 1_000_000  # bytes

    def __init__(self, file_path: Path | str, clean=True) -> None:
        self.data_path = file_path
        self.clean = False
        self.parquet_path = file_path.with_suffix(".parquet")
        self.buffer = StringIO()
        self.file_handle = open(file_path, "w")  # Open file in truncate mode

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.flush()
        self.file_handle.close()
        self.convert_to_parquet()

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

    def convert_to_parquet(self) -> None:
        from varpile.allele_counts import data_tsv_to_parquet

        data_tsv_to_parquet(self.data_path, self.data_path.with_suffix(".parquet"))
        if self.clean:
            self.data_path.unlink()


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
