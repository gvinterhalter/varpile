import logging
import math
from pathlib import Path
import contextlib
from typing import Iterator

import pysam

from varpile.utils import Region1, is_test

logger = logging.getLogger(__name__)


class VariantFile(contextlib.AbstractContextManager):

    def __init__(self, file_path: Path | str):
        self.file_path: Path = Path(file_path)
        self._handle = pysam.VariantFile(file_path)
        self.header = self._handle.header

    def __exit__(self, exc_type, exc_value, traceback, /):
        self._handle.close()

    def fetch(self, region: Region1 | str) -> Iterator[pysam.VariantRecord]:
        """Return an iterator over the records."""

        if isinstance(region, str):
            region = Region1.from_string(region)

        if region.contig not in self.header.contigs:
            logger.warning(f"Skiping, contig '{region.contig}' not found in the header")
            return iter([])

        if is_test():
            # Inside  the test we can fake that the file is indexed
            if not self._handle.index:
                print("No index")
                # Index is not present
                contig = region.contig
                begin = region.begin or 0
                end = region.end or math.inf
                return (
                    record
                    for record in self._handle.fetch()
                    if record.contig == contig and begin < record.stop and record.start < end
                )

        # Otherwise try to fetch (this will fail if no index was detected
        # This is ok, we want things to fail, force the user to have an index
        return self._handle.fetch(*region.to_pysam_tuple())
