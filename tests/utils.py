import textwrap
from pathlib import Path
from typing import Optional

import pysam


# This function is not currently used but is provided for cases
#  where a VariantHeader needs to be created from a dictionary.
def header_from_dict(d: dict[str, str | list]):
    # NOTE: pysam automatically chooses the format version
    header = pysam.VariantHeader()
    for key, value in d.items():
        if isinstance(value, str):
            header.add_meta(key, value)
        elif isinstance(value, list):
            if key == "samples":
                # special handling
                header.add_samples(value)
            else:
                for item in value:
                    header.add_meta(key, items=item.items())
        else:
            raise ValueError(f"Invalid value {value!r}, only str|list supported")
    return header


def write_vcf(path: Path, content: str, header: Optional[str] = None):
    """Write a proper VCF file given text that is using spaces.

    NOTE: When content is one string, it's not possible to use the ' ' as the value in the columns.
    The header part of the VCF we keep as is given.
    The columns need to be separated by \t instead of spaces.
    Number of expected columns is dictated by the #CHROM line (because there can be several samples)
    """
    content = textwrap.dedent(content)
    with open(path, "w") as f:
        if header:
            # In case we want to reuse a header
            header = textwrap.dedent(header)
            assert content.startswith("#CHROM")
            f.write(header.rstrip() + "\n")

        expected_number_of_columns = 0
        for i, line in enumerate(content.split("\n")):
            if not (line := line.strip()):
                continue  # skip empty lines

            if line.startswith("##"):
                f.write(line + "\n")
            else:
                columns = line.split()
                if line.startswith("#CHROM"):
                    expected_number_of_columns = len(columns)
                else:
                    assert len(columns) == expected_number_of_columns, f"Line {i}: {columns}"
                f.write("\t".join(line.split()) + "\n")
