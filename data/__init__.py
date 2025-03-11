"""These datasets where used during development, they are not commited to the repo."""

from pathlib import Path

DATA_DIR = Path(__file__).parent

xy_vcf = DATA_DIR / "source" / "xy.vcf.gz"
xx_vcf = DATA_DIR / "source" / "xx.vcf.gz"
merged_vcf = DATA_DIR / "source" / "merged.vcf.gz"

vcfs = DATA_DIR / "vcfs"
