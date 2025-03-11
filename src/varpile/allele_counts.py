import json
import shutil
from pathlib import Path
from typing import Final

import duckdb
import pysam

from varpile.infer_sex import SamplesSex, in_non_par_Y, in_non_par_X, PAR1_END_1, PAR2_X_BEGIN_1
from varpile.utils import OutFile, Region1

CON = duckdb.connect(":memory:")
CON.query("set threads to 1")

VARIANT_PILE_COLUMNS = {
    "pos": "INT",
    "ref": "VARCHAR",
    "alt": "VARCHAR",
    "XX_AC": "INT",
    "XX_AC_hom": "INT",
    "XX_AC_hemi": "INT",
    "XX_n_DP_discarded": "INT",
    "XY_AC": "INT",
    "XY_AC_hom": "INT",
    "XY_AC_hemi": "INT",
    "XY_n_DP_discarded": "INT",
    "DP": "INT",
}


def process_chromosome(
    vcf_path: Path, region: Region1, sex_info: SamplesSex, out_dir: Path, min_DP: int, debug: bool = False
):
    # define the location where we will save the chromosome data (out_path is treated as directory)
    variant_pile_path = out_dir / "data.parquet"

    out_file = OutFile(variant_pile_path, columns=VARIANT_PILE_COLUMNS)
    vcf = pysam.VariantFile(str(vcf_path))
    empty_values = "0\t0\t0\t0"  # used when there are no values (example sex=XX and we need to fill XY values)
    zero_counts = "0\t0\t0\t1"  # used when counts are zero
    with out_file, vcf:
        for (PASS, rec, sex, sample, dp), alt, (ac, ac_hom, ac_hemi) in iter_alleles(vcf, region, sex_info):
            if alt == "*":  # TODO: is this correct
                continue

            if dp >= min_DP:
                if PASS:
                    str_counts = f"{ac}\t{ac_hom}\t{ac_hemi}\t0"
                else:
                    str_counts = empty_values  # AC is 0 but we don't decrease AN
            else:
                str_counts = zero_counts

            if sex == "XX":
                XX_counts, XY_counts = (str_counts, empty_values)
            else:
                XX_counts, XY_counts = (empty_values, str_counts)

            line = f"{rec.pos}\t{rec.ref}\t{alt}\t{XX_counts}\t{XY_counts}\t{dp}\n"
            out_file.write_line(line)


def iter_alleles(vcf_file: pysam.VariantFile, region: Region1, sex_info: SamplesSex):

    if region.contig not in vcf_file.header.contigs:
        return

    vcf_records = vcf_file.fetch(*region.to_pysam_tuple())

    sex_list = list(sex_info.values())

    # constants for counts (AC, AC_hom, AC_hemi)
    zero_counts: Final = (0, 0, 0)  # count in case of no variants
    _hom_alt_counts: Final = (2, 1, 0)
    _het_counts: Final = (1, 0, 0)
    hemi_counts: Final = (1, 0, 1)

    is_chrX = region.contig in ("chrX", "X")
    is_chrY = region.contig in ("chrY", "Y")

    for record in vcf_records:

        alts = record.alts

        for sex, sample in zip(sex_list, record.samples.values()):

            # Depending on the sex and chromosome we change what these 2 counts are
            hom_alt_counts = _hom_alt_counts
            het_counts = _het_counts

            if is_chrY:
                if sex == "XX":
                    # TODO: Why do we have chrY variants for XX samples
                    continue
                if sex == "XY":
                    if in_non_par_Y(record.pos):
                        # we want to keep the variant, but it should be counted as hemizygote
                        hom_alt_counts = het_counts = hemi_counts
                    else:
                        continue  # ignore variants in PAR regions

            if is_chrX and sex == "XY":
                # for non-PAR region the variant is hemizygote, but in PAR region it's hom/het
                if in_non_par_X(record.pos):
                    hom_alt_counts = het_counts = hemi_counts

            try:
                # DP can be missing for various reason and in various ways:
                # - DP not present in format specification  (leads to KeyError)
                # - DP is '.' or sometimes None. (leads to ValueError or TypeError)
                #   Sometimes pysam returns None when in the vcf it's .
                dp = int(sample["DP"])
                # Exceptions can be that DP is "." o doesn't exist
            except Exception:  # noqa
                # If DP is missing we must ignore this variant for this sample
                # Most likely the variate does not exist for this sample in multisample vcf
                # or this is a structural variant or something which does not have DP
                continue

            # At this DP is defined for the variant, so we can be pretty sure GT is as well.
            gt = sample["GT"]

            try:
                GQ = int(sample["GQ"])
            except Exception:
                GQ = 0

            # custom filtering
            PASS = GQ >= 20

            common = [PASS, record, sex, sample, dp]

            match gt:
                # Diploid
                case (0, 0) | (None, None):
                    yield common, record.ref, zero_counts  # HOM_ref
                case (0, a) | (a, 0):  # HET
                    # Allelic balance filtering
                    try:
                        AD = sample["AD"]
                        AB = AD[0] / sum(AD)  #  Allele balance
                        if AB <= 0.2:
                            common[0] = False
                    except Exception:
                        common[0] = False  # PASS = false

                    yield common, alts[a - 1], het_counts
                case (a1, a2):
                    if a1 == a2:
                        # HOM_alt
                        yield common, alts[a1 - 1], hom_alt_counts
                    else:
                        # Multi allelic
                        yield common, alts[a1 - 1], het_counts
                        yield common, alts[a2 - 1], het_counts

                # Haploid
                case (0,) | (None,):
                    yield common, record.ref, zero_counts  # HEMI ref
                case (a,):
                    yield common, alts[a - 1], hemi_counts  # HEMI

                # TODO triploid or Multi-ploid... ?


def merge_piles(dir_path: Path, debug: bool = False) -> None:
    """Combine parquet files (piles of variants) into a single file containing counts.

    This is the first merge operation done to produce count datasets in a single center.
    """
    file_glob: str = str(Path(dir_path) / "*" / "data.parquet")
    out_path: str = str(dir_path / "data.parquet")
    rel = CON.query(
        f"""
        select 
        pos, ref, alt,
        XX_AC: sum(XX_AC)::int,
        XX_AC_hom: sum(XX_AC_hom)::int,
        XX_AC_hemi: sum(XX_AC_hemi)::int,
        XY_AC: sum(XY_AC)::int,
        XY_AC_hom: sum(XY_AC_hom)::int,
        XY_AC_hemi: sum(XY_AC_hemi)::int,
        -- DP stat counts
        XX_n_DP_discarded: sum(XX_n_DP_discarded)::int, -- number of samples that are DP discarded
        XY_n_DP_discarded: sum(XY_n_DP_discarded)::int, -- number of samples that are DP discarded
        n_samples: count(*)::int,  -- total number of samples (this is for DP statistics)
        DP_sum: sum(DP)::double,
        DP2_sum: sum(DP**2)::double,
        -- array_agg(DP) DPs, -- for debug
        from read_parquet('{file_glob}', hive_partitioning = false)
        group by pos, ref, alt
        order by pos, ref, alt
    """
    )

    rel.write_parquet(out_path, compression="ZSTD")

    if not debug:
        for file in dir_path.iterdir():
            if file.is_dir():
                shutil.rmtree(file)


def finalize_contig(in_dir: Path, contig: str, out_dir: Path):
    # get info from info.json
    info = json.loads((in_dir / "info.json").read_text())
    XX_sample_number = info["sample_number"]["XX"]
    XY_sample_number = info["sample_number"]["XY"]

    path: Path = in_dir / contig / "data.parquet"
    out_path: Path = out_dir / contig / "result.parquet"
    rel = CON.read_parquet(str(path))

    rel = rel.filter("XX_AC > 0 or XY_AC > 0")  # discard

    XX_n = f"({XX_sample_number} - XX_n_DP_discarded)"
    XY_n = f"({XY_sample_number} - XY_n_DP_discarded)"

    XX_multiplier = XY_multiplier = 2
    if contig == "chrM":
        XX_multiplier = XY_multiplier = 1
    elif contig == "chrY":
        # it's 1 because we only keep non-autosomal variants on XY samples
        XY_multiplier = 1
        XX_multiplier = 0
    elif contig == "chrX":
        # if not in PAR then 1 else 2
        XY_multiplier = f"if({PAR1_END_1} < pos and pos < {PAR2_X_BEGIN_1}, 1, 2)"

    XX_AN_computation = f"{XX_n} * {XX_multiplier}"
    XY_AN_computation = f"{XY_n} * {XY_multiplier}"

    rel = rel.select(
        f"""pos, ref, alt,
        XX_AN: {XX_AN_computation}, 
        XX_AC, XX_AC_hom, XX_AC_hemi,
        XY_AN: {XY_AN_computation},
        XY_AC, XY_AC_hom, XY_AC_hemi,
        DP_mean: (DP_sum/n_samples),
        DP_std: sqrt((DP2_sum - 2*DP_mean*DP_sum + n_samples*(DP_mean**2)) / n_samples),
        """
    )

    out_path.parent.mkdir(exist_ok=True)
    # rel.write_csv(str(out_path), sep="\t")
    rel.write_parquet(str(out_path), compression="ZSTD")
