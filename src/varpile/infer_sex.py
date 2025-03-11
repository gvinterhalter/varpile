"""
Determine the sex of the sample based on the homozygous and heterozygous ratio on non-autosomal
regions of chrX.

The rationale for ratio based calculation is based on the genetic differences
between males (XY genotype) and females (XX genotype) in the non-PAR regions of the X chromosome:
- Males, having only one copy of chromosome X (hemizygous for X), will typically show very few
  heterozygous genotypes in the non-PAR region but will have a higher count of homozygous
  or single-copy variants.
- Females, having two copies of chromosome X, will have a higher proportion of heterozygous
  genotypes in this region as they inherit one copy from each parent.

By calculating the fraction of heterozygous events relative to the total (heterozygous + homozygous),
the function determines if the sample corresponds to a male or female:
- A low heterozygous fraction (e.g., < 0.2) suggests a male (XY).
- A higher heterozygous fraction suggests a female (XX).

From ensembl GRCH38 (https://www.ensembl.org/info/genome/genebuild/human_PARS.html)
The Y chromosome has 5 regions (3rd and 5th are unique to Y chromosome)
chromosome:GRCh38:Y:1 - 10,000 is unique to Y but is a string of 10000 Ns
chromosome:GRCh38:Y:10,001 - 2,781,479 is shared with X: 10,001 - 2,781,479 (PAR1)
chromosome:GRCh38:Y:2,781,480 - 56,887,902 is unique to Y
chromosome:GRCh38:Y:56,887,903 - 57,217,415 is shared with X: 155,701,383 - 156,030,895 (PAR2)
chromosome:GRCh38:Y:57,217,416 - 57,227,415 is unique to Y

156,040,895
"""

from pathlib import Path
from typing import Literal
import pysam

# GRCh38 PAR regions

# PAR1 has the same coordinates for X and Y (coordinates are 1-based)
PAR1_BEGIN_1 = 10_001
PAR1_END_1 = 2_781_479

PAR2_X_BEGIN_1 = 155_701_383
PAR2_X_END_1 = 156_030_895

PAR2_Y_BEGIN_1 = 56_887_903
PAR2_Y_END_1 = 57_217_415


def in_non_par_X(pos_1: int):
    return PAR1_END_1 < pos_1 < PAR2_X_BEGIN_1


def in_non_par_Y(pos_1: int):
    """True if the position is outside of the PAR region.

    Args:
        pos_1 (int): The genomic position on the Y chromosome to evaluate (TODO: do we need begin/end).
    """
    # TODO: from ensembl description it would appear we need to look at > PAR2_Y_END_1
    #       Double check this, also, what about
    return PAR1_END_1 < pos_1 < PAR2_Y_BEGIN_1 or pos_1 > PAR2_Y_END_1


Sex = Literal["XX", "XY"]  # Sex is either "XX" or "XY"
# I understand gnomeAD does something more complex, but I haven't explored it yet.
# TODO: https://gnomad.broadinstitute.org/news/2023-11-gnomad-v4-0/#sex-inference
SamplesSex = dict[str, Sex]


FRACTION_LIMIT = 0.2  # 20 percent


def infer_sex(input_file: Path | str, sample_rank: int = 0) -> Sex:
    """Infer the sex ('XX' or 'XY') of the sample based on genotype data.

    This function analyzes the number of homozygous and heterozygous events
    in the non-PAR region of chromosome X to determine the sample's sex.

    Args:
        input_file (Path | str): Path to the VCF input file.
        sample_rank (int): The sample index within the VCF file to analyze. Defaults to 0.

    Returns:
        Sex: 'XX' if the sample is inferred to be female, or 'XY' if inferred to be male.
    """

    with pysam.VariantFile(str(input_file)) as f:

        # This is the large non-PAR
        X_non_par_region_iter = f.fetch("chrX", PAR1_END_1 - 1, PAR2_X_BEGIN_1)

        hom_event = 0
        het_event = 0
        for r in X_non_par_region_iter:
            s = r.samples[sample_rank]

            match s["GT"]:
                case (None,) | (None, None):
                    continue  # Ignore records with no genotype
                case (a, b) if a == b:
                    if a > 0:
                        hom_event += 1
                    elif a == 0:
                        # Count homozygous reference genotypes
                        hom_event += 1
                case (a, b):  # Treat other cases as heterozygous
                    het_event += 1
                case (a,) if a > 1:  # Treat hemizygous cases as heterozygous
                    het_event += 1
                case wtf:
                    raise ValueError("Unexpected genotype format:", wtf)

        # TODO: what if there are no variants (0 divided by 0)
        het_fraction = het_event / (hom_event + het_event)
        # print(hom_event, het_event,   "het_fraction:", het_fraction)

        if het_fraction < FRACTION_LIMIT:
            return "XY"  # Male
        else:
            return "XX"  # Female


def infer_samples_sex(input_file: Path | str) -> SamplesSex:
    with pysam.VariantFile(str(input_file)) as vcf:
        # get the samples in the file (we don't really need the names)
        samples: list[str] = list(vcf.header.samples)
        samples_sex = {sample: infer_sex(input_file, rank) for rank, sample in enumerate(samples)}
        return samples_sex
