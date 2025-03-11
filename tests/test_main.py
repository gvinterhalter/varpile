"""
Input:
- folder or a list of files, vcf files.


Result:
 - vid
 - AN_male
 - AN_female
 - AC_male
 - AC_female


- 10 people (= 20 AN)
- 5 people homo
- 2 people het
AF = (5*2 + 2) / 20

pysam (but check pyvcf)


VCF files:
 - some will be gvcf
 - version of the file:  should be irrelevant (most will be 4.2)
    - troublesome spanning deletions                <-----------------
 - genotype determines hom/het

                                            Allele count
                                            hom        het      |  AC_f   N_hom_f  AC_m  N_hom_m
 (A, G)   ->  0/1        het                  0          1      |  1      0         0       0
 (A, G)   ->  1/1        hom                  2          0      |  2      1         0       0
 (A, G, C)   ->  1/2     het for both                           |
 (A, G)                                       0          1      |  1      0
 (A, C)                                       0          1      |  1      0

 chrX  (hemi)                                                   |
 chrY                                                           |
 chrM

Males hemi (X and Y)

Phasing
 (T, G)   ->  0|1        het
 (C, G)   ->  1|0        hom

#
1. # of heterozygous variants outside of pseudo-autosomal regions  on chrX
       < 5000 it's a male
       < 20%

2. Filtering on coverage 10x
"""

import duckdb
import pysam
import pysam.bcftools
import pytest

from tests.utils import write_vcf
from varpile.allele_counts import iter_alleles


# variants = [
#     ["chr1", 100, ".", "A", ("G",), None, "PASS", {"DP": 100}, {"GT": "0/1"}],
#     ["chr2", 101, ".", "A", ("G",), None, "PASS", {"DP": 100}, {"GT": "0/1"}],
# ]

# def test_something(tmp_path):
#     test_path = tmp_path / "test.vcf"
#     print(test_path)
#
#     # Dictionary to define header fields
#     header = {
#         "reference": "GRCh38",
#         "contig": [
#             {"ID": "chr1", "length": 248956422},
#         ],
#         "INFO": [
#             {"ID": "DP", "Number": "1", "Type": "Integer", "Description": "Total Depth"},
#             {"ID": "AF", "Number": "A", "Type": "Float", "Description": "Allele Frequency"},
#         ],
#         "FORMAT": [
#             {"ID": "GT", "Number": "1", "Type": "String", "Description": "Genotype"},
#             {"ID": "DP", "Number": "1", "Type": "Integer", "Description": "Read Depth"},
#         ],
#         "samples": ["SAMPLE1"],
#     }
#
#     header = header_from_dict(header)
#     with pysam.VariantFile(str(test_path), "w", header=header) as f:
#         header.new_record(
#             start=12345,
#             stop=12346,
#             alleles=('A', 'T'),
#             id='rs123',
#             qual=50,
#             filter='PASS',
#             info={'DP': 100},
#             samples=[{'GT': (0, 1)}],  # Heterozygous genotype for SAMPLE1
#         )


VCF_HEADER = """\
    ##fileformat=VCFv4.3    
    ##fileDate=20110413
    ##source=VCFtools
    ##reference=file:///refs/human_NCBI36.fasta
    ##contig=<ID=1,length=249250621,md5=1b22b98cdeb4a9304cb5d48026a85128,species="Homo Sapiens">
    ##contig=<ID=X,length=155270560,md5=1b22b98cdeb4a9304cb5d48026a85128,species="Homo Sapiens">
    ##INFO=<ID=AA,Number=1,Type=String,Description="Ancestral Allele">
    ##INFO=<ID=H2,Number=0,Type=Flag,Description="HapMap2 membership">
    ##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
    ##FORMAT=<ID=GQ,Number=1,Type=String,Description="Genotype Quality">
    ##FORMAT=<ID=DP,Number=1,Type=String,Description="Read Depth">
    ##ALT=<ID=DEL,Description="Deletion">
    ##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">
    ##INFO=<ID=END,Number=1,Type=Integer,Description="End position of the variant">
"""


@pytest.fixture(scope="module")
def vcf_directory(tmp_path_factory):
    dir_path = tmp_path_factory.mktemp("vcf_files")
    print(dir_path)
    return dir_path


@pytest.fixture(scope="module")
def example_vcf(vcf_directory):
    tmp_path = vcf_directory / "example.vcf"
    # Example from: The variant call format and VCFtools
    content = """\
    #CHROM   POS   ID     REF   ALT     QUAL   FILTER   INFO                 FORMAT     SAMPLE1   SAMPLE2
    1        1     .      ACG   A,AT    40     PASS     .                    GT:DP      1/1:13    2/2:29
    1        2     .      C     T,CT    .      PASS     H2;AA=T              GT         0/1       2/2
    1        5     rs12   A     G       67     PASS     .                    GT:DP      1|0:16    2/2:20
    1         6    .      C     A       255.7  .        .                    GT:DP:GQ    0/1:35:99    ./.:.:.
    X        100   .      T     <DEL>   .      PASS     SVTYPE=DEL;END=299   GT:GQ:DP   1:12:.    0/0:20:36
    """
    write_vcf(tmp_path, content, header=VCF_HEADER)
    return tmp_path


@pytest.fixture(scope="module")
def example_bcf(example_vcf):
    bcf_path = example_vcf.with_suffix(".bcf")
    with pysam.VariantFile(str(example_vcf)) as f:
        with pysam.VariantFile(str(bcf_path), "wb", header=f.header) as out:
            for r in f:
                out.write(r)

    pysam.bcftools.index("--csi", str(bcf_path))
    return bcf_path


def test_example_vcf_parsing(example_vcf):
    """Test that the example vcf  can be parsed."""
    with pysam.VariantFile(str(example_vcf)) as f:
        for _ in f:
            pass


# def test_allele_count_example(example_vcf):
#     """Test that the genotype can be parsed."""
#     values = []
#     with pysam.VariantFile(str(example_vcf)) as f:
#         for r, allele, counts, _dp in iter_alleles(f):
#             values.append((r.chrom, r.pos, r.ref, allele, counts))
#
#     assert values == [
#         #(AC AC_hom AC_hemi)
#         ('1', 1, 'ACG', 'A', (2, 1, 0)),
#         ('1', 1, 'ACG', 'AT', (2, 1, 0)),
#         ('1', 2, 'C', 'T', (1, 0, 0)),
#         ('1', 2, 'C', 'CT', (2, 1, 0)),
#         ('1', 5, 'A', 'G', (1, 0, 0)),
#         ('1', 5, 'A', 'G', (1, 0, 0)),
#         ('X', 100, 'T', '<DEL>', (1, 0, 1)),
#     ]


# def test_variants(tmp_path):
#     tmp_path /= "example.vcf"
#     content = """\
#         #CHROM   POS   ID     REF   ALT     QUAL   FILTER   INFO FORMAT     SAMPLE1   SAMPLE2
#         1        1     .      ACG   A,AT    40     .        .    GT:DP      0/1:13    2/0:29
#         1        2     .      C     T,CT    .      PASS     .    GT:DP      1/1:20    2/2:15
#         1        3     .      A     C,G     .      PASS     .    GT:DP      0/0:20    1/2:15
#         1        5     .      A     G       67     .        .    GT:DP      1|0:16    2/2:20
#         X        3     .      GGGCG G       .      .        .    GT:DP      0/1:13    ./.:29
#         X        5     .      G     *,G     .      .        .    GT:DP      0/1:13    ./.:29
#         X        100   .      T     <DEL>   .      .        .    GT:GQ:DP   1:12:.    0/0:20:36
#     """
#     write_vcf(tmp_path, content, header=VCF_HEADER)
#
#     values = []
#     with pysam.VariantFile(str(tmp_path)) as f:
#         for r, allele, counts in iter_alleles(f):
#             values.append(" ".join(map(str, [allele, *counts])))
#     print(values)
#
#     assert values == [
#         # Allele AC AC_hom AC_hemi
#         "A 1 0 0",
#         "AT 1 0 0",
#         "AT 1 0 0",
#         "T 2 1 0",
#         "CT 2 1 0",
#         "C 1 0 0",
#         "G 1 0 0",
#         "G 1 0 0",
#         "G 1 0 0",
#         "G 1 0 0",
#         "* 1 0 0",
#         "<DEL> 1 0 1",
#     ]


# @pytest.mark.parametrize(
#     "vcf_line, expected",
#     [  #  CHROM POS ID REF     ALT     QUAL FILTER INFO FORMAT   SAMPLE1  SAMPLE2
#         ("1     1   .  ACG     A,AT    .    .      .    GT:DP    0/1:13   2/0:29", ["A 1 0 0 13", "AT 1 0 0 29"]),
#         ("1     2   .  C       T,CT    .    PASS   .    GT:DP    1/1:20   2/2:15", ["T 2 1 0 20", "CT 2 1 0 15"]),
#         ("1     3   .  A       C,G     .    PASS   .    GT:DP    0/0:20   1/2:15", ["C 1 0 0 15", "G 1 0 0 15"]),
#         # 2/2:2 for SAMPLE2 is ignored as alt 2 does not exist
#         ("1     5   .  A       G       .    .      .    GT:DP    1|0:16   2/2:20", ["G 1 0 0 16"]),
#         ("1     6   .  C       G       .    .      .    GT:DP    1|0:16   2/2:20", ["G 1 0 0 16"]),
#         (  # spanning deletions
#             "X  10  .  GGGCG   G       .    .      .    GT:DP    1/1:13   ./.:29\n"
#             "X  12  .  G       *,G     .    .      .    GT:DP    1/2:13   0/0:29\n",
#             ["G 2 1 0 13", "* 1 0 0 13", "G 1 0 0 13"],  #                                            <------ TODO
#         ),
#         # ("X     100 .  T       <DEL>   .    .      .    GT:GQ:DP 1:12:.   0/0:20:36", []),
#     ],
# )
# def test_iter_alleles(tmp_path, vcf_line, expected):
#     tmp_path /= "test.vcf"
#     content = "#CHROM POS ID REF ALT QUAL FILTER INFO FORMAT SAMPLE1 SAMPLE2\n" + vcf_line
#     write_vcf(tmp_path, content, header=VCF_HEADER)
#
#     values = []
#     with pysam.VariantFile(str(tmp_path)) as f:
#         for r, allele, counts, dp in iter_alleles(f):
#             values.append(" ".join(map(str, [allele, *counts, dp])))
#
#     assert values == expected
#
#
# def test_vcf_to_allele_counts(example_bcf, tmp_path):
#     out_file = tmp_path / "out"
#     convert_vcf_to_allele_counts(example_bcf, out_file)
#     df = duckdb.read_parquet(str(out_file) + "/**/*.parquet").pl()
#     # parse the
#     # TODO: how do we assert this?
