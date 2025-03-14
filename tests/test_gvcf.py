import data
from varpile.VariantFile import VariantFile


def test_gvcf_parssing():
    with VariantFile(data.gvcf_broad_example) as vf:
        records = list(vf.fetch("20"))
        for r in records:

            # discard GVCF blocks
            if r.alts[0] == "<NON_REF>":
                continue
            pass
            # TODO: figure out a way to write easy tests
