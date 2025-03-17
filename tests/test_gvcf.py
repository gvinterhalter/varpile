import data
from varpile.VariantFile import VariantFile


def test_gvcf_parssing():
    with VariantFile(data.gvcf_broad_example) as vf:
        records_raw = list(vf.fetch("chr20"))
        assert len(records_raw) == 11  # in total 11 records

        records = []
        for r in records_raw:

            # discard GVCF blocks
            if r.alts[0] == "<NON_REF>":
                continue

            records.append(r)
            # TODO: figure out a way to write easy tests

        assert len(records) == 4
