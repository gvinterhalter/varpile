import json
import shutil
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import duckdb
from tqdm import tqdm

from varpile.infer_sex import PAR1_END_1, PAR2_X_BEGIN_1
from varpile.utils import Region1


def finalize(in_path: Path, out_path: Path, threads: int):

    if out_path.exists():
        shutil.rmtree(out_path)

    out_path.mkdir(exist_ok=True)
    info_path = in_path / "info.json"
    shutil.copyfile(info_path, out_path / "info.json")

    # Get a list of all directories under in_path
    directories = [path for path in in_path.iterdir() if path.is_dir()]
    regions = [Region1.from_string(x.name) for x in directories]

    # Process the directories in parallel using ProcessPoolExecutor
    with ProcessPoolExecutor(threads) as executor:
        futures = [executor.submit(finalize_region, in_path, region, out_path) for region in regions]
        for future in tqdm(futures, desc="Finalizing dataset"):
            future.result()


def finalize_region(in_dir: Path, region: Region1, out_dir: Path):
    # get info from info.json
    info = json.loads((in_dir / "info.json").read_text())
    XX_sample_number = info["sample_number"]["XX"]
    XY_sample_number = info["sample_number"]["XY"]

    path: Path = in_dir / str(region) / "data.parquet"
    out_path: Path = out_dir / str(region) / "result.parquet"
    contig = region.contig

    con = duckdb.connect()
    con.query("set threads=1")

    rel = con.read_parquet(str(path))

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
