import json
import logging
import shutil
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
from typing import TypedDict, Optional, Final, final

from tqdm import tqdm

import varpile
from varpile.allele_counts import process_chromosome, merge_piles, finalize_contig
from varpile.infer_sex import infer_samples_sex, SamplesSex

logger = logging.getLogger(__name__)
info = logger.info

# default chromosomes to use (if regions are not specified)
CHROMOSOMES: final = [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY", "chrM"]

# We only support indexed vcf so only gziped vcf or bcf file should be supported
SUPPORTED_EXTENSIONS: Final = [".vcf.gz", ".vcf.bgz", ".bcf"]


def is_vcf(path: Path) -> bool:
    return any(path.name.endswith(ext) for ext in SUPPORTED_EXTENSIONS)


def find_input_files(paths: list[Path]) -> list[Path]:
    """Find and return a list of VCF (Variant Call Format) input files."""
    input_files = []
    for path in paths:
        if path.is_dir():
            input_files.extend(file for file in path.iterdir() if is_vcf(file))
        else:
            input_files.append(path)

    return input_files


def get_vcf_file_name(input_file: Path | str) -> str:
    """
    Extracts the base file name from the given file path by removing any recognized
    VCF-related file extensions.
    Args:
        input_file: The file path (as Path or str) to extract the name from

    Returns: The base file name without the recognized VCF-related extensions.
    """
    name = Path(input_file).name

    for ext in SUPPORTED_EXTENSIONS:
        if name.endswith(ext):
            return name.removesuffix(ext)

    # Return the name as-is if no matching extension is found (try why not)
    logger.warning("Unknown VCF extension '%s'", name)
    return name


class IOptions(TypedDict):
    """Interface for Arguments."""

    paths: list[Path]
    output: Path
    regions: Optional[list]
    min_DP: int
    threads: int
    debug: bool


def count(opt: IOptions) -> None:

    input_files: list[Path] = find_input_files(opt["paths"])
    output: Path = opt["output"]
    threads = opt["threads"]
    regions = opt["regions"] or CHROMOSOMES
    min_DP = opt["min_DP"]
    debug = opt["debug"]

    with ProcessPoolExecutor(threads) as executor:

        info(f"Infer sex of input files")
        futures = {executor.submit(infer_samples_sex, input_file): input_file for input_file in input_files}
        vcf_sex_info: dict[Path, SamplesSex]
        vcf_sex_info = {futures[future]: future.result() for future in tqdm(futures, desc="Inferring sex")}

        sample_number = defaultdict(int)  # number of XX, and XY samples
        for sex_info in vcf_sex_info.values():
            for sex in sex_info.values():
                sample_number[sex] += 1

        info(f"Identified {sample_number["XX"]} XX samples and {sample_number["XY"]} XY.")

        # Clean output directory (in case it already exists so we can cleanly overwrite data)
        if output.exists():
            if output.is_dir():
                shutil.rmtree(output)
            else:
                output.unlink()

        output.mkdir()  # create the Output directory

        output_info = output / "info.json"
        output_info.write_text(
            json.dumps(
                {
                    "version": varpile.__VERSION__,
                    "min_DP": min_DP,
                    "sample_number": sample_number,
                }
            )
        )

        info(f"Processing chromosomes/regions:")
        for region in regions:
            region_dir = output / region
            futures = []
            for input_file in input_files:

                # infer sex for each sample in the vcf.
                sex_info = vcf_sex_info[input_file]

                # determine the output directory
                file_name: str = get_vcf_file_name(input_file)
                file_output = region_dir / file_name
                file_output.mkdir(parents=True, exist_ok=True)

                # Submit the task to the process pool
                futures.append(
                    executor.submit(
                        process_chromosome,
                        input_file,
                        region,
                        sex_info,
                        file_output,
                        min_DP,
                        debug=debug,
                    )
                )

            # Wait for all tasks to be completed
            for future in tqdm(futures, desc=region):
                future.result()

            info(f"Merging files for {region}")
            merge_piles(region_dir, debug=debug)


def finalize(in_path: Path, out_path: Path, threads: int):

    if out_path.exists():
        shutil.rmtree(out_path)

    out_path.mkdir(exist_ok=True)
    info_path = in_path / "info.json"
    shutil.copyfile(info_path, out_path / "info.json")

    # Get a list of all directories under in_path
    directories = [path for path in in_path.iterdir() if path.is_dir()]

    # Process the directories in parallel using ProcessPoolExecutor
    with ProcessPoolExecutor(threads) as executor:
        futures = [executor.submit(finalize_contig, in_path, directory.name, out_path) for directory in directories]
        for future in tqdm(futures, desc="Finalizing dataset"):
            future.result()
