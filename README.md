
# Install and Running

Install **uv**, [instructions](https://docs.astral.sh/uv/getting-started/installation/#standalone-installer).

Enter the `varpile` root directory, then run `uv run varpile --help`.
This will do several things:
  - it will install correct version of Python if not already present
  - create a `.venv` virtual environment in the directory
  - install necessary packages and varpile in venv

To be able to run `varpile` without `uv run` and outside the `varpile` root directory you can
activate the virtual environment with varpile, run `source .venv/bin/activate`.

# Workflow

There are three actions (sub-commands) of varpile: `count`, `merge`, `finlize`.
- Every datacenter runs `varpile count` on its vcf files. This will produce a directory with parquet
files that contain per-chromosome variant counts. This is not the final file.
- The count files (directory) from all data-centers are merged using action `merge`.
- Finally, action `finalize` calculates AN, mean DP and std DP for every variant, and it filters out
  variants with 0 allele counts (they have been kept for DP statistics).
  Each center can run this on their datasets, but this is not the result that is feed to `merge`.

In the output directory `info.json` file will keep track of the sample number (per sex) and other meta data.

# varpile count

`varpile count a.vcf.gz b.bcf directory_a  directory_b/**/*.vcf.gz -o <output_path>  -@ 4`

- All input files need to be indexed.
- `bcf`, `vcf.gz`, `vcf.bgz` files accepted.
- Above we provide a list of 2 files a directory and then another list of files (globing is done by shell)
- In case of directory all accepted files in the directory will be processed.
- `-@` is specifying the number of threads to use. Multithreading is done per file.
- Filtering:
  - We do filtering on Depth (DP), Allelic balance (AB) and Genotype Quality(GQ)
  - If the variants fall any of these criteria it's not present in the allele count (we still keep track of DP values)
  - The default values are `--min-DP 10 --min-GQ 20  --min-AB  0.2` 

Processing can be limited to chromosomes or regions using the `-r`, `--regions`.
Format is familiar comma separated `'chromosome[:start[-stop]]'` (1-based position).


# varpile merge

Not implemented


# varpile finalize

`varpile finalize in_dir -o out_dir -@ 4`


# Inspecting the parquet files
To view the parquet file add this helper method to your `.bash_profile` or `.bashrc`.

NOTE: this requires duckdb (which is present in the virtual environment).

```bash
function pcat() {
    if [ -n "$2" ]; then
      filter_expr="WHERE $2"
    else
      filter_expr=""
    fi

    duckdb -c "
        SET threads TO 1;
        COPY (
            SELECT *
            FROM read_parquet('$1')
            $filter_expr
        ) TO '/dev/stdout' (FORMAT csv, DELIMITER '\t', HEADER);
    "
}
```

The `pcat` function (parquet cat) prints out a parquet file in `tsv` format. 

Examples:
  - `pcat <in.parquet> | less`
  - `pcat <in.parquet> "XX_AC > 5" | less`  additional filtering

Example output:\
Note that for each variant we split the counts for XX and XY samples.
```tsv
pos	ref	alt	XX_AN	XX_AC	XX_AC_hom	XX_AC_hemi	XY_AN	XY_AC	XY_AC_hom	XY_AC_hemi	DP_mean	DP_std
270787	T	C	8	2	1	0	0	0	0	0	10.0	0.0
270906	G	C	8	1	0	0	0	0	0	0	56.0	0.0
270907	G	A	8	1	0	0	0	0	0	0	55.0	0.0
271075	G	A	8	4	1	0	0	0	0	0	45.666666666666664	6.599663291074451
```


``