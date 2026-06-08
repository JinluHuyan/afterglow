# TTDB CSV placeholder

This directory is a placeholder for local TTDB CSV files. Development checkouts may contain local CSV snapshots, but public source tarballs intentionally exclude the large TTDB CSV tables through `.Rbuildignore`. For manuscript reproduction, use `写作/resources/manuscript_input_manifest.tsv` to check the expected TTDB file names and configure a local TTDB directory with `set_ttdb_dir()` or `AFTERGLOW_TTDB_DIR`.

For public release, the large TTDB CSV tables are intentionally excluded from the source tarball by `.Rbuildignore`:

```text
^inst/extdata/ttdb/.*\.csv$
```

Download TTDB from:

- https://sysbio.gzzoc.com/ttdb/download.html

Expected file names:

- `computed_feature.csv`
- `species_stability_0927.csv`
- `species_stability_combined_0927_2_no_threshold.csv`
- `species_top_100_gene_correlations.csv`
- `stuty_info_0913.csv`

Recommended use after installing afterglow:

```r
library(afterglow)
set_ttdb_dir("~/data/ttdb")
ttdb <- read_ttdb_downloads()
```

or pass the directory explicitly:

```r
ttdb <- read_ttdb_downloads("~/data/ttdb")
```

Non-interactive pipelines can set:

```bash
export AFTERGLOW_TTDB_DIR=/path/to/ttdb
```

Then call:

```r
ttdb <- read_ttdb_downloads()
```
