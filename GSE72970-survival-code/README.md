# GSE72970 Lnc-DZANK1 survival analysis

These scripts reproduce the overall survival (OS) and progression-endpoint analyses for Affymetrix probe `1561307_at` in GSE72970. The analysis computes the cohort median and assigns expression values above the median to High and values at or below the median to Low.

The progression endpoint uses the GEO fields `pfs` and `pfs censored`. Some output filenames and figures use the label PFI for this source endpoint.

## Run the analysis

From this folder, install the Python dependencies and download the series matrix from [NCBI GEO GSE72970](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE72970).

```sh
python3 -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements.txt
mkdir -p work/geo
curl -fL https://ftp.ncbi.nlm.nih.gov/geo/series/GSE72nnn/GSE72970/matrix/GSE72970_series_matrix.txt.gz -o work/geo/GSE72970_series_matrix.txt.gz
MPLBACKEND=Agg python work/analyze_gse72970.py
python work/analyze_gse72970_pfi.py
```

The scripts write Kaplan–Meier figures, log-rank and Cox results, and patient-level tables to `outputs/`. The input data and generated outputs are excluded from this package and ignored by Git.

## Render the Cox forest plots

The Python and R plotting scripts use the fixed summary estimates from the original analysis. They reproduce those figures independently of the cohort input. If you rerun the analysis with a different cohort or cutoff, review the estimates embedded in the plotting scripts before reusing them.

```sh
python work/create_forest_plot.py
python work/create_pfi_forest_plot.py
```

Alternatively, install the R packages `ggplot2`, `dplyr`, and `svglite`, then run the R plot scripts. Run either plotting implementation for a given endpoint to keep its own output files.

```sh
mkdir -p outputs
(cd outputs && Rscript ../work/Cox_OS_forest_plot.R)
(cd outputs && Rscript ../work/Cox_PFI_forest_plot.R)
```
