# TCGA MYC / Lnc-DZANK1 analysis

This repository contains reproducible Python analyses of TCGA colorectal cancer (COAD + READ, reported as CRC) and TCGA esophageal cancer (ESCA).

## Analyses

1. `scripts/analyze_tcga_myc_lncdzank1.py` compares Lnc-DZANK1 expression between MYC RNA-expression high and low primary-tumor groups and calculates Spearman correlations.
2. `scripts/analyze_tcga_myc_cnv_lncdzank1.py` compares Lnc-DZANK1 expression by MYC GISTIC discrete copy-number status, including +2 amplification, +1/+2 gain sensitivity analysis, +2 versus diploid sensitivity analysis, MYC RNA validation, and copy-number dosage correlations.

Both analyses use the same GDC/Xena STAR-TPM matrices. The matrices are distributed as `log2(TPM + 1)` values; no raw counts, edgeR, or TMM normalization is used.

## Data sources

RNA expression matrices are downloaded from the GDC/Xena STAR-TPM resources:

- https://gdc.xenahubs.net/download/TCGA-COAD.star_tpm.tsv.gz
- https://gdc.xenahubs.net/download/TCGA-READ.star_tpm.tsv.gz
- https://gdc.xenahubs.net/download/TCGA-ESCA.star_tpm.tsv.gz

MYC CNV is downloaded at runtime from the cBioPortal TCGA PanCancer Atlas GISTIC discrete profiles:

- `coadread_tcga_pan_can_atlas_2018_gistic`
- `esca_tcga_pan_can_atlas_2018_gistic`

The discrete encoding is `-2` deep deletion, `-1` shallow deletion, `0` diploid, `+1` gain, and `+2` high-level amplification. MYC is queried using Entrez Gene ID 4609. See `results/TCGA_MYC_CNV_source.json` for the recorded source and encoding.

## Installation

```bash
python -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
```

## Running

The scripts accept an input directory containing `coad.star_tpm.tsv.gz`, `read.star_tpm.tsv.gz`, `esca.star_tpm.tsv.gz`, and `probemap.tsv`.

```bash
python scripts/analyze_tcga_myc_lncdzank1.py \
  --input-dir /path/to/tcga_input \
  --output-dir results

python scripts/analyze_tcga_myc_cnv_lncdzank1.py \
  --input-dir /path/to/tcga_input \
  --output-dir results
```

The CNV script requires internet access to query cBioPortal. It keeps only sample type `01` primary tumors, normalizes TCGA barcodes, matches RNA and CNV at the patient level, and retains one deterministic RNA sample per patient.

## Main outputs

- `results/TCGA_MYC_LncDZANK1_summary.csv`: RNA MYC high/low comparison.
- `results/TCGA_MYC_CNV_LncDZANK1_summary.csv`: CNV comparisons and dosage correlations.
- `results/TCGA_CRC_MYC_CNV_LncDZANK1_analysis.csv` and `results/TCGA_ESCA_MYC_CNV_LncDZANK1_analysis.csv`: matched sample-level data.
- `results/TCGA_MYC_CNV_QC.json`: sample counts and QC checks.
- `figures/`: PDF and 300-dpi PNG figures.

The bundled `results/` files are the outputs generated for the current input data and are provided for inspection. Raw TCGA expression matrices are not committed to this repository; download them before rerunning.
