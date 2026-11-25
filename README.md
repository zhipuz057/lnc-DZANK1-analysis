# lnc-DZANK1 — Bioinformatics Analysis

## Overview

This repository contains bioinformatics analysis scripts, processed results, and computational figures associated with our study of **lnc-DZANK1, PURA, BRCA1 transcription, and homologous recombination (HR) in colorectal cancer**.

The repository is organized by individual figures or analysis modules to facilitate data inspection, methodological transparency, and reproducibility.

Each analysis folder is maintained independently and may include the scripts, input data references, processed outputs, and figures associated with a specific analysis.

## Repository Organization

Bioinformatics analyses are organized into separate folders according to their corresponding figures in the manuscript.

Example directory structure:

```text
lnc-DZANK1-analysis/
│
├── README.md
│
├── Figure_01/
│   ├── scripts/
│   ├── results/
│   └── figures/
│
├── Figure_02/
│   ├── scripts/
│   ├── results/
│   └── figures/
│
├── Figure_S01/
│   ├── scripts/
│   ├── results/
│   └── figures/
│
└── ...
```

Folders will be added as the corresponding analyses are prepared and verified. Figure numbering and directory contents will reflect the final manuscript organization.

### Folder Contents

Each figure-specific folder may contain:

- **scripts/** — R, Python, or other scripts used for computational analysis and figure generation.
- **results/** — Processed datasets, statistical summaries, and intermediate analysis outputs.
- **figures/** — Generated plots and visualizations.
- **README.md** — Analysis-specific documentation, including data sources, software dependencies, parameters, and execution instructions.

Not all folders require every subdirectory. Their contents depend on the analyses performed.

## Data Availability

The following sequencing and proteomics datasets associated with this study have been deposited in public scientific repositories:

| Dataset | Repository | Accession |
|---|---|---|
| lnc-DZANK1 knockdown RNA-seq | NCBI GEO | GSE355541 |
| c-MYC knockdown RNA-seq | NCBI GEO | GSE355542 |
| PURA ChIP-seq | NCBI GEO | GSE355832 |
| RNA pull-down LC-MS/MS | ProteomeXchange / PRIDE | PXD087936 |

The datasets may remain private during manuscript peer review. Reviewer access information is provided separately to the journal editorial office.

## Reproducibility

Analysis-specific documentation will be provided alongside the corresponding scripts, including relevant software versions, input data sources, computational parameters, and execution instructions where applicable.

Large raw sequencing and mass spectrometry files are maintained in their respective scientific repositories rather than duplicated on GitHub.

## Repository Status

This repository is being populated with the computational analyses associated with the manuscript. Additional figure-specific folders and analysis documentation will be added as they are prepared and verified.

## Contact

For questions regarding this repository or the associated analyses, please contact the corresponding authors of the manuscript.
