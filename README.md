# Allostery is a widespread cause of loss-of-function variant pathogenicity

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
[![R](https://img.shields.io/badge/R-%3E%3D4.3-276DC3?logo=r)](https://www.r-project.org/)
[![bioRxiv](https://img.shields.io/badge/bioRxiv-2025.06.20.660737-b31b1b)](https://doi.org/10.1101/2025.06.20.660737)

Analysis code and metadata for:

> **Liao & Lehner (2025).** Allostery is a widespread cause of loss-of-function variant pathogenicity. *bioRxiv*. https://doi.org/10.1101/2025.06.20.660737

---

## Overview

Many disease-causing mutations occur far from protein active sites, but the mechanisms remain poorly understood. This repository provides the analysis code and metadata to reproduce all findings from our accompanying preprint. Here we introduce a framework to decouple stability-driven effects from functional mutational effects, allowing for the systematic detection of allostery in experimental or computational variant effect maps.

---

## Repository structure

```
.
├── src/                          # R Markdown analysis scripts (one per figure panel)
│   ├── revision_fig1_*.Rmd       # Fig 1 — ESM-1v R² across Domainome & VAMP-seq datasets
│   ├── revision_fig2_*.Rmd       # Fig 2 — GCK & PTEN dual-phenotype residual analysis
│   ├── revision_fig3*.Rmd        # Fig 3 — BRCA1, BAP1, RAD51C SGE benchmarks
│   ├── revision_fig4.Rmd         # Fig 4 — PDZ3 & SH3 residual decay analysis
│   ├── revision_fig5_*.Rmd       # Fig 5 — Proteome-wide decay (CASP1/CHK1/GCK/IDH/PDK1/PTP1B)
│   ├── supp_fig4.Rmd             # Supp Fig 4 — Extended PDZ3 & SH3 allosteric analysis
│   ├── revision_supp_*.Rmd       # Supplementary figures (AM, MSA, DCA, PhyloP comparisons)
│   ├── revision_*_esm2.Rmd       # ESM-2 vs ESM-1v residual comparisons
│   ├── revision_*_phylop.Rmd     # PhyloP conservation benchmarks
│   ├── revision_*_comp.Rmd       # Model comparison panels
│   └── revision_hypermutant.Rmd  # Hypermutant position analysis
│
├── lib/
│   ├── helpers.R                 # Shared helper functions (PDB mapping, LOESS fitting,
│   │                             #   exponential decay bootstrap, scatter/violin plots)
│   └── globals.R                 # Project-wide constants (paths, colour palettes, AA codes)
│
├── data/
│   ├── cleaned_ddg/              # Per-protein MoCHI folding & binding ddG tables
│   │   └── (pdz3, sh3, …)
│   ├── paper_supplements/        # Externally downloaded supplement tables
│   │   ├── domainome/            # Domainome 1.0 (Beltran et al. 2025)
│   │   ├── megascale/            # Mega-scale stability (Tsuboyama et al. 2023)
│   ├── proteome_meta/            # ClinVar variant table and UniProt metadata
│   ├── vampseq/                  # VAMP-seq datasets (7 full-length human proteins)
│   ├── scores/                   # ESM-1v, ESM-2, ThermoMPNN, AlphaMissense scores
│   ├── decay_pdb/ & residual_pdb/ # PDB files 
│   └── sasa/                     # Solvent-accessible surface area annotations
│
├── figs/                         # Output figure panels (PDF/PNG)
├── munge/                        # ProjectTemplate data-preprocessing scripts
├── cache/                        # ProjectTemplate object cache
└── src_hpc/                      # HPC batch scripts for large-scale computations
```

---

## Data access

Large input datasets for proteome-wide analyses are hosted on Zenodo:

| Dataset | Description | DOI |
|---|---|---|
| Proteome-wide raw data | ThermoMPNN predictions | [10.5281/zenodo.18381534](https://zenodo.org/records/18381534) |
| Proteome-wide meta data | Per-protein ESM-1v fitness, ThermoMPNN ddGf scores, and ClinVar annotations  | [10.5281/zenodo.18386427](https://zenodo.org/records/18386427) |

---

## Requirements

All analyses run in **R ≥ 4.3** using the [ProjectTemplate](http://projecttemplate.net/) framework, which handles library loading and data munging automatically via `load.project()`.

---

## Reproducing the analysis

1. **Clone this repository**

   ```bash
   git clone <repo-url>
   ```

2. **Open the R project** (`01.protein-seq-evo-v1.Rproj`) in RStudio or set the working directory to the project root, then run any analysis script:

   ```r
   library(ProjectTemplate)
   load.project()          # loads lib/, munge/, and caches data
   rmarkdown::render("src/revision_fig4.Rmd")
   ```

   Each `.Rmd` in `src/` is self-contained and documents its inputs, model steps, and outputs with markdown section headers.

---

## Citation

```bibtex
@article{liao2025allostery,
  title   = {Allostery is a widespread cause of loss-of-function variant pathogenicity},
  author  = {Liao, Xiaotian and Lehner, Ben},
  journal = {bioRxiv},
  year    = {2025},
  doi     = {10.1101/2025.06.20.660737}
}
```

---

## Contact

For questions, bug reports, or data requests please open a GitHub issue or contact **xl7@sanger.ac.uk**.
