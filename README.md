# Allostery is a widespread cause of loss-of-function variant pathogenicity

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
[![R](https://img.shields.io/badge/R-%3E%3D4.3-276DC3?logo=r)](https://www.r-project.org/)
[![bioRxiv](https://img.shields.io/badge/bioRxiv-2025.06.20.660737-b31b1b)](https://doi.org/10.1101/2025.06.20.660737)

Analysis code and metadata for:

> **Liao & Lehner (2025).** Allostery is a widespread cause of loss-of-function variant pathogenicity. *bioRxiv*. https://doi.org/10.1101/2025.06.20.660737

---

## Overview

This repository reproduces all main and supplementary figures from the above preprint. The central framework decomposes single-amino-acid variant effects from deep mutational scanning (DMS) experiments into a **stability component** and a **residual functional component** that is unexplained by stability alone. The residual signal is enriched at ligand-contacting orthosteric sites and known allosteric sites. 

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
│   ├── cleaned_ddg/              # Per-protein MoCHI-refit folding & binding ddG tables
│   │   └── (kras, pdz3, sh3, src, …)
│   ├── paper_supplements/        # Externally downloaded supplement tables
│   │   ├── domainome/            # Domainome 1.0 (Hoefler et al.)
│   │   ├── megascale/            # Mega-scale stability (Tsuboyama et al. 2023)
│   │   ├── kras_chenchun/        # KRAS DMS (Chenchun et al.)
│   │   └── src_toni/             # SRC DMS
│   ├── proteome_meta/            # ClinVar variant table and UniProt metadata
│   ├── vampseq/                  # VAMP-seq datasets (7 full-length human proteins)
│   ├── scores/                   # ESM-1v, ESM-2, ThermoMPNN, AlphaMissense scores
│   ├── decay_pdb/ & residual_pdb/ # PDB files with B-factors replaced by model values
│   └── sasa/                     # Solvent-accessible surface area annotations
│
├── figs/                         # Output figure panels (PDF/PNG)
├── munge/                        # ProjectTemplate data-preprocessing scripts
├── cache/                        # ProjectTemplate object cache
└── src_hpc/                      # HPC batch scripts for large-scale computations
```

---

## Data access

Large input datasets are hosted on Zenodo and must be downloaded separately before running the analysis scripts:

| Dataset | Description | DOI |
|---|---|---|
| Proteome-wide raw data | Per-protein DMS tables, ESM/AM scores | [10.5281/zenodo.18381534](https://zenodo.org/records/18381534) |
| Proteome-wide meta data | ClinVar annotations, UniProt mappings | [10.5281/zenodo.18386427](https://zenodo.org/records/18386427) |

---

## Requirements

All analyses run in **R ≥ 4.3** using the [ProjectTemplate](http://projecttemplate.net/) framework, which handles library loading and data munging automatically via `load.project()`.

Key R packages:

| Package | Purpose |
|---|---|
| `tidyverse` | Data wrangling and ggplot2 visualisation |
| `data.table` | Fast in-memory operations on large tables |
| `bio3d` | PDB I/O and B-factor mapping |
| `ggExtra` | Marginal density overlays on scatter plots |
| `minpack.lm` | Non-linear least-squares exponential decay fitting (`nlsLM`) |
| `mgcv` | GAM models for spatial decay baselines |
| `lme4` | Linear mixed models for per-residue random effects |
| `broom` | Tidy model summaries |
| `Metrics` | R² and other regression metrics |

Install all packages with:

```r
install.packages(c(
  "ProjectTemplate", "tidyverse", "data.table", "bio3d",
  "ggExtra", "minpack.lm", "mgcv", "lme4", "broom",
  "Metrics", "ggrepel", "pROC", "PRROC"
))
```

---

## Reproducing the analysis

1. **Clone this repository**

   ```bash
   git clone <repo-url>
   cd 01.protein-seq-evo-v1
   ```

2. **Download datasets** from the Zenodo links above and place them under `data/` according to the paths referenced in each script (see `lib/globals.R` for `BASE_DIR` and `SUPPLEMENTS_DIR`).

3. **Open the R project** (`01.protein-seq-evo-v1.Rproj`) in RStudio or set the working directory to the project root, then run any analysis script:

   ```r
   library(ProjectTemplate)
   load.project()          # loads lib/, munge/, and caches data
   rmarkdown::render("src/revision_fig4.Rmd")
   ```

   Each `.Rmd` in `src/` is self-contained and documents its inputs, model steps, and outputs with markdown section headers.

---

## Key methods

| Step | Method | Scripts |
|---|---|---|
| Stability residual computation | LOESS regression of DMS score on ESM-1v/ThermoMPNN ddGf; residual = observed − predicted | All `revision_fig*.Rmd` |
| Spatial decay analysis | Exponential decay fit (NLS + bootstrap CI) of median residual vs. minimum heavy-atom distance to ligand | `revision_fig4.Rmd`, `revision_fig5_*.Rmd` |
| ClinVar pathogenicity enrichment | Wilcoxon test + stability class contingency across benign/pathogenic ClinVar variants | `revision_fig2_*.Rmd`, `revision_fig3*.Rmd` |
| 3D structure mapping | Per-residue residuals written to PDB B-factor column for ChimeraX visualisation | `lib/helpers.R → map_residuals_to_pdb()` |
| Model comparison | ANOVA likelihood-ratio tests and bootstrap-adjusted R² for nested linear/mixed models | `supp_fig4.Rmd` |
| Conservation benchmarks | PhyloP scores and MSA-derived Neff compared against LOESS residuals | `revision_*_phylop.Rmd`, `revision_supp_msa.Rmd` |

---

## Citation

```bibtex
@article{liao2025allostery,
  title   = {Allostery is a widespread cause of loss-of-function variant pathogenicity},
  author  = {Liao, Xiang and Lehner, Ben},
  journal = {bioRxiv},
  year    = {2025},
  doi     = {10.1101/2025.06.20.660737}
}
```

---

## Contact

For questions, bug reports, or data requests please open a GitHub issue or contact **xl7@sanger.ac.uk**.
