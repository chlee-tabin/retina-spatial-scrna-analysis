# Retina Spatial scRNA-seq Analysis

Spatial gene expression analysis for developing retina across chick, human, and mouse.

> "Integration of in situ hybridization and scRNA-seq data provides a 2D topographical map of the developing retina across species"

## Citation

If you use this code, please cite:

> Joisher HNV, Lee C, Prabhakara C, van der Weide I, Si Y, Lonfat N, Cepko C. Integration of in situ hybridization and scRNA-seq data provides a 2D topographical map of the developing retina across species. *bioRxiv* (2026). doi: [10.64898/2026.01.04.697548](https://doi.org/10.64898/2026.01.04.697548)

## Data Availability

- **Interactive viewer**: explore the cross-species 2D topographic gene-expression maps in your browser — [Retina scRNA-seq Pattern Viewer (Hugging Face Spaces)](https://huggingface.co/spaces/chlee-tabin/retina-scrnaseq-viewer)
- **Raw and processed data**: [GEO accession GSE322831](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE322831)
- **Public reference datasets**: GSE138002, GSE234963, GSE246169 (human); GSE118614, GSE139904, GSE149040, GSE122466 (mouse)
- **Mouse re-alignment (Cell Ranger 9.0.1 / GRCm39)**: the mouse arm has been re-aligned with Cell Ranger 9.0.1 to the GRCm39-2024-A reference and integrated across seven libraries spanning the four mouse GEO series above (E13.5–E16). This supersedes an earlier four-library alignment (GSE139904, GSE118614; older reference) that pooled wild-type and Fgfr1/2-mutant cells from a mislabeled source library (GSE139904); the re-alignment restricts that library to control cells only. The prior version is retained in the interactive viewer as a labeled `legacy` dataset for reproducibility only.

## Overview

This repository contains:

1. **`retina_spatial_scrna/`** — Installable Python package for spatial expression analysis: 2D topographic maps, gene similarity, and spatial clustering
2. **`scripts/preprocessing/`** — R pipeline for QC, integration, and spatial axis scoring
3. **`scripts/figures/`** — Standalone scripts reproducing Figures 5-8 and Supplementary Figures 12-23

## Installation

```bash
git clone https://github.com/chlee-tabin/retina-spatial-scrna-analysis.git
cd retina-spatial-scrna-analysis

# Install Python package
pip install -e .

# R dependencies (run in R)
# install.packages(c("Seurat", "tidyverse", "patchwork", "svglite", "harmony",
#                     "glmGamPoi", "viridis", "ggforce", "ggh4x", "tictoc", "glue",
#                     "here"))
# BiocManager::install(c("scDblFinder", "SingleCellExperiment"))
```

### Requirements

- **Python** >= 3.10: scanpy, anndata, numpy, pandas, scipy, matplotlib, mygene
- **R** >= 4.4.0: Seurat v5, tidyverse, patchwork, harmony, scDblFinder, here

See `requirements.txt` for full Python dependencies.

## Repository Structure

```
retina-spatial-scrna-analysis/
├── README.md
├── LICENSE
├── config.sh.example              # Documents the env vars run_preprocess.sh / 05_export reads
├── requirements.txt
├── setup.py
│
├── retina_spatial_scrna/           # Installable Python package
│   ├── __init__.py
│   ├── spatial_expression_analysis.py
│   └── control_genes.yaml         # Species-specific spatial anchor genes
│
├── spatial_expression_analysis.py  # Convenience module (same as package)
├── control_genes.yaml
│
├── data/
│   ├── chick_W_genes.tsv          # W chromosome gene list
│   └── chick_Z_genes.tsv          # Z chromosome gene list
│
├── scripts/
│   ├── README.md                  # Detailed execution guide
│   ├── preprocessing/
│   │   ├── 00_utils.R             # Shared R utilities
│   │   ├── 01_chick_preprocessing.R
│   │   ├── 02_chick_dv_nt_scoring.R
│   │   ├── 03_human_preprocessing.R
│   │   ├── 04_mouse_preprocessing.R
│   │   ├── 05_export_h5ad.R
│   │   └── run_preprocess.sh      # SLURM batch wrapper
│   └── figures/
│       ├── fig5_dv_nt_scores.R
│       ├── fig6ag_chick_topographic.py
│       ├── fig6h_sfig23_area_deg.R
│       ├── fig7_cross_species.py
│       ├── fig8_human_haa.py
│       ├── sfig12-14 (DV/NT scores).R
│       ├── sfig15_grid_sensitivity.py
│       ├── sfig16-17 (chick clusters/signaling).py
│       ├── sfig18_mouse_human_scores.R
│       ├── sfig19_mouse_human_clusters.py
│       └── sfig20_22_pathway_maps.py
│
└── notebooks/figures/             # (figure legends + methods text added at publication)
    ├── figure_legends.md
    └── methods_section.md
```

## Quick Start

### Reproducing Figures

The deposited [GEO GSE322831](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE322831) objects are the supported entry point — no preprocessing run is required.

1. **Download the processed objects** from GEO and strip the `GSE322831_` prefix into `data/`:

   ```bash
   # from the repository root, with the GSE322831_* downloads in the current directory:
   mkdir -p data
   for f in GSE322831_*.h5ad GSE322831_*.rds; do [ -e "$f" ] && mv "$f" "data/${f#GSE322831_}"; done
   ```

   | GEO supplementary file | becomes | used by |
   |---|---|---|
   | `GSE322831_20250604_chick_RPC.h5ad` | `data/20250604_chick_RPC.h5ad` | fig6ag, fig7, sfig15–17 |
   | `GSE322831_20250604_human_RPC.h5ad` | `data/20250604_human_RPC.h5ad` | fig7, fig8, sfig19–22 |
   | `GSE322831_20260528_mouse_RPC_cr9_e13e16.h5ad` | `data/20260528_mouse_RPC_cr9_e13e16.h5ad` | fig7, sfig19 |
   | `GSE322831_20250604_02_fabp7.rds` | `data/20250604_02_fabp7.rds` | fig5, sfig13, sfig14, fig6h |
   | `GSE322831_20250604_01_retina.rds` | `data/20250604_01_retina.rds` | sfig12 |

   The remaining `GSE322831_*` supplementary files are not needed for the figures: `..._mouse_RPC_cr9.h5ad` is the full seven-library mouse object behind the interactive viewer (the figures use the E13.5–E16 subset), and the GTF/genome files support re-alignment from raw reads. GEO's bundled `readme.txt` predates the mouse re-alignment; this table is current.

2. **Generate figures** — each script is independent and runnable from anywhere inside the repository:
   ```bash
   python scripts/figures/fig6ag_chick_topographic.py    # Figure 6A-G
   Rscript scripts/figures/fig5_dv_nt_scores.R           # Figure 5B-E
   ```

   12 of the 14 figure scripts run from the GEO objects alone. The two exceptions consume *in-session* Seurat objects — objects that exist only as variables inside a running R session, never written to disk, so `Rscript` (a fresh session per script) cannot chain them:

   - **SF23 half of `fig6h_sfig23`** (the Figure 6H half runs from the GEO objects and exits cleanly): needs `human`. In one R session: `source("scripts/preprocessing/03_human_preprocessing.R", chdir = TRUE)`, then source the figure script. Note 03 itself reads an archived pre-publication intermediate (internal; available on request) — see *Reproducibility scope*.
   - **`sfig18`** (all panels): needs `human` (as above) **and** the CR9 re-aligned `mouse` object, which is produced outside this repo — not reproducible from public data alone.

See [`scripts/README.md`](scripts/README.md) for the complete figure-to-script mapping and data flow.

### Using the Python Package

```python
from retina_spatial_scrna import SpatialExpressionAnalyzer, SpatialAnalysisParams

# Heuristic starting-point parameters (the published panels set parameters
# explicitly in scripts/figures/ — see get_species_defaults docstring)
params = SpatialAnalysisParams.get_species_defaults('chick')

# Run analysis
analyzer = SpatialExpressionAnalyzer(params)
results = analyzer.run_full_analysis('data/20250604_chick_RPC.h5ad')

# Visualize genes
analyzer.show_gene_image('FGF8')

# Get gene correlations
correlations = analyzer.get_gene_correlations('FGF8', top_n=10)
```

## Figure-to-Script Mapping

| Figure | Script | Language |
|--------|--------|----------|
| F5B-E | `fig5_dv_nt_scores.R` | R |
| F6A-G | `fig6ag_chick_topographic.py` | Python |
| F6H + SF23 | `fig6h_sfig23_area_deg.R` | R |
| F7A-G | `fig7_cross_species.py` | Python |
| F8B-C | `fig8_human_haa.py` | Python |
| SF12A-E | `sfig12_scrnaseq_qc.R` | R |
| SF13A-E | `sfig13_dv_score.R` | R |
| SF14A-E | `sfig14_nt_score.R` | R |
| SF15 | `sfig15_grid_sensitivity.py` | Python |
| SF16 + Table S1 | `sfig16_chick_spatial_clusters.py` | Python |
| SF17A-D | `sfig17_chick_signaling.py` | Python |
| SF18A-D | `sfig18_mouse_human_scores.R` | R |
| SF19 + Tables S2-S3 | `sfig19_mouse_human_clusters.py` | Python |
| SF20-22 | `sfig20_22_pathway_maps.py` | Python |

### Figures Not in Scope

Figures generated by RNA-FISH imaging (MATLAB pipeline): F1-F3, SF1-SF11, F4, F5A, F8A.

## Reproducibility scope

**Supported entry point: the GEO GSE322831 processed objects** (table above). 12 of the 14 figure scripts run directly from them; the two exceptions are listed under *Reproducing Figures*.

**`scripts/preprocessing/` is the provenance record** of how those objects were made. It is published for transparency and is not fully re-runnable from public data alone:

- `01_chick_preprocessing.R` reads the raw Cell Ranger/cellsnp/vireo tree (raw reads are available from GEO; the aligned tree is not distributed), plus a 2024-08 merged chick object used for barcode cross-checks.
- `03_human_preprocessing.R` and `04_mouse_preprocessing.R` read pre-publication intermediates (the merged human and integrated mouse Seurat objects). Their own derivation — from the public GSE138002/GSE234963/GSE246169 (human) and GSE118614/GSE139904 (mouse) data — is retained in the lab's archived analysis notebook (internal; available on request).
- The deposited mouse h5ads come from a later Cell Ranger 9.0.1 / GRCm39 re-alignment (see *Data Availability*), not from `04` — `04` documents the superseded pre-CR9 mouse arm.

The chain was re-executed end-to-end from the archived inputs (2026-03):

- **Human**: re-run 23,031 cells = original 23,031 at the R-export stage — exact. (The deposited `20250604_human_RPC.h5ad` contains 21,793 cells after additional filtering during Python-side assembly.)
- **Mouse**: re-run 25,202 = original 25,202 — exact, for the superseded pre-CR9 mouse arm that `04` documents. The deposited CR9 `..._e13e16` object (26,505 cells) is validated separately by the re-alignment pipeline.
- **Chick**: re-run 83,915 vs original 85,135 cells at the integrated-retina stage (~1.4%), cascading to 29,987 vs 29,025 in the RPC subset — the drift traces to unseeded `scDblFinder` doublet calls propagating through Harmony/Leiden.

Even after a re-run, make figures from the deposited GEO objects — they are the objects the published panels were built from.

## Data Format

Input `.h5ad` files should contain:
- **Expression matrix**: `adata.X` (cells x genes)
- **Spatial coordinates**: `adata.obs['DV.Score']`, `adata.obs['NT.Score']`
- **Gene names**: `adata.var_names`

## Species and Gene Naming

| Species | Gene case | Example |
|---------|-----------|---------|
| Human (*Homo sapiens*) | UPPERCASE | FGF8, TBX5 |
| Chick (*Gallus gallus*) | UPPERCASE | FGF8, TBX5 |
| Mouse (*Mus musculus*) | Sentence case | Fgf8, Tbx5 |

## Contributing

Contributions welcome. In particular, RNA-FISH reproduction scripts (MATLAB) are planned for inclusion under `scripts/figures/` or `scripts/fish/`.

## License

This project is released under the [MIT License](LICENSE).
