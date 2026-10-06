# Retina Spatial scRNA-seq Analysis

Spatial gene expression analysis for developing retina across chick, human, and mouse.

> "Integration of in situ hybridization and scRNA-seq data provides a 2D topographical map of the developing retina across species"

## Citation

If you use this code, please cite:

> Joisher HNV, Lee C, Prabhakara C, van der Weide I, Si Y, Lonfat N, Cepko C. Integration of in situ hybridization and scRNA-seq data provides a 2D topographical map of the developing retina across species. *bioRxiv* (2026). doi: [10.64898/2026.01.04.697548](https://doi.org/10.64898/2026.01.04.697548)

## Data Availability

- **Interactive viewer**: explore the cross-species 2D topographic gene-expression maps in your browser — [Retina scRNA-seq Pattern Viewer (Hugging Face Spaces)](https://huggingface.co/spaces/chlee-tabin/retina-scrnaseq-viewer) ([source code](https://github.com/chlee-tabin/retina-scrnaseq-viewer), MIT)
- **Raw and processed data**: [GEO accession GSE322831](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE322831)
- **Public reference datasets**: GSE138002, GSE234963, GSE246169 (human); GSE118614, GSE139904, GSE149040, GSE122466 (mouse)
- **Mouse re-alignment (Cell Ranger 9.0.1 / GRCm39)**: the mouse arm has been re-aligned with Cell Ranger 9.0.1 to the GRCm39-2024-A reference and integrated across ten libraries spanning the four mouse GEO series above (E13.5–P0; 43,991 RPCs). The E13.5–E16 subset used in figures comprises seven libraries (26,505 cells). This supersedes an earlier four-library alignment (GSE139904, GSE118614; older reference) that pooled wild-type and Fgfr1/2-mutant cells from a mislabeled source library (GSE139904); the re-alignment restricts that library to control cells only. The prior version is retained in the interactive viewer as a labeled `legacy` dataset for reproducibility only. Re-alignment producers and inputs are documented in [`scripts/realign_mouse/README.md`](scripts/realign_mouse/README.md).

## Overview

This repository contains:

1. **`retina_spatial_scrna/`** — Installable Python package for spatial expression analysis: 2D topographic maps, gene similarity, and spatial clustering
2. **`scripts/preprocessing/`** — R pipeline for QC, integration, and spatial axis scoring
3. **`scripts/figures/`** — Standalone scripts reproducing the scRNA-seq panels of Figures 5-8, Supplementary Figures 12-24 and Supplementary Tables 1-2
4. **`scripts/realign_mouse/`** — Mouse CR9/GRCm39 re-alignment, atlas build, RPC export and h5ad assembly provenance

## Installation

```bash
git clone https://github.com/chlee-tabin/retina-spatial-scrna-analysis.git
cd retina-spatial-scrna-analysis

# Install Python package
pip install -e .

# R dependencies (run in R)
# install.packages(c("Seurat", "tidyverse", "patchwork", "svglite", "harmony",
#                     "glmGamPoi", "viridis", "ggforce", "ggh4x", "tictoc", "glue",
#                     "here", "reticulate"))   # reticulate: R scripts that read the human h5ad
# BiocManager::install(c("scDblFinder", "SingleCellExperiment", "glmGamPoi"))
# remotes::install_github("immunogenomics/presto")   # used by 01_chick_preprocessing.R
```

### Requirements

- **Python** >= 3.10: scanpy, anndata, numpy, pandas, scipy, matplotlib, mygene
- **R** >= 4.4.0: Seurat v5, tidyverse, patchwork, harmony, scDblFinder, glmGamPoi, here, reticulate (reads the human h5ad via Python anndata)

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
│   ├── chick_Z_genes.tsv          # Z chromosome gene list
│   └── fig8d_pathway_modules/     # per-bin module maps rendered as Figure 8D
│
├── tests/                         # Smoke tests for the Python package
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
│   ├── realign_mouse/             # CR9/GRCm39 provenance; see its README for env vars and steps
│   └── figures/
│       ├── fig5_dv_nt_scores.R
│       ├── fig6ag_chick_topographic.py
│       ├── fig7_cross_species.py
│       ├── fig8d_pathway_module_data.R / fig8d_pathway_module_topography.R
│       ├── sfig12-14 (DV/NT scores).R
│       ├── sfig15_grid_sensitivity.py
│       ├── sfig15c_flower_reprojection.py / sfig15c_extract_dvnt.py
│       ├── sfig16-17 (chick clusters/signaling).py
│       ├── sfig18_mouse_human_scores.R
│       ├── sfig19_mouse_human_clusters.py
│       ├── sfig20_22_pathway_maps.py
│       ├── sfig20_22_pathway_maps_mouse.py
│       ├── sfig23_human_cyp26_correlations.py
│       ├── sfig24_area_deg.R
│       └── supptables1_2_area_deg.R
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
   | `GSE322831_20250604_human_RPC.h5ad` | `data/20250604_human_RPC.h5ad` | fig7, sfig19–23 |
   | `GSE322831_20260528_mouse_RPC_cr9_e13e16.h5ad` | `data/20260528_mouse_RPC_cr9_e13e16.h5ad` | fig7, sfig19, sfig20–22 (mouse) |
   | `GSE322831_20250604_02_fabp7.rds` | `data/20250604_02_fabp7.rds` | fig5, sfig13, sfig14, fig8d (data), sfig24, supptables1_2 |
   | `GSE322831_20250604_01_retina.rds` | `data/20250604_01_retina.rds` | sfig12 |

   **Human R-export.** The human area-DEG (Supplementary Table 2, Fig. S24B) and the Fig. 8D human module maps were computed on the R-export-stage human RPC object (23,031 cells), before the blood-contamination filter that produced the GEO h5ad (21,793 cells; cells with ≥1% of reads from hemoglobin genes removed, almost all from GSE234963 whole-eye libraries). On the GEO object the effect sizes are the same (log2FC r = 0.97; 0.999 within the fovea gate), but the fovea gate loses its 67 blood-contaminated cells and two whole-eye libraries fall below the 50-cell pseudobulk minimum, so fewer fovea genes reach significance and Table 2 does not reproduce exactly. The R-export (MEX: raw counts, barcodes, features, metadata) is archived at [Zenodo DOI to be added]; unpack it to `data/20250604human.RPC/`. It is read by `read_human_rexport()` in `scripts/preprocessing/00_utils.R` and used by fig8d (data), sfig24 (human) and supptables1_2 (Table 2).

   The remaining `GSE322831_*` supplementary files are not needed for the figures: `..._mouse_RPC_cr9.h5ad` is the full ten-library E13.5–P0 mouse object (43,991 RPCs) behind the interactive viewer (the figures use the seven-library E13.5–E16 subset, 26,505 cells), and the GTF/genome files support re-alignment from raw reads. GEO's bundled `readme.txt` predates the mouse re-alignment; this table is current.

2. **Generate figures** — each script is independent and runnable from anywhere inside the repository:
   ```bash
   python scripts/figures/fig6ag_chick_topographic.py    # Figure 6A-G
   Rscript scripts/figures/fig5_dv_nt_scores.R           # Figure 5B-E
   python scripts/figures/sfig15c_flower_reprojection.py # SF15C; locked parameters are the defaults
   python scripts/figures/sfig20_22_pathway_maps_mouse.py # SF20A-SF22A; builds the mouse cache automatically
   ```

   Every figure/table script except `sfig18` runs from the GEO objects alone (`fig8d_pathway_module_topography.R` needs only the committed CSVs in `data/fig8d_pathway_modules/`). R scripts that read the human object (`fig8d_pathway_module_data.R`, `sfig24_area_deg.R`, `supptables1_2_area_deg.R`) load the h5ad through `read_h5ad_as_seurat()` in `00_utils.R`, which needs the R package `reticulate` and the Python package `anndata`.

   The exception, **`sfig18`** (all panels), consumes *in-session* Seurat objects, so `Rscript` (a fresh session per script) cannot chain them: `human` from `source("scripts/preprocessing/03_human_preprocessing.R", chdir = TRUE)` (03 itself reads an archived pre-publication intermediate; internal, available on request — see *Reproducibility scope*) **and** the CR9 re-aligned `mouse` object. The mouse provenance is now in `scripts/realign_mouse/`; its deposited E13.5–E16 h5ad can be loaded in the same R session with `mouse <- NormalizeData(read_h5ad_as_seurat(here::here("data", "20260528_mouse_RPC_cr9_e13e16.h5ad")))` after sourcing `00_utils.R`. The archived human input still prevents reproducing all SF18 panels from public data alone.

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

Numbering follows the revised manuscript (v15).

| Figure / table | Script | Language |
|--------|--------|----------|
| F5B-E | `fig5_dv_nt_scores.R` | R |
| F6A-G | `fig6ag_chick_topographic.py` | Python |
| F7A-G | `fig7_cross_species.py` | Python |
| F8A (2D maps) | `fig6ag_chick_topographic.py`: DUSP6 / SPRY1 (`Figure_SF17/spatial_fgf8_downstream/`) and MYOF / NPY (`Figure6/F6A-G_spatial_maps/`) | Python |
| F8B (2D map) | the human BAMBI per-gene map is emitted by `sfig20_22_pathway_maps.py` (`Figure_SF22/individual_genes/`) | Python |
| F8C (2D maps) | the human CYP26A1 / CYP26C1 maps are the target-gene panels of `sfig23_human_cyp26_correlations.py` | Python |
| F8D | `fig8d_pathway_module_data.R` (per-bin module maps from the GEO objects) → `fig8d_pathway_module_topography.R` (render) | R |
| SF12A-E | `sfig12_scrnaseq_qc.R` | R |
| SF13A-D | `sfig13_dv_score.R` | R |
| SF14A-D | `sfig14_nt_score.R` | R |
| SF15A-B | `sfig15_grid_sensitivity.py` | Python |
| SF15C | `sfig15c_flower_reprojection.py` (+ `sfig15c_extract_dvnt.py`; locked defaults, GEO chick object) | Python |
| SF16A | `sfig16_chick_spatial_clusters.py` (+ cluster-member table) | Python |
| SF16B-D | `sfig17_chick_signaling.py` (FGF8 / CYP26C1 / BMP2 correlations) | Python |
| SF17A-C | `fig6ag_chick_topographic.py` (FGF family / FGF8 downstream / BMP maps) | Python |
| SF17D | `sfig17_chick_signaling.py` (CYP1B1 correlations) | Python |
| SF18A, B; A′, B′; C | `sfig18_mouse_human_scores.R` (human DV / NT; human validation; mouse DV / NT) | R |
| SF19A-B | `sfig19_mouse_human_clusters.py` (+ cluster-member tables) | Python |
| SF20-22 B (human) | `sfig20_22_pathway_maps.py` | Python |
| SF20-22 A (mouse) | `sfig20_22_pathway_maps_mouse.py` (+ `scripts/realign_mouse/build_mouse_pathway_pickle.py`; CR9 E13.5–E16 object) | Python |
| SF23A-B | `sfig23_human_cyp26_correlations.py` | Python |
| SF24A-B | `sfig24_area_deg.R` | R |
| Supplementary Tables 1-2 | `supptables1_2_area_deg.R` | R |
| Mouse CR9 provenance | `scripts/realign_mouse/` ([steps and inputs](scripts/realign_mouse/README.md)) | Bash / R / Python |

Supplementary Table 3 (FGF/BMP pathway gene survey in human and mouse) and SF13E/SF14E (scRNA-seq correlation heatmaps) have no producer in this repository.

### Figures Not in Scope

RNA-FISH imaging and quantification (MATLAB pipeline): F1-F4, SF1-SF11, and the RNA-FISH halves of F8A-C. Schematics: F5A (BioRender), F6H.

## Reproducibility scope

**Supported entry point: the GEO GSE322831 processed objects** (table above). All figure/table scripts except `sfig18` run directly from them; the exception is described under *Reproducing Figures*.

**`scripts/preprocessing/` is the provenance record** of how those objects were made. It is published for transparency and is not fully re-runnable from public data alone:

- `01_chick_preprocessing.R` reads the raw Cell Ranger/cellsnp/vireo tree (raw reads are available from GEO; the aligned tree is not distributed), plus a 2024-08 merged chick object used for barcode cross-checks.
- `03_human_preprocessing.R` and `04_mouse_preprocessing.R` read pre-publication intermediates (the merged human and integrated mouse Seurat objects). Their own derivation — from the public GSE138002/GSE234963/GSE246169 (human) and GSE118614/GSE139904 (mouse) data — is retained in the lab's archived analysis notebook (internal; available on request).
- The deposited mouse h5ads come from the later Cell Ranger 9.0.1 / GRCm39 re-alignment documented in [`scripts/realign_mouse/`](scripts/realign_mouse/README.md), not from `04` — `04` documents the superseded pre-CR9 mouse arm.

The chain was re-executed end-to-end from the archived inputs (2026-03):

- **Human**: re-run 23,031 cells = original 23,031 at the R-export stage — exact. (The deposited `20250604_human_RPC.h5ad` contains 21,793 cells after the blood-contamination filter, `percent.rbc < 0.01`, applied during Python-side assembly.)
- **Mouse**: re-run 25,202 = original 25,202 — exact, for the superseded pre-CR9 mouse arm that `04` documents. The deposited CR9 `..._e13e16` object (26,505 cells) is validated separately by the re-alignment pipeline.
- **Chick**: re-run 83,915 vs original 85,135 cells at the integrated-retina stage (~1.4%), cascading to 29,987 vs 29,025 in the RPC subset — the drift traces to unseeded `scDblFinder` doublet calls propagating through Harmony/Leiden.

Even after a re-run, make figures from the deposited GEO objects — they are the objects the published panels were built from.

## Data Format

Input `.h5ad` files should contain:
- **Expression matrix**: `adata.X` (cells x genes)
- **Spatial coordinates**: `adata.obs['DV.Score']`, `adata.obs['NT.Score']`
- **Gene names**: `adata.var_names`

## Species and Gene Naming

Data objects use each species' native gene symbols (UPPERCASE chick/human, Sentence case mouse); the manuscript displays them in italicised Sentence case.

| Species | Gene case | Example |
|---------|-----------|---------|
| Human (*Homo sapiens*) | UPPERCASE | FGF8, TBX5 |
| Chick (*Gallus gallus*) | UPPERCASE | FGF8, TBX5 |
| Mouse (*Mus musculus*) | Sentence case | Fgf8, Tbx5 |

## Contributing

Contributions welcome. In particular, RNA-FISH reproduction scripts (MATLAB) are planned for inclusion under `scripts/figures/` or `scripts/fish/`.

## License

This project is released under the [MIT License](LICENSE).
