# Reproducible Analysis Scripts

Standalone scripts in [jupytext percent format](https://jupytext.readthedocs.io/) for reproducing the scRNA-seq computational analysis from:

> "Integration of in situ hybridization and scRNA-seq data provides a 2D topographical map of the developing retina across species"

These scripts generate Figures 5–8 and Supplementary Figures 12–23 of the manuscript.

## Prerequisites

### R (preprocessing + R figure scripts)
- R >= 4.4.0
- Seurat v5, tidyverse, patchwork, svglite, viridis, ggforce, ggh4x
- harmony, scDblFinder, presto, glmGamPoi, SingleCellExperiment
- glue, tictoc

### Python (Python figure scripts)
- Python >= 3.10
- scanpy, anndata, numpy, pandas, scipy, matplotlib
- mygene (for cross-species ortholog lookup)
- Project modules: `spatial_expression_analysis.py`, `control_genes.yaml` (in repository root)

## Data Requirements

Download from GEO (accession TBD) or provide your own data in the expected layout:

```
data/
  # Chick (from GEO submission)
  20250604_chick_RPC.h5ad         # Chick RPC with DV/NT scores
  20250604chick/                  # Full chick MEX export
  20250604chick.RPC/              # Chick RPC MEX export

  # Human (public datasets: GSE138002, GSE234963, GSE246169)
  20250604_human_RPC.h5ad         # Human RPC with DV/NT scores

  # Mouse (public datasets: GSE139904, GSE118614)
  20250604_mouse_RPC.h5ad         # Mouse RPC with DV/NT scores

  # Intermediate R objects (produced by preprocessing scripts)
  20250604_01_gex.rds             # Checkpoint 1: raw Seurat list
  20250604_01_retina.rds          # Checkpoint 2: integrated Seurat
  20250604_02_fabp7.rds           # Checkpoint 3: scored RPC subset
  chick_W_genes.tsv               # W chromosome gene list
  chick_Z_genes.tsv               # Z chromosome gene list
```

Set `RETINA_DATA_DIR` environment variable to override the default `../data` path:
```bash
export RETINA_DATA_DIR=/path/to/your/data
```

## Execution Order

### 1. Preprocessing (R)

Run in order — each script depends on outputs from the previous:

```bash
cd scripts/preprocessing/
Rscript 01_chick_preprocessing.R      # ~30 min, produces checkpoints 1 & 2
Rscript 02_chick_dv_nt_scoring.R      # ~5 min, produces checkpoint 3 + h5ad
Rscript 03_human_preprocessing.R      # ~5 min
Rscript 04_mouse_preprocessing.R      # ~5 min
Rscript 05_export_h5ad.R              # Optional: re-export h5ad files
```

### 2. Figure scripts (R and Python)

Run in any order after preprocessing. Each script is independent:

| Script | Figure | Language |
|--------|--------|----------|
| `figures/fig5_dv_nt_scores.R` | F5B-E | R |
| `figures/sfig12_scrnaseq_qc.R` | SF12A-E | R |
| `figures/sfig13_dv_score.R` | SF13A-E | R |
| `figures/sfig14_nt_score.R` | SF14A-E | R |
| `figures/fig6h_sfig23_area_deg.R` | F6H + SF23 | R |
| `figures/sfig18_mouse_human_scores.R` | SF18A-D | R |
| `figures/fig6ag_chick_topographic.py` | F6A-G | Python |
| `figures/fig7_cross_species.py` | F7A-G | Python |
| `figures/fig8_human_haa.py` | F8B-C | Python |
| `figures/sfig15_grid_sensitivity.py` | SF15 | Python |
| `figures/sfig16_chick_spatial_clusters.py` | SF16 + Table S1 | Python |
| `figures/sfig17_chick_signaling.py` | SF17A-D | Python |
| `figures/sfig19_mouse_human_clusters.py` | SF19A-B + Tables S2/S3 | Python |
| `figures/sfig20_22_pathway_maps.py` | SF20-22 | Python |

### Running as Jupyter notebooks

These scripts can be opened as notebooks using jupytext:

```bash
# Open in Jupyter
jupytext --to notebook scripts/figures/fig6ag_chick_topographic.py
jupyter notebook scripts/figures/fig6ag_chick_topographic.ipynb

# Or pair for live sync
jupytext --set-formats py:percent,ipynb scripts/figures/fig6ag_chick_topographic.py
```

## Data Flow

```
Cell Ranger outputs (10X)
  │
  ├── 01_chick_preprocessing.R
  │   └── 20250604_01_retina.rds (integrated Seurat)
  │       └── h5ad export: 20250604chick/
  │
  ├── 02_chick_dv_nt_scoring.R
  │   └── 20250604_02_fabp7.rds (RPC + DV/NT scores)
  │       └── h5ad export: 20250604_chick_RPC.h5ad
  │
  ├── 03_human_preprocessing.R → human RPC h5ad
  └── 04_mouse_preprocessing.R → mouse RPC h5ad
       │
       ├── Python figure scripts (consume h5ad files)
       │   ├── fig6ag: chick topos
       │   ├── fig7: cross-species
       │   ├── fig8: human HAA
       │   ├── sfig15-17: chick analysis
       │   ├── sfig19: mouse/human clusters
       │   └── sfig20-22: pathway maps
       │
       └── R figure scripts (consume .rds checkpoints)
           ├── fig5: DV/NT scores
           ├── sfig12: QC panels
           ├── sfig13-14: score details
           ├── fig6h+sfig23: area DEG
           └── sfig18: mouse/human scores
```

## Figures NOT in scope

- F1–F3, SF1–SF11: RNA-FISH images (MATLAB pipeline)
- F4: RNA-FISH quantification
- F5A: Schematic workflow (BioRender)
- F8A: RNA-FISH cross-section

## Notes

- Scripts use `source("00_utils.R")` for shared R functions; ensure working directory is `scripts/preprocessing/` or adjust the source path.
- The `plot_axial_expression()` function (in 00_utils.R) is used by fig5, sfig13, sfig14, and sfig18 scripts.
- Python scripts import from `spatial_expression_analysis.py` and `control_genes.yaml` in the repository root. Ensure these are on your Python path.
