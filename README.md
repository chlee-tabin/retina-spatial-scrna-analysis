# Retina Spatial scRNA-seq Analysis

Spatial gene expression analysis for developing retina across chick, human, and mouse.

> "Integration of in situ hybridization and scRNA-seq data provides a 2D topographical map of the developing retina across species"

## Citation

If you use this code, please cite:

> Joisher HNV, Lee C, Prabhakara C, van der Weide I, Si Y, Lonfat N, Cepko C. Integration of in situ hybridization and scRNA-seq data provides a 2D topographical map of the developing retina across species. *bioRxiv* (2026). doi: [10.64898/2026.01.04.697548](https://doi.org/10.64898/2026.01.04.697548)

## Data Availability

- **Interactive viewer**: explore the cross-species 2D topographic gene-expression maps in your browser — [Retina scRNA-seq Pattern Viewer (Hugging Face Spaces)](https://huggingface.co/spaces/chlee-tabin/retina-scrnaseq-viewer)
- **Raw and processed data**: [GEO accession GSE322831](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE322831)
- **Public reference datasets**: GSE138002, GSE234963, GSE246169 (human); GSE139904, GSE118614 (mouse)

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
#                     "glmGamPoi", "viridis", "ggforce", "ggh4x", "tictoc", "glue"))
# BiocManager::install(c("scDblFinder", "SingleCellExperiment"))
```

### Requirements

- **Python** >= 3.10: scanpy, anndata, numpy, pandas, scipy, matplotlib, mygene
- **R** >= 4.4.0: Seurat v5, tidyverse, patchwork, harmony, scDblFinder

See `requirements.txt` for full Python dependencies.

## Repository Structure

```
retina-spatial-scrna-analysis/
├── README.md
├── LICENSE
├── config.sh.example              # Template for user paths (copy to config.sh)
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
└── notebooks/figures/
    ├── figure_legends.md
    └── methods_section.md
```

## Quick Start

### Reproducing Figures

1. **Download data** from GEO (accession above)
2. **Configure paths**: `cp config.sh.example config.sh` and edit
3. **Run preprocessing** (R, ~45 min total):
   ```bash
   cd scripts/preprocessing/
   Rscript 01_chick_preprocessing.R
   Rscript 02_chick_dv_nt_scoring.R
   Rscript 03_human_preprocessing.R
   Rscript 04_mouse_preprocessing.R
   ```
4. **Generate figures** (each script is independent):
   ```bash
   cd scripts/figures/
   python fig6ag_chick_topographic.py    # Figure 6A-G
   Rscript fig5_dv_nt_scores.R           # Figure 5B-E
   ```

See [`scripts/README.md`](scripts/README.md) for the complete figure-to-script mapping and data flow.

### Using the Python Package

```python
from retina_spatial_scrna import SpatialExpressionAnalyzer, SpatialAnalysisParams

# Get species-optimized parameters
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
