# Mouse Cell Ranger 9 re-alignment

Provenance for the mouse retinal atlas aligned with **Cell Ranger 9.0.1** to
**refdata-gex-GRCm39-2024-A**, followed by QC, Harmony integration, annotation,
RPC extraction and DV/NT scoring. Figures use the seven-library E13.5–E16 RPC
subset (26,505 cells); the full E13.5–P0 RPC object has ten libraries (43,991 cells).
Download the deposited objects from GEO GSE322831 to reproduce figures without
running this resource-intensive pipeline.

## Inputs and configuration

[manifest.tsv](manifest.tsv) records the accessions, stages, chemistry and download
routes. The atlas includes Clark GSE118614 (SRR7699772–SRR7699777),
Balasubramanian GSE139904 (SRR10398966), Wu GSE149040 (SRR11582194), and
Lo Giudice GSE122466 (SRR8181428–SRR8181429). The two Georges E-MTAB-9395
candidates are retained in the manifest for provenance but excluded from the
atlas because their BAM headers are incompatible with `cellranger bamtofastq`.

Install Cell Ranger 9.0.1, SRA Toolkit 3.2.0, wget, pigz/gzip and samtools on
`PATH`, and use Bash >= 4 and SLURM for the submission wrappers. Activate an
environment with R (built with 4.4.3), Seurat v5, SingleCellExperiment,
glmGamPoi, Matrix, tidyverse, harmony, presto, scDblFinder, scater, reticulate,
and Python with leidenalg, anndata, numpy, pandas and scipy before submitting.
The batch files inherit this environment; adjust resource requests to your system.

Set these variables from the repository root (paths may point to external storage):

```bash
export RETINA_REPO="$PWD"
export WORKDIR="$PWD/data/mouse_cr9"
export CR_REF="$PWD/data/refdata-gex-GRCm39-2024-A"
export RETICULATE_PYTHON="$(command -v python3)"
mkdir -p "$WORKDIR/logs"
```

Obtain the following author metadata before the atlas build. Each path can be
overridden by the named environment variable; defaults are under `WORKDIR`:

| Variable | Default path relative to `WORKDIR` | Input |
|---|---|---|
| `SCP_CTRL` | `SCP1618/SCP1618_control_bcs.txt` | SCP1618 Control (WT) cells, one original 10x barcode per line including `-1`; 4,642 barcodes in the published run |
| `CLARK_META` | `GSE118614_meta/GSE118614_barcodes.tsv.gz` | Clark GSE118614 barcode metadata with `umap2_CellType` annotations |
| `BALA_GRAPH` | `GSE139904/analysis/clustering/graphclust/clusters.csv` | Balasubramanian author graph-cluster CSV with `Barcode` and `Cluster` columns |

The Control barcode list restricts the pooled Balasubramanian library to WT cells;
mutant cells are excluded before QC. Download temporary files use `TMPDIR`,
defaulting to `WORKDIR/tmp/<sample>` for the SRA route.

## Order of steps

1. Run `bash scripts/realign_mouse/submit_downloads.sh all` and wait for downloads
   to finish. BAM-submitted runs use `download_bam_single.sbatch` and
   `cellranger bamtofastq`; FASTQ-submitted runs use `download_fastq.sbatch`.
   FASTQs are written to `WORKDIR/<series>/<sample>/fastq/`.
2. Run `bash scripts/realign_mouse/submit_counts.sh` and wait for all ten counts.
   `cellranger_count.sbatch` pins v2/v3 chemistry from the source libraries and
   writes `WORKDIR/counts/<sample>/outs/`, including BAMs and filtered matrices.
3. Submit the atlas build and wait for completion:

   ```bash
   sbatch --output="$WORKDIR/logs/build_%j.out" --error="$WORKDIR/logs/build_%j.err" \
     scripts/realign_mouse/mouse_cr9_01_build.sbatch
   ```

   `mouse_cr9_01_build.R` preserves the seed (1234), QC thresholds, clustering,
   scoring genes and curated annotation checks. It writes
   `WORKDIR/atlas/mouse_cr9_integrated.rds`, diagnostic TSVs and full-atlas MEX files.
4. Submit the RPC export and h5ad assembly:

   ```bash
   sbatch --output="$WORKDIR/logs/export_%j.out" --error="$WORKDIR/logs/export_%j.err" \
     scripts/realign_mouse/mouse_cr9_02_export.sbatch
   ```

   `mouse_cr9_02_export_rpc.R` writes `atlas/export_rpc/` and
   `atlas/export_rpc_e13e16/` with the inherited DV/NT scores. The wrapper then
   calls `mouse_cr9_03_assemble_h5ad.py` for both RPC exports and the full atlas.

Final outputs under `WORKDIR/atlas/h5ad/` are `20260528_mouse_RPC_cr9.h5ad`,
`20260528_mouse_RPC_cr9_e13e16.h5ad` and `20260528_mouse_atlas_cr9_full.h5ad`.
Their `.X` holds log-normalized expression, `.raw.X` raw counts, `.obs` the
annotations and DV/NT scores, and `.obsm` the exported embeddings. The assembler
can also run directly with `--dir`, `--prefix` and `--out`.

[build_mouse_pathway_pickle.py](build_mouse_pathway_pickle.py) is the lightweight
SF20A–SF22A analyzer helper. It reads the deposited E13.5–E16 h5ad directly and
does not require re-alignment; the mouse pathway figure script calls it automatically.
