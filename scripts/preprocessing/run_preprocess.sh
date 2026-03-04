#!/bin/bash
#SBATCH -p short
#SBATCH -t 0-04:00
#SBATCH -c 4
#SBATCH --mem=96G
#SBATCH -o logs/preprocess_%j.out
#SBATCH -e logs/preprocess_%j.err
#SBATCH --job-name=retina_preprocess

# ---- User configuration ----
# Set REPO_DIR to point to your local clone of this repository.
# Example:
#   export REPO_DIR=/path/to/retina-spatial-scrna-analysis
# Or edit the line below directly.
REPO_DIR="${REPO_DIR:?Set REPO_DIR to your local clone of retina-spatial-scrna-analysis}"

set -uo pipefail

# Adapt the module/conda lines to your HPC environment
module load conda/miniforge3/24.11.3-0
eval "$(conda shell.bash hook)"
conda activate jupyter

cd "${REPO_DIR}/scripts/preprocessing"
mkdir -p ../../logs

Rscript --no-save --no-restore - <<'REOF'
options(warn = 1)
library(tictoc)

cat("=== Preprocessing Pipeline ===\n")
cat("Start:", format(Sys.time()), "\n")
cat("Working dir:", getwd(), "\n")
cat("here::here():", here::here(), "\n\n")

# --- Script 01: Chick preprocessing ---
cat("================================================================\n")
cat("RUNNING: 01_chick_preprocessing.R\n")
tic("01_chick_preprocessing")
source("01_chick_preprocessing.R")
toc()

# --- Script 02: DV/NT scoring ---
cat("\n================================================================\n")
cat("RUNNING: 02_chick_dv_nt_scoring.R\n")
tic("02_chick_dv_nt_scoring")
source("02_chick_dv_nt_scoring.R")
toc()

# --- Script 03: Human preprocessing ---
cat("\n================================================================\n")
cat("RUNNING: 03_human_preprocessing.R\n")
tic("03_human_preprocessing")
source("03_human_preprocessing.R")
toc()

# --- Script 04: Mouse preprocessing ---
cat("\n================================================================\n")
cat("RUNNING: 04_mouse_preprocessing.R\n")
tic("04_mouse_preprocessing")
source("04_mouse_preprocessing.R")
toc()

# --- Script 05: Export h5ad ---
cat("\n================================================================\n")
cat("RUNNING: 05_export_h5ad.R\n")
tic("05_export_h5ad")
source("05_export_h5ad.R")
toc()

# --- Validation ---
cat("\n================================================================\n")
cat("VALIDATION\n")
cat("  retina:", ncol(retina), "cells,", nrow(retina), "genes\n")
cat("  fabp7:", ncol(fabp7), "cells,", nrow(fabp7), "genes\n")
cat("  human:", ncol(human), "cells,", nrow(human), "genes\n")
cat("  mouse:", ncol(mouse), "cells,", nrow(mouse), "genes\n")
cat("  Expected: retina=85135, fabp7=29025, human~21793, mouse~25202\n")
cat("================================================================\n")
cat("End:", format(Sys.time()), "\n")
sessionInfo()
REOF
