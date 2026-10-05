#!/bin/bash
#SBATCH -J mcr9_build
#SBATCH -t 0-8:00
#SBATCH -c 16
#SBATCH --mem 128G
#SBATCH -o build_%j.out
#SBATCH -e build_%j.err
set -euo pipefail
WORKDIR=${WORKDIR:?Set WORKDIR to the mouse re-alignment working directory}
REPO=${RETINA_REPO:?Set RETINA_REPO to the repository root}
export RETICULATE_PYTHON="${RETICULATE_PYTHON:-$(command -v python3)}"
export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-16}" OPENBLAS_NUM_THREADS="${SLURM_CPUS_PER_TASK:-16}"
mkdir -p "${WORKDIR}/atlas"
cd "${WORKDIR}"
echo "host=$(hostname) start=$(date)"
Rscript "$REPO/scripts/realign_mouse/mouse_cr9_01_build.R"
echo "end=$(date)"
