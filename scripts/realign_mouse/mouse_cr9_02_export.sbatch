#!/bin/bash
#SBATCH -J mcr9_export
#SBATCH -t 0-3:00
#SBATCH -c 8
#SBATCH --mem 128G
#SBATCH -o export_%j.out
#SBATCH -e export_%j.err
set -euo pipefail
WORKDIR=${WORKDIR:?Set WORKDIR to the mouse re-alignment working directory}
REPO=${RETINA_REPO:?Set RETINA_REPO to the repository root}
export RETICULATE_PYTHON="${RETICULATE_PYTHON:-$(command -v python3)}"
export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-8}" OPENBLAS_NUM_THREADS="${SLURM_CPUS_PER_TASK:-8}"
mkdir -p "${WORKDIR}/atlas"
cd "${WORKDIR}"
echo "host=$(hostname) start=$(date)"

Rscript "$REPO/scripts/realign_mouse/mouse_cr9_02_export_rpc.R"

mkdir -p "$WORKDIR/atlas/h5ad"
echo "=== assemble RPC (all-stage) ==="
python3 "$REPO/scripts/realign_mouse/mouse_cr9_03_assemble_h5ad.py" \
  --dir "$WORKDIR/atlas/export_rpc" --prefix "20260528_mouse_RPC_cr9_" \
  --out "$WORKDIR/atlas/h5ad/20260528_mouse_RPC_cr9.h5ad"
echo "=== assemble RPC (E13-E16) ==="
python3 "$REPO/scripts/realign_mouse/mouse_cr9_03_assemble_h5ad.py" \
  --dir "$WORKDIR/atlas/export_rpc_e13e16" --prefix "20260528_mouse_RPC_cr9_e13e16_" \
  --out "$WORKDIR/atlas/h5ad/20260528_mouse_RPC_cr9_e13e16.h5ad"
echo "=== assemble full atlas ==="
python3 "$REPO/scripts/realign_mouse/mouse_cr9_03_assemble_h5ad.py" \
  --dir "$WORKDIR/atlas/export_full" --prefix "mouse_cr9_full_" \
  --out "$WORKDIR/atlas/h5ad/20260528_mouse_atlas_cr9_full.h5ad"

echo "end=$(date)"
ls -la "$WORKDIR/atlas/h5ad/"
