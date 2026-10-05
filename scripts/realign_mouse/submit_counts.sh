#!/bin/bash
# Submit `cellranger count` for every library whose FASTQs are ready.
# Idempotent: skips libraries already counted or still downloading.
# Counts emit BAMs (cellranger_count.sbatch sets --create-bam=true; CR9 default
# is false). Run from scripts/realign_mouse/ on a SLURM submit host.
#
# Usage:
#   bash submit_counts.sh            # submit all ready libraries
#   bash submit_counts.sh <crid>     # submit one specific library
set -euo pipefail
cd "$(dirname "$0")"
export WORKDIR=${WORKDIR:?Set WORKDIR to the mouse re-alignment working directory}
: "${CR_REF:?Set CR_REF to refdata-gex-GRCm39-2024-A}"
ROOT="${WORKDIR}"
mkdir -p "${ROOT}/logs"
ONLY="${1:-}"

# Chemistry per library (all 10x 3' v2 except bala = v3). cellranger also
# auto-detects, but we pin known values for reproducibility.
declare -A CHEM=(
  [clark_E14_rep1]=SC3Pv2 [clark_E14_rep2]=SC3Pv2 [clark_E16]=SC3Pv2
  [clark_E18_rep2]=SC3Pv2 [clark_E18_rep3]=SC3Pv2 [clark_P0]=SC3Pv2
  [wu_E13p5_atoh7het]=SC3Pv2
  [logiudice_E15p5_rep1]=SC3Pv2 [logiudice_E15p5_rep2]=SC3Pv2
  [bala_E13p5_WT]=SC3Pv3
)

for fqdir in "$ROOT"/*/*/fastq; do
    crid=$(basename "$(dirname "$fqdir")")
    [ -n "$ONLY" ] && [ "$crid" != "$ONLY" ] && continue
    case "$crid" in georges_*) continue;; esac           # deferred
    ls "$fqdir"/*_R1_*.fastq.gz >/dev/null 2>&1 || continue
    ls "$fqdir"/*_R2_*.fastq.gz >/dev/null 2>&1 || continue
    if [ -s "$ROOT/counts/$crid/outs/filtered_feature_bc_matrix.h5" ]; then echo "skip (counted): $crid"; continue; fi
    running=$(squeue --me -h -o '%j' 2>/dev/null)
    # a job named exactly $crid = its download is still running
    if grep -qx "$crid" <<<"$running"; then echo "skip (download running): $crid"; continue; fi
    # a job named ct_$crid = its count is already queued/running (avoid duplicate submit)
    if grep -qx "ct_$crid" <<<"$running"; then echo "skip (count already queued/running): $crid"; continue; fi
    chem=${CHEM[$crid]:-auto}
    echo "submit count: $crid (chem=$chem)"
    sbatch --job-name="ct_${crid}" \
        --output="${ROOT}/logs/count_%x_%j.out" \
        --error="${ROOT}/logs/count_%x_%j.err" \
        --export=ALL,CRID="$crid",FASTQ_DIR="$fqdir",CHEM="$chem" \
        cellranger_count.sbatch
done
echo "done. monitor: squeue --me"
