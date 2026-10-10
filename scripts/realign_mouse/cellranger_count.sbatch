#!/bin/bash
#SBATCH --job-name=cr9_count
#SBATCH --time=0-12:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=96G
# cellranger uses --localcores=$SLURM_CPUS_PER_TASK with the resources above.
#SBATCH --output=cr9_count_%x_%j.out
#SBATCH --error=cr9_count_%x_%j.err
#
# Submit with:
#   sbatch --job-name=clark_E14_rep1 \
#          --export=ALL,CRID=clark_E14_rep1,FASTQ_DIR=${WORKDIR}/GSE118614/clark_E14_rep1/fastq,CHEM=SC3Pv2 \
#          cellranger_count.sbatch
#
# Required env vars:
#   CRID       - cellranger sample id (matches fastq filename prefix)
#   FASTQ_DIR  - directory with <CRID>_S1_L00X_R[12]_001.fastq.gz files
#   CHEM       - chemistry (SC3Pv2, SC3Pv3, or auto)

set -euo pipefail
WORKDIR=${WORKDIR:?Set WORKDIR to the mouse re-alignment working directory}

: "${CRID:?CRID required}"
: "${FASTQ_DIR:?FASTQ_DIR required}"
: "${CHEM:=auto}"

REF=${CR_REF:?Set CR_REF to refdata-gex-GRCm39-2024-A}
OUTROOT="${WORKDIR}/counts"
mkdir -p "${OUTROOT}"
cd "${OUTROOT}"


echo "[$(date)] cellranger count starting"
echo "  CRID=${CRID}"
echo "  FASTQ_DIR=${FASTQ_DIR}"
echo "  CHEM=${CHEM}"
echo "  REF=${REF}"
echo "  hostname=$(hostname)"
echo "  cores=${SLURM_CPUS_PER_TASK}"

# remove a stale partial run if present so cellranger doesn't refuse to start
if [ -d "${CRID}" ]; then
    echo "WARN: existing ${OUTROOT}/${CRID} -- moving aside"
    mv "${CRID}" "${CRID}.old.$(date +%s)"
fi

cellranger count \
    --id="${CRID}" \
    --transcriptome="${REF}" \
    --fastqs="${FASTQ_DIR}" \
    --sample="${CRID}" \
    --chemistry="${CHEM}" \
    --create-bam=true \
    --localcores="${SLURM_CPUS_PER_TASK}" \
    --localmem=90

echo "[$(date)] cellranger count done for ${CRID}"
ls -lh "${OUTROOT}/${CRID}/outs/" || true
