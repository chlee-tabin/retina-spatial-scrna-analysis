#!/bin/bash
#SBATCH --job-name=sra_dl
#SBATCH --time=0-12:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --output=sra_dl_%x_%j.out
#SBATCH --error=sra_dl_%x_%j.err
#
# Submit:
#   sbatch --job-name=clark_E14_rep1 \
#     --export=ALL,SRR=SRR7699772,CRID=clark_E14_rep1,GSE=GSE118614 \
#     download_fastq.sbatch
#
# Required env vars:
#   SRR  - SRA run accession
#   CRID - cellranger sample id (used as fastq file prefix)
#   GSE  - GEO series (used for output subdir)

set -euo pipefail
WORKDIR=${WORKDIR:?Set WORKDIR to the mouse re-alignment working directory}

: "${SRR:?SRR required}"
: "${CRID:?CRID required}"
: "${GSE:?GSE required}"


OUTROOT="${WORKDIR}"
TMPDIR="${TMPDIR:-${WORKDIR}/tmp/${CRID}}"
mkdir -p "${TMPDIR}"
DEST="${OUTROOT}/${GSE}/${CRID}/fastq"
mkdir -p "${DEST}"
cd "${DEST}"

echo "[$(date)] === ${CRID} (${SRR}) ==="

# prefetch uses the explicit output directory below for this run.
echo "[$(date)] prefetch starting"
prefetch --max-size 100g --output-directory "${DEST}" "${SRR}"

echo "[$(date)] fasterq-dump starting"
# 10x scRNA-seq: --split-files gives R1/R2/I1 separately
fasterq-dump --split-files --include-technical \
    --threads "${SLURM_CPUS_PER_TASK}" \
    --temp "${TMPDIR}" \
    --outdir "${DEST}" \
    "${DEST}/${SRR}/${SRR}.sra"

echo "[$(date)] post-process: rename + gzip"
# Inspect what came out
ls -lh "${DEST}/${SRR}"_*.fastq 2>/dev/null || ls -lh "${DEST}/"

# Standardize naming for cellranger.
# 10x v2/v3 outputs from SRA dump come as <SRR>_1.fastq (typically R1 = bc+UMI, 26bp),
# <SRR>_2.fastq (R2 = cDNA, 98bp), and sometimes <SRR>_3.fastq (I1 = 8bp index).
# But occasionally the dump labels them differently. We'll detect by length:
#   - shortest = R1 (cell barcode + UMI)
#   - longest  = R2 (cDNA insert)
#   - middle (if 3 files) = I1 (sample index)
# To keep it simple here, we use the natural numerical order (1=R1, 2=R2, 3=I1)
# which holds for the SRA dumps we expect from Clark, Wu, and Balasubramanian.
for f in "${DEST}/${SRR}"_*.fastq; do
    [ -f "$f" ] || continue
    base=$(basename "$f" .fastq)
    suffix=${base##*_}
    case "$suffix" in
        1) newname="${CRID}_S1_L001_R1_001.fastq" ;;
        2) newname="${CRID}_S1_L001_R2_001.fastq" ;;
        3) newname="${CRID}_S1_L001_I1_001.fastq" ;;
        *) echo "WARN: unknown suffix ${suffix} on ${f}" >&2; continue ;;
    esac
    mv -v "${f}" "${DEST}/${newname}"
    pigz -p "${SLURM_CPUS_PER_TASK}" "${DEST}/${newname}" || gzip "${DEST}/${newname}"
done

# Clean up the .sra cache directory to save space
rm -rf "${DEST}/${SRR}"

echo "[$(date)] done ${CRID}"
ls -lh "${DEST}/"
