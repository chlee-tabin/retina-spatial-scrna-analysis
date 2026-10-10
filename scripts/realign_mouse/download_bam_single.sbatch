#!/bin/bash
#SBATCH --job-name=ena_bam1
#SBATCH --time=0-08:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --output=ena_bam1_%x_%j.out
#SBATCH --error=ena_bam1_%x_%j.err
#
# Download a SINGLE original 10x BAM from ENA and reconstruct FASTQs via
# cellranger bamtofastq. For BAM-submitted SRA runs (Clark GSE118614, Bala
# GSE139904) where fasterq-dump yields cDNA-only (no cell barcode).
#
# Submit:
#   sbatch --job-name=clark_E16 \
#     --export=ALL,BAM_URL=ftp.sra.ebi.ac.uk/vol1/run/SRR769/SRR7699774/E16.bam.1,CRID=clark_E16,SUBDIR=GSE118614 \
#     download_bam_single.sbatch
#
# Env vars:
#   BAM_URL - ENA ftp path (without scheme) to the .bam.1 file
#   CRID    - cellranger sample id (fastq prefix)
#   SUBDIR  - working-directory grouping (GSE118614 / GSE139904)

set -euo pipefail
WORKDIR=${WORKDIR:?Set WORKDIR to the mouse re-alignment working directory}
: "${BAM_URL:?BAM_URL required}"
: "${CRID:?CRID required}"
: "${SUBDIR:?SUBDIR required}"

OUTROOT="${WORKDIR}"
BAMDIR="${OUTROOT}/${SUBDIR}/${CRID}/bam"
FQDIR="${OUTROOT}/${SUBDIR}/${CRID}/fastq"
mkdir -p "${BAMDIR}" "${FQDIR}"

echo "[$(date)] === ${CRID} ==="
BAM_LOCAL="${BAMDIR}/${CRID}.bam"
if [ ! -s "${BAM_LOCAL}" ]; then
    echo "[$(date)] downloading ${BAM_URL}"
    wget -q -O "${BAM_LOCAL}.tmp" "https://${BAM_URL}" && mv "${BAM_LOCAL}.tmp" "${BAM_LOCAL}"
fi
echo "[$(date)] BAM size: $(du -h "${BAM_LOCAL}" | cut -f1)"


# cellranger bamtofastq reconstructs the original 10x FASTQs (incl. barcode read)
# from the cellranger-produced BAM, using the @CO "10x_bam_to_fastq" header to
# detect the read structure. We do NOT auto-fallback to --cr11: that flag is for
# SC3' v1 (14 bp barcode) and would silently mis-parse v2/v3 barcodes. If the
# default fails (e.g. a re-sorted BAM with the @CO header stripped, as seen for
# the Georges E-MTAB-9395 BAMs), we dump the header for manual diagnosis instead.
echo "[$(date)] cellranger bamtofastq"
rm -rf "${BAMDIR}/bf"
if ! cellranger bamtofastq --nthreads="${SLURM_CPUS_PER_TASK}" \
        "${BAM_LOCAL}" "${BAMDIR}/bf"; then
    echo "ERROR: bamtofastq failed for ${CRID}. BAM @CO/header dump follows:" >&2
    samtools view -H "${BAM_LOCAL}" | grep -E '^@CO|^@PG|^@RG' >&2 || true
    echo "ERROR: leaving BAM in place at ${BAM_LOCAL} for diagnosis; not deleting." >&2
    exit 2
fi

# Collect FASTQs; renumber lanes sequentially so cellranger count sees one sample.
echo "[$(date)] collecting FASTQs into ${FQDIR}"
LANE_IDX=1
while IFS= read -r r1; do
    r2=${r1/_R1_/_R2_}
    i1=${r1/_R1_/_I1_}
    lane=$(printf "%03d" "${LANE_IDX}")
    cp -v "${r1}" "${FQDIR}/${CRID}_S1_L${lane}_R1_001.fastq.gz"
    [ -f "${r2}" ] && cp -v "${r2}" "${FQDIR}/${CRID}_S1_L${lane}_R2_001.fastq.gz"
    [ -f "${i1}" ] && cp -v "${i1}" "${FQDIR}/${CRID}_S1_L${lane}_I1_001.fastq.gz"
    LANE_IDX=$((LANE_IDX+1))
done < <(find "${BAMDIR}/bf" -name "*_R1_*.fastq.gz" | sort)

# Sanity: report R1/R2 read lengths so we can confirm barcode recovery.
echo "[$(date)] read-length sanity:"
for r in R1 R2; do
    f=$(ls "${FQDIR}/${CRID}_S1_L001_${r}_001.fastq.gz" 2>/dev/null | head -1)
    [ -n "${f}" ] && echo "  ${r}: $(zcat "${f}" | head -400 | awk 'NR%4==2{print length($0)}' | sort -n | uniq -c | tr '\n' ' ')"
done

# Free space: drop the BAM + intermediate now that FASTQs are extracted.
rm -rf "${BAMDIR}"
echo "[$(date)] done ${CRID}"
ls -lh "${FQDIR}/"
