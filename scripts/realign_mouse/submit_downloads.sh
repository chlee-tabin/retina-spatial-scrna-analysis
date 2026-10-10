#!/bin/bash
# Submit all download jobs for the 10-library atlas realignment set, routed by how each
# dataset was submitted to the archive (decisive for whether barcodes survive):
#
#   bam_single (download_bam_single.sbatch): BAM-submitted SRA runs where
#       fasterq-dump yields cDNA-only. Recover via cellranger bamtofastq.
#       -> Clark x6 (CR 1.2.1 BAMs) + Bala (CR possorted BAM)
#   fasterq    (download_fastq.sbatch): FASTQ-submitted runs; fasterq-dump
#       --split-files --include-technical preserves R1(bc)+R2(cDNA)+I1.
#       -> Wu, Lo Giudice x2
# The two E-MTAB-9395 candidates in manifest.tsv are excluded because their
# BAM headers cannot be parsed by cellranger bamtofastq.
#
# Usage: bash submit_downloads.sh [canary|clark|fasterq|all]
set -euo pipefail
export WORKDIR=${WORKDIR:?Set WORKDIR to the mouse re-alignment working directory}
mkdir -p "${WORKDIR}/logs"
SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
cd "${SCRIPT_DIR}"
MODE="${1:-all}"

submit_bam_single () {  # crid  ena_bam_path  subdir
    sbatch --job-name="$1" \
        --output="${WORKDIR}/logs/download_%x_%j.out" \
        --error="${WORKDIR}/logs/download_%x_%j.err" \
        --export=ALL,BAM_URL="$2",CRID="$1",SUBDIR="$3" \
        download_bam_single.sbatch
}
submit_fasterq () {     # crid  srr  gse
    sbatch --job-name="$1" \
        --output="${WORKDIR}/logs/download_%x_%j.out" \
        --error="${WORKDIR}/logs/download_%x_%j.err" \
        --export=ALL,SRR="$2",CRID="$1",GSE="$3" \
        download_fastq.sbatch
}

case "$MODE" in
  canary)
    # Prove each route on the smallest sample before committing the rest.
    submit_bam_single clark_E16 ftp.sra.ebi.ac.uk/vol1/run/SRR769/SRR7699774/E16.bam.1 GSE118614
    submit_fasterq    logiudice_E15p5_rep1 SRR8181428 GSE122466
    ;;
  clark)
    submit_bam_single clark_E14_rep1 ftp.sra.ebi.ac.uk/vol1/run/SRR769/SRR7699772/E14_rep1.bam.1 GSE118614
    submit_bam_single clark_E14_rep2 ftp.sra.ebi.ac.uk/vol1/run/SRR769/SRR7699773/E14_rep2.bam.1 GSE118614
    submit_bam_single clark_E18_rep2 ftp.sra.ebi.ac.uk/vol1/run/SRR769/SRR7699775/E18_rep2.bam.1 GSE118614
    submit_bam_single clark_E18_rep3 ftp.sra.ebi.ac.uk/vol1/run/SRR769/SRR7699776/E18_rep3.bam.1 GSE118614
    submit_bam_single clark_P0       ftp.sra.ebi.ac.uk/vol1/run/SRR769/SRR7699777/P0.bam.1       GSE118614
    submit_bam_single bala_E13p5_WT  ftp.sra.ebi.ac.uk/vol1/run/SRR103/SRR10398966/possorted_genome_bam.bam.1 GSE139904
    ;;
  fasterq)
    submit_fasterq logiudice_E15p5_rep2 SRR8181429 GSE122466
    submit_fasterq wu_E13p5_atoh7het SRR11582194 GSE149040
    ;;
  all)
    bash "${SCRIPT_DIR}/submit_downloads.sh" canary
    bash "${SCRIPT_DIR}/submit_downloads.sh" clark
    bash "${SCRIPT_DIR}/submit_downloads.sh" fasterq
    ;;
  *) echo "unknown mode: $MODE"; exit 1 ;;
esac
echo "Submitted ($MODE). Check: squeue --me"
