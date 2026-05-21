#!/usr/bin/env bash
# Standalone fastp QC wrapper — run outside the pipeline.
# Usage: 01_fastp_qc.sh <sample_id> <R1.fastq.gz> <R2.fastq.gz> <threads> <outdir>
set -euo pipefail

sample_id=$1
r1=$2
r2=$3
threads=${4:-4}
outdir=${5:-.}
mkdir -p "$outdir"

fastp \
    --in1 "$r1" \
    --in2 "$r2" \
    --disable_adapter_trimming \
    --disable_quality_filtering \
    --disable_length_filtering \
    --thread "$threads" \
    --json "$outdir/${sample_id}_fastp.json" \
    --html "$outdir/${sample_id}_fastp.html"
