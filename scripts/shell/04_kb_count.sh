#!/usr/bin/env bash
# Standalone kb count.
# Usage: 04_kb_count.sh <index.idx> <t2g.txt> <chemistry> <threads> <mem_gb> <outdir> <R1> <R2>
set -euo pipefail
index=$1; t2g=$2; chem=$3; threads=$4; mem=$5; outdir=$6; r1=$7; r2=$8
mkdir -p "$outdir"
kb count \
    --kallisto "$(command -v kallisto)" \
    --bustools "$(command -v bustools)" \
    -i "$index" \
    -g "$t2g" \
    -x "$chem" \
    -t "$threads" \
    -m "${mem}G" \
    -o "$outdir" \
    "$r1" "$r2"
