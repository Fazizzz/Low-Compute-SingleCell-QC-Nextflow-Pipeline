#!/usr/bin/env bash
set -euo pipefail
sample_id=$1; matrix_dir=$2; model=$3; outdir=${4:-.}
script_dir=$(cd "$(dirname "$0")/../python" && pwd)
mkdir -p "$outdir"
python3 "$script_dir/celltypist_annotate.py" \
    --sample_id "$sample_id" --matrix_dir "$matrix_dir" --model "$model" \
    --out_csv "$outdir/${sample_id}_celltypist_predictions.csv" \
    --out_meta "$outdir/${sample_id}_celltypist_meta.json"
