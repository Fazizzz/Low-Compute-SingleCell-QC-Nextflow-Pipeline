#!/usr/bin/env bash
set -euo pipefail
sample_id=$1; matrix_dir=$2; doublet_csv=$3; species=$4; outdir=${5:-.}
script_dir=$(cd "$(dirname "$0")/../python" && pwd)
mkdir -p "$outdir"
python3 "$script_dir/qc_metrics.py" \
    --sample_id "$sample_id" --matrix_dir "$matrix_dir" \
    --doublet_csv "$doublet_csv" --species "$species" \
    --out_csv "$outdir/${sample_id}_qc_metrics.csv" \
    --out_summary "$outdir/${sample_id}_qc_summary.json"
