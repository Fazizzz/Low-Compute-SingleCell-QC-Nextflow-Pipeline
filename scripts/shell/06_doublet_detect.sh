#!/usr/bin/env bash
set -euo pipefail
sample_id=$1; matrix_dir=$2; outdir=${3:-.}
script_dir=$(cd "$(dirname "$0")/../python" && pwd)
mkdir -p "$outdir"
python3 "$script_dir/doublet_detect.py" \
    --sample_id "$sample_id" --matrix_dir "$matrix_dir" \
    --out_scores "$outdir/${sample_id}_doublet_scores.csv" \
    --out_meta   "$outdir/${sample_id}_doublet_meta.json"
