#!/usr/bin/env bash
# Standalone cell calling.
# Usage: 05_cell_calling.sh <sample_id> <kb_out_dir> <threshold> <outdir>
set -euo pipefail
sample_id=$1; kb_dir=$2; threshold=$3; outdir=${4:-.}
script_dir=$(cd "$(dirname "$0")/../python" && pwd)
mkdir -p "$outdir"
python3 "$script_dir/cell_calling.py" \
    --sample_id "$sample_id" --kb_dir "$kb_dir" --threshold "$threshold" --outdir "$outdir"
