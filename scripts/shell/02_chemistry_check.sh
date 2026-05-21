#!/usr/bin/env bash
# Standalone chemistry check.
# Usage: 02_chemistry_check.sh <sample_id> <R1.fastq.gz> <chemistry> <out.json>
set -euo pipefail
sample_id=$1
r1=$2
chemistry=$3
out=$4
script_dir=$(cd "$(dirname "$0")/../python" && pwd)
python3 "$script_dir/chemistry_check.py" \
    --sample_id "$sample_id" \
    --r1 "$r1" \
    --chemistry "$chemistry" \
    --out "$out"
