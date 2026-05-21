#!/usr/bin/env bash
# Standalone sample report.
# Args: sample_id species fastp.json chemistry.json kb_out cell_calling.json knee_data.json doublet_meta.json qc.csv qc_summary.json out.html
set -euo pipefail
script_dir=$(cd "$(dirname "$0")/../python" && pwd)
template="$script_dir/templates/sample_report.html.j2"
python3 "$script_dir/sample_report.py" \
    --sample_id "$1" --species "$2" \
    --fastp_json "$3" --chemistry_json "$4" \
    --kb_dir "$5" --cell_calling_json "$6" --knee_data_json "$7" \
    --doublet_meta "$8" --qc_csv "$9" --qc_summary "${10}" \
    --template "$template" --out_html "${11}"
