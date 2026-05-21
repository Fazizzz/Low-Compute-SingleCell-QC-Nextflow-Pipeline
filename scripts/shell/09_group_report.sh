#!/usr/bin/env bash
# Standalone group report — pass space-separated lists per category.
# Usage: 09_group_report.sh <out.html> <software_versions.yml> <trace.txt|-> <chem_jsons> <cc_jsons> <dbl_metas> <qc_summaries> <run_infos>
set -euo pipefail
out=$1; sw=$2; trace=$3; chem=$4; cc=$5; dbl=$6; qc=$7; runs=$8
script_dir=$(cd "$(dirname "$0")/../python" && pwd)
template="$script_dir/templates/group_report.html.j2"
[[ "$trace" == "-" ]] && trace=""
python3 "$script_dir/group_report.py" \
    --chemistry_jsons $chem \
    --cell_calling_jsons $cc \
    --doublet_metas $dbl \
    --qc_summaries $qc \
    --run_infos $runs \
    ${trace:+--trace "$trace"} \
    --software_versions "$sw" \
    --template "$template" \
    --out_html "$out"
