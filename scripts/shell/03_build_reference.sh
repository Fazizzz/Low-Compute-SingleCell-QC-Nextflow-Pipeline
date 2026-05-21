#!/usr/bin/env bash
# Build a kallisto index from genome FASTA + GTF.
#
# Usage:
#   03_build_reference.sh <genome.fa(.gz)> <genome.gtf(.gz)> <outdir> [--full]
#
# Without --full (default): builds a cDNA-only index for `--workflow quick`.
#   Outputs: index.idx, t2g.txt, cdna.fa
#
# With --full: builds a nac (spliced + nascent) index for `--workflow full`.
#   Outputs: index_nac.idx, t2g_nac.txt, cdna.fa, nascent.fa,
#            cdna_t2c.txt, nascent_t2c.txt
#
# IMPORTANT: --kallisto / --bustools are forced to the binaries on PATH (=
# conda's pinned kallisto 0.51.1). Without these flags, kb-python silently
# picks its bundled 0.52.0 binary, producing an index that crashes the
# pipeline's count step with std::bad_array_new_length.
set -euo pipefail

genome_fa=$1
genome_gtf=$2
outdir=${3:-.}
mode=${4:-quick}
case "$mode" in
    --full|full|nac) mode=full ;;
    *) mode=quick ;;
esac

mkdir -p "$outdir"
cd "$outdir"

if [ "$mode" = "full" ]; then
    kb ref --workflow nac \
        --kallisto "$(command -v kallisto)" \
        --bustools "$(command -v bustools)" \
        -i index_nac.idx \
        -g t2g_nac.txt \
        -f1 cdna.fa \
        -f2 nascent.fa \
        -c1 cdna_t2c.txt \
        -c2 nascent_t2c.txt \
        "$genome_fa" \
        "$genome_gtf"
else
    kb ref \
        --kallisto "$(command -v kallisto)" \
        --bustools "$(command -v bustools)" \
        -i index.idx \
        -g t2g.txt \
        -f1 cdna.fa \
        "$genome_fa" \
        "$genome_gtf"
fi
