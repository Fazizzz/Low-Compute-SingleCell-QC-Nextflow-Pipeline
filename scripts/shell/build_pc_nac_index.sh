#!/usr/bin/env bash
# Build a protein-coding-only nac index for `--workflow full` on memory-constrained hosts.
#
# Why: the standard human nac index (kb ref --workflow nac on the full GTF) is
# ~6 GB and needs ~10 to 12 GB of free RAM to load via kallisto bus. Filtering
# the GTF to protein-coding genes drops the combined cdna + nascent FASTA from
# ~3.2 GB to ~2.0 GB, producing an index in the ~3 to 4 GB range.
#
# Usage:
#   ./build_pc_nac_index.sh <genome.fa(.gz)> <annotation.gtf(.gz)> <outdir>
#
# Outputs (placed in <outdir>):
#   index_nac.idx       kallisto nac index (protein-coding only)
#   t2g_nac.txt         transcript-to-gene mapping for nac mode
#   cdna_t2c.txt        cdna transcript-to-capture file
#   nascent_t2c.txt     nascent transcript-to-capture file
#
# Pass these to the pipeline as:
#   --workflow full
#   --prebuilt_index <outdir>/index_nac.idx
#   --t2g <outdir>/t2g_nac.txt
#   --cdna_t2c <outdir>/cdna_t2c.txt
#   --nascent_t2c <outdir>/nascent_t2c.txt

set -euo pipefail

if [ $# -lt 3 ]; then
    echo "Usage: $0 <genome.fa(.gz)> <annotation.gtf(.gz)> <outdir>" >&2
    exit 1
fi

# Resolve to absolute paths BEFORE cd-ing, because relative args become invalid
# after we change into $OUTDIR.
GENOME=$(cd "$(dirname "$1")" && pwd)/$(basename "$1")
GTF=$(cd "$(dirname "$2")" && pwd)/$(basename "$2")
mkdir -p "$3"
OUTDIR=$(cd "$3" && pwd)
cd "$OUTDIR"

KALLISTO=$(command -v kallisto)
BUSTOOLS=$(command -v bustools)
KB=$(command -v kb)

if [ -z "$KALLISTO" ] || [ -z "$BUSTOOLS" ] || [ -z "$KB" ]; then
    echo "kallisto / bustools / kb must be on PATH (use the scq conda env)" >&2
    exit 1
fi

# Step 1: filter GTF to protein-coding only
PC_GTF="${OUTDIR}/$(basename "$GTF" .gz).pc.gz"
if [ ! -s "$PC_GTF" ]; then
    echo "[$(date +%H:%M:%S)] filtering GTF to protein-coding..."
    { gunzip -c "$GTF" | grep "^#"; \
      gunzip -c "$GTF" | grep 'gene_biotype "protein_coding"'; } \
        | gzip > "$PC_GTF"
    echo "[$(date +%H:%M:%S)] filtered GTF: $(ls -lh "$PC_GTF" | awk '{print $5}')"
fi

# Step 2: run kb ref --workflow nac to produce the intermediate FASTAs and t2c files
echo "[$(date +%H:%M:%S)] running kb ref --workflow nac (this will start a slow d-list step;"
echo "                    we will kill it and rebuild the index without d-list)..."

"$KB" ref --workflow nac \
    --kallisto "$KALLISTO" --bustools "$BUSTOOLS" \
    --tmp "${OUTDIR}/tmp_kbref" \
    -i index_nac.idx -g t2g_nac.txt \
    -f1 cdna_pc.fa -f2 nascent_pc.fa \
    -c1 cdna_t2c.txt -c2 nascent_t2c.txt \
    "$GENOME" "$PC_GTF" &
KB_PID=$!

# Wait until kb_ref has produced the t2c files + concatenated combined FASTA, then kill it
echo "[$(date +%H:%M:%S)] waiting for intermediate FASTAs and t2c files..."
while ! { [ -s cdna_pc.fa ] && [ -s nascent_pc.fa ] && [ -s cdna_t2c.txt ] && [ -s nascent_t2c.txt ] && [ -s t2g_nac.txt ]; }; do
    sleep 10
    if ! kill -0 "$KB_PID" 2>/dev/null; then
        echo "kb ref exited before intermediates were complete" >&2
        exit 1
    fi
done
echo "[$(date +%H:%M:%S)] intermediates ready; killing slow d-list step..."
pkill -P "$KB_PID" -f "kallisto index" 2>/dev/null || true
kill "$KB_PID" 2>/dev/null || true
wait "$KB_PID" 2>/dev/null || true

# Step 3: build the kallisto index without -d (d-list) on the combined FASTA
# kb_python concatenates cdna + nascent into the combined FASTA before
# passing to kallisto. We re-do the concatenation explicitly.
echo "[$(date +%H:%M:%S)] concatenating cdna_pc.fa + nascent_pc.fa for indexing..."
cat cdna_pc.fa nascent_pc.fa > combined_pc.fa

echo "[$(date +%H:%M:%S)] building kallisto index (no d-list)..."
"$KALLISTO" index -i index_nac.idx -k 31 -t 4 combined_pc.fa
rm -f combined_pc.fa

echo "[$(date +%H:%M:%S)] done."
ls -lh index_nac.idx t2g_nac.txt cdna_t2c.txt nascent_t2c.txt
