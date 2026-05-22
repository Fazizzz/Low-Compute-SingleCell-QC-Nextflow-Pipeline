#!/usr/bin/env python3
"""Sanity-check the user-declared chemistry against R1 read length.

This is a cheap pre-flight check: we read up to N reads from R1, compute the
modal read length, and compare to what the declared chemistry expects. It is
intentionally lightweight; kb-python's own --x parser does the authoritative
layout work later inside KB_COUNT. The output JSON is consumed by the
per-sample HTML report so the user sees whether the declared chemistry looks
plausible before they trust the cell-level metrics.

Supports every chemistry recognised by kb-python 0.29.5 / kallisto 0.51.1.
End-to-end QC validation has only been performed on 10xv3 so far; the
remaining chemistries pass through the same code path (kb count handles
them natively).
"""
import argparse
import gzip
import json
import sys
from collections import Counter

# Expected R1 length (in bp) per chemistry, derived from kb-python's
# technology table (`kb --list`). The R1 file is the one carrying barcode +
# UMI for most droplet chemistries; the value here is the minimum length
# of the barcode + UMI block on that file.
EXPECTED_R1_LEN = {
    '10xv1':        14,   # barcode-only on R1; UMI on R2
    '10xv2':        26,   # 16 BC + 10 UMI
    '10xv3':        28,   # 16 BC + 12 UMI
    '10xv3.1':      28,   # multiome GEX, same layout as v3
    '10xv3_ultima': 50,   # Ultima-sequenced 10x v3
    '10xv4':        28,   # 16 BC + 12 UMI
    'bdwta':        60,   # BD Rhapsody WTA: three BC blocks plus 8 nt UMI
    'celseq':       12,   # 8 BC + 4 UMI
    'celseq2':      12,   # 6 BC + 6 UMI
    'dropseq':      20,   # 12 BC + 8 UMI
    'indropsv1':    48,   # 11 BC + linker + 11 BC + UMI
    'indropsv3':    8,    # BC only on file 0 (UMI on file 1)
    'scrubseq':     16,   # 6 BC + 10 UMI
    'smartseq3':    19,   # UMI on R1 (no per-cell barcode)
    'surecell':     59,   # three BC blocks plus 8 nt UMI
    'visium':       28,   # 16 BC + 12 UMI (spatial)
}

# Chemistries where R1 length is NOT a deterministic indicator
# (no per-cell barcode on R1, or barcode spans multiple input files).
# For these we still emit a JSON record but skip the length comparison.
LAYOUT_VARIABLE = {
    'bulk',         # no barcode at all
    'smartseq2',    # plate-based, no per-cell barcode
    'indropsv2',    # barcode on file 1, not R1
    'split-seq',    # SPLiT-seq uses file 1 for barcodes
    'stormseq',     # UMI on file 1
}

# Length tolerance (bp). FASTP and other trimmers can trim a base or two
# off R1; we do not want to flag a sample for a 1 bp deviation.
LEN_TOLERANCE = 2

# Chemistries that share an R1 length (e.g. 10xv3 / 10xv3.1 / 10xv4 / visium
# all expect 28). These cannot be told apart from read length alone; we
# trust the user-declared one as long as the modal length matches.


def sample_lengths(r1_path, n):
    counts = Counter()
    opener = gzip.open if r1_path.endswith('.gz') else open
    with opener(r1_path, 'rt') as fh:
        line_no = 0
        seen = 0
        for line in fh:
            line_no += 1
            # FASTQ sequence lines are line numbers 2, 6, 10, ... (1-indexed)
            if line_no % 4 == 2:
                counts[len(line.rstrip('\n'))] += 1
                seen += 1
                if seen >= n:
                    break
    return counts


def evaluate(specified, modal_len):
    """Return (status, note). status is PASS / WARN / INFO / UNKNOWN."""
    if modal_len is None:
        return 'UNKNOWN', 'no reads sampled from R1'
    if specified in LAYOUT_VARIABLE:
        return 'INFO', (
            "R1 length not used as a check for chemistry "
            "'{spec}' (barcode / UMI layout is not on R1 alone); "
            "modal R1 length is {n} bp"
        ).format(spec=specified, n=modal_len)
    expected = EXPECTED_R1_LEN.get(specified)
    if expected is None:
        # Chemistry is in the validator's allow-list but we have no length
        # entry for it; let it through with an informational note.
        return 'INFO', (
            "no R1-length expectation registered for chemistry "
            "'{spec}'; modal R1 length is {n} bp"
        ).format(spec=specified, n=modal_len)
    if abs(modal_len - expected) <= LEN_TOLERANCE:
        return 'PASS', (
            "modal R1 length {n} bp matches expected {e} bp"
            .format(n=modal_len, e=expected)
        )
    return 'WARN', (
        "modal R1 length {n} bp does not match expected {e} bp "
        "for chemistry '{spec}'"
    ).format(n=modal_len, e=expected, spec=specified)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--sample_id', required=True)
    ap.add_argument('--r1', required=True)
    ap.add_argument('--chemistry', required=True, help='user-specified chemistry')
    ap.add_argument('--out', required=True)
    ap.add_argument('--n', type=int, default=10000)
    args = ap.parse_args()

    counts = sample_lengths(args.r1, n=args.n)
    modal_len = counts.most_common(1)[0][0] if counts else None
    status, note = evaluate(args.chemistry, modal_len)

    result = {
        'sample_id': args.sample_id,
        'specified_chemistry': args.chemistry,
        'modal_r1_length': modal_len,
        'expected_r1_length': EXPECTED_R1_LEN.get(args.chemistry),
        'reads_sampled': sum(counts.values()),
        'length_histogram': dict(counts),
        'status': status,
        'note': note,
    }
    with open(args.out, 'w') as fh:
        json.dump(result, fh, indent=2)
    print(
        "chemistry_check {sid}: status={st} ({note})".format(
            sid=args.sample_id, st=status, note=note
        ),
        file=sys.stderr,
    )


if __name__ == '__main__':
    main()
