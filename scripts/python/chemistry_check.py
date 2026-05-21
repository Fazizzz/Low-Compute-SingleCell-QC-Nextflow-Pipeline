#!/usr/bin/env python3
"""Detect chemistry from R1 read length distribution."""
import argparse, gzip, json, sys
from collections import Counter

LEN_TO_CHEM = {
    26: "10xv2",
    28: "10xv3-class",  # cannot distinguish v3 / v3.1 / Multiome
    20: "dropseq",
}

V3_CLASS = {"10xv3", "10xv3.1"}


def sample_lengths(r1_path: str, n: int = 10000) -> Counter:
    counts: Counter = Counter()
    opener = gzip.open if r1_path.endswith(".gz") else open
    with opener(r1_path, "rt") as fh:
        line_no = 0
        seen = 0
        for line in fh:
            line_no += 1
            # FASTQ: lines 2, 6, 10, ... are sequences (1-indexed) → (line_no % 4) == 2
            if line_no % 4 == 2:
                counts[len(line.rstrip("\n"))] += 1
                seen += 1
                if seen >= n:
                    break
    return counts


def classify(counts: Counter) -> str:
    if not counts:
        return "unknown"
    modal_len, _ = counts.most_common(1)[0]
    return LEN_TO_CHEM.get(modal_len, "unknown")


def compare(detected: str, specified: str) -> str:
    if detected == "unknown":
        return "UNKNOWN"
    if detected == "10xv3-class" and specified in V3_CLASS:
        return "PASS"
    if detected == specified:
        return "PASS"
    return "WARN"


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--sample_id", required=True)
    ap.add_argument("--r1", required=True)
    ap.add_argument("--chemistry", required=True, help="user-specified chemistry")
    ap.add_argument("--out", required=True)
    ap.add_argument("--n", type=int, default=10000)
    args = ap.parse_args()

    counts = sample_lengths(args.r1, n=args.n)
    detected = classify(counts)
    status = compare(detected, args.chemistry)

    result = {
        "sample_id": args.sample_id,
        "specified_chemistry": args.chemistry,
        "detected_chemistry": detected,
        "modal_length": counts.most_common(1)[0][0] if counts else None,
        "reads_sampled": sum(counts.values()),
        "length_histogram": dict(counts),
        "status": status,
    }
    with open(args.out, "w") as fh:
        json.dump(result, fh, indent=2)
    print(f"chemistry_check {args.sample_id}: detected={detected} status={status}", file=sys.stderr)


if __name__ == "__main__":
    main()
