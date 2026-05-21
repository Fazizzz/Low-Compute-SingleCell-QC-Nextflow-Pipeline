#!/usr/bin/env python3
"""Per-cell QC metrics: UMI, genes detected, MT%, doublet score."""
import argparse, csv, gzip, json, os, sys
import numpy as np
import scipy.io as sio


MT_PREFIX = {"human": "MT-", "mouse": "mt-"}


def load_filtered(matrix_dir: str):
    """Return (cells_x_genes, barcodes, features)."""
    with gzip.open(os.path.join(matrix_dir, "matrix.mtx.gz"), "rt") as fh:
        mat = sio.mmread(fh).tocsr()  # genes x cells
    with gzip.open(os.path.join(matrix_dir, "barcodes.tsv.gz"), "rt") as fh:
        barcodes = [l.strip() for l in fh if l.strip()]
    features = []
    with gzip.open(os.path.join(matrix_dir, "features.tsv.gz"), "rt") as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            features.append((parts[0], parts[1] if len(parts) > 1 else parts[0]))
    return mat.T.tocsr(), barcodes, features  # cells x genes


def load_doublet_scores(csv_path: str) -> dict:
    scores = {}
    with open(csv_path) as fh:
        reader = csv.DictReader(fh)
        for row in reader:
            try:
                scores[row["barcode"]] = (float(row["score"]), int(row["predicted_doublet"]))
            except (ValueError, KeyError):
                continue
    return scores


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--sample_id", required=True)
    ap.add_argument("--matrix_dir", required=True)
    ap.add_argument("--doublet_csv", required=True)
    ap.add_argument("--species", required=True)
    ap.add_argument("--out_csv", required=True)
    ap.add_argument("--out_summary", required=True)
    args = ap.parse_args()

    cells_x_genes, barcodes, features = load_filtered(args.matrix_dir)
    gene_names = [g[1] for g in features]
    prefix = MT_PREFIX.get(args.species, "MT-")
    mt_mask = np.array([n.startswith(prefix) for n in gene_names], dtype=bool)
    n_cells, n_genes = cells_x_genes.shape

    umi_per_cell = np.asarray(cells_x_genes.sum(axis=1)).ravel().astype(np.int64)
    genes_per_cell = np.asarray((cells_x_genes > 0).sum(axis=1)).ravel().astype(np.int64)
    if mt_mask.any():
        mt_umi = np.asarray(cells_x_genes[:, mt_mask].sum(axis=1)).ravel().astype(np.int64)
    else:
        mt_umi = np.zeros(n_cells, dtype=np.int64)
    with np.errstate(divide="ignore", invalid="ignore"):
        pct_mt = np.where(umi_per_cell > 0, 100.0 * mt_umi / umi_per_cell, 0.0)

    doublet_map = load_doublet_scores(args.doublet_csv)

    with open(args.out_csv, "w", newline="") as fh:
        writer = csv.writer(fh)
        writer.writerow(["barcode", "umi", "genes", "pct_mt", "doublet_score", "predicted_doublet"])
        for i, bc in enumerate(barcodes):
            ds, pd_flag = doublet_map.get(bc, (None, None))
            writer.writerow([
                bc,
                int(umi_per_cell[i]),
                int(genes_per_cell[i]),
                round(float(pct_mt[i]), 4),
                "" if ds is None else round(ds, 6),
                "" if pd_flag is None else pd_flag,
            ])

    def q(arr, pct):
        return float(np.percentile(arr, pct)) if arr.size else 0.0

    mt_flag_likely_failed = bool(
        args.species in {"human", "mouse"} and mt_mask.any() == False
    ) or bool(args.species in {"human", "mouse"} and pct_mt.size and (pct_mt == 0).all())

    summary = {
        "sample_id": args.sample_id,
        "species": args.species,
        "n_cells": int(n_cells),
        "n_genes_detected": int((cells_x_genes.sum(axis=0) > 0).sum()),
        "median_umi": int(np.median(umi_per_cell)) if n_cells else 0,
        "mean_umi": float(np.mean(umi_per_cell)) if n_cells else 0.0,
        "p25_umi": q(umi_per_cell, 25),
        "p75_umi": q(umi_per_cell, 75),
        "median_genes": int(np.median(genes_per_cell)) if n_cells else 0,
        "mean_pct_mt": float(np.mean(pct_mt)) if n_cells else 0.0,
        "median_pct_mt": float(np.median(pct_mt)) if n_cells else 0.0,
        "mt_genes_in_features": int(mt_mask.sum()),
        "mt_detection_warning": mt_flag_likely_failed,
    }
    with open(args.out_summary, "w") as fh:
        json.dump(summary, fh, indent=2)
    print(f"qc_metrics {args.sample_id}: n_cells={n_cells} median_umi={summary['median_umi']}", file=sys.stderr)


if __name__ == "__main__":
    main()
