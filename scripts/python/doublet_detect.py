#!/usr/bin/env python3
"""Doublet detection via scanpy's scrublet wrapper. Graceful fail on any error."""
import argparse, gzip, json, os, sys
import numpy as np
import scipy.io as sio


def load_filtered(matrix_dir: str):
    with gzip.open(os.path.join(matrix_dir, "matrix.mtx.gz"), "rt") as fh:
        mat = sio.mmread(fh).tocsr()  # genes x cells
    with gzip.open(os.path.join(matrix_dir, "barcodes.tsv.gz"), "rt") as fh:
        barcodes = [l.strip() for l in fh if l.strip()]
    features = []
    with gzip.open(os.path.join(matrix_dir, "features.tsv.gz"), "rt") as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            features.append(parts[1] if len(parts) > 1 else parts[0])
    return mat.T.tocsr(), barcodes, features  # cells x genes


def write_skip(out_meta, out_scores, sample_id, status, reason, n_cells):
    with open(out_meta, "w") as fh:
        json.dump({
            "sample_id": sample_id, "status": status, "reason": reason,
            "n_cells": n_cells, "doublet_rate": None,
            "expected_doublet_rate": None, "threshold": None,
            "backend": "scanpy.pp.scrublet",
        }, fh, indent=2)
    with open(out_scores, "w") as fh:
        fh.write("barcode,score,predicted_doublet\n")


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--sample_id", required=True)
    ap.add_argument("--matrix_dir", required=True)
    ap.add_argument("--out_scores", required=True)
    ap.add_argument("--out_meta", required=True)
    args = ap.parse_args()

    try:
        cells_x_genes, barcodes, features = load_filtered(args.matrix_dir)
    except Exception as e:
        write_skip(args.out_meta, args.out_scores, args.sample_id,
                   "skipped", f"load_failed:{e}", 0)
        return

    n_cells = cells_x_genes.shape[0]
    if n_cells < 100:
        write_skip(args.out_meta, args.out_scores, args.sample_id,
                   "skipped", f"n_cells={n_cells} below threshold (100)", n_cells)
        return

    try:
        import anndata as ad
        import scanpy as sc
        adata = ad.AnnData(X=cells_x_genes)
        adata.obs_names = barcodes
        adata.var_names = features
        # silence verbosity; scanpy logs to stderr
        sc.settings.verbosity = 1
        with np.errstate(invalid="ignore", divide="ignore"):
            sc.pp.scrublet(
                adata,
                expected_doublet_rate=0.06,
                random_state=42,
                verbose=False,
            )
        scores = adata.obs["doublet_score"].to_numpy()
        predicted = adata.obs["predicted_doublet"].to_numpy()
        if scores is None or not np.isfinite(scores).any():
            raise RuntimeError("scrublet returned non-finite scores")
        meta = adata.uns.get("scrublet", {})
        params = meta.get("parameters", {}) if isinstance(meta, dict) else {}
        threshold = float(meta.get("threshold", float("nan"))) if isinstance(meta, dict) else float("nan")
        if np.isnan(threshold):
            threshold = float(params.get("threshold", float("nan")))
    except Exception as e:
        write_skip(args.out_meta, args.out_scores, args.sample_id,
                   "skipped", f"scrublet_failed:{type(e).__name__}:{e}", n_cells)
        return

    doublet_rate = float(np.mean(predicted)) if predicted is not None else None

    # Sanity check: real biology rarely exceeds ~10-15% doublets. Anything above
    # 25% almost always means scrublet's bimodal threshold finder picked a bad
    # cutoff (poor data quality, low alignment, weak structure).
    status = "ok"
    quality_reason = None
    if doublet_rate is not None and doublet_rate > 0.25:
        status = "unreliable"
        quality_reason = (
            f"doublet rate {doublet_rate * 100:.1f}% implausibly high; "
            "scrublet threshold likely picked noise floor (poor bimodal structure)"
        )

    with open(args.out_scores, "w") as fh:
        fh.write("barcode,score,predicted_doublet\n")
        for bc, sc_, pred in zip(barcodes, scores, predicted):
            fh.write(f"{bc},{float(sc_):.6f},{int(bool(pred))}\n")

    with open(args.out_meta, "w") as fh:
        json.dump({
            "sample_id": args.sample_id,
            "status": status,
            "reason": quality_reason,
            "n_cells": int(n_cells),
            "doublet_rate": doublet_rate,
            "expected_doublet_rate": 0.06,
            "threshold": None if np.isnan(threshold) else threshold,
            "backend": "scanpy.pp.scrublet",
        }, fh, indent=2)
    print(f"doublet_detect {args.sample_id}: rate={doublet_rate}", file=sys.stderr)


if __name__ == "__main__":
    main()
