#!/usr/bin/env python3
"""Optional CellTypist annotation. Graceful fail on <50 cells or any exception."""
import argparse, csv, gzip, json, os, sys


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--sample_id", required=True)
    ap.add_argument("--matrix_dir", required=True)
    ap.add_argument("--model", required=True)
    ap.add_argument("--out_csv", required=True)
    ap.add_argument("--out_meta", required=True)
    args = ap.parse_args()

    def write_skip(reason: str, n_cells: int = 0) -> None:
        with open(args.out_csv, "w") as fh:
            fh.write("barcode,predicted_labels,majority_voting\n")
        with open(args.out_meta, "w") as fh:
            json.dump(
                {"sample_id": args.sample_id, "status": "skipped",
                 "reason": reason, "n_cells": n_cells, "model": args.model},
                fh, indent=2,
            )

    try:
        import anndata as ad
        import scipy.io as sio
        import numpy as np

        with gzip.open(os.path.join(args.matrix_dir, "matrix.mtx.gz"), "rt") as fh:
            mat = sio.mmread(fh).tocsr()
        with gzip.open(os.path.join(args.matrix_dir, "barcodes.tsv.gz"), "rt") as fh:
            barcodes = [l.strip() for l in fh if l.strip()]
        features = []
        with gzip.open(os.path.join(args.matrix_dir, "features.tsv.gz"), "rt") as fh:
            for line in fh:
                parts = line.rstrip("\n").split("\t")
                features.append(parts[1] if len(parts) > 1 else parts[0])

        adata = ad.AnnData(X=mat.T.tocsr())  # cells x genes
        adata.obs_names = barcodes
        adata.var_names = features

        n_cells = adata.n_obs
        if n_cells < 50:
            write_skip(f"n_cells={n_cells} below threshold (50)", n_cells)
            return

        import scanpy as sc
        # CellTypist requires log1p-normalized expression to 1e4 counts/cell
        sc.pp.normalize_total(adata, target_sum=1e4)
        sc.pp.log1p(adata)

        import celltypist
        result = celltypist.annotate(adata, model=args.model, majority_voting=True)
        df = result.predicted_labels

    except Exception as e:
        write_skip(f"celltypist_failed:{e}")
        return

    with open(args.out_csv, "w", newline="") as fh:
        writer = csv.writer(fh)
        cols = ["barcode"] + list(df.columns)
        writer.writerow(cols)
        for bc, row in zip(df.index, df.itertuples(index=False)):
            writer.writerow([bc] + list(row))

    with open(args.out_meta, "w") as fh:
        json.dump(
            {"sample_id": args.sample_id, "status": "ok", "n_cells": int(adata.n_obs),
             "model": args.model, "labels_n_unique": int(df.iloc[:, 0].nunique())},
            fh, indent=2,
        )
    print(f"celltypist {args.sample_id}: ok", file=sys.stderr)


if __name__ == "__main__":
    main()
