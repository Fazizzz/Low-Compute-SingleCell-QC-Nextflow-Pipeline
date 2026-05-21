#!/usr/bin/env python3
"""Cell calling: auto knee for display + hard UMI threshold + Cell Ranger MTX output."""
import argparse, gzip, json, os, sys
import numpy as np
import scipy.io as sio
import scipy.sparse as sp


def detect_knee(umis_sorted: np.ndarray) -> dict:
    """Run KneeLocator on log-log rank curve; sanity-check the result."""
    from kneed import KneeLocator

    nonzero = umis_sorted[umis_sorted > 0]
    if nonzero.size < 10:
        return {"auto_knee_umi": None, "auto_knee_rank": None, "knee_status": "insufficient_data"}

    ranks = np.arange(1, nonzero.size + 1)
    log_ranks = np.log10(ranks)
    log_umis = np.log10(nonzero)

    try:
        kl = KneeLocator(
            log_ranks,
            log_umis,
            curve="convex",
            direction="decreasing",
        )
    except Exception as e:
        return {"auto_knee_umi": None, "auto_knee_rank": None, "knee_status": f"failed:{e}"}

    if kl.knee is None:
        return {"auto_knee_umi": None, "auto_knee_rank": None, "knee_status": "no_knee"}

    knee_rank = int(round(10 ** kl.knee))
    knee_rank = max(1, min(knee_rank, nonzero.size))
    knee_umi = int(nonzero[knee_rank - 1])

    p5, p95 = np.percentile(nonzero, [5, 95])
    status = "reliable"
    if knee_umi < p5 or knee_umi > p95:
        status = "unreliable"

    return {
        "auto_knee_umi": knee_umi,
        "auto_knee_rank": knee_rank,
        "knee_status": status,
    }


def rank_data_for_plot(umis_sorted: np.ndarray, threshold: int, max_points: int = 5000) -> dict:
    """Downsample log-log rank curve for plotting."""
    nonzero = umis_sorted[umis_sorted > 0]
    n = nonzero.size
    if n == 0:
        return {"ranks": [], "umis": [], "called_mask": []}
    if n > max_points:
        idx = np.unique(np.geomspace(1, n, max_points).astype(int)) - 1
    else:
        idx = np.arange(n)
    ranks = (idx + 1).tolist()
    umis = nonzero[idx].astype(int).tolist()
    called_mask = (nonzero[idx] >= threshold).tolist()
    return {"ranks": ranks, "umis": umis, "called_mask": called_mask}


def write_cellranger_mtx(matrix: sp.spmatrix, barcodes: list, features: list, outdir: str) -> None:
    """Write Cell Ranger style MTX: genes x cells, gzipped."""
    os.makedirs(outdir, exist_ok=True)
    mtx_path = os.path.join(outdir, "matrix.mtx")
    sio.mmwrite(mtx_path, matrix, field="integer")
    with open(mtx_path, "rb") as src, gzip.open(mtx_path + ".gz", "wb") as dst:
        dst.write(src.read())
    os.remove(mtx_path)

    with gzip.open(os.path.join(outdir, "barcodes.tsv.gz"), "wt") as fh:
        for bc in barcodes:
            fh.write(bc + "\n")

    with gzip.open(os.path.join(outdir, "features.tsv.gz"), "wt") as fh:
        for gid, gname in features:
            fh.write(f"{gid}\t{gname}\tGene Expression\n")


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--sample_id", required=True)
    ap.add_argument("--kb_dir", required=True, help="kb_out directory")
    ap.add_argument("--threshold", type=int, required=True)
    ap.add_argument("--outdir", default=".")
    args = ap.parse_args()

    unfiltered = os.path.join(args.kb_dir, "counts_unfiltered")
    # In nac (--workflow full) mode kb_python emits cells_x_genes.total.* (sum of
    # mature + nascent + ambiguous). In quick mode it emits cells_x_genes.* only.
    if os.path.exists(os.path.join(unfiltered, "cells_x_genes.total.mtx")):
        prefix = "cells_x_genes.total"
    else:
        prefix = "cells_x_genes"
    mtx_path = os.path.join(unfiltered, f"{prefix}.mtx")
    barcodes_path = os.path.join(unfiltered, f"{prefix}.barcodes.txt")
    genes_path = os.path.join(unfiltered, f"{prefix}.genes.txt")
    names_path = os.path.join(unfiltered, f"{prefix}.genes.names.txt")
    if not os.path.exists(barcodes_path):
        barcodes_path = os.path.join(unfiltered, "cells_x_genes.barcodes.txt")
    if not os.path.exists(genes_path):
        genes_path = os.path.join(unfiltered, "cells_x_genes.genes.txt")
    if not os.path.exists(names_path):
        names_path = os.path.join(unfiltered, "cells_x_genes.genes.names.txt")

    matrix = sio.mmread(mtx_path).tocsr()  # cells x genes
    with open(barcodes_path) as fh:
        barcodes = [l.strip() for l in fh if l.strip()]
    with open(genes_path) as fh:
        gene_ids = [l.strip() for l in fh if l.strip()]
    if os.path.exists(names_path):
        with open(names_path) as fh:
            gene_names = [l.strip() for l in fh if l.strip()]
    else:
        gene_names = gene_ids
    if len(gene_names) != len(gene_ids):
        gene_names = gene_ids
    features = list(zip(gene_ids, gene_names))

    umis_per_bc = np.asarray(matrix.sum(axis=1)).ravel().astype(np.int64)
    order = np.argsort(-umis_per_bc)
    umis_sorted = umis_per_bc[order]

    knee = detect_knee(umis_sorted)
    rank_curve = rank_data_for_plot(umis_sorted, args.threshold)

    keep_idx = np.where(umis_per_bc >= args.threshold)[0]
    cells_called = int(keep_idx.size)

    filtered = matrix[keep_idx, :].T.tocsc()  # transpose to genes x cells
    kept_barcodes = [barcodes[i] for i in keep_idx]

    filtered_dir = os.path.join(args.outdir, "filtered_matrix")
    write_cellranger_mtx(filtered, kept_barcodes, features, filtered_dir)

    summary = {
        "sample_id": args.sample_id,
        "threshold_used": args.threshold,
        "cells_called": cells_called,
        "n_barcodes_total": int(umis_per_bc.size),
        "median_umi_called": int(np.median(umis_per_bc[keep_idx])) if cells_called else 0,
        "total_umi_called": int(umis_per_bc[keep_idx].sum()) if cells_called else 0,
        **knee,
    }
    with open(os.path.join(args.outdir, "cell_calling.json"), "w") as fh:
        json.dump(summary, fh, indent=2)
    with open(os.path.join(args.outdir, "knee_data.json"), "w") as fh:
        json.dump({"threshold": args.threshold, **knee, **rank_curve}, fh)
    print(f"cell_calling {args.sample_id}: cells={cells_called} knee={knee}", file=sys.stderr)


if __name__ == "__main__":
    main()
