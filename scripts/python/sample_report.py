#!/usr/bin/env python3
"""Render per-sample interactive HTML report."""
import argparse, csv, datetime, json, os, sys
import numpy as np
import plotly.graph_objects as go
import plotly.offline as pyo
from jinja2 import Environment, FileSystemLoader, select_autoescape


def fig_to_dict(fig: go.Figure) -> dict:
    return json.loads(fig.to_json())


def build_knee_plot(knee_data: dict) -> dict:
    ranks = np.asarray(knee_data.get("ranks") or [])
    umis  = np.asarray(knee_data.get("umis")  or [])
    mask  = np.asarray(knee_data.get("called_mask") or [], dtype=bool)
    if ranks.size == 0:
        return {"data": [], "layout": {"title": "No barcode data"}}

    fig = go.Figure()
    if (~mask).any():
        fig.add_trace(go.Scattergl(x=ranks[~mask], y=umis[~mask], mode="markers",
                                   marker=dict(size=4, color="lightgrey"), name="background"))
    if mask.any():
        fig.add_trace(go.Scattergl(x=ranks[mask], y=umis[mask], mode="markers",
                                   marker=dict(size=4, color="#1f77b4"), name="cells called"))

    threshold = knee_data.get("threshold")
    if threshold:
        fig.add_hline(y=threshold, line_dash="dash", line_color="#c0392b",
                      annotation_text=f"filter threshold ({threshold} UMI)",
                      annotation_position="top right",
                      annotation_font=dict(size=11, color="#c0392b"))

    auto_umi = knee_data.get("auto_knee_umi")
    knee_status = knee_data.get("knee_status", "")
    if auto_umi and knee_status == "reliable":
        fig.add_hline(y=auto_umi, line_dash="dot", line_color="#e67e22",
                      annotation_text=f"predicted knee (advisory): {auto_umi} UMI",
                      annotation_position="bottom right",
                      annotation_font=dict(size=11, color="#e67e22"))
    elif auto_umi:
        fig.add_hline(y=auto_umi, line_dash="dot", line_color="#bdc3c7",
                      annotation_text=f"predicted knee (unreliable, advisory only): {auto_umi} UMI",
                      annotation_position="bottom right",
                      annotation_font=dict(size=10, color="#7f8c8d"))

    cells_called = int(mask.sum())
    if cells_called > 0:
        x_max = min(int(ranks.max()), max(cells_called * 10, 1000))
    else:
        x_max = int(ranks.max())
    x_min = 1
    log_x_range = [np.log10(x_min), np.log10(x_max)]

    visible_umis = umis[ranks <= x_max]
    y_min = max(1, int(visible_umis.min())) if visible_umis.size else 1
    y_max = max(int(visible_umis.max()), threshold or 1) if visible_umis.size else 1
    log_y_range = [np.log10(y_min) - 0.1, np.log10(y_max) + 0.1]

    log_axis = lambda title, rng: dict(
        title=title,
        type="log",
        range=rng,
        dtick=1,                     # one major tick per decade
        exponentformat="power",      # render as 10^n
        showexponent="all",
        minor=dict(showgrid=False, ticks=""),  # suppress 1/2/5 subticks
        gridcolor="#e6e9ef",
        showline=True, linecolor="black", linewidth=1, mirror=True,
        ticks="outside",
    )

    fig.update_layout(
        title="Barcode rank plot",
        xaxis=log_axis("Barcode rank", log_x_range),
        yaxis=log_axis("UMI count", log_y_range),
        margin=dict(l=60, r=20, t=40, b=90),
        height=440,
        plot_bgcolor="white",
        showlegend=True,
        legend=dict(orientation="h", yanchor="top", y=-0.20, xanchor="right", x=1.0),
    )
    return fig_to_dict(fig)


# 12-color categorical palette anchored on sea green (#2E8B57).
# Ordered so adjacent colors never collide: greens (anchor + analogous),
# warm accents (yellow/amber/orange/rose), cool tones (purple/indigo/blue),
# and a neutral cap. None of the colors collide with the warning red (#c0392b)
# used for threshold lines. Cycling only repeats above 12 categories.
PALETTE = [
    "#2E8B57",  # sea green (primary)
    "#1F6E44",  # forest
    "#3CB4A0",  # teal mint
    "#A4B92F",  # chartreuse
    "#F4D03F",  # yellow
    "#F39C12",  # amber
    "#E67E22",  # orange
    "#D45D79",  # rose
    "#9B59B6",  # purple
    "#5B5EA6",  # indigo
    "#5DADE2",  # sky blue
    "#7F8C8D",  # graphite
]
BAR_OUTLINE = dict(line=dict(color="black", width=1))


def build_histogram(values: np.ndarray, title: str, xlabel: str, log_x: bool = False, color: str = "#2E8B57") -> dict:
    if values.size == 0:
        return {"data": [], "layout": {"title": title}}
    if log_x:
        values = np.log10(np.maximum(values, 1))
        xlabel = f"log10({xlabel})"
    fig = go.Figure(data=[go.Histogram(
        x=values, nbinsx=60, marker_color=color,
        marker_line_color="black", marker_line_width=1,
    )])
    fig.update_layout(title=title, xaxis_title=xlabel, yaxis_title="cells",
                      margin=dict(l=50, r=20, t=40, b=40), height=320, bargap=0.05,
                      plot_bgcolor="white")
    fig.update_xaxes(showline=True, linecolor="black", linewidth=1, mirror=True, ticks="outside")
    fig.update_yaxes(showline=True, linecolor="black", linewidth=1, mirror=True, ticks="outside")
    return fig_to_dict(fig)


def build_mt_violin(pct_mt: np.ndarray) -> dict:
    if pct_mt.size == 0:
        return {"data": [], "layout": {"title": "% MT"}}
    fig = go.Figure(data=[go.Violin(
        y=pct_mt, box_visible=True, line_color="black",
        fillcolor="#2E8B57", opacity=0.7, name="% MT",
        marker=dict(line=dict(color="black", width=1)),
    )])
    fig.add_hline(y=10, line_dash="dash", line_color="#F39C12", annotation_text="10%")
    fig.add_hline(y=20, line_dash="dash", line_color="#c0392b", annotation_text="20%")
    fig.update_layout(title="Mitochondrial %", yaxis_title="% MT",
                      margin=dict(l=50, r=20, t=40, b=40), height=320, showlegend=False,
                      plot_bgcolor="white")
    fig.update_xaxes(showline=True, linecolor="black", linewidth=1, mirror=True)
    fig.update_yaxes(showline=True, linecolor="black", linewidth=1, mirror=True, ticks="outside")
    return fig_to_dict(fig)


def build_doublet_plot(scores: np.ndarray, threshold: float) -> dict:
    if scores.size == 0:
        return {"data": [], "layout": {"title": "Doublet scores"}}
    fig = go.Figure(data=[go.Histogram(
        x=scores, nbinsx=50, marker_color="#F39C12",
        marker_line_color="black", marker_line_width=1,
    )])
    if threshold is not None and not np.isnan(threshold):
        fig.add_vline(x=threshold, line_dash="dash", line_color="#c0392b",
                      annotation_text=f"threshold={threshold:.3f}")
    fig.update_layout(title="Scrublet doublet scores", xaxis_title="score",
                      yaxis_title="cells", margin=dict(l=50, r=20, t=40, b=40),
                      height=320, bargap=0.05, plot_bgcolor="white")
    fig.update_xaxes(showline=True, linecolor="black", linewidth=1, mirror=True, ticks="outside")
    fig.update_yaxes(showline=True, linecolor="black", linewidth=1, mirror=True, ticks="outside")
    return fig_to_dict(fig)


def load_qc_csv(path: str):
    umi, genes, mt, dscore = [], [], [], []
    with open(path) as fh:
        reader = csv.DictReader(fh)
        for row in reader:
            try:
                umi.append(int(row["umi"]))
                genes.append(int(row["genes"]))
                mt.append(float(row["pct_mt"]))
                if row.get("doublet_score"):
                    dscore.append(float(row["doublet_score"]))
            except (ValueError, KeyError):
                continue
    return np.array(umi), np.array(genes), np.array(mt), np.array(dscore)


def load_celltypist_counts(csv_path: str) -> "OrderedDict[str, int]":
    """Return cell counts per predicted label, sorted descending."""
    from collections import Counter, OrderedDict
    counter: Counter = Counter()
    with open(csv_path) as fh:
        reader = csv.DictReader(fh)
        if not reader.fieldnames:
            return OrderedDict()
        # Prefer majority_voting if present (cleaner), else predicted_labels
        label_col = "majority_voting" if "majority_voting" in reader.fieldnames else "predicted_labels"
        if label_col not in reader.fieldnames:
            return OrderedDict()
        for row in reader:
            label = (row.get(label_col) or "").strip()
            if label:
                counter[label] += 1
    return OrderedDict(counter.most_common())


def build_celltypist_plot(counts) -> dict:
    if not counts:
        return {"data": [], "layout": {"title": "No cell-type predictions"}}
    labels = list(counts.keys())
    values = list(counts.values())
    colors = [PALETTE[i % len(PALETTE)] for i in range(len(labels))]
    fig = go.Figure(data=[go.Bar(
        x=labels, y=values, marker_color=colors,
        marker_line_color="black", marker_line_width=1,
    )])
    fig.update_layout(
        title="Cell-type composition (CellTypist)",
        xaxis_title="predicted label",
        yaxis_title="cells",
        margin=dict(l=60, r=20, t=40, b=140),
        height=400,
        xaxis=dict(tickangle=-30, showline=True, linecolor="black", linewidth=1, mirror=True),
        yaxis=dict(showline=True, linecolor="black", linewidth=1, mirror=True, ticks="outside"),
        plot_bgcolor="white",
    )
    return fig_to_dict(fig)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--sample_id", required=True)
    ap.add_argument("--species", required=True)
    ap.add_argument("--fastp_json", required=True)
    ap.add_argument("--chemistry_json", required=True)
    ap.add_argument("--kb_dir", required=True)
    ap.add_argument("--cell_calling_json", required=True)
    ap.add_argument("--knee_data_json", required=True)
    ap.add_argument("--doublet_meta", required=True)
    ap.add_argument("--qc_csv", required=True)
    ap.add_argument("--qc_summary", required=True)
    ap.add_argument("--celltypist_csv", required=True)
    ap.add_argument("--celltypist_meta", required=True)
    ap.add_argument("--template", required=True)
    ap.add_argument("--out_html", required=True)
    ap.add_argument("--pipeline_version", default="0.1.0")
    ap.add_argument("--scq_workflow", default="quick", choices=["quick", "full"])
    ap.add_argument("--trace", default=None,
                    help="optional pipeline trace.txt to extract per-step times for this sample")
    args = ap.parse_args()

    with open(args.fastp_json) as fh: fastp_raw = json.load(fh)
    with open(args.chemistry_json) as fh: chemistry = json.load(fh)
    with open(args.cell_calling_json) as fh: cell_calling = json.load(fh)
    with open(args.knee_data_json) as fh: knee_data = json.load(fh)
    with open(args.doublet_meta) as fh: doublet = json.load(fh)
    with open(args.qc_summary) as fh: qc = json.load(fh)
    with open(args.celltypist_meta) as fh: celltypist_meta = json.load(fh)
    celltypist_counts = (
        load_celltypist_counts(args.celltypist_csv)
        if celltypist_meta.get("status") == "ok"
        else {}
    )

    run_info_path = os.path.join(args.kb_dir, "run_info.json")
    if os.path.exists(run_info_path):
        with open(run_info_path) as fh: run_info = json.load(fh)
    else:
        run_info = {}
    alignment = {
        "n_processed": run_info.get("n_processed", 0),
        "n_pseudoaligned": run_info.get("n_pseudoaligned", 0),
        "p_pseudoaligned": run_info.get("p_pseudoaligned", 0.0),
    }

    summary = fastp_raw.get("summary", {})
    before = summary.get("before_filtering", {})
    fastp_view = {
        "before_total_reads": before.get("total_reads", 0),
        "before_q30": before.get("q30_rate", 0.0),
        "duplication": fastp_raw.get("duplication", {}).get("rate", 0.0),
    }

    umi_arr, genes_arr, mt_arr, dscore_arr = load_qc_csv(args.qc_csv)

    knee_json = build_knee_plot(knee_data)
    umi_hist_json = build_histogram(umi_arr, "UMI per cell (log10)", "UMI", log_x=True)
    genes_hist_json = build_histogram(genes_arr, "Genes per cell", "genes", log_x=False)
    mt_violin_json = build_mt_violin(mt_arr)
    doublet_plot_json = build_doublet_plot(dscore_arr, doublet.get("threshold"))
    celltypist_plot_json = build_celltypist_plot(celltypist_counts)

    # Optional: per-sample step times harvested from Nextflow trace.txt
    step_times = []
    if args.trace and os.path.exists(args.trace):
        try:
            with open(args.trace) as fh:
                reader = csv.DictReader(fh, delimiter="\t")
                for r in reader:
                    name = r.get("name", "")
                    if args.sample_id in name or (
                        "(" not in name and r.get("status") == "COMPLETED"
                    ):
                        step_times.append({
                            "name": name,
                            "realtime": r.get("realtime", ""),
                            "cpus": r.get("cpus", ""),
                            "memory": r.get("memory", ""),
                        })
        except Exception as e:
            print(f"WARN: failed to parse trace {args.trace}: {e}", file=sys.stderr)

    template_dir = os.path.dirname(os.path.abspath(args.template))
    template_name = os.path.basename(args.template)
    env = Environment(loader=FileSystemLoader(template_dir),
                      autoescape=select_autoescape(["html"]))
    template = env.get_template(template_name)

    plotly_js = pyo.get_plotlyjs()
    html = template.render(
        sample_id=args.sample_id,
        species=args.species,
        generated_at=datetime.datetime.now().strftime("%Y-%m-%d %H:%M"),
        pipeline_version=args.pipeline_version,
        fastp=fastp_view,
        chemistry=chemistry,
        cell_calling=cell_calling,
        alignment=alignment,
        qc=qc,
        doublet=doublet,
        celltypist=celltypist_meta,
        celltypist_n_labels=len(celltypist_counts),
        celltypist_total_cells=sum(celltypist_counts.values()),
        scq_workflow=args.scq_workflow,
        step_times=step_times,
        plotly_js=plotly_js,
        knee_json=json.dumps(knee_json),
        umi_hist_json=json.dumps(umi_hist_json),
        genes_hist_json=json.dumps(genes_hist_json),
        mt_violin_json=json.dumps(mt_violin_json),
        doublet_plot_json=json.dumps(doublet_plot_json),
        celltypist_plot_json=json.dumps(celltypist_plot_json),
    )
    with open(args.out_html, "w") as fh:
        fh.write(html)
    print(f"sample_report {args.sample_id}: wrote {args.out_html}", file=sys.stderr)


if __name__ == "__main__":
    main()
