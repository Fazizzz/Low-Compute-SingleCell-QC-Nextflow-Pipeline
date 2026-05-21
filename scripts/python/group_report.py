#!/usr/bin/env python3
"""Render group HTML report from collected per-sample JSONs + trace + software_versions."""
import argparse, csv, datetime, json, os, re, sys
import plotly.graph_objects as go
import plotly.offline as pyo
from jinja2 import Environment, FileSystemLoader, select_autoescape


CHART_SAMPLE_LIMIT = 20

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


def load_json(path: str) -> dict:
    with open(path) as fh:
        return json.load(fh)


def index_by_sample(paths, key: str = "sample_id") -> dict:
    out = {}
    for p in paths:
        d = load_json(p)
        if key in d:
            out[d[key]] = d
        else:
            out[os.path.basename(p)] = d
    return out


def derive_status(cells_called, p_aligned, doublet_rate, chemistry_status, doublet_status,
                  alignment_warn_threshold: float) -> tuple:
    warnings = []
    status = "PASS"
    if p_aligned is not None and p_aligned < alignment_warn_threshold:
        warnings.append(f"low alignment rate (<{alignment_warn_threshold:.0f}%)")
        status = "WARN"
    if cells_called is not None and cells_called < 100:
        warnings.append("fewer than 100 cells called")
        status = "WARN"
    if doublet_rate is not None and doublet_rate > 0.20:
        warnings.append("doublet rate above 20%")
        status = "WARN"
    if chemistry_status == "WARN":
        warnings.append("chemistry mismatch")
        status = "WARN"
    elif chemistry_status == "UNKNOWN":
        warnings.append("chemistry could not be detected")
    if doublet_status == "skipped":
        warnings.append("doublet detection skipped/failed")
    elif doublet_status == "unreliable":
        warnings.append("doublet rate flagged unreliable (implausibly high)")
    return status, warnings


def bar_chart(samples, key, title, ytitle, scale=1.0, hline=None, base_color="#2E8B57"):
    xs = [s["sample_id"] for s in samples]
    ys = [(s[key] * scale) if s.get(key) is not None else 0 for s in samples]
    colors = [PALETTE[i % len(PALETTE)] for i in range(len(xs))]
    fig = go.Figure(data=[go.Bar(
        x=xs, y=ys, marker_color=colors,
        marker_line_color="black", marker_line_width=1,
    )])
    if hline is not None:
        fig.add_hline(y=hline, line_dash="dash", line_color="#c0392b",
                      annotation_text=str(hline), annotation_position="top right")
    fig.update_layout(title=title, yaxis_title=ytitle,
                      margin=dict(l=50, r=20, t=40, b=60), height=340,
                      plot_bgcolor="white")
    fig.update_xaxes(showline=True, linecolor="black", linewidth=1, mirror=True,
                     tickangle=-30)
    fig.update_yaxes(showline=True, linecolor="black", linewidth=1, mirror=True,
                     ticks="outside")
    return json.loads(fig.to_json())


def parse_trace(trace_path: str):
    if not trace_path or not os.path.exists(trace_path):
        return []
    rows = []
    with open(trace_path) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for r in reader:
            rows.append({
                "name": r.get("name", ""),
                "realtime": r.get("realtime", ""),
                "cpu": r.get("%cpu", ""),
                "peak_rss": r.get("peak_rss", ""),
                "cpus": r.get("cpus", ""),
                "memory": r.get("memory", ""),
                "status": r.get("status", ""),
            })
    return rows


_DUR_UNITS = {"ms": 1e-3, "s": 1, "m": 60, "h": 3600, "d": 86400}


def parse_duration(s: str) -> float:
    """Parse a nextflow trace duration like '8m 20s' or '1h 30m 5s' or '300ms' into seconds."""
    if not s:
        return 0.0
    total = 0.0
    for num, unit in re.findall(r"(\d+(?:\.\d+)?)\s*(ms|s|m|h|d)", s):
        total += float(num) * _DUR_UNITS.get(unit, 0)
    return total


def parse_memory_gb(s: str) -> float:
    if not s:
        return 0.0
    m = re.search(r"(\d+(?:\.\d+)?)\s*([KMGT]?B?)", s, re.IGNORECASE)
    if not m:
        return 0.0
    val = float(m.group(1))
    unit = (m.group(2) or "").upper().rstrip("B")
    return {"": val / 1024**3, "K": val / 1024**2, "M": val / 1024,
            "G": val, "T": val * 1024}.get(unit, val)


def format_duration(seconds: float) -> str:
    if seconds < 60:
        return f"{seconds:.0f}s"
    if seconds < 3600:
        return f"{seconds/60:.1f}m"
    return f"{seconds/3600:.2f}h"


def compute_summary_from_trace(trace_rows: list) -> dict:
    if not trace_rows:
        return {}
    total_secs = 0.0
    cpu_hours = 0.0
    max_mem = 0.0
    n_tasks = 0
    for r in trace_rows:
        if r.get("status") != "COMPLETED":
            continue
        n_tasks += 1
        secs = parse_duration(r.get("realtime", ""))
        total_secs += secs
        try:
            cpus = int(r.get("cpus", "0") or 0)
        except ValueError:
            cpus = 0
        cpu_hours += (secs / 3600.0) * cpus
        mem_gb = parse_memory_gb(r.get("memory", ""))
        if mem_gb > max_mem:
            max_mem = mem_gb
    return {
        "n_tasks": n_tasks,
        "total_realtime": format_duration(total_secs),
        "cpu_hours": f"{cpu_hours:.2f}",
        "max_memory": f"{max_mem:.0f} GB" if max_mem >= 1 else f"{max_mem*1024:.0f} MB",
    }


def summary_stats(samples_view: list) -> list:
    """For compact mode: min/median/max for the key metrics."""
    import statistics

    def stats_for(key, scale=1.0, fmt="{:.0f}"):
        vals = [
            (s[key] * scale)
            for s in samples_view
            if s.get(key) is not None
        ]
        if not vals:
            return None
        return {
            "min": fmt.format(min(vals)),
            "median": fmt.format(statistics.median(vals)),
            "max": fmt.format(max(vals)),
        }

    rows = []
    for label, key, scale, fmt in [
        ("Cells called",       "cells_called",     1.0, "{:.0f}"),
        ("Median UMI",         "median_umi",       1.0, "{:.0f}"),
        ("Pseudoalignment (%)",    "p_pseudoaligned",  1.0, "{:.1f}%"),
        ("Doublet rate",       "doublet_rate",     100.0, "{:.2f}%"),
    ]:
        s = stats_for(key, scale, fmt)
        if s:
            rows.append({"metric": label, **s})
    return rows


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--chemistry_jsons", nargs="*", default=[])
    ap.add_argument("--cell_calling_jsons", nargs="*", default=[])
    ap.add_argument("--doublet_metas", nargs="*", default=[])
    ap.add_argument("--qc_summaries", nargs="*", default=[])
    ap.add_argument("--run_infos", nargs="*", default=[],
                    help="kb_count run_info.json files; filename prefix encodes sample_id")
    ap.add_argument("--trace", default=None)
    ap.add_argument("--software_versions", required=True)
    ap.add_argument("--template", required=True)
    ap.add_argument("--out_html", required=True)
    ap.add_argument("--scq_workflow", default="quick", choices=["quick", "full"])
    args = ap.parse_args()

    alignment_warn_threshold = 40.0 if args.scq_workflow == "quick" else 60.0

    chem  = index_by_sample(args.chemistry_jsons)
    cc    = index_by_sample(args.cell_calling_jsons)
    dbl   = index_by_sample(args.doublet_metas)
    qc    = index_by_sample(args.qc_summaries)

    run_infos = {}
    for p in args.run_infos:
        base = os.path.basename(p)
        if not base.endswith("_run_info.json"):
            print(f"WARN: run_info file {p!r} not named <sample>_run_info.json; skipping", file=sys.stderr)
            continue
        sid = base[: -len("_run_info.json")]
        run_infos[sid] = load_json(p)

    sample_ids = sorted(set(chem) | set(cc) | set(qc))
    samples_view = []
    warnings_all = []

    for sid in sample_ids:
        ch = chem.get(sid, {})
        c  = cc.get(sid, {})
        q  = qc.get(sid, {})
        d  = dbl.get(sid, {})
        ri = run_infos.get(sid) or {}
        if not ri:
            print(f"WARN: no run_info for {sid}", file=sys.stderr)

        cells = c.get("cells_called")
        median_umi = q.get("median_umi", 0)
        p_aligned = ri.get("p_pseudoaligned", 0.0)
        doublet_rate = d.get("doublet_rate")
        chemistry_status = ch.get("status", "UNKNOWN")
        doublet_status = d.get("status", "ok")
        status, warns = derive_status(
            cells, p_aligned, doublet_rate, chemistry_status, doublet_status,
            alignment_warn_threshold,
        )
        for w in warns:
            warnings_all.append(f"{sid}: {w}")

        samples_view.append({
            "sample_id": sid,
            "status": status,
            "cells_called": cells,
            "median_umi": median_umi,
            "p_pseudoaligned": p_aligned or 0.0,
            "doublet_rate": doublet_rate,
            "chemistry_status": chemistry_status,
        })

    n_samples = len(samples_view)
    compact_mode = n_samples > CHART_SAMPLE_LIMIT

    if compact_mode:
        bar_charts = {}
        stats = summary_stats(samples_view)
    else:
        bar_charts = {
            "cells":  bar_chart(samples_view, "cells_called",    "Cells called", "cells"),
            "umi":    bar_chart(samples_view, "median_umi",      "Median UMI per cell", "UMI"),
            "aligned": bar_chart(samples_view, "p_pseudoaligned", "Pseudoalignment (%)", "%",
                                 hline=alignment_warn_threshold),
            "doublet": bar_chart(
                [{**s, "_dr": (s["doublet_rate"] * 100 if s["doublet_rate"] is not None else 0)}
                 for s in samples_view],
                "_dr", "Doublet rate", "%",
            ),
        }
        stats = []

    trace_rows = parse_trace(args.trace)
    compute_summary = compute_summary_from_trace(trace_rows)

    with open(args.software_versions) as fh:
        software_versions = fh.read()

    template_dir = os.path.dirname(os.path.abspath(args.template))
    env = Environment(loader=FileSystemLoader(template_dir),
                      autoescape=select_autoescape(["html"]))
    template = env.get_template(os.path.basename(args.template))
    html = template.render(
        n_samples=n_samples,
        samples=samples_view,
        warnings=warnings_all,
        compute_summary=compute_summary,
        compact_mode=compact_mode,
        chart_sample_limit=CHART_SAMPLE_LIMIT,
        summary_stats=stats,
        software_versions=software_versions,
        scq_workflow=args.scq_workflow,
        generated_at=datetime.datetime.now().strftime("%Y-%m-%d %H:%M"),
        plotly_js=pyo.get_plotlyjs(),
        bar_charts_json=json.dumps(bar_charts),
    )
    with open(args.out_html, "w") as fh:
        fh.write(html)
    print(f"group_report: wrote {args.out_html} ({n_samples} samples, "
          f"{len(warnings_all)} warnings, compact={compact_mode})", file=sys.stderr)


if __name__ == "__main__":
    main()
