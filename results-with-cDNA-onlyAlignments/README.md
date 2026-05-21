# Quick mode demo results

These are the `--workflow quick` outputs from the four-sample validation set:
SRR21064278, SRR21064279, SRR33398036, SRR33398039 (see the project README for the dataset overview).

## What is preserved here

- **Per-sample HTML reports** (`<sample>/<sample>_report.html`) — the primary deliverable.
- **Group report** (`group_report.html`).
- **JSON / CSV metadata** — cell calling, doublet detection, QC summary, CellTypist predictions, chemistry detection, fastp.
- **One filtered_matrix as a layout example** — `SRR21064278/filtered_matrix/` (14 MB). Contains the Cell Ranger style `matrix.mtx.gz`, `barcodes.tsv.gz`, `features.tsv.gz` you would feed straight into scanpy or Seurat.
- **Pipeline diagnostics** — `pipeline_info/trace.txt`, `execution_report.html`, `software_versions.yml`.

## What is NOT preserved (regenerable from a re-run)

- `kb_count/kb_out/` — the unfiltered BUS file, equivalence classes, and unfiltered MTX matrices (these are 2-13 GB per sample).
- `filtered_matrix/` for the other three samples (each is 52-70 MB; same layout as the preserved SRR21064278 example).

## How to regenerate everything

```bash
nextflow run main.nf -profile local --workflow quick \
    --samplesheet your_samples.csv \
    --prebuilt_index your/cdna.idx --t2g your/t2g.txt
```

Quick mode runs in 10-30 min per sample on a laptop and reproduces every artefact in this directory.
