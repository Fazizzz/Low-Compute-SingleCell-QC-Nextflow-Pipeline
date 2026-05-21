# scq — Single Cell Quick QC Pipeline: Build Plan

A Nextflow DSL2 pipeline for fast, low-compute QC of single-cell RNA-seq samples
using kallisto|bustools. Generates Cell Ranger-style per-sample HTML reports and
a multi-sample group summary report. Targets local compute environments and is
easy to install via conda.

---

## 1. Project Context

### Why this pipeline

Cell Ranger is the de facto QC tool for scRNA-seq but requires a 10x license, is
compute-heavy (16+ cores, 64GB RAM), and takes 2–6 hours per sample. Kallisto|bustools
runs on a laptop in ~10 minutes with 8GB RAM. This pipeline fills a real gap:
fast, open-source, license-free single cell QC that produces Cell Ranger-style
outputs (reports, plots, filtered MTX matrices) without needing Cell Ranger.

### Test data already on disk

- `test_data/SRR21064278_R1.fastq.gz` — 10x Chromium v3 chemistry
- `test_data/SRR21064278_R2.fastq.gz`
- `ref/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz` — genome FASTA
- `ref/Homo_sapiens.GRCh38.115.gtf.gz` — annotation GTF
- `ref/index.idx` — pre-built kallisto index
- `ref/t2g.txt` — transcript-to-gene map
- `ref/cdna.fa` — cDNA FASTA from kb ref

### Reference environment

User's `SingleCell` conda environment already has: kallisto 0.51.1, bustools 0.44.1,
kb-python 0.29.5, nextflow 24.10.6, scipy 1.15.3, numpy 2.2.6, pandas 2.3.3,
plotly 6.7.0, jinja2 3.1.6, matplotlib 3.8.4, scikit-learn 1.7.2, scanpy 1.11.5,
anndata 0.11.4, umap-learn 0.5.12.

**Missing packages** (must be added to pipeline env): `fastp`, `scrublet`,
`celltypist`, `kneed`.

---

## 2. Final Architecture Decisions

| Decision | Choice |
|---|---|
| Workflow language | Nextflow DSL2 |
| Containers | Docker + Singularity; single base `sc_qc_base` + extending `sc_qc_annotate` |
| fastp role | QC-only, NO trimming/filtering (preserves barcode+UMI in R1) |
| MultiQC | Not used — fastp HTML + custom group report covers everything |
| Reference handling | Path-only; user provides FASTA+GTF or pre-built index. No downloads. |
| Multi-species | Not supported. `species` column is informational (for MT prefix). |
| RNA velocity | Removed — would slow runs significantly |
| Cell calling | Auto-detected knee shown as annotation line only. Hard threshold (`--min_umi_threshold`, default 200) does actual filtering. |
| Sample filtering | All samples processed regardless of cell count. Failures + warnings collected in group report. |
| Doublet detection | Scrublet, graceful fail on <100 cells |
| Cell type annotation | CellTypist (optional, requires user-provided model path) |
| Per-sample reports | Full detail, self-contained interactive HTML |
| Group report | Slim: pass/fail table, bar charts, warnings, resource usage from trace |
| MTX output format | Cell Ranger compatible (barcodes.tsv.gz, features.tsv.gz, matrix.mtx.gz) |
| Stub testing | Every module has a `stub:` block for CI |
| Shell scripts | Standalone `scripts/shell/*.sh` wrappers for running each step outside the pipeline |
| Software versions | Captured in `pipeline_info/software_versions.yml` and embedded in group report |
| Cost metrics | Resource usage (CPU-hours, wall time, peak RAM) only — no pricing assumptions |

---

## 3. Pipeline Flow

```
Samplesheet (CSV: sample_id, fastq_r1, fastq_r2, chemistry, species)
    │
    ├─► VALIDATE_SAMPLESHEET     fail-fast: columns, files, chemistry enum
    │
    ├─► FASTP                    QC-only, no trimming
    │                            outputs: {sample}_fastp.json, .html
    │
    ├─► CHEMISTRY_CHECK          sample 10K reads from R1, measure length
    │                            outputs: {sample}_chemistry.json
    │
    ├─► KB_REF                   build kallisto index (storeDir cached)
    │                            skipped if --prebuilt_index provided
    │
    ├─► KB_COUNT                 pseudoalignment → unfiltered count matrix
    │
    ├─► CELL_CALLING             1. auto-detect knee (KneeLocator) → for display
    │                            2. apply --min_umi_threshold for filter
    │                            3. write Cell Ranger MTX
    │                            outputs: filtered_matrix/, cell_calling.json
    │
    ├─► DOUBLET_DETECT           Scrublet, random_state=42
    │                            graceful fail if n_cells < 100
    │
    ├─► QC_METRICS               per-cell: UMI, genes, MT%, doublet score
    │
    ├─► CELLTYPIST (optional)    requires --run_celltypist + --celltypist_model_path
    │                            graceful fail if n_cells < 50
    │
    ├─► SAMPLE_REPORT            per-sample interactive HTML (full detail)
    │
    └─► GROUP_REPORT             collect all samples:
                                 - pass/fail table
                                 - bar charts
                                 - warnings panel
                                 - resource usage from trace.txt
                                 - software versions
```

---

## 4. Samplesheet Format

```csv
sample_id,fastq_r1,fastq_r2,chemistry,species
PBMC_01,/data/PBMC_01_R1.fastq.gz,/data/PBMC_01_R2.fastq.gz,10xv3,human
```

- `chemistry`: `10xv2`, `10xv3`, `10xv3.1`, `dropseq`
- `species`: `human`, `mouse` (used for MT prefix: `MT-` vs `mt-`)

---

## 5. Pipeline Parameters

```
Required (one mode):
  --genome_fasta      Path to genome FASTA (gzipped)
  --genome_gtf        Path to GTF annotation
OR:
  --prebuilt_index    Path to pre-built kallisto index (.idx)
  --t2g               Path to transcript-to-gene file (required with prebuilt_index)

Always:
  --samplesheet       Path to input CSV

Optional:
  --outdir            Output directory (default: results)
  --min_umi_threshold Hard UMI filter (default: 200)
  --threads           CPUs per kb count (default: 4)
  --ref_cache_dir     Where KB_REF stores cached index (default: ${HOME}/.scq_cache)
  --run_celltypist    Enable annotation (default: false)
  --celltypist_model_path  Path to local .pkl model file (required if run_celltypist)
```

Validation enforces: not both genome_fasta and prebuilt_index; t2g required with
prebuilt_index; all paths must exist; chemistry+species enum check on each row.

---

## 6. Project File Structure

```
scq/
├── main.nf
├── nextflow.config
├── nextflow_schema.json
├── CHANGELOG.md
├── README.md
├── .gitignore
│
├── params/
│   └── default.yaml
│
├── workflows/
│   └── scq.nf
│
├── modules/
│   └── local/
│       ├── validate_samplesheet.nf
│       ├── fastp.nf
│       ├── chemistry_check.nf
│       ├── kb_ref.nf
│       ├── kb_count.nf
│       ├── cell_calling.nf
│       ├── doublet_detect.nf
│       ├── qc_metrics.nf
│       ├── celltypist_annotate.nf
│       ├── sample_report.nf
│       ├── group_report.nf
│       └── software_versions.nf
│
├── scripts/
│   ├── python/
│   │   ├── chemistry_check.py
│   │   ├── cell_calling.py
│   │   ├── doublet_detect.py
│   │   ├── qc_metrics.py
│   │   ├── sample_report.py
│   │   ├── group_report.py
│   │   ├── celltypist_annotate.py
│   │   └── templates/
│   │       ├── sample_report.html.j2
│   │       └── group_report.html.j2
│   └── shell/
│       ├── 01_fastp_qc.sh
│       ├── 02_chemistry_check.sh
│       ├── 03_build_reference.sh
│       ├── 04_kb_count.sh
│       ├── 05_cell_calling.sh
│       ├── 06_doublet_detect.sh
│       ├── 07_qc_metrics.sh
│       ├── 08_sample_report.sh
│       ├── 09_group_report.sh
│       └── 10_celltypist_annotate.sh
│
├── containers/
│   ├── sc_qc_base/
│   │   ├── Dockerfile
│   │   └── requirements.txt
│   └── sc_qc_annotate/
│       ├── Dockerfile
│       └── requirements.txt
│
├── conf/
│   ├── base.config
│   ├── docker.config
│   └── singularity.config
│
├── envs/
│   └── scq.yml
│
├── test/
│   ├── samplesheet.csv          # points to test_data/SRR21064278
│   └── stub_samplesheet.csv
│
├── docs/
│   ├── genome_download.md       # wget commands for GRCh38 + GRCm39
│   └── celltypist_models.md     # how to list + download models
│
├── test_data/                   # (gitignored) SRR21064278 fastqs
├── ref/                         # (gitignored) genome FASTA, GTF, index, t2g
└── results/                     # (gitignored) pipeline outputs
```

---

## 7. Module Specifications

### VALIDATE_SAMPLESHEET
- Python script reads CSV, checks headers, validates file paths, chemistry enum, species enum
- Fails fast with per-row errors
- Stub: emit dummy channel

### FASTP
- `fastp --in1 R1 --in2 R2 --disable_adapter_trimming --disable_quality_filtering --disable_length_filtering --thread N --json out.json --html out.html`
- Disable flags are HARDCODED — must not be user-configurable (preserves barcode+UMI integrity)
- Output: `{sample}_fastp.json`, `{sample}_fastp.html`
- Resource label: `process_low`

### CHEMISTRY_CHECK
- Python: stream first 10K reads from R1.fastq.gz, measure length distribution
- Map length to chemistry class:
  - 26bp → `10xv2`
  - 28bp → `10xv3-class` (cannot distinguish v3/v3.1/Multiome)
  - 20bp → `dropseq`
  - other → `unknown`
- Compare to samplesheet-specified chemistry → PASS / WARN / UNKNOWN
- Output: `{sample}_chemistry.json` (advisory only, doesn't block pipeline)
- Resource label: `process_low`

### KB_REF
- `kb ref -i index.idx -g t2g.txt -f1 cdna.fa FASTA GTF`
- `storeDir` keyed on FASTA filename + GTF filename
- Skipped entirely if `--prebuilt_index` provided
- Output: `index.idx`, `t2g.txt`, `cdna.fa`
- Resource label: `process_high` (memory-intensive)

### KB_COUNT
- `kb count -i index.idx -g t2g.txt -x {chemistry} -t {threads} -m {memory} -o out R1 R2`
- Output: `counts_unfiltered/`, `run_info.json`, `kb_info.json`, `inspect.json`
- Resource label: `process_high`

### CELL_CALLING
- Python script:
  1. Load `counts_unfiltered/cells_x_genes.mtx`
  2. Compute UMI per barcode, sort descending
  3. Run `KneeLocator(ranks, log10(umis), curve='convex', direction='decreasing')` on log-log curve → record auto-detected UMI value
  4. Filter at `--min_umi_threshold` (actual filter)
  5. Transpose matrix (kb outputs cells×genes; Cell Ranger expects genes×cells)
  6. Write Cell Ranger MTX: `barcodes.tsv.gz`, `features.tsv.gz`, `matrix.mtx.gz`
  7. Write `cell_calling.json` with: auto_knee_umi, auto_knee_rank, threshold_used, cells_called, knee_status (reliable/unreliable based on 5th-95th percentile sanity check)
  8. Write `knee_data.json` with rank curve data for interactive Plotly rendering in report
- Resource label: `process_medium`

### DOUBLET_DETECT
- Python script using Scrublet
- `random_state=42` for reproducibility
- Wrap in try/except: if `n_cells < 100` or any exception → status JSON with `{"status": "skipped", "reason": "..."}`
- Output: `doublet_scores.csv` (barcode, score, predicted_doublet), `doublet_meta.json`
- Resource label: `process_medium`

### QC_METRICS
- Load filtered MTX + doublet scores
- Per-cell: total UMI, genes detected, MT%, doublet score
- MT prefix: `MT-` (human), `mt-` (mouse) from species column
- Output: `qc_metrics.csv`, `qc_summary.json` (median/mean/p25/p75)
- Guard: if MT% is 0% for all cells AND species is human/mouse, flag warning (likely MT detection failed)
- Resource label: `process_low`

### CELLTYPIST (optional)
- Python: load filtered MTX as AnnData, run `celltypist.annotate(model=path, majority_voting=True)`
- Wrap in try/except: if `n_cells < 50` → skip with status JSON
- Output: `celltypist_predictions.csv`, `celltypist_meta.json`
- Uses `sc_qc_annotate` container (extends `sc_qc_base` with celltypist+anndata)
- Resource label: `process_medium`

### SAMPLE_REPORT
- Jinja2 + Plotly (embedded, self-contained HTML)
- Sections:
  - Header: sample ID, date, pipeline version, software versions
  - Read QC: fastp Q30%, mean quality, duplication
  - Chemistry: detected vs specified, badge
  - Alignment: % pseudoaligned, total reads, runtime
  - Barcode rank plot: log-log interactive Plotly with grey=all, blue=cells, red dashed=threshold, orange annotated=auto-detected knee
  - UMI distribution histogram
  - Genes per cell histogram
  - MT% violin with 10%/20% reference lines
  - Doublet score histogram with call threshold (or skip-reason box)
  - Cell type composition bar chart (only if CellTypist ran)
  - QC stats table
- Output: `{sample}_report.html`
- Resource label: `process_low`

### GROUP_REPORT
- Runs last, collects all per-sample JSONs + `trace.txt` + `software_versions.yml`
- Sections:
  - Pass/Fail table: per-sample status with color coding (green/amber/red)
    - Warning triggers: alignment <60%, cells <100, doublet rate >20%, chemistry WARN, Scrublet failed
  - Bar charts (Plotly): cells, median UMI, % aligned (60% line), doublet rate
  - Warnings panel: explicit bulleted list of all warnings across samples
  - Resource usage table: per-step wall time, CPU-hours, peak RAM from trace
  - Software versions section
  - Links to per-sample reports
- Gracefully degrades if trace.txt unavailable
- Output: `group_report.html`
- Resource label: `process_low`

### SOFTWARE_VERSIONS
- Runs in `sc_qc_base` container
- Captures: kallisto, bustools, kb, fastp, scrublet, celltypist (if available), kneed, scanpy, scipy, numpy versions
- Output: `pipeline_info/software_versions.yml`
- Resource label: `process_low`

---

## 8. Conda Environment (`envs/scq.yml`)

```yaml
name: scq
channels:
  - bioconda
  - conda-forge
  - defaults
dependencies:
  - python=3.10
  - kallisto=0.51.1
  - bustools=0.44.1
  - fastp>=1.3
  - nextflow>=24.10.6
  - scipy>=1.15
  - numpy>=2.2
  - pandas>=2.3
  - matplotlib>=3.8
  - plotly>=6.7
  - jinja2>=3.1
  - scikit-learn>=1.7
  - scanpy>=1.11
  - anndata>=0.11
  - pip
  - pip:
      - kb-python==0.29.5
      - scrublet>=0.2.3
      - celltypist>=1.6
      - kneed>=0.8
```

---

## 9. Containers

| Container | Base | Adds |
|---|---|---|
| `sc_qc_base` | `mambaorg/micromamba:1.5` | Full conda env from `scq.yml` minus celltypist |
| `sc_qc_annotate` | `sc_qc_base` | `pip install celltypist anndata` |

Dockerfile pattern for `sc_qc_base`:
```dockerfile
FROM mambaorg/micromamba:1.5
COPY envs/scq.yml /tmp/scq.yml
RUN micromamba install -y -n base -f /tmp/scq.yml && micromamba clean --all --yes
ENV PATH=/opt/conda/bin:$PATH
```

---

## 10. Resource Profile (`conf/base.config`)

```groovy
process {
    errorStrategy = 'terminate'
    maxRetries    = 1

    withLabel: process_low {
        cpus   = 1
        memory = '2 GB'
        time   = '30m'
    }
    withLabel: process_medium {
        cpus   = 4
        memory = '8 GB'
        time   = '2h'
    }
    withLabel: process_high {
        cpus   = 8
        memory = '16 GB'
        time   = '4h'
    }
}

trace {
    enabled = true
    file = "${params.outdir}/pipeline_info/trace.txt"
    fields = 'task_id,name,status,exit,realtime,%cpu,peak_rss,peak_vmem,cpus,memory'
}
```

---

## 11. Build Task List (Execute in Order)

### Phase 0 — Setup (DONE)
- [x] Move `test_data/` and `ref/` into subdirs
- [x] Delete legacy notebook + scripts

### Phase 1 — Scaffold
- [ ] Create directory structure (workflows/, modules/local/, scripts/python/, scripts/shell/, containers/, conf/, envs/, params/, test/, docs/)
- [ ] `main.nf` — entry point with parameter validation, samplesheet parsing, workflow invocation
- [ ] `nextflow.config` — params block, process resource labels, trace config, profile placeholders
- [ ] `nextflow_schema.json` — parameter schema for `--help`
- [ ] `params/default.yaml`
- [ ] `conf/base.config`, `conf/docker.config`, `conf/singularity.config`
- [ ] `envs/scq.yml`
- [ ] `workflows/scq.nf` — main workflow with channel wiring (modules empty for now)
- [ ] `modules/local/validate_samplesheet.nf` + stub
- [ ] `.gitignore` (test_data/, ref/, results/, work/, .nextflow/)
- [ ] `test/samplesheet.csv` pointing to test_data/SRR21064278

### Phase 2 — FASTP
- [ ] `modules/local/fastp.nf` with stub
- [ ] `scripts/shell/01_fastp_qc.sh`
- [ ] Wire into `workflows/scq.nf`

### Phase 3 — CHEMISTRY_CHECK
- [ ] `modules/local/chemistry_check.nf` with stub
- [ ] `scripts/python/chemistry_check.py`
- [ ] `scripts/shell/02_chemistry_check.sh`
- [ ] Wire into workflow

### Phase 4 — KB_REF
- [ ] `modules/local/kb_ref.nf` with storeDir, stub
- [ ] `scripts/shell/03_build_reference.sh`
- [ ] Conditional skip logic when `--prebuilt_index` provided

### Phase 5 — KB_COUNT
- [ ] `modules/local/kb_count.nf` with stub
- [ ] `scripts/shell/04_kb_count.sh`

### Phase 6 — CELL_CALLING
- [ ] `modules/local/cell_calling.nf` with stub
- [ ] `scripts/python/cell_calling.py` (knee detection + filter + Cell Ranger MTX output + knee data JSON)
- [ ] `scripts/shell/05_cell_calling.sh`

### Phase 7 — DOUBLET_DETECT
- [ ] `modules/local/doublet_detect.nf` with stub
- [ ] `scripts/python/doublet_detect.py` (Scrublet + graceful fail)
- [ ] `scripts/shell/06_doublet_detect.sh`

### Phase 8 — QC_METRICS
- [ ] `modules/local/qc_metrics.nf` with stub
- [ ] `scripts/python/qc_metrics.py`
- [ ] `scripts/shell/07_qc_metrics.sh`

### Phase 9 — SAMPLE_REPORT
- [ ] `modules/local/sample_report.nf` with stub
- [ ] `scripts/python/sample_report.py`
- [ ] `scripts/python/templates/sample_report.html.j2`
- [ ] `scripts/shell/08_sample_report.sh`

### Phase 10 — GROUP_REPORT
- [ ] `modules/local/group_report.nf` with stub
- [ ] `scripts/python/group_report.py` (collect, charts, trace parsing, warnings)
- [ ] `scripts/python/templates/group_report.html.j2`
- [ ] `scripts/shell/09_group_report.sh`

### Phase 11 — CELLTYPIST (optional, wired but disabled by default)
- [ ] `modules/local/celltypist_annotate.nf` with stub
- [ ] `scripts/python/celltypist_annotate.py`
- [ ] `scripts/shell/10_celltypist_annotate.sh`
- [ ] `containers/sc_qc_annotate/Dockerfile`
- [ ] Conditional wiring (only runs if `--run_celltypist` true and model path provided)

### Phase 12 — Finalization
- [ ] `modules/local/software_versions.nf`
- [ ] `containers/sc_qc_base/Dockerfile` + `requirements.txt`
- [ ] `docs/genome_download.md` (Ensembl wget commands for GRCh38, GRCm39)
- [ ] `docs/celltypist_models.md` (how to browse and download models)
- [ ] `README.md` (quick start, full usage, requirements, examples)
- [ ] `CHANGELOG.md`

### Validation
- [ ] Run stub: `nextflow run main.nf -stub-run --samplesheet test/samplesheet.csv --prebuilt_index ref/index.idx --t2g ref/t2g.txt`
- [ ] Run full: `nextflow run main.nf -profile local --samplesheet test/samplesheet.csv --prebuilt_index ref/index.idx --t2g ref/t2g.txt`
- [ ] Inspect `results/{sample_id}/{sample_id}_report.html`
- [ ] Inspect `results/group_report.html`
- [ ] User OK → init new git repo, push to GitHub

---

## 12. Issues to Watch For

1. **fastp R1 trim flags** — `--disable_adapter_trimming --disable_quality_filtering --disable_length_filtering` are MANDATORY. Hardcode them; do not expose as user params.

2. **MTX orientation** — kb count outputs `cells_x_genes.mtx` (cells as rows). Cell Ranger MTX format is genes × barcodes. The CELL_CALLING Python script MUST transpose with `.T` before writing.

3. **storeDir for KB_REF** — must be an absolute path outside `work/` (which `nextflow clean` wipes). Default: `${HOME}/.scq_cache`. User-configurable via `--ref_cache_dir`. Key on FASTA+GTF filename+size.

4. **Scrublet failure modes** — fails on <100 cells (numpy decomposition) and on bad bimodal distributions. Wrap in try/except, emit status JSON, surface in group report.

5. **CellTypist on low-cell samples** — same: wrap and skip if <50 cells.

6. **KneeLocator on flat curves** — can return extreme values. Sanity check: if auto-detected UMI is below 5th or above 95th percentile of barcode UMI distribution, mark as "unreliable" in the report.

7. **MT prefix detection** — `MT-` (human), `mt-` (mouse). If MT% is 0% for all cells, warn that MT detection may have failed.

8. **trace.txt may be missing** — group report must check existence before parsing, degrade gracefully.

9. **Single-sample group report** — must work with N=1 sample (single-bar charts, one-row table).

10. **Parameter conflict** — main.nf must error if both `--genome_fasta` and `--prebuilt_index` are provided.

11. **Plotly self-contained HTML size** — embedded JS makes each report ~3MB but ensures portability. Acceptable for local use with handful of samples.

12. **Chemistry check is advisory** — never blocks the pipeline. User-specified chemistry always wins. Mismatches surface as warnings.

---

## 13. Reference Commands (After Build)

### Stub test (~10 seconds)
```bash
cd /Users/m.faizankh/Coding/GitHub/Low-compute-SingleCell-QC-Nextflow-Pipeline
nextflow run main.nf -stub-run \
  --samplesheet test/samplesheet.csv \
  --prebuilt_index ref/index.idx \
  --t2g ref/t2g.txt
```

### Full test on SRR21064278 (~10-15 minutes)
```bash
nextflow run main.nf -profile local \
  --samplesheet test/samplesheet.csv \
  --prebuilt_index ref/index.idx \
  --t2g ref/t2g.txt \
  --outdir results
```

### Outputs to inspect
- `results/SRR21064278/SRR21064278_report.html` — per-sample report
- `results/group_report.html` — group summary
- `results/SRR21064278/filtered_matrix/` — Cell Ranger MTX
- `results/pipeline_info/trace.txt` — resource usage
- `results/pipeline_info/software_versions.yml` — tool versions

### After user approval
```bash
git init
git add .
git commit -m "Initial commit: scq single cell QC pipeline"
gh repo create scq --public --source=. --push
```

---

## 14. Test Samplesheet Content

`test/samplesheet.csv`:
```csv
sample_id,fastq_r1,fastq_r2,chemistry,species
SRR21064278,test_data/SRR21064278_R1.fastq.gz,test_data/SRR21064278_R2.fastq.gz,10xv3,human
```
