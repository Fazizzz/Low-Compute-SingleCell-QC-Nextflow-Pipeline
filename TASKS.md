# scq Build Task List

Compact phase-by-phase checklist. See `BUILD_PLAN.md` for full architecture.

## Phase 0 — Setup
- [x] `test_data/` contains SRR21064278 fastqs
- [x] `ref/` contains genome FASTA, GTF, prebuilt index, t2g, cdna
- [x] Legacy notebook + raw scripts deleted

## Phase 1 — Scaffold
- [x] Create directories: `workflows/ modules/local/ scripts/python/templates/ scripts/shell/ containers/sc_qc_base/ containers/sc_qc_annotate/ conf/ envs/ params/ test/ docs/`
- [x] `.gitignore` (test_data/, ref/, results/, work/, .nextflow/, *.html outputs)
- [x] `envs/scq.yml` — conda env definition
- [x] `nextflow.config` — params block, trace, resource labels
- [x] `nextflow_schema.json` — parameter schema
- [x] `params/default.yaml`
- [x] `conf/base.config`, `conf/docker.config`, `conf/singularity.config`
- [x] `main.nf` — validation, samplesheet parsing, channel wiring
- [x] `workflows/scq.nf` — module orchestration
- [x] `modules/local/validate_samplesheet.nf` (+ stub)
- [x] `test/samplesheet.csv`

## Phase 2 — FASTP
- [x] `modules/local/fastp.nf` (+ stub) — hardcoded `--disable_*` flags
- [x] `scripts/shell/01_fastp_qc.sh`
- [x] Wire into `workflows/scq.nf`

## Phase 3 — CHEMISTRY_CHECK
- [x] `modules/local/chemistry_check.nf` (+ stub)
- [x] `scripts/python/chemistry_check.py`
- [x] `scripts/shell/02_chemistry_check.sh`

## Phase 4 — KB_REF
- [x] `modules/local/kb_ref.nf` (+ stub, storeDir)
- [x] `scripts/shell/03_build_reference.sh`
- [x] Skip logic when `--prebuilt_index` set (in main.nf or workflows/scq.nf)

## Phase 5 — KB_COUNT
- [x] `modules/local/kb_count.nf` (+ stub)
- [x] `scripts/shell/04_kb_count.sh`

## Phase 6 — CELL_CALLING
- [x] `modules/local/cell_calling.nf` (+ stub)
- [x] `scripts/python/cell_calling.py` (knee + filter + MTX transpose + Cell Ranger format)
- [x] `scripts/shell/05_cell_calling.sh`

## Phase 7 — DOUBLET_DETECT
- [x] `modules/local/doublet_detect.nf` (+ stub)
- [x] `scripts/python/doublet_detect.py` (Scrublet + try/except + status JSON)
- [x] `scripts/shell/06_doublet_detect.sh`

## Phase 8 — QC_METRICS
- [x] `modules/local/qc_metrics.nf` (+ stub)
- [x] `scripts/python/qc_metrics.py` (MT% with species-aware prefix)
- [x] `scripts/shell/07_qc_metrics.sh`

## Phase 9 — SAMPLE_REPORT
- [x] `modules/local/sample_report.nf` (+ stub)
- [x] `scripts/python/sample_report.py`
- [x] `scripts/python/templates/sample_report.html.j2` (embedded Plotly)
- [x] `scripts/shell/08_sample_report.sh`

## Phase 10 — GROUP_REPORT
- [x] `modules/local/group_report.nf` (+ stub)
- [x] `scripts/python/group_report.py` (pass/fail, bar charts, trace parsing, warnings)
- [x] `scripts/python/templates/group_report.html.j2`
- [x] `scripts/shell/09_group_report.sh`

## Phase 11 — CELLTYPIST (optional, off by default)
- [x] `modules/local/celltypist_annotate.nf` (+ stub)
- [x] `scripts/python/celltypist_annotate.py` (try/except, status JSON)
- [x] `scripts/shell/10_celltypist_annotate.sh`
- [x] `containers/sc_qc_annotate/Dockerfile` (FROM sc_qc_base + celltypist)
- [x] Conditional wiring in workflows/scq.nf

## Phase 12 — Finalization
- [x] `modules/local/software_versions.nf`
- [x] `containers/sc_qc_base/Dockerfile` + `requirements.txt`
- [x] `docs/genome_download.md` (GRCh38 + GRCm39 wget commands)
- [x] `docs/celltypist_models.md`
- [x] `README.md` (overview, install, quick start, full usage, parameters table)
- [x] `CHANGELOG.md`

## Validation
- [x] `nextflow run main.nf -stub-run --samplesheet test/samplesheet.csv --prebuilt_index ref/index.idx --t2g ref/t2g.txt`
- [x] `nextflow run main.nf -profile local --samplesheet test/samplesheet.csv --prebuilt_index ref/index.idx --t2g ref/t2g.txt`
- [ ] Visual inspect `results/SRR21064278/SRR21064278_report.html`  ← **you**
- [ ] Visual inspect `results/group_report.html`                      ← **you**
- [x] Check `results/pipeline_info/trace.txt` exists
- [ ] Check `results/SRR21064278/filtered_matrix/` loads via scanpy/Seurat ← **you**

## Post-Approval
- [ ] `git init` in this folder
- [ ] Initial commit with full pipeline
- [ ] `gh repo create scq --public --source=. --push`
