# AGENTS.md — operating guide for AI agents running this pipeline

This file orients an autonomous or semi-autonomous agent that has been asked to run, troubleshoot, or extend the **SingleCell Quick QC (scq) pipeline**. Read it once before acting.

## 1. What this pipeline does (1-paragraph mental model)

`scq` is a Nextflow DSL2 pipeline that takes paired-end single-cell RNA-seq fastqs (10x Chromium 3' v3, Element AVITI, or other supported chemistries — chemistry is auto-detected), pseudoaligns them with **kallisto | bustools** via `kb-python`, calls cells from the barcode-rank knee, runs **Scrublet** doublet detection, computes per-cell QC metrics, optionally annotates with **CellTypist**, and emits per-sample + group HTML reports. There are two operating modes:

- `--workflow quick` — cDNA-only index. Small (≈1 GB), runs on a laptop. Captures ~50-60% of reads because intronic reads are missed. Use for **rapid integrity checks and compute-cost prediction**.
- `--workflow full` — `nac` index (spliced + nascent). Large (~7 GB for human). Captures ~80-85% of reads. Produces 10x-Cell-Ranger-comparable count matrices suitable for downstream Seurat / scanpy. Requires HPC, cloud, or workstation with ≥32 GB RAM.

The pipeline writes a Cell-Ranger-style count matrix per sample and an HTML report that is the **primary deliverable** for end users.

## 2. Pick a profile before invoking

Profiles in `nextflow.config`:

| Compute profile | Target host | process_high envelope | When to use |
|---|---|---|---|
| `local` (default) | Laptop / desktop, ≤16 GB RAM | 4 cpu / 14 GB / 4h | Quick mode only. Full mode here OOMs the nac index. |
| `hpc` | HPC node or workstation, 32+ GB RAM, 16+ cores | 16 cpu / 64 GB / 6h | Full mode; allows parallel KB_COUNT if the scheduler permits. |
| `aws` | EC2 m6i.4xlarge (16 vCPU / 64 GB RAM, 300 GB gp3) | 8 cpu / 48 GB / 6h | Validated reference instance; forces serial KB_COUNT to avoid co-tenant OOM. |

Container profiles (combine with one compute profile using a comma):
- `docker` — uses local docker daemon and the `sc_qc_base` / `sc_qc_annotate` images.
- `singularity` — same containers as `.sif`.

**Examples:**

```bash
# Laptop, quick mode, conda env on PATH
nextflow run main.nf -profile local \
    --samplesheet test/samplesheet_single.csv \
    --prebuilt_index path/to/cdna.idx --t2g path/to/t2g.txt

# HPC, full mode with singularity
nextflow run main.nf -profile hpc,singularity \
    --samplesheet samples.csv --workflow full \
    --prebuilt_index nac.idx --t2g t2g.txt \
    --cdna_t2c cdna_t2c.txt --nascent_t2c nascent_t2c.txt

# AWS docker
nextflow run main.nf -profile aws,docker \
    --samplesheet samples.csv --workflow full \
    --prebuilt_index nac.idx --t2g t2g.txt \
    --cdna_t2c cdna_t2c.txt --nascent_t2c nascent_t2c.txt \
    --run_celltypist --celltypist_model_path Immune_All_Low.pkl
```

## 3. Reference index — user's responsibility

**The pipeline does not build the kallisto index for you.** A reference build is a one-shot, can take 30-60 min on 16 cores, and is identical across runs against the same genome — keeping it outside the pipeline lets users version-control and share a single index. Build it with:

```bash
# Quick (cDNA-only)
scripts/shell/03_build_reference.sh genome.fa.gz genome.gtf.gz out/

# Full (nac, spliced + nascent) — pass --full to the script
scripts/shell/03_build_reference.sh genome.fa.gz genome.gtf.gz out/ --full

# Memory-constrained host: build a protein-coding-only nac index
scripts/shell/build_pc_nac_index.sh genome.fa.gz genome.gtf.gz out/
```

**CRITICAL — kallisto version pinning:** The index MUST be built with the same kallisto major version as KB_COUNT will use. The container ships kallisto 0.51.1. kb-python 0.29.5 also ships a bundled 0.52.0 binary at `/opt/conda/lib/python3.10/site-packages/kb_python/bins/linux/kallisto/kallisto`. The pipeline forces `--kallisto $(command -v kallisto)` (= conda 0.51.1) everywhere, but if you build the index outside the container without that flag, kb-python will silently pick the bundled binary and produce a 0.52.0-format index that **crashes kallisto 0.51.1 with `std::bad_array_new_length` at read time**. Always pass `--kallisto /opt/conda/bin/kallisto --bustools /opt/conda/bin/bustools` if invoking `kb ref` inside the container manually.

## 4. Standalone shell scripts

Each pipeline stage has a `scripts/shell/0N_*.sh` wrapper that runs that stage outside Nextflow. Useful for debugging, custom reruns, or pipelining into a different orchestrator. They take CLI args (no AWS/profile assumptions) and find python helpers via `$script_dir/../python/`.

| Stage | Script | Outputs |
|---|---|---|
| 1 | `01_fastp_qc.sh` | `*_fastp.{json,html}` |
| 2 | `02_chemistry_check.sh` | `*_chemistry.json` |
| 3 | `03_build_reference.sh` | `index.idx`, `t2g.txt`, optional `*_t2c.txt` |
| 4 | `04_kb_count.sh` | `kb_out/` (matrices, BUS file, `run_info.json`) |
| 5 | `05_cell_calling.sh` | `*_cell_calling.json`, `*_knee_data.json` |
| 6 | `06_doublet_detect.sh` | `*_doublet_scores.csv`, `*_doublet_meta.json` |
| 7 | `07_qc_metrics.sh` | `*_qc_metrics.csv`, `*_qc_summary.json` |
| 8 | `08_sample_report.sh` | `*_report.html` |
| 9 | `09_group_report.sh` | `group_report.html` |
| 10 | `10_celltypist_annotate.sh` | `*_celltypist_predictions.csv`, `*_celltypist_meta.json` |

All scripts work with both `--workflow quick` and `--workflow full` outputs — they consume the same matrix layout.

## 5. Failure modes you will encounter (and how to recover)

This pipeline has been beaten on. The validated recovery patterns are below — when you see these signatures, jump to the fix.

| Symptom | Cause | Fix |
|---|---|---|
| `kallisto bus ... died with std::bad_array_new_length` | Index built with kallisto 0.52.x, count run with 0.51.x (or vice versa) | Rebuild index with the same kallisto the container ships. Force `--kallisto /opt/conda/bin/kallisto` in kb ref. |
| `kallisto bus ... died with SIGKILL` | Docker OOM-killed the container (memory cap < index + working memory) | Switch to `-profile hpc` or `-profile aws`. Local profile cannot run full mode on the human nac index. |
| FASTP hangs at <1% CPU for >20 min on one sample | gzip-block pathology, EBS throughput throttle, or random sra-tools artifact | Already handled: `withName: 'FASTP'` has a 20-min cap, retries once, then `errorStrategy='ignore'`. KB_COUNT doesn't depend on FASTP, so other samples and downstream proceed. That sample's HTML report won't be produced. |
| `fasterq-dump quit with error code 3` (during ad-hoc fastq prep) | sra-tools 2.11 flakiness; not a pipeline failure but blocks data prep | Pre-prep fastqs once locally with sra-tools, upload to S3 (or shared storage), have downstream agents `aws s3 cp` them instead of re-running `fasterq-dump`. |
| Nextflow compile error `No such variable: sample_id` | Running on Nextflow 26.04.1+ (publishDir tightening) | Pin Nextflow to 24.10.x: `export NXF_VER=24.10.6` before `nextflow run`. |
| Nextflow error `Invalid include source: workflows/scq.nf` | Tarball/clone missing the `workflows/` dir | Re-fetch the full repo. The pipeline needs `main.nf`, `workflows/`, `modules/`, `conf/`, `scripts/`, `assets/`, `envs/`, `containers/`. |
| Celltypist permission denied writing model | Docker container ran as mambauser (uid 1000) into a root-owned mount | When invoking `kb`-style docker commands manually outside Nextflow, pass `-u 0`. Inside Nextflow it's handled by `runOptions = '-u $(id -u):$(id -g)'`. |

## 6. Pre-built artifact reuse (the "warm start" pattern)

For repeat runs against the same reference + samples, S3 (or any shared object store) lets you skip the slow steps:

```
s3://<bucket>/built_index/   # nac index — ~7 GB, saves ~45 min kb ref
  index_nac.idx
  t2g_nac.txt
  cdna_t2c.txt
  nascent_t2c.txt

s3://<bucket>/fastqs/        # pre-prepped fastqs — saves fasterq-dump + upload
  <SRR>_R1.fastq.gz
  <SRR>_R2.fastq.gz

s3://<bucket>/scq_images.tar.gz  # docker image, saves rebuild + push
```

**Agent automation rule:** before triggering a fresh run, check whether the index for the requested genome+gtf already exists in the shared store. If yes, `--prebuilt_index` it. If not, build with the SAME container image that will count, then push to the store before kicking off downstream.

## 7. Validation dataset (the four-sample demo)

The pipeline is validated against four public samples that exercise the chemistry-detection paths and the doublet-detection signal:

| Accession | Platform / chemistry | Cells | Source |
|---|---|---|---|
| SRR21064278 | Element AVITI, 10x 3' v3 | ~1K | PBMC, demonstrator |
| SRR21064279 | Element AVITI, 10x 3' v3 | ~10K | PBMC, demonstrator |
| SRR33398036 | Illumina NextSeq 2000, 10x 3' v3 | ~13K | "scCLEAN" Jumpcode-depleted (from *A CRISPR/Cas9-based enhancement of high-throughput single-cell transcriptomics*) |
| SRR33398039 | Illumina NextSeq 2000, 10x 3' v3 | ~13K | Control matched to SRR33398036 |

**What this dataset proves:**
- The chemistry detector works on both Element AVITI and Illumina output (chemistry is identical 10x 3' v3, but read-length and quality profiles differ).
- The pseudoalignment + cell-calling pipeline produces 10x-Cell-Ranger-comparable cell counts on both platforms.
- The pseudo-experiment (scCLEAN vs. matched control) reproduces published depletion signal: scCLEAN sample shows higher pseudoalignment rate (84% vs. 81%) and dramatically lower mitochondrial fraction (0.02% vs. 5.1%).
- Doublet rate stability: quick-mode runs of the smallest sample (~1K cells) flag the doublet rate as "unreliable" due to low per-cell UMI counts; full-mode runs collapse it to a realistic value, demonstrating the **integrity-check use case** of running quick mode first.

Group + per-sample HTML reports for both modes ship under `results-with-cDNA-onlyAlignments/` (quick) and `results-with-splicedAlignments/` (full) so users can compare side-by-side without running the pipeline themselves.

## 8. Cost prediction (rule of thumb)

The pipeline is designed so a user can predict full-mode cost from a quick-mode dry run.

| Mode | Per-sample wall on m6i.4xlarge | Per-sample wall on laptop (M1/Intel, 16 GB) | $/sample on AWS |
|---|---|---|---|
| `quick` (cDNA only) | 5-15 min | 10-30 min | ≈ $0.10 |
| `full` (nac, w/ prebuilt index) | 30-90 min | not feasible | ≈ $0.50 - $1.20 |
| `full` (nac, building index inline) | + 45 min one-time | not feasible | + ≈ $0.60 one-time |

A user processing 10 samples should expect ≈ $5-12 on AWS at validated configuration. Local quick-mode is free of compute charges but only captures ~50-60% of the read signal.

## 9. Automation surface (what an agent should be able to do)

Concrete autonomous actions, in priority order:

1. **Triage a samplesheet.** Inspect each row's chemistry hint, read length, and fastq size. If any sample is >20 GB R2 and the host has <32 GB RAM, switch the planned mode to `quick` and warn the user that downstream estimates will be partial.
2. **Reference selection.** Resolve genome species → check shared store → if absent, build with `scripts/shell/03_build_reference.sh` and publish to the shared store before triggering the pipeline.
3. **Profile selection.** Match the host's `nproc` and `free -g` against the profile table in §2; default to `local` if uncertain.
4. **Recovery actions.** On finalize-trap failure, parse `scq_bootstrap.log` for any of the symptoms in §5 and apply the matching fix in §5 without operator intervention. Re-launch with `-resume` from the existing work directory, OR with the same `--outdir` if the run was on a fresh ephemeral instance.
5. **Result delivery.** On pipeline success, copy `<outdir>/group_report.html`, all `<outdir>/<sample>/<sample>_report.html`, and `<outdir>/pipeline_info/` to the user's destination. The HTML reports are self-contained (Plotly inlined) — no further hosting needed.
6. **Cost reconciliation.** On AWS runs, read the trace.txt for actual realtime per process, compute total CPU-hours, and report against the predicted cost in §8. Flag drifts >50%.

## 10. Repo invariants — do not change these without coordination

- `kb_count.nf` and `kb_ref.nf` both pass `--kallisto $(command -v kallisto)` explicitly. **Never remove this.** It's the only thing preventing kb-python's bundled 0.52.0 binary from silently producing an incompatible index.
- `conf/base.config` `withName: 'FASTP'` retry+ignore is intentional. Removing it makes the pipeline brittle to gzip-block pathologies on any single fastq.
- The per-sample HTML report uses `FASTP.out.json.join(KB_COUNT.out.counts)` (inner join). A sample with a failed FASTP will be silently absent from per-sample reports but present in the group report (which uses KB_COUNT output, not FASTP). This is the intended fallback — preserve it.
- The Nextflow version is pinned via `export NXF_VER=24.10.6` in the AWS bootstrap and the user-facing README. Bumping requires testing publishDir input-var compatibility (broken in 26.04.1).
