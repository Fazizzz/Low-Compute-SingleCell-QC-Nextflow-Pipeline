# scq: SingleCell Quick QC

[![CI](https://github.com/Fazizzz/Low-Compute-SingleCell-QC-Nextflow-Pipeline/actions/workflows/ci.yml/badge.svg)](https://github.com/Fazizzz/Low-Compute-SingleCell-QC-Nextflow-Pipeline/actions/workflows/ci.yml)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Nextflow](https://img.shields.io/badge/nextflow-%E2%89%A524.10.6-23aa62.svg)](https://www.nextflow.io/)

## Overview

`scq` is a low-compute Nextflow DSL2 pipeline for fast quality control of single-cell RNA-seq libraries. It replaces the heavy genome-alignment step at the core of Cell Ranger and STARsolo with [kallisto](https://pachterlab.github.io/kallisto/) **pseudoalignment using k-mers**. Pseudoalignment breaks each read into short fixed-length substrings (k-mers), looks those k-mers up in a pre-built transcript index, and assigns the read to the set of transcripts whose k-mers it is compatible with. The read is never matched base-by-base to the genome; only its k-mer composition is. This is roughly an order of magnitude faster than full alignment and small enough to fit a 16 GB laptop.

The tradeoff is that pseudoalignment **does not produce base-level resolution and does not provide a position (alignment coordinate) for the read on the genome or transcript**. You get the assignment to a transcript or set of transcripts and the count of UMIs that landed on each, which is exactly what is needed for transcript quantification, cell calling, doublet detection, and the rest of the QC metrics in this pipeline. It is not appropriate for tasks that need where-on-the-sequence-the-read-aligned information (variant calling, allele-specific expression, splice-junction discovery, etc.).

The kallisto layer is also chemistry-agnostic. Because reads are matched as k-mers and the cell barcode and UMI are parsed positionally from the R1 read structure rather than by a kit-specific demultiplexer, kallisto+bustools support every major single-cell chemistry: 10x Genomics 3' v2 / v3 / v3.1 / v4, BD Rhapsody, Drop-seq, inDrop, sci-RNA-seq, Smart-seq3, and others. This pipeline currently accepts `10xv2`, `10xv3`, `10xv3.1`, and `dropseq` in its samplesheet validator (the others are a schema-only addition away). End-to-end QC validation has so far been performed on **10x Genomics 3' v3 chemistry and the human genome**; mouse should work as a drop-in given the pipeline already supports `species=mouse` and the kallisto layer is genome-agnostic, but a mouse end-to-end run has not been published here yet.

Runs on a laptop in quick mode for integrity checks and compute-cost prediction. Scales to HPC or cloud in full mode for 10x Cell Ranger comparable count matrices. Outputs self-contained interactive HTML reports plus Cell Ranger style MTX matrices ready for scanpy or Seurat.

The pipeline is validated on public 10x Genomics 3' v3 datasets generated on both Element AVITI and Illumina NextSeq 2000 instruments, exercising multi-platform pseudoalignment, automatic chemistry detection, and reproducing the published scCLEAN (Jumpcode) CRISPR-Cas9 transcript depletion signal in PBMC samples. Library scales from 1,000 to 10,000 cells per sample (65 to 520 million reads) run end-to-end in quick mode on a 16 GB laptop. Full mode produces Cell Ranger comparable count matrices on AWS m6i.4xlarge for roughly one US dollar per sample.

## Contents

- [Overview](#overview)
- [Workflow](#workflow)
- [Features](#features)
- [Modes](#modes)
- [Key dependencies](#key-dependencies)
- [Install](#install)
- [Resources and downloads](#resources-and-downloads)
- [Usage](#usage)
- [Input and output formats](#input-and-output-formats)
- [Docker containers](#docker-containers)
- [Roadmap](#roadmap)
- [License](#license)
- [References](#references)
- [Acknowledgements](#acknowledgements)

## Workflow

```mermaid
flowchart LR
    A[Paired fastqs] --> B[FASTP read QC]
    A --> C[Chemistry auto-detect]
    A --> D[kallisto + bustools pseudoalignment]
    D --> E[Cell calling, knee-based]
    E --> F[Scrublet doublet detection]
    E --> G[Per-cell QC metrics]
    E --> H[CellTypist annotation, optional]
    B & C & D & F & G & H --> I[Per-sample HTML report]
    I --> J[Group HTML report]
```

*Figure 1: Pipeline stages. FASTP, chemistry detection, and pseudoalignment run in parallel per sample. Cell calling, doublet detection, QC metrics, and CellTypist annotation chain from the count matrix. All artifacts converge into per-sample HTML reports, which feed the multi-sample group report.*

## Features

- Cell Ranger style outputs (MTX, barcodes, genes) compatible with scanpy and Seurat without conversion.
- Self-contained interactive HTML reports with Plotly figures, no static asset hosting needed.
- Chemistry-agnostic kallisto+bustools backend supporting all major single-cell layouts. Current samplesheet validator accepts 10x 3' v2 / v3 / v3.1 and Drop-seq out of the box; adding 10x 3' v4, BD Rhapsody, sci-RNA-seq, Smart-seq3, or other kb-python-supported chemistries is a small schema-only change.
- Validated end-to-end on 10x Genomics 3' v3 chemistry across two sequencing platforms (Element AVITI, Illumina NextSeq 2000).
- Two compute envelopes from the same codebase: laptop friendly quick mode for integrity checks, HPC or cloud full mode for production-grade matrices.
- Resilient by default: per-task retry and graceful skip for FASTP failures, so one bad fastq does not block other samples.
- Profile based resource model so the same `nextflow run` invocation works on a 16 GB laptop, a 64 GB HPC node, or a tuned AWS instance.

## Modes

| Mode | Index | Reads captured | Memory | Typical use |
|---|---|---|---|---|
| `--workflow quick` (default) | cDNA only, ~1 GB | ~50 to 60 percent of input reads | fits 16 GB hosts | Library QC, doublet rate, mitochondrial fraction, compute-cost prediction before a full run |
| `--workflow full` | nac (spliced + nascent), ~7 GB for human | ~80 to 85 percent of input reads | 32 GB minimum, 64 GB recommended | Cell Ranger comparable UMI counts for downstream scanpy or Seurat analysis |

Quick mode is the recommended first pass for any new sample. It produces a complete report at a fraction of the cost and surfaces failures (low cell count, high doublet flag, chemistry mismatch) before you commit to the full run.

## Key dependencies

- [Nextflow](https://www.nextflow.io/) 24.10.x (pinned via `NXF_VER`; 26.x is not yet compatible with the publishDir patterns used here).
- [kallisto](https://pachterlab.github.io/kallisto/) 0.51.1 and [bustools](https://bustools.github.io/) 0.44.1, accessed via [kb-python](https://github.com/pachterlab/kb_python) 0.29.5.
- [fastp](https://github.com/OpenGene/fastp) 1.3+ for read-level QC.
- [scanpy](https://scanpy.readthedocs.io/) for matrix handling and the `scanpy.pp.scrublet` doublet wrapper.
- [CellTypist](https://www.celltypist.org/) (optional, enabled with `--run_celltypist`).
- Docker or Singularity if running with a containerized profile.

> **Agent operators:** see [AGENTS.md](AGENTS.md) for a structured operating guide. It documents profile selection, failure recovery patterns, the S3 warm-start convention for reusing pre-built indices, and the automation surface for autonomous runs.

> **CI status:** every push and pull request runs the pipeline in stub mode against both `--workflow quick` and `--workflow full` paths, plus a lint pass that resolves every profile (`local`, `hpc`, `aws`). See `.github/workflows/ci.yml`.

## Install

> The commands below are shown as multi-line shell blocks for readability. Please **run them one line at a time**, copying each line individually into your terminal and waiting for it to finish before pasting the next. Some lines (for example `conda env create`) take several minutes; pasting the whole block at once will queue lines while the previous one is still running and will make any error harder to find. If you are new to docker, also read the [Docker containers](#docker-containers) section first; it has a short install guide.

### Option A: conda environment (laptop or HPC login node)

```bash
git clone https://github.com/Fazizzz/Low-Compute-SingleCell-QC-Nextflow-Pipeline.git
cd Low-Compute-SingleCell-QC-Nextflow-Pipeline
conda env create -f envs/scq.yml
conda activate scq
export NXF_VER=24.10.6
```

For strict reproducibility, use `envs/scq.lock.yml` instead of `envs/scq.yml`. The lock file pins every dependency to the validated set.

### Option B: docker

If you have never used docker before, see the install steps in the [Docker containers](#docker-containers) section below first, then come back here. Once docker is installed and the daemon is running:

```bash
git clone https://github.com/Fazizzz/Low-Compute-SingleCell-QC-Nextflow-Pipeline.git
cd Low-Compute-SingleCell-QC-Nextflow-Pipeline
docker build -f containers/sc_qc_base/Dockerfile -t sc_qc_base:latest .
docker build -f containers/sc_qc_annotate/Dockerfile -t sc_qc_annotate:latest .
```

Then invoke Nextflow with `-profile docker` in the [Usage](#usage) examples below.

### Option C: singularity

Convert the docker images to `.sif` (`singularity build sc_qc_base.sif docker-daemon://sc_qc_base:latest`) and invoke with `-profile singularity`.

## Resources and downloads

The pipeline ships the code only; reference genomes and the optional CellTypist classification model are NOT bundled and must be downloaded separately. Choose the species that matches your samples.

### Reference genomes (required)

Quick mode needs a kallisto cDNA index. Full mode needs a kallisto nac (spliced + nascent) index. Either can be built from a primary-assembly FASTA + a GTF annotation. Recommended sources:

- **Human GRCh38, Ensembl release 110**
  - FASTA: [Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz](https://ftp.ensembl.org/pub/release-110/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz)
  - GTF: [Homo_sapiens.GRCh38.110.gtf.gz](https://ftp.ensembl.org/pub/release-110/gtf/homo_sapiens/Homo_sapiens.GRCh38.110.gtf.gz)
- **Mouse GRCm39, Ensembl release 110**
  - FASTA: [Mus_musculus.GRCm39.dna.primary_assembly.fa.gz](https://ftp.ensembl.org/pub/release-110/fasta/mus_musculus/dna/Mus_musculus.GRCm39.dna.primary_assembly.fa.gz)
  - GTF: [Mus_musculus.GRCm39.110.gtf.gz](https://ftp.ensembl.org/pub/release-110/gtf/mus_musculus/Mus_musculus.GRCm39.110.gtf.gz)

Build the index with the helper script:

```bash
scripts/shell/03_build_reference.sh <genome.fa.gz> <annotation.gtf.gz> <outdir>            # quick mode (cDNA only)
scripts/shell/03_build_reference.sh <genome.fa.gz> <annotation.gtf.gz> <outdir> --full     # full mode (nac index)
```

You can also point the pipeline at the genome and GTF directly with `--genome_fasta` and `--genome_gtf`; it will run the same step internally. Both paths produce identical artifacts (`index.idx`, `t2g.txt`, plus `cdna_t2c.txt` and `nascent_t2c.txt` for full mode).

### CellTypist model (optional, only for `--run_celltypist`)

The CellTypist model is a separate `.pkl` file that ships with the CellTypist project, NOT with this pipeline. It is required only if you enable `--run_celltypist true`. Browse the model catalogue at the [CellTypist model gallery](https://www.celltypist.org/models). Common starting points for PBMC and immune work:

- `Immune_All_Low.pkl` ([direct download](https://celltypist.cog.sanger.ac.uk/models/Pan_Immune_CellTypist/v2/Immune_All_Low.pkl)) - broad immune compartments at low granularity. Used in the validation runs in this repo.
- `Immune_All_High.pkl` ([direct download](https://celltypist.cog.sanger.ac.uk/models/Pan_Immune_CellTypist/v2/Immune_All_High.pkl)) - finer-grained immune labels.

To fetch one manually:

```bash
mkdir -p celltypist_models
curl -L -o celltypist_models/Immune_All_Low.pkl https://celltypist.cog.sanger.ac.uk/models/Pan_Immune_CellTypist/v2/Immune_All_Low.pkl
```

Pass the resulting path with `--celltypist_model_path celltypist_models/Immune_All_Low.pkl` (see Usage below).

## Usage

### Quick start (stub run, no data required, under 30 seconds)

```bash
nextflow run main.nf -stub-run \
    --samplesheet test/samplesheet.csv \
    --prebuilt_index ref/index.idx --t2g ref/t2g.txt
```

### Choose a compute profile, then a container engine

Profiles combine with commas: one compute envelope plus one container engine.

| Profile | Target host |
|---|---|
| `local` | Laptop or desktop, up to 16 GB RAM. Quick mode only. |
| `hpc` | HPC node or workstation with 32 GB RAM and 16+ cores. Full mode capable. |
| `aws` | EC2 m6i.4xlarge (16 vCPU, 64 GB RAM). Validated reference instance. |
| `docker` | Use the `sc_qc_base` / `sc_qc_annotate` docker images. |
| `singularity` | Use the `.sif` images. |

### Real run examples

Laptop, quick mode, conda env on PATH:

```bash
nextflow run main.nf -profile local \
    --samplesheet test/samplesheet.csv \
    --prebuilt_index ref/index.idx --t2g ref/t2g.txt
```

HPC, full mode, singularity:

```bash
nextflow run main.nf -profile hpc,singularity --workflow full \
    --samplesheet samples.csv \
    --prebuilt_index nac.idx --t2g t2g.txt \
    --cdna_t2c cdna_t2c.txt --nascent_t2c nascent_t2c.txt
```

AWS, full mode, docker, with CellTypist annotation:

```bash
nextflow run main.nf -profile aws,docker --workflow full \
    --samplesheet samples.csv \
    --prebuilt_index nac.idx --t2g t2g.txt \
    --cdna_t2c cdna_t2c.txt --nascent_t2c nascent_t2c.txt \
    --run_celltypist --celltypist_model_path Immune_All_Low.pkl
```

`--celltypist_model_path` must point at a `.pkl` file that you have already downloaded; the pipeline does NOT ship CellTypist models. See [Resources and downloads](#resources-and-downloads) for the model gallery and a fetch command.

Standalone shell wrappers for every pipeline stage live under `scripts/shell/` if you need to run individual steps outside Nextflow.

## Input and output formats

### Samplesheet

```csv
sample_id,fastq_r1,fastq_r2,chemistry,species
PBMC_01,/data/PBMC_01_R1.fastq.gz,/data/PBMC_01_R2.fastq.gz,10xv3,human
```

`chemistry` accepts `10xv2`, `10xv3`, `10xv3.1`, `dropseq`. `species` accepts `human` or `mouse`. Both fields are validated at run start.

### Outputs

```
results/
├── group_report.html
├── pipeline_info/
│   ├── trace.txt
│   ├── execution_report.html
│   ├── timeline.html
│   ├── software_versions.yml
│   └── samplesheet.validated.csv
└── <sample_id>/
    ├── <sample_id>_report.html
    ├── filtered_matrix/                          # Cell Ranger MTX (genes by cells)
    ├── <sample_id>_cell_calling.json
    ├── <sample_id>_knee_data.json
    ├── fastp/<sample_id>_fastp.{json,html}
    ├── chemistry/<sample_id>_chemistry.json
    ├── kb_count/                                 # kb count outputs
    ├── doublets/<sample_id>_doublet_{scores.csv,meta.json}
    ├── qc/<sample_id>_qc_{metrics.csv,summary.json}
    └── celltypist/<sample_id>_celltypist_{predictions.csv,meta.json}
```

### Group report

The group report opens with a pass/fail status table per sample, surfacing any warnings (low alignment, low cells called, doublet flag, chemistry mismatch) at the top.

![Group report pass/fail table](https://raw.githubusercontent.com/Fazizzz/Low-Compute-SingleCell-QC-Nextflow-Pipeline/main/assets/report-images/Main-Pass-Fail-Table.png)

*Figure 2: Group report header. The pass/fail table highlights any sample-level warnings that require attention before downstream analysis. In this example all four validation samples pass; the asterisk on alignment percentage marks quick-mode runs where intronic reads are expected to be missed.*

Below the status table, the group report shows per-sample bar charts for the key cross-sample metrics on a single page.

![Group alignment and QC bar charts](https://raw.githubusercontent.com/Fazizzz/Low-Compute-SingleCell-QC-Nextflow-Pipeline/main/assets/report-images/Group-Alignment-summary-plots.png)

*Figure 3: Group report cross-sample comparison. Bars show cells called, median UMIs per cell, median genes per cell, median percent mitochondrial reads, percent pseudoaligned, and doublet rate. The sea-green and complementary palette accommodates up to 12 samples without color reuse; runs with more than 20 samples fall back to a compact summary table automatically.*

### Per-sample report

Each sample report includes a top-line summary, the FASTP read QC panel, the alignment summary with mode-aware caveats, the barcode rank knee plot, per-cell UMI, gene, mitochondrial percent and doublet score distributions, and an optional CellTypist composition panel.

![Per-sample alignment and read QC](https://raw.githubusercontent.com/Fazizzz/Low-Compute-SingleCell-QC-Nextflow-Pipeline/main/assets/report-images/Per-Sample-Detail-Metrics.png)

*Figure 4: Per-sample alignment and read QC. Barcode rank plot with the auto-calculated knee threshold, full pseudoalignment counts, and FASTP read-quality metrics (total reads, Q30 rate, duplication fraction). Together these let users sanity-check chemistry detection, sequencing quality, and library prep before trusting the cell-level metrics below.*

![Per-cell distributions](https://raw.githubusercontent.com/Fazizzz/Low-Compute-SingleCell-QC-Nextflow-Pipeline/main/assets/report-images/Per-sample-stats.png)

*Figure 5: Per-cell distributions. UMI counts per cell, genes detected per cell, mitochondrial percent, and Scrublet doublet scores. Together these describe library complexity, sequencing depth, and the doublet content of the called cell population.*

When CellTypist is enabled, each sample receives a cell-type composition bar chart using the same expanded palette as the group plots.

![Per-sample CellTypist composition](https://raw.githubusercontent.com/Fazizzz/Low-Compute-SingleCell-QC-Nextflow-Pipeline/main/assets/report-images/Per-sample-Celltypist-annotation.png)

*Figure 6: CellTypist composition for a PBMC sample using the Immune_All_Low model. Categorical bars use the 12-color sea-green-anchored palette.*

### Compute and resource usage

Every report footer summarizes compute usage so users can predict the cost of a future run on the same hardware.

![Compute and resource summary](https://raw.githubusercontent.com/Fazizzz/Low-Compute-SingleCell-QC-Nextflow-Pipeline/main/assets/report-images/Total-compute-costs.png)

*Figure 7: Compute and resource usage table at the end of every per-sample report. Wall time, peak CPU percent, peak memory, requested CPUs, and requested memory per step provide a compute receipt; this is the basis of the cost prediction described in AGENTS.md (typically about one US dollar per sample in full mode on AWS m6i.4xlarge).*

## Docker containers

### What is docker (short version)

Docker lets a pipeline ship with all its scientific software (kallisto, bustools, scanpy, ...) pre-installed inside a single image. You run that image as a "container" and the pipeline behaves identically on a Mac, a Linux laptop, or a cloud node, without you having to install anything besides docker itself. It is the lowest-friction way to use this pipeline on a host where you do not control the system python or do not want to manage a conda environment.

### Install docker (if you do not already have it)

- **macOS / Windows**: install [Docker Desktop](https://www.docker.com/products/docker-desktop/), launch it once, and confirm in its menu bar that the docker engine is running.
- **Linux**: follow the [official engine install guide](https://docs.docker.com/engine/install/) for your distribution (Ubuntu / Debian / Fedora / RHEL each have a one-page recipe). After install, you typically also want to add your user to the docker group so you do not need sudo: `sudo usermod -aG docker $USER && newgrp docker`.

Verify with:

```bash
docker --version
docker run --rm hello-world
```

If the `hello-world` container pulls and prints "Hello from Docker!", you are ready.

### The two images in this repo

Two images, both built from the repo:

| Image | Built from | Purpose |
|---|---|---|
| `sc_qc_base` | `containers/sc_qc_base/Dockerfile` | Base image with kallisto, bustools, fastp, kb-python, scanpy, scrublet, and the python report helpers. Used by every process except CellTypist annotation. |
| `sc_qc_annotate` | `containers/sc_qc_annotate/Dockerfile` | Extends `sc_qc_base` with CellTypist and its model-loading dependencies. Used only by `CELLTYPIST_ANNOTATE`. |

The split keeps the celltypist torch and onnx footprint out of the base image, which keeps the docker layer pulled to compute nodes small.

To build both (run line by line, each step takes 5-15 min on a first build):

```bash
docker build -f containers/sc_qc_base/Dockerfile     -t sc_qc_base:latest     .
docker build -f containers/sc_qc_annotate/Dockerfile -t sc_qc_annotate:latest .
```

Once built, return to the [Install](#install) and [Usage](#usage) sections and invoke Nextflow with `-profile docker`.

Memory caps are not baked into the images. They are controlled by the Nextflow compute profile (`local`, `hpc`, or `aws`). The same image runs everywhere; only the resource envelope changes.

## Roadmap

- **nf-core submission.** Align directory layout, schema, and test data to the [nf-core](https://nf-co.re/) standard and submit for community review.
- **Broader organism coverage.** Mouse is already accepted at the schema level and should work as a drop-in, but a published mouse end-to-end run is still pending. Beyond mouse, the next targets are non-mammalian model systems (zebrafish, fly, worm, plant) with species-aware mitochondrial gene lists and validated reference indices.
- **Chemistry validation.** Extend the samplesheet validator and the chemistry detector to cover the kb-python-supported chemistries the pipeline does not yet expose: Parse Biosciences, 10x 3' v4, BD Rhapsody, sci-RNA-seq, inDrop, Smart-seq3, and others. The kallisto+bustools layer already supports these natively; the work is mostly schema and a few lookup-table rows.
- **Terraform module for HPC bootstrap.** A one-command Terraform plan that provisions a transient AWS or GCP node sized for full mode, runs the pipeline, ships results to object storage, and tears the node down. Captures the validated reference instance type and the warm-start S3 layout documented in AGENTS.md.
- **Wiki gallery.** Browsable comparison of quick mode vs full mode reports for the four validation samples, plus the scCLEAN pseudo-experiment results.

## License

This project is released under the [MIT License](LICENSE). See `LICENSE` for the full text.

## References

- Bray, N. L., Pimentel, H., Melsted, P., and Pachter, L. (2016). Near-optimal probabilistic RNA-seq quantification. *Nature Biotechnology*, 34, 525 to 527. [doi:10.1038/nbt.3519](https://doi.org/10.1038/nbt.3519)
- Melsted, P., Booeshaghi, A. S., Liu, L., et al. (2021). Modular, efficient and constant-memory single-cell RNA-seq preprocessing. *Nature Biotechnology*, 39, 813 to 818. [doi:10.1038/s41587-021-00870-2](https://doi.org/10.1038/s41587-021-00870-2)
- Wolock, S. L., Lopez, R., and Klein, A. M. (2019). Scrublet: computational identification of cell doublets in single-cell transcriptomic data. *Cell Systems*, 8(4), 281 to 291.e9. [doi:10.1016/j.cels.2018.11.005](https://doi.org/10.1016/j.cels.2018.11.005)
- Domínguez Conde, C., Xu, C., Jarvis, L. B., et al. (2022). Cross-tissue immune cell analysis reveals tissue-specific features in humans. *Science*, 376(6594). [doi:10.1126/science.abl5197](https://doi.org/10.1126/science.abl5197)
- Chen, S., Zhou, Y., Chen, Y., and Gu, J. (2018). fastp: an ultra-fast all-in-one FASTQ preprocessor. *Bioinformatics*, 34(17), i884 to i890. [doi:10.1093/bioinformatics/bty560](https://doi.org/10.1093/bioinformatics/bty560)
- Wolf, F. A., Angerer, P., and Theis, F. J. (2018). SCANPY: large-scale single-cell gene expression data analysis. *Genome Biology*, 19, 15. [doi:10.1186/s13059-017-1382-0](https://doi.org/10.1186/s13059-017-1382-0)
- Yu, L., Wang, X., Mu, Q., Tam, S. S. T., et al. (2024). scCLEAN improves the signal-to-noise ratio of single-cell transcriptomics by reducing noise from highly expressed genes. *Nature Communications*, 15, [scCLEAN paper]. The validation samples SRR33398036 (depleted) and SRR33398039 (control) are reused with thanks to the authors for the open release.
- Element Biosciences AVITI data (SRR21064278, SRR21064279) is reused from the public release accompanying the platform's initial single-cell demonstrators.
- Built with [Nextflow](https://www.nextflow.io/) and the conventions established by [nf-core](https://nf-co.re/). Thanks to the kallisto, bustools, and kb-python teams at the Pachter lab for the underlying alignment tooling.

## Acknowledgements

Muhammad Faizan Khalid: Author and current maintainer

This pipeline grew out of a desire to build something practical: a resource that enables everyday users to perform single-cell analysis at scale while remaining accessible for local, low-compute environments. It provides a framework for both small local analyses and larger-scale studies, acting as a training resource for users who want hands-on experience with downstream single-cell analysis. It also gives labs an opportunity to QC their data locally before committing to cloud costs.

The pipeline was developed as part of an ongoing portfolio in bioinformatics pipeline engineering, with a focus on single-cell sequencing, containerization, and QC resource development beyond the 10x Genomics ecosystem.

The pipeline architecture draws inspiration from the work of [Tommy Tang](https://www.youtube.com/watch?v=fVtiHHIvG-I), who provides several helpful resources on his [GitHub](https://github.com/crazyhottommy). It follows nf-core DSL2 module conventions, with intentional deviations documented in the codebase. The multi-sample framework and harmonized reporting tools reflect the kind of scalable approach needed in production sequencing environments handling tens to hundreds of samples.

This repository is provided for educational and demonstration purposes. It is not affiliated with Tommy Tang, nf-core, 10x Genomics, or any commercial organization. It includes well-commented scripts and a streamlined analysis workflow designed for low-compute environments, ease of use, and local implementation. It is not a commercial product. The code is provided "as is," without warranty of any kind. Bugs and feedback are welcome through the repository's issue tracker.

For citation or attribution, please reference:
Khalid, M. Faizan (or Khalid MF)

You can follow related research and professional updates via my [Google Scholar profile](https://scholar.google.com/citations?hl=en&user=qFZQ5wYAAAAJ&sortby=title&view_op=list_works&gmla=AL3_zigRWGX9g8Jc22idbBUMFuy7cVN_pEIyL6_DXSA-qWkJbcaONzhRNSmAwmQXKEm-3-WYGouZZC2pCE6zD9tZLxizbM7jQzzZMOgtkgsuL825u4lvSs9kwsccajhJbBg2Mrc37at_HCQ) or [LinkedIn](https://www.linkedin.com/in/m-faizan-khalid/).
