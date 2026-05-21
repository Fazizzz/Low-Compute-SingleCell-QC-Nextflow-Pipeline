# Full mode demo results

These are the `--workflow full` (nac index, spliced + nascent) outputs from the four-sample validation set: SRR21064278, SRR21064279, SRR33398036, SRR33398039. Generated on AWS EC2 m6i.4xlarge with `-profile aws,docker`. Total pipeline runtime ~7 hours; per-sample compute cost ~$1.00 USD.

Compare these reports against the quick-mode counterparts in `../results-with-cDNA-onlyAlignments/` to see how full (nac) alignment lifts pseudoalignment from ~55% to ~80%, recovers ~40-60% more median genes per cell, and resolves doublet detection on small samples (SRR21064278 drops from a flagged 90% rate to a realistic 3.7%).

## What is preserved here

- **Per-sample HTML reports** (`<sample>/<sample>_report.html`) — the primary deliverable.
- **Group report** (`group_report.html`).
- **JSON / CSV metadata** — cell calling, doublet detection, QC summary, CellTypist predictions, chemistry detection, fastp.
- **Pipeline diagnostics** — `pipeline_info/trace.txt`, `execution_report.html`, `software_versions.yml`.

## What is NOT preserved (regenerable from a re-run)

- `kb_count/kb_out/` — the nac-mode BUS file and equivalence classes (`matrix.ec` is 64-200 MB per sample).
- `filtered_matrix/` — the Cell Ranger style MTX outputs (50-99 MB per sample). For a layout example see `../results-with-cDNA-onlyAlignments/SRR21064278/filtered_matrix/`; the full-mode equivalents have the same files, with denser counts because intronic reads are included.

## How to regenerate everything

```bash
nextflow run main.nf -profile aws,docker --workflow full \
    --samplesheet your_samples.csv \
    --prebuilt_index your/index_nac.idx --t2g your/t2g_nac.txt \
    --cdna_t2c your/cdna_t2c.txt --nascent_t2c your/nascent_t2c.txt
```

Full mode requires the nac index built with kallisto 0.51.1; build it with `scripts/shell/03_build_reference.sh genome.fa.gz genome.gtf.gz outdir --full`. See [AGENTS.md](../AGENTS.md) for the warm-start S3 layout if you are running this multiple times.
