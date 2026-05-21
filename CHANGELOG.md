# Changelog

## 0.1.0 (2026-05-16)

Initial release.

- Nextflow DSL2 pipeline scaffold
- FASTP (QC-only), chemistry detection, kallisto|bustools alignment
- Knee-detected + hard-threshold cell calling, Cell Ranger MTX output
- Scrublet doublet detection with graceful fail
- Per-cell QC metrics including MT% with species-aware prefix
- Per-sample interactive HTML report (Plotly, self-contained)
- Group HTML report: pass/fail table, bar charts, warnings, trace usage
- Optional CellTypist annotation
- Stub blocks on every module for CI
- Standalone shell wrappers for every Python step
