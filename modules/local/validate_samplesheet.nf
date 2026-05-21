process VALIDATE_SAMPLESHEET {
    label 'process_low'
    publishDir "${params.outdir}/pipeline_info", mode: 'copy'

    input:
    path samplesheet

    output:
    path 'samplesheet.validated.csv', emit: csv

    script:
    """
    python3 <<'PY'
import csv, sys, os

ALLOWED_CHEM = {'10xv2','10xv3','10xv3.1','dropseq'}
ALLOWED_SPECIES = {'human','mouse'}
REQUIRED = ['sample_id','fastq_r1','fastq_r2','chemistry','species']

errors = []
with open('${samplesheet}', newline='') as fh:
    reader = csv.DictReader(fh)
    missing = [c for c in REQUIRED if c not in reader.fieldnames]
    if missing:
        errors.append(f"Missing columns: {missing}")
    rows = list(reader)

for i, row in enumerate(rows, start=2):
    if row.get('chemistry') not in ALLOWED_CHEM:
        errors.append(f"row {i}: chemistry '{row.get('chemistry')}' not in {ALLOWED_CHEM}")
    if row.get('species') not in ALLOWED_SPECIES:
        errors.append(f"row {i}: species '{row.get('species')}' not in {ALLOWED_SPECIES}")
    for col in ('fastq_r1','fastq_r2'):
        if not row.get(col):
            errors.append(f"row {i}: missing {col}")

if errors:
    for e in errors:
        print('ERROR:', e, file=sys.stderr)
    sys.exit(1)

with open('samplesheet.validated.csv', 'w', newline='') as fh:
    writer = csv.DictWriter(fh, fieldnames=REQUIRED)
    writer.writeheader()
    for row in rows:
        writer.writerow({k: row[k] for k in REQUIRED})
PY
    """

    stub:
    """
    cp ${samplesheet} samplesheet.validated.csv
    """
}
