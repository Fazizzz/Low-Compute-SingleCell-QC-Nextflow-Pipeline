process DOUBLET_DETECT {
    tag "${sample_id}"
    label 'process_medium'
    publishDir "${params.outdir}/${sample_id}/doublets", mode: 'copy'

    input:
    tuple val(sample_id), path(filtered_dir), val(species)
    path script

    output:
    tuple val(sample_id), path("${sample_id}_doublet_scores.csv"), emit: scores
    tuple val(sample_id), path("${sample_id}_doublet_meta.json"),  emit: meta

    script:
    """
    set +e
    python3 ${script} \\
        --sample_id ${sample_id} \\
        --matrix_dir ${filtered_dir} \\
        --out_scores ${sample_id}_doublet_scores.csv \\
        --out_meta ${sample_id}_doublet_meta.json
    rc=\$?
    if [ \$rc -ne 0 ] || [ ! -s ${sample_id}_doublet_meta.json ]; then
        echo "scrublet exited with status \$rc; emitting skipped status" >&2
        echo "barcode,score,predicted_doublet" > ${sample_id}_doublet_scores.csv
        cat > ${sample_id}_doublet_meta.json <<JSON
    {"sample_id":"${sample_id}","status":"skipped","reason":"scrublet_crashed_exit_\$rc","n_cells":null,"doublet_rate":null,"expected_doublet_rate":null,"threshold":null}
    JSON
    fi
    exit 0
    """

    stub:
    """
    echo "barcode,score,predicted_doublet" > ${sample_id}_doublet_scores.csv
    echo "BC_STUB_1,0.12,0" >> ${sample_id}_doublet_scores.csv
    cat > ${sample_id}_doublet_meta.json <<JSON
    {"sample_id":"${sample_id}","status":"ok","reason":null,"n_cells":1,"doublet_rate":0.0,"expected_doublet_rate":0.06,"threshold":0.25}
    JSON
    """
}
