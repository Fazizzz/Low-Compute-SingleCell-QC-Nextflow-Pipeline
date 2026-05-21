process CELLTYPIST_ANNOTATE {
    tag "${sample_id}"
    label 'process_medium'
    publishDir "${params.outdir}/${sample_id}/celltypist", mode: 'copy'

    input:
    tuple val(sample_id), path(filtered_dir), val(species)
    path model
    path script

    output:
    tuple val(sample_id), path("${sample_id}_celltypist_predictions.csv"),
                           path("${sample_id}_celltypist_meta.json"), emit: results

    script:
    """
    python3 ${script} \\
        --sample_id ${sample_id} \\
        --matrix_dir ${filtered_dir} \\
        --model ${model} \\
        --out_csv ${sample_id}_celltypist_predictions.csv \\
        --out_meta ${sample_id}_celltypist_meta.json
    """

    stub:
    """
    echo "barcode,predicted_labels,majority_voting" > ${sample_id}_celltypist_predictions.csv
    cat > ${sample_id}_celltypist_meta.json <<JSON
    {"sample_id":"${sample_id}","status":"skipped","reason":"stub","n_cells":0,"model":"${model}"}
    JSON
    """
}
