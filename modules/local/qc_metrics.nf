process QC_METRICS {
    tag "${sample_id}"
    label 'process_low'
    publishDir "${params.outdir}/${sample_id}/qc", mode: 'copy'

    input:
    tuple val(sample_id), path(filtered_dir), val(species), path(doublet_csv)
    path script

    output:
    tuple val(sample_id), path("${sample_id}_qc_metrics.csv"), emit: csv
    tuple val(sample_id), path("${sample_id}_qc_summary.json"), emit: summary

    script:
    """
    python3 ${script} \\
        --sample_id ${sample_id} \\
        --matrix_dir ${filtered_dir} \\
        --doublet_csv ${doublet_csv} \\
        --species ${species} \\
        --out_csv ${sample_id}_qc_metrics.csv \\
        --out_summary ${sample_id}_qc_summary.json
    """

    stub:
    """
    echo "barcode,umi,genes,pct_mt,doublet_score,predicted_doublet" > ${sample_id}_qc_metrics.csv
    echo "BC_STUB_1,250,42,3.5,0.12,0" >> ${sample_id}_qc_metrics.csv
    cat > ${sample_id}_qc_summary.json <<JSON
    {"sample_id":"${sample_id}","species":"${species}","n_cells":1,"n_genes_detected":42,"median_umi":250,"mean_umi":250.0,"p25_umi":250.0,"p75_umi":250.0,"median_genes":42,"mean_pct_mt":3.5,"median_pct_mt":3.5,"mt_genes_in_features":13,"mt_detection_warning":false}
    JSON
    """
}
