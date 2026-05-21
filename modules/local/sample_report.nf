process SAMPLE_REPORT {
    tag "${sample_id}"
    label 'process_low'
    publishDir "${params.outdir}/${sample_id}", mode: 'copy'

    input:
    tuple val(sample_id),
          path(fastp_json),
          path(chemistry_json),
          path(kb_dir),
          val(species),
          path(cell_calling_json),
          path(knee_data_json),
          path(doublet_meta),
          path(qc_csv),
          path(qc_summary),
          path(celltypist_csv),
          path(celltypist_meta)
    path script
    path template

    output:
    tuple val(sample_id), path("${sample_id}_report.html"), emit: html

    script:
    def species_str = species ?: 'human'
    def outdir_abs = file(params.outdir).toAbsolutePath()
    def trace_path = "${outdir_abs}/pipeline_info/trace.txt"
    """
    python3 ${script} \\
        --sample_id ${sample_id} \\
        --species ${species_str} \\
        --fastp_json ${fastp_json} \\
        --chemistry_json ${chemistry_json} \\
        --kb_dir ${kb_dir} \\
        --cell_calling_json ${cell_calling_json} \\
        --knee_data_json ${knee_data_json} \\
        --doublet_meta ${doublet_meta} \\
        --qc_csv ${qc_csv} \\
        --qc_summary ${qc_summary} \\
        --celltypist_csv ${celltypist_csv} \\
        --celltypist_meta ${celltypist_meta} \\
        --template ${template} \\
        --out_html ${sample_id}_report.html \\
        --pipeline_version ${workflow.manifest.version} \\
        --scq_workflow ${params.workflow} \\
        --trace ${trace_path}
    """

    stub:
    """
    echo "<html><body>stub sample report ${sample_id}</body></html>" > ${sample_id}_report.html
    """
}
