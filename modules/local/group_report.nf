process GROUP_REPORT {
    label 'process_low'
    publishDir "${params.outdir}", mode: 'copy'

    input:
    path chemistry_jsons,    stageAs: 'chem/*'
    path cell_calling_jsons, stageAs: 'cc/*'
    path doublet_metas,      stageAs: 'dbl/*'
    path qc_summaries,       stageAs: 'qc/*'
    path run_infos,          stageAs: 'ri/*'
    path sample_htmls,       stageAs: 'rpt/*'
    path software_versions
    path script
    path template

    output:
    path 'group_report.html', emit: html

    script:
    def outdir_abs = file(params.outdir).toAbsolutePath()
    def trace_path = "${outdir_abs}/pipeline_info/trace.txt"
    """
    python3 ${script} \\
        --chemistry_jsons ${chemistry_jsons} \\
        --cell_calling_jsons ${cell_calling_jsons} \\
        --doublet_metas ${doublet_metas} \\
        --qc_summaries ${qc_summaries} \\
        --run_infos ${run_infos} \\
        --trace ${trace_path} \\
        --software_versions ${software_versions} \\
        --template ${template} \\
        --scq_workflow ${params.workflow} \\
        --out_html group_report.html
    """

    stub:
    """
    echo "<html><body>stub group report</body></html>" > group_report.html
    """
}
