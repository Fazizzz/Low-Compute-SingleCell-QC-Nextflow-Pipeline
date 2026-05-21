process FASTP {
    tag "${sample_id}"
    label 'process_medium'
    publishDir "${params.outdir}/${sample_id}/fastp", mode: 'copy'

    input:
    tuple val(sample_id), path(r1), path(r2), val(chemistry), val(species)

    output:
    tuple val(sample_id), path("${sample_id}_fastp.json"), emit: json
    tuple val(sample_id), path("${sample_id}_fastp.html"), emit: html

    script:
    """
    fastp \\
        --in1 ${r1} \\
        --in2 ${r2} \\
        --disable_adapter_trimming \\
        --disable_quality_filtering \\
        --disable_length_filtering \\
        --thread ${task.cpus} \\
        --json ${sample_id}_fastp.json \\
        --html ${sample_id}_fastp.html
    """

    stub:
    """
    echo '{"summary":{"before_filtering":{"total_reads":1000,"total_bases":100000,"q30_rate":0.95},"after_filtering":{"total_reads":1000,"total_bases":100000,"q30_rate":0.95}},"duplication":{"rate":0.05}}' > ${sample_id}_fastp.json
    echo '<html><body>stub fastp report</body></html>' > ${sample_id}_fastp.html
    """
}
