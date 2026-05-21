process CHEMISTRY_CHECK {
    tag "${sample_id}"
    label 'process_low'
    publishDir "${params.outdir}/${sample_id}/chemistry", mode: 'copy'

    input:
    tuple val(sample_id), path(r1), path(r2), val(chemistry), val(species)
    path script

    output:
    tuple val(sample_id), path("${sample_id}_chemistry.json"), emit: json

    script:
    """
    python3 ${script} \\
        --sample_id ${sample_id} \\
        --r1 ${r1} \\
        --chemistry ${chemistry} \\
        --out ${sample_id}_chemistry.json
    """

    stub:
    """
    cat > ${sample_id}_chemistry.json <<JSON
    {"sample_id":"${sample_id}","specified_chemistry":"${chemistry}","detected_chemistry":"10xv3-class","modal_length":28,"reads_sampled":10000,"length_histogram":{"28":10000},"status":"PASS"}
    JSON
    """
}
