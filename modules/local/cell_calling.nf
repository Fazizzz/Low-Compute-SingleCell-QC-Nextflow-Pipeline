process CELL_CALLING {
    tag "${sample_id}"
    label 'process_medium'
    publishDir "${params.outdir}/${sample_id}", mode: 'copy'

    input:
    tuple val(sample_id), path(kb_dir), val(species)
    path script

    output:
    tuple val(sample_id), path('filtered_matrix'), val(species),    emit: filtered
    tuple val(sample_id), path("${sample_id}_cell_calling.json"),    emit: json
    tuple val(sample_id), path("${sample_id}_knee_data.json"),       emit: knee_data

    script:
    """
    python3 ${script} \\
        --sample_id ${sample_id} \\
        --kb_dir ${kb_dir} \\
        --threshold ${params.min_umi_threshold} \\
        --outdir .
    mv cell_calling.json ${sample_id}_cell_calling.json
    mv knee_data.json    ${sample_id}_knee_data.json
    """

    stub:
    """
    mkdir -p filtered_matrix
    echo "%%MatrixMarket matrix coordinate integer general" > filtered_matrix/matrix.mtx
    echo "%" >> filtered_matrix/matrix.mtx
    echo "2 1 2" >> filtered_matrix/matrix.mtx
    echo "1 1 250" >> filtered_matrix/matrix.mtx
    echo "2 1 300" >> filtered_matrix/matrix.mtx
    gzip filtered_matrix/matrix.mtx
    echo "BC_STUB_1" | gzip > filtered_matrix/barcodes.tsv.gz
    printf "ENSG1\\tGENE1\\tGene Expression\\nENSG2\\tGENE2\\tGene Expression\\n" | gzip > filtered_matrix/features.tsv.gz
    cat > ${sample_id}_cell_calling.json <<JSON
    {"sample_id":"${sample_id}","threshold_used":${params.min_umi_threshold},"cells_called":1,"n_barcodes_total":3,"median_umi_called":250,"total_umi_called":250,"auto_knee_umi":150,"auto_knee_rank":2,"knee_status":"reliable"}
    JSON
    cat > ${sample_id}_knee_data.json <<JSON
    {"threshold":${params.min_umi_threshold},"auto_knee_umi":150,"auto_knee_rank":2,"knee_status":"reliable","ranks":[1,2,3],"umis":[300,250,5],"called_mask":[true,true,false]}
    JSON
    """
}
