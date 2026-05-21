process KB_COUNT {
    tag "${sample_id}"
    label 'process_high'
    publishDir "${params.outdir}/${sample_id}/kb_count", mode: 'copy'

    input:
    tuple val(sample_id), path(r1), path(r2), val(chemistry), val(species)
    path index
    path t2g
    path cdna_t2c     // empty file in quick mode
    path nascent_t2c  // empty file in quick mode

    output:
    tuple val(sample_id), path("kb_out"), val(species), emit: counts
    tuple val(sample_id), path("${sample_id}_run_info.json"), emit: run_info

    script:
    def mem_gb = (task.memory ? task.memory.toGiga() : 8)
    // bustools sort -m: min of the task memory and the profile cap.
    // Profile defaults: local=10, aws=32, hpc=48. Override: --bustools_mem_cap_gb=N.
    def cap_gb = (params.bustools_mem_cap_gb ?: 12) as int
    def kb_mem = Math.min(mem_gb as int, cap_gb)
    // kallisto bus thread count. Full (nac) mode keeps the ~7 GB index
    // resident and each thread adds working memory, so per-profile defaults
    // are smaller. Override: --kb_threads_full / --kb_threads_quick.
    def threads = (params.workflow == 'full'
                       ? (params.kb_threads_full  ?: 2)
                       : (params.kb_threads_quick ?: task.cpus)) as int
    def kallisto_bin = "\$(command -v kallisto)"
    def bustools_bin = "\$(command -v bustools)"
    if (params.workflow == 'full') {
        """
        mkdir -p kb_out
        kb count \\
            --workflow nac \\
            --kallisto ${kallisto_bin} \\
            --bustools ${bustools_bin} \\
            -i ${index} \\
            -g ${t2g} \\
            -c1 ${cdna_t2c} \\
            -c2 ${nascent_t2c} \\
            -x ${chemistry} \\
            -t ${threads} \\
            -m ${kb_mem}G \\
            -o kb_out \\
            --sum=total \\
            ${r1} ${r2}
        cp kb_out/run_info.json ${sample_id}_run_info.json
        """
    } else {
        """
        mkdir -p kb_out
        kb count \\
            --kallisto ${kallisto_bin} \\
            --bustools ${bustools_bin} \\
            -i ${index} \\
            -g ${t2g} \\
            -x ${chemistry} \\
            -t ${threads} \\
            -m ${kb_mem}G \\
            -o kb_out \\
            ${r1} ${r2}
        cp kb_out/run_info.json ${sample_id}_run_info.json
        """
    }

    stub:
    """
    mkdir -p kb_out/counts_unfiltered
    echo "%%MatrixMarket matrix coordinate integer general" > kb_out/counts_unfiltered/cells_x_genes.mtx
    echo "%" >> kb_out/counts_unfiltered/cells_x_genes.mtx
    echo "3 2 4" >> kb_out/counts_unfiltered/cells_x_genes.mtx
    echo "1 1 5" >> kb_out/counts_unfiltered/cells_x_genes.mtx
    echo "1 2 1" >> kb_out/counts_unfiltered/cells_x_genes.mtx
    echo "2 1 250" >> kb_out/counts_unfiltered/cells_x_genes.mtx
    echo "3 2 300" >> kb_out/counts_unfiltered/cells_x_genes.mtx
    printf "BC1\\nBC2\\nBC3\\n" > kb_out/counts_unfiltered/cells_x_genes.barcodes.txt
    printf "ENSG1\\nENSG2\\n" > kb_out/counts_unfiltered/cells_x_genes.genes.txt
    printf "GENE1\\nGENE2\\n" > kb_out/counts_unfiltered/cells_x_genes.genes.names.txt
    echo '{"n_processed":1000,"n_pseudoaligned":850,"p_pseudoaligned":85.0}' > kb_out/run_info.json
    echo '{"version":"0.51.1"}' > kb_out/kb_info.json
    echo '{"numBarcodes":3,"numReads":1000}' > kb_out/inspect.json
    cp kb_out/run_info.json ${sample_id}_run_info.json
    """
}
