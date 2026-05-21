process KB_REF {
    label 'process_high'
    storeDir { params.ref_cache_dir }

    input:
    path genome_fasta
    path genome_gtf

    output:
    path 'index.idx',         emit: index
    path 't2g.txt',           emit: t2g
    path 'cdna.fa',           emit: cdna
    path 'cdna_t2c.txt',      emit: cdna_t2c,    optional: true
    path 'nascent_t2c.txt',   emit: nascent_t2c, optional: true

    script:
    def kallisto_bin = "\$(command -v kallisto)"
    def bustools_bin = "\$(command -v bustools)"
    if (params.workflow == 'full') {
        """
        kb ref \\
            --workflow nac \\
            --kallisto ${kallisto_bin} \\
            --bustools ${bustools_bin} \\
            -i index.idx \\
            -g t2g.txt \\
            -f1 cdna.fa \\
            -f2 nascent.fa \\
            -c1 cdna_t2c.txt \\
            -c2 nascent_t2c.txt \\
            ${genome_fasta} \\
            ${genome_gtf}
        """
    } else {
        """
        kb ref \\
            --kallisto ${kallisto_bin} \\
            --bustools ${bustools_bin} \\
            -i index.idx \\
            -g t2g.txt \\
            -f1 cdna.fa \\
            ${genome_fasta} \\
            ${genome_gtf}
        """
    }

    stub:
    if (params.workflow == 'full') {
        """
        touch index.idx t2g.txt cdna.fa cdna_t2c.txt nascent_t2c.txt
        """
    } else {
        """
        touch index.idx t2g.txt cdna.fa
        """
    }
}
