#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { SCQ } from './workflows/scq.nf'

def validateParams() {
    def errors = []

    if (!params.samplesheet) {
        errors << "Missing required parameter: --samplesheet"
    } else if (!file(params.samplesheet).exists()) {
        errors << "Samplesheet not found: ${params.samplesheet}"
    }

    def hasFasta = params.genome_fasta != null
    def hasIndex = params.prebuilt_index != null

    if (hasFasta && hasIndex) {
        errors << "Provide either --genome_fasta/--genome_gtf OR --prebuilt_index/--t2g, not both"
    }
    if (!hasFasta && !hasIndex) {
        errors << "Provide either --genome_fasta/--genome_gtf OR --prebuilt_index/--t2g"
    }
    if (hasFasta && !params.genome_gtf) {
        errors << "--genome_gtf is required when --genome_fasta is provided"
    }
    if (hasIndex && !params.t2g) {
        errors << "--t2g is required when --prebuilt_index is provided"
    }
    if (hasFasta && !file(params.genome_fasta).exists()) {
        errors << "Genome FASTA not found: ${params.genome_fasta}"
    }
    if (params.genome_gtf && !file(params.genome_gtf).exists()) {
        errors << "GTF not found: ${params.genome_gtf}"
    }
    if (hasIndex && !file(params.prebuilt_index).exists()) {
        errors << "Prebuilt index not found: ${params.prebuilt_index}"
    }
    if (params.t2g && !file(params.t2g).exists()) {
        errors << "t2g not found: ${params.t2g}"
    }

    if (!(params.workflow in ['quick', 'full'])) {
        errors << "--workflow must be 'quick' or 'full' (got '${params.workflow}')"
    }
    if (params.workflow == 'full' && hasIndex) {
        if (!params.cdna_t2c || !params.nascent_t2c) {
            errors << "--workflow full with --prebuilt_index requires --cdna_t2c and --nascent_t2c"
        } else {
            if (!file(params.cdna_t2c).exists()) errors << "cdna_t2c not found: ${params.cdna_t2c}"
            if (!file(params.nascent_t2c).exists()) errors << "nascent_t2c not found: ${params.nascent_t2c}"
        }
    }
    if (params.run_celltypist && !params.celltypist_model_path) {
        errors << "--celltypist_model_path is required when --run_celltypist is enabled"
    }
    if (params.celltypist_model_path && !file(params.celltypist_model_path).exists()) {
        errors << "CellTypist model not found: ${params.celltypist_model_path}"
    }

    if (errors) {
        log.error "Parameter validation failed:\n  - " + errors.join('\n  - ')
        System.exit(1)
    }
}

workflow {
    validateParams()

    samples = Channel
        .fromPath(params.samplesheet)
        .splitCsv(header: true)
        .map { row ->
            def r1 = file(row.fastq_r1)
            def r2 = file(row.fastq_r2)
            tuple(row.sample_id, r1, r2, row.chemistry, row.species)
        }

    SCQ(samples)
}
