include { VALIDATE_SAMPLESHEET } from '../modules/local/validate_samplesheet.nf'
include { FASTP }                from '../modules/local/fastp.nf'
include { CHEMISTRY_CHECK }      from '../modules/local/chemistry_check.nf'
include { KB_REF }               from '../modules/local/kb_ref.nf'
include { KB_COUNT }             from '../modules/local/kb_count.nf'
include { CELL_CALLING }         from '../modules/local/cell_calling.nf'
include { DOUBLET_DETECT }       from '../modules/local/doublet_detect.nf'
include { QC_METRICS }           from '../modules/local/qc_metrics.nf'
include { CELLTYPIST_ANNOTATE }  from '../modules/local/celltypist_annotate.nf'
include { SAMPLE_REPORT }        from '../modules/local/sample_report.nf'
include { GROUP_REPORT }         from '../modules/local/group_report.nf'
include { SOFTWARE_VERSIONS }    from '../modules/local/software_versions.nf'

workflow SCQ {
    take:
    samples  // tuple(sample_id, r1, r2, chemistry, species)

    main:
    VALIDATE_SAMPLESHEET(file(params.samplesheet))

    SOFTWARE_VERSIONS()

    FASTP(samples)

    chemistry_script = file("${projectDir}/scripts/python/chemistry_check.py")
    CHEMISTRY_CHECK(samples, chemistry_script)

    // Placeholders used when --workflow quick (KB_COUNT ignores the t2c inputs).
    // Two distinct files to avoid Nextflow input-name collisions.
    placeholder_cdna_t2c    = file("${projectDir}/assets/cdna_t2c_placeholder.txt")
    placeholder_nascent_t2c = file("${projectDir}/assets/nascent_t2c_placeholder.txt")

    if (params.prebuilt_index) {
        index_ch = Channel.value(file(params.prebuilt_index))
        t2g_ch   = Channel.value(file(params.t2g))
        if (params.workflow == 'full') {
            cdna_t2c_ch    = Channel.value(file(params.cdna_t2c))
            nascent_t2c_ch = Channel.value(file(params.nascent_t2c))
        } else {
            cdna_t2c_ch    = Channel.value(placeholder_cdna_t2c)
            nascent_t2c_ch = Channel.value(placeholder_nascent_t2c)
        }
    } else {
        KB_REF(file(params.genome_fasta), file(params.genome_gtf))
        index_ch = KB_REF.out.index
        t2g_ch   = KB_REF.out.t2g
        if (params.workflow == 'full') {
            cdna_t2c_ch    = KB_REF.out.cdna_t2c
            nascent_t2c_ch = KB_REF.out.nascent_t2c
        } else {
            cdna_t2c_ch    = Channel.value(placeholder_cdna_t2c)
            nascent_t2c_ch = Channel.value(placeholder_nascent_t2c)
        }
    }

    KB_COUNT(samples, index_ch, t2g_ch, cdna_t2c_ch, nascent_t2c_ch)

    cell_calling_script = file("${projectDir}/scripts/python/cell_calling.py")
    CELL_CALLING(KB_COUNT.out.counts, cell_calling_script)

    doublet_script = file("${projectDir}/scripts/python/doublet_detect.py")
    DOUBLET_DETECT(CELL_CALLING.out.filtered, doublet_script)

    qc_script = file("${projectDir}/scripts/python/qc_metrics.py")
    qc_input = CELL_CALLING.out.filtered.join(DOUBLET_DETECT.out.scores)
    QC_METRICS(qc_input, qc_script)

    // optional celltypist: always feed a (sample_id, csv, meta) tuple into
    // SAMPLE_REPORT; use placeholder files when disabled.
    if (params.run_celltypist) {
        ct_script = file("${projectDir}/scripts/python/celltypist_annotate.py")
        ct_model  = Channel.value(file(params.celltypist_model_path))
        CELLTYPIST_ANNOTATE(CELL_CALLING.out.filtered, ct_model, ct_script)
        ct_out = CELLTYPIST_ANNOTATE.out.results
    } else {
        no_ct_csv  = file("${projectDir}/assets/no_celltypist.csv")
        no_ct_meta = file("${projectDir}/assets/no_celltypist_meta.json")
        ct_out = samples.map { sid, _r1, _r2, _ch, _sp -> tuple(sid, no_ct_csv, no_ct_meta) }
    }

    // per-sample report inputs: collect all per-sample outputs by sample_id
    sample_report_script = file("${projectDir}/scripts/python/sample_report.py")
    sample_report_template = file("${projectDir}/scripts/python/templates/sample_report.html.j2")

    // join all per-sample channels
    report_input = FASTP.out.json
        .join(CHEMISTRY_CHECK.out.json)
        .join(KB_COUNT.out.counts)
        .join(CELL_CALLING.out.json)
        .join(CELL_CALLING.out.knee_data)
        .join(DOUBLET_DETECT.out.meta)
        .join(QC_METRICS.out.csv)
        .join(QC_METRICS.out.summary)
        .join(ct_out)

    SAMPLE_REPORT(report_input, sample_report_script, sample_report_template)

    // group report: collect everything
    group_report_script = file("${projectDir}/scripts/python/group_report.py")
    group_report_template = file("${projectDir}/scripts/python/templates/group_report.html.j2")

    GROUP_REPORT(
        CHEMISTRY_CHECK.out.json.map{ it[1] }.collect(),
        CELL_CALLING.out.json.map{ it[1] }.collect(),
        DOUBLET_DETECT.out.meta.map{ it[1] }.collect(),
        QC_METRICS.out.summary.map{ it[1] }.collect(),
        KB_COUNT.out.run_info.map{ it[1] }.collect(),
        SAMPLE_REPORT.out.html.map{ it[1] }.collect(),
        SOFTWARE_VERSIONS.out.yml,
        group_report_script,
        group_report_template
    )

    emit:
    group_html = GROUP_REPORT.out.html
}
