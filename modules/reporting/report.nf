process GENERATE_REPORT {

    tag "report"
    label 'low_compute'

    publishDir "${params.outdir}/reports", mode: 'copy'
    publishDir "${params.outdir}/logs", mode: 'copy', pattern: "report_log.txt"

    input:
    path qc_table
    path tree
    path provenance_tsv
    path figure_manifest
    path variants_tsv
    path variant_qc_tsv
    path metadata_validated
    path metadata_issues

    output:
    path "zika_summary.csv", emit: summary
    path "zika_report.html", emit: html
    path "report_log.txt", emit: log

    script:

    def metadata_arg = metadata_validated \
        ? "--metadata ${metadata_validated}" \
        : "--metadata /dev/null"

    def metadata_issues_arg = metadata_issues \
        ? "--metadata-issues ${metadata_issues}" \
        : "--metadata-issues /dev/null"

    """
    #!/usr/bin/env bash
    set -euo pipefail

    LOG="report_log.txt"

    echo "Stage: Report" > "\$LOG"
    echo "Started: \$(date -u +"%Y-%m-%dT%H:%M:%SZ")" >> "\$LOG"

    python3 ${projectDir}/bin/generate_report.py \\
        --qc ${qc_table} \\
        --tree ${tree} \\
        --provenance ${provenance_tsv} \\
        --manifest ${figure_manifest} \\
        --variants ${variants_tsv} \\
        --variant-qc ${variant_qc_tsv} \\
        ${metadata_arg} \\
        ${metadata_issues_arg} \\
        --pipeline-version ${params.pipeline_version} \\
        >> "\$LOG" 2>&1

    echo "Completed" >> "\$LOG"
    """
}
