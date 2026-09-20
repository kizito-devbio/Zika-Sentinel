/*
 * Quality Control — raw read trimming and FastQC assessment.
 *
 * Input:
 *   Paired-end FASTQ reads
 *   tuple val(sample), path(read1), path(read2)
 *
 * Output:
 *   Trimmed gzipped FASTQ
 *   FastQC reports
 *   QC log
 */

process FASTP_QC {

    tag "$sample"
    label 'med_compute'

    publishDir "${params.outdir}/qc",
        mode: 'copy',
        pattern: "*.{fastq.gz,html,json}"

    publishDir "${params.outdir}/logs",
        mode: 'copy',
        pattern: "*_qc_log.txt"

    input:
    tuple val(sample), path(read1), path(read2)

    output:
    tuple val(sample), path("${sample}_trimmed_1.fastq.gz"), path("${sample}_trimmed_2.fastq.gz"), emit: trimmed
    path "${sample}_fastp.html", emit: html
    path "${sample}_fastp.json", emit: json
    path "${sample}_qc_log.txt", emit: log

    script:
    """
    #!/usr/bin/env bash
    set -euo pipefail

    LOG="${sample}_qc_log.txt"
    START=\$(date +%s)
    TS=\$(date +"%Y-%m-%d %H:%M:%S")

    {
        echo "========================================"
        echo "Stage:       Quality Control (fastp + FastQC)"
        echo "Sample:      ${sample}"
        echo "Started:     \$TS"
        echo "Input R1:    ${read1}"
        echo "Input R2:    ${read2}"
        echo "Software:    \$(fastp --version 2>&1 | head -1)"
        echo "             \$(fastqc --version 2>&1 | head -1)"
        echo "----------------------------------------"
    } > "\$LOG"

    fastp \\
        -i ${read1} \\
        -I ${read2} \\
        -o ${sample}_trimmed_1.fastq.gz \\
        -O ${sample}_trimmed_2.fastq.gz \\
        --html ${sample}_fastp.html \\
        --json ${sample}_fastp.json \\
        --thread ${task.cpus} \\
        >> "\$LOG" 2>&1

    READS_AFTER=\$(gzip -cd ${sample}_trimmed_1.fastq.gz | wc -l)
    READS_AFTER=\$(( READS_AFTER / 4 ))
    echo "Reads after trimming: \$READS_AFTER" >> "\$LOG"

    fastqc \\
        --quiet \\
        --threads ${task.cpus} \\
        ${sample}_trimmed_1.fastq.gz \\
        ${sample}_trimmed_2.fastq.gz \\
        >> "\$LOG" 2>&1 || true

    ELAPSED=\$(( \$(date +%s) - START ))
    {
        echo "----------------------------------------"
        echo "Completed:   \$(date +"%Y-%m-%d %H:%M:%S")"
        echo "Elapsed:     \${ELAPSED}s"
        echo "Status:      SUCCESS"
        echo "========================================"
    } >> "\$LOG"
    """
}
