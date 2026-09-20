/*
 * Reference mapping with BWA-MEM + samtools sort/index.
 *
 * Input:
 *   tuple val(sample), path(r1), path(r2)
 *   path reference
 *
 * Output:
 *   tuple val(sample), path(bam), path(bai)
 */

process MAP_BWA {

    tag "$sample"
    label 'high_compute'

    publishDir "${params.outdir}/mapping",
        mode: 'copy',
        pattern: "*.{bam,bai}"

    publishDir "${params.outdir}/logs",
        mode: 'copy',
        pattern: "*_map_log.txt"

    input:
    tuple val(sample), path(r1), path(r2)
    path reference

    output:
    tuple val(sample), path("${sample}.sorted.bam"), path("${sample}.sorted.bam.bai"), emit: bam
    path "${sample}_map_log.txt", emit: log

    script:
    """
    #!/usr/bin/env bash
    set -euo pipefail

    LOG="${sample}_map_log.txt"
    START=\$(date +%s)
    TS=\$(date +"%Y-%m-%d %H:%M:%S")

    {
        echo "========================================"
        echo "Stage:       BWA-MEM mapping + sort/index"
        echo "Sample:      ${sample}"
        echo "Started:     \$TS"
        echo "Reference:   ${reference}"
        echo "Software:    bwa + samtools"
        echo "----------------------------------------"
    } > "\$LOG"

    # Index reference if needed
    if [ ! -f "${reference}.bwt" ]; then
        bwa index ${reference} >> "\$LOG" 2>&1
    fi

    bwa mem -t ${task.cpus} ${reference} ${r1} ${r2} 2>> "\$LOG" \\
        | samtools sort -@ ${task.cpus} -o ${sample}.sorted.bam -

    samtools index -@ ${task.cpus} ${sample}.sorted.bam

    samtools flagstat ${sample}.sorted.bam >> "\$LOG"
    samtools coverage ${sample}.sorted.bam >> "\$LOG" 2>&1 || true

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
