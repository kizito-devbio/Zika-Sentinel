
process IVAR_CONSENSUS {

    tag "$sample"
    label 'med_compute'

    publishDir "${params.outdir}/consensus",
        mode: 'copy',
        pattern: "*.{fasta,tsv}"

    publishDir "${params.outdir}/logs",
        mode: 'copy',
        pattern: "*_consensus_log.txt"

    input:
    tuple val(sample), path(bam), path(bai)
    path reference

    output:
    tuple val(sample), path("${sample}.consensus.fasta"), emit: consensus
    path "${sample}.variants.tsv", optional: true, emit: variants
    path "${sample}.depth.tsv", emit: depth
    path "${sample}_consensus_log.txt", emit: log

    script:
    def primer_bed = params.primer_bed ? file(params.primer_bed) : null
    def use_primers = primer_bed != null
    def primer_bed_path = primer_bed ? primer_bed.toString() : ''

    """
    #!/usr/bin/env bash
    set -euo pipefail

    LOG="${sample}_consensus_log.txt"
    START=\$(date +%s)
    TS=\$(date +"%Y-%m-%d %H:%M:%S")

    {
        echo "========================================"
        echo "Stage:       iVar consensus (ZIKV)"
        echo "Sample:      ${sample}"
        echo "Started:     \$TS"
        echo "Reference:   ${reference}"
        echo "Min depth:   ${params.min_depth}"
        echo "Min freq:    ${params.min_freq}"
        echo "Min qual:    ${params.min_qual}"
        echo "Max pileup depth: 10000"

        if [ "${use_primers}" = "true" ]; then
            echo "Primer BED:  ${primer_bed_path}"
        else
            echo "Primer BED:  none (no primer trimming)"
        fi

        echo "Software:    \$(ivar version 2>&1 | head -1 || echo 'ivar')"
        echo "             \$(samtools --version 2>&1 | head -1)"
        echo "----------------------------------------"
    } > "\$LOG"

    WORK_BAM="${bam}"

    if [ "${use_primers}" = "true" ]; then
        echo "Primer trimming enabled for sample ${sample}" >> "\$LOG"

        ivar trim \\
            -i ${bam} \\
            -b ${primer_bed_path} \\
            -p ${sample}.trimmed \\
            -e \\
            >> "\$LOG" 2>&1

        samtools sort \\
            -@ ${task.cpus} \\
            -o ${sample}.trimmed.sorted.bam \\
            ${sample}.trimmed.bam

        samtools index ${sample}.trimmed.sorted.bam

        WORK_BAM="${sample}.trimmed.sorted.bam"
    else
        echo "No primer BED supplied — skipping primer trim" >> "\$LOG"
    fi

    # Generate consensus.
    # Positions below the minimum depth are represented as N by iVar.
    # Maximum pileup depth is capped at 10000 to prevent excessive memory use.
    samtools mpileup \\
        -aa \\
        -A \\
        -d 10000 \\
        -Q 0 \\
        -f ${reference} \\
        "\$WORK_BAM" \\
        | ivar consensus \\
            -p ${sample}.consensus \\
            -m ${params.min_depth} \\
            -t ${params.min_freq} \\
            -q ${params.min_qual} \\
            >> "\$LOG" 2>&1

    # Preserve the biological sample identifier as the FASTA header.
    if [ -f ${sample}.consensus.fa ]; then
        sed -i "1s/^>.*/>${sample}/" ${sample}.consensus.fa
        mv ${sample}.consensus.fa ${sample}.consensus.fasta
    else
        echo "ERROR: iVar did not produce consensus FASTA" >> "\$LOG"
        exit 1
    fi

    # Variants relative to the reference.
    samtools mpileup \\
        -aa \\
        -A \\
        -d 10000 \\
        -Q 0 \\
        -f ${reference} \\
        "\$WORK_BAM" \\
        | ivar variants \\
            -p ${sample}.variants \\
            -q ${params.min_qual} \\
            -t ${params.min_freq} \\
            -r ${reference} \\
            >> "\$LOG" 2>&1 || true

    if [ -f ${sample}.variants.tsv ]; then
        echo "Variant table written for ${sample}" >> "\$LOG"
    fi

    # Per-position depth.
    samtools depth -a "\$WORK_BAM" > ${sample}.depth.tsv

    TOTAL=\$(wc -l < ${sample}.depth.tsv || echo 0)
    COVERED=\$(awk -v d=${params.min_depth} '\$3 >= d {c++} END {print c+0}' ${sample}.depth.tsv)

    echo "Positions with depth >= ${params.min_depth}: \$COVERED / \$TOTAL" >> "\$LOG"

    # Record sequence length and N content.
    SEQ_LEN=\$(grep -v '^>' ${sample}.consensus.fasta | tr -d '\\n' | wc -c)
    N_COUNT=\$(grep -v '^>' ${sample}.consensus.fasta | tr -d '\\n' | tr -cd 'Nn' | wc -c)

    echo "Consensus length: \$SEQ_LEN" >> "\$LOG"
    echo "N bases: \$N_COUNT" >> "\$LOG"

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

