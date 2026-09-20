/*
 * Multiple sequence alignment with MAFFT.
 *
 * Input:
 *   path sequences (multi-FASTA of consensus sequences)
 *
 * Output:
 *   path aligned.fasta
 */

process MAFFT_ALIGN {

    tag "mafft"
    label 'high_compute'

    publishDir "${params.outdir}/alignment",
        mode: 'copy'

    publishDir "${params.outdir}/logs",
        mode: 'copy',
        pattern: "mafft_log.txt"

    input:
    path sequences

    output:
    path "zika_aligned.fasta", emit: alignment
    path "mafft_log.txt", emit: log

    script:
    """
    #!/usr/bin/env bash
    set -euo pipefail

    LOG="mafft_log.txt"
    START=\$(date +%s)
    TS=\$(date +"%Y-%m-%d %H:%M:%S")

    {
        echo "========================================"
        echo "Stage:       Multiple sequence alignment (MAFFT)"
        echo "Started:     \$TS"
        echo "Input:       ${sequences}"
        echo "Software:    \$(mafft --version 2>&1 | head -1)"
        echo "Options:     ${params.mafft_opts}"
        echo "----------------------------------------"
    } > "\$LOG"

    NSEQ=\$(grep -c '^>' ${sequences} || echo 0)
    echo "Number of sequences: \$NSEQ" >> "\$LOG"

    if [ "\$NSEQ" -lt 2 ]; then
        echo "WARNING: fewer than 2 sequences; copying input as alignment" >> "\$LOG"
        cp ${sequences} zika_aligned.fasta
    else
        mafft ${params.mafft_opts} --thread ${task.cpus} ${sequences} \\
            > zika_aligned.fasta 2>> "\$LOG"
    fi

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
