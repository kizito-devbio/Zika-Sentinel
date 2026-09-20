/*
 * Collect per-sample consensus FASTA files into a single multi-FASTA
 * while preserving every original sample identifier as the FASTA header.
 *
 * No renaming. Headers remain exactly the sample IDs supplied upstream.
 */

process MERGE_CONSENSUS {

    tag "merge_consensus"
    label 'low_compute'

    publishDir "${params.outdir}/consensus",
        mode: 'copy',
        pattern: "all_consensus.fasta"

    publishDir "${params.outdir}/logs",
        mode: 'copy',
        pattern: "merge_consensus_log.txt"

    input:
    path consensus_files

    output:
    path "all_consensus.fasta", emit: merged
    path "merge_consensus_log.txt", emit: log

    script:
    """
    #!/usr/bin/env bash
    set -euo pipefail

    LOG="merge_consensus_log.txt"
    START=\$(date +%s)

    {
        echo "========================================"
        echo "Stage:       Merge consensus FASTA"
        echo "Started:     \$(date +"%Y-%m-%d %H:%M:%S")"
        echo "----------------------------------------"
    } > "\$LOG"

    : > all_consensus.fasta

    for f in ${consensus_files}; do
        if [ -s "\$f" ]; then
            cat "\$f" >> all_consensus.fasta
            echo "" >> all_consensus.fasta
            HDR=\$(head -1 "\$f" | sed 's/^>//')
            echo "Included: \$HDR  (from \$f)" >> "\$LOG"
        else
            echo "WARNING: empty or missing file \$f — skipped" >> "\$LOG"
        fi
    done

    NSEQ=\$(grep -c '^>' all_consensus.fasta || echo 0)
    echo "Total sequences in merged FASTA: \$NSEQ" >> "\$LOG"

    if [ "\$NSEQ" -eq 0 ]; then
        echo "ERROR: no consensus sequences to merge" >> "\$LOG"
        exit 1
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
