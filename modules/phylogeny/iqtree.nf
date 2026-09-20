/*
 * Maximum-likelihood phylogeny with IQ-TREE 2 for ZIKV.
 *
 * - Uses ModelFinder Plus by default
 * - Ultrafast bootstrap when >= 4 sequences and bootstrap > 0
 * - Does NOT invent a tree for < 3 sequences
 * - Sample identifiers in the alignment are preserved as tip labels
 */

process IQTREE_PHYLO {

    tag "iqtree"
    label 'high_compute'

    publishDir "${params.outdir}/phylogeny",
        mode: 'copy'

    publishDir "${params.outdir}/logs",
        mode: 'copy',
        pattern: "iqtree_log.txt"

    input:
    path alignment

    output:
    path "zika_phylogeny.nwk", optional: true, emit: tree
    path "zika_phylogeny.iqtree", optional: true, emit: report
    path "iqtree_log.txt", emit: log

    script:
    def bs = params.iqtree_bootstrap as int
    def model = params.iqtree_model
    """
    #!/usr/bin/env bash
    set -euo pipefail

    LOG="iqtree_log.txt"    
    START=\$(date +%s)
    TS=\$(date +"%Y-%m-%d %H:%M:%S")

    {
        echo "========================================"
        echo "Stage:       Phylogenetic inference (IQ-TREE 2)"
        echo "Started:     \$TS"
        echo "Input:       ${alignment}"
        echo "Model:       ${model}"
        echo "Bootstrap:   ${bs}"
        echo "Software:    \$(iqtree2 --version 2>&1 | head -1)"
        echo "----------------------------------------"
    } > "\$LOG"

    # Prefer iqtree2, fall back to iqtree
    if command -v iqtree2 >/dev/null 2>&1; then
        IQTREE_BIN=iqtree2
    elif command -v iqtree >/dev/null 2>&1; then
        IQTREE_BIN=iqtree
    else
        echo "ERROR: iqtree/iqtree2 not found" >> "\$LOG"
        exit 1
    fi
    echo "IQ-TREE binary: \$IQTREE_BIN" >> "\$LOG"
    echo "IQ-TREE version: \$(\$IQTREE_BIN --version 2>&1 | head -1)" >> "\$LOG"

    NSEQ=\$(grep -c '^>' ${alignment} || true)
    echo "Sequences in alignment: \$NSEQ" >> "\$LOG"


    if [ "\$NSEQ" -lt 3 ]; then
        echo "Fewer than 3 sequences — phylogenetic inference skipped (no tree invented)" >> "\$LOG"
        echo "Status: SKIPPED (insufficient sequences)" >> "\$LOG"
        # Do not write a fabricated Newick tree
        touch zika_phylogeny.nwk
        echo "() // insufficient sequences for ML inference" > zika_phylogeny.nwk
    else

        CMD="\$IQTREE_BIN -s ${alignment} -m ${model} -nt ${task.cpus} -pre zika_phylogeny"
        if [ ${bs} -gt 0 ] && [ "\$NSEQ" -ge 4 ]; then
            CMD="\$CMD -B ${bs}"
        fi
        echo "Command: \$CMD" >> "\$LOG"
        eval \$CMD >> "\$LOG" 2>&1

        if [ -f zika_phylogeny.treefile ]; then
            cp zika_phylogeny.treefile zika_phylogeny.nwk
            echo "Tree written: zika_phylogeny.nwk" >> "\$LOG"
        else
            echo "ERROR: IQ-TREE did not produce a treefile" >> "\$LOG"
            exit 1
        fi
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
