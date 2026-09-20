#!/usr/bin/env nextflow

/*
 * Zika-Sentinel v1.0
 *
 * Modular Nextflow DSL2 workflow for genomic surveillance of Zika virus (ZIKV).
 * Refactored from GBS-Sentinel for viral pathways.
 *
 * Sample identifiers from the input dataset are preserved end-to-end.
 * The software name "Zika-Sentinel" is never used as a sample label.
 *
 * Input pathways (exactly one required):
 *   --raw_dir       Directory of paired-end Illumina FASTQ
 *   --sequences     Directory or multi-FASTA of consensus / public sequences
 *                   (alias: --curated_dir)
 */

nextflow.enable.dsl = 2

// ============================================================================
// MODULE IMPORTS
// ============================================================================

include { FASTP_QC }               from './modules/qc.nf'
include { MAP_BWA }                from './modules/mapping/map_bwa.nf'
include { IVAR_CONSENSUS }         from './modules/consensus/ivar_consensus.nf'
include { MERGE_CONSENSUS }        from './modules/merge/merge_consensus.nf'
include { GENOME_QC }              from './modules/genome_qc/genome_qc.nf'
include { MAFFT_ALIGN }             from './modules/align/mafft.nf'
include { IQTREE_PHYLO }           from './modules/phylogeny/iqtree.nf'
include { FETCH_NCBI_METADATA }    from './modules/metadata/fetch_ncbi_metadata.nf'
include { VALIDATE_METADATA }      from './modules/metadata/validate_metadata.nf'
include { AGGREGATE_IVAR_VARIANTS;
          VARIANTS_FROM_CONSENSUS } from './modules/variants/aggregate_variants.nf'
include { GENERATE_FIGURES }       from './modules/visualization/generate_figures.nf'
include { WRITE_PROVENANCE }       from './modules/reporting/provenance.nf'
include { GENERATE_REPORT }        from './modules/reporting/report.nf'


// ============================================================================
// HELP
// ============================================================================

def printHelp() {

    log.info """
    ======================================================================
    Zika-Sentinel v${params.pipeline_version} — Zika virus genomic surveillance
    ======================================================================

    USAGE

      nextflow run pipeline.nf -profile docker [options]

    REQUIRED (exactly one input pathway)

      --raw_dir <dir>       Paired-end FASTQ directory

      --sequences <path>    Consensus / public sequence FASTA(s) or directory
                            (alias: --curated_dir)

    OPTIONAL

      --metadata <file>     Sample metadata TSV/CSV (sample_id primary key).
                            For public FASTA without --metadata, NCBI metadata
                            is retrieved automatically from sequence accessions.

      --outdir <dir>        Output directory (default: ./results)

      --reference <fasta>   ZIKV reference genome

      --primer_bed <bed>    Primer BED for amplicon schemes only

      --min_depth <int>     Min depth for consensus (default: 10)

      --min_freq <float>    Majority frequency (default: 0.5)

      --qc_min_length_frac  Min genome length as fraction of reference

      --qc_max_n_frac       FAIL threshold for N fraction

      --qc_warn_n_frac      WARN threshold for N fraction

      --max_cpus <int>      CPU ceiling

      --max_memory <mem>    Memory ceiling

      --help                Print this help

    PROFILES

      -profile docker | singularity | conda | cluster | test | local

    EXAMPLES

      nextflow run pipeline.nf -profile docker \
          --sequences test_data/zika_public.fasta \
          --outdir results_public

      nextflow run pipeline.nf -profile docker \
          --raw_dir path/to/fastq \
          --reference assets/reference/ZIKV_NC_012532.1.fasta \
          --outdir results_raw

    ======================================================================
    """
}


// ============================================================================
// PARAMETER VALIDATION
// ============================================================================

def validateParams() {

    if (params.help) {
        printHelp()
        exit 0
    }

    def has_raw = params.raw_dir && params.raw_dir != false

    def has_seq = (
        (params.sequences && params.sequences != false) ||
        (params.curated_dir && params.curated_dir != false)
    )

    if (!has_raw && !has_seq) {
        error "No input supplied. Provide exactly one of --raw_dir or --sequences (alias --curated_dir)."
    }

    if (has_raw && has_seq) {
        error "Provide only one input pathway (--raw_dir OR --sequences)."
    }

    if (has_raw && !file(params.raw_dir).exists()) {
        error "Raw-read directory does not exist: ${params.raw_dir}"
    }

    if (has_seq) {
        def p = params.sequences ?: params.curated_dir

        if (!file(p).exists()) {
            error "Sequences path does not exist: ${p}"
        }
    }

    if (!file(params.reference).exists()) {
        error "Reference genome does not exist: ${params.reference}"
    }
}


// ============================================================================
// CHANNEL HELPERS
// ============================================================================

def createRawChannel(dir) {

    /*
     * Discover paired-end FASTQ files in a generic way.
     *
     * Supported naming:
     *   sample_1.fastq.gz / sample_2.fastq.gz
     *   sample_1.fastq    / sample_2.fastq
     *   sample_1.fq.gz    / sample_2.fq.gz
     *   sample_R1.fastq.gz / sample_R2.fastq.gz
     *   sample_R1.fastq    / sample_R2.fastq
     *   sample_R1.fq.gz    / sample_R2.fq.gz
     *   sample_R1.fq      / sample_R2.fq
     *
     * The function validates that every sample has exactly one R1
     * and one R2 before creating the Nextflow channel.
     */

    def supported = files("${dir}/*", checkIfExists: true)
        .findAll { f ->
            f.name =~ /.+(_R?[12])\.(fastq|fq)(\.gz)?$/
        }

    if (!supported) {
        error "No supported paired-end FASTQ files found in: ${dir}"
    }

    def pairs = [:].withDefault { [:] }

    supported.each { f ->

        def m = (f.name =~ /^(.+)(_R?([12]))\.(fastq|fq)(\.gz)?$/)

        if (!m.matches()) {
            return
        }

        def sample = m.group(1)
        def read = m.group(3) == '1' ? 'R1' : 'R2'

        if (pairs[sample][read]) {

            error """
Duplicate ${read} detected for sample '${sample}'.

Files:

  ${pairs[sample][read].name}
  ${f.name}

Each sample must have exactly one R1 and one R2 file.
"""
        }

        pairs[sample][read] = f
    }

    def incomplete = pairs.findAll { sample, reads ->
        !reads.R1 || !reads.R2
    }

    if (!incomplete.isEmpty()) {

        def details = incomplete.collect { sample, reads ->

            def r1 = reads.R1 ? reads.R1.name : "MISSING"
            def r2 = reads.R2 ? reads.R2.name : "MISSING"

            "  ${sample}: R1=${r1}, R2=${r2}"

        }.join("\n")

        error """
Incomplete paired-end FASTQ data detected.

Every sample must have both R1 and R2.

${details}
"""
    }

    def tuples = pairs.keySet().sort().collect { sample ->

        tuple(
            sample,
            pairs[sample].R1,
            pairs[sample].R2
        )
    }

    Channel.fromList(tuples)
}


def createSequenceChannel(path) {

    def f = file(path)

    if (f.isDirectory()) {

        Channel
            .fromPath("${path}/*.{fa,fna,fasta}")
            .ifEmpty {
                error "No FASTA files found in ${path}"
            }

    } else {

        Channel.fromPath(path)
    }
}


// ============================================================================
// MAIN WORKFLOW
// ============================================================================

workflow {

    validateParams()

    log.info """
    ╔══════════════════════════════════════════════════════════════════════╗
    ║                    Zika-Sentinel v${params.pipeline_version}          ║
    ║          Zika virus (ZIKV) genomic surveillance                      ║
    ╚══════════════════════════════════════════════════════════════════════╝
    """

    ch_reference = Channel.value(file(params.reference))

    // Primer BED only when a real file is supplied

    ch_primer = (params.primer_bed && params.primer_bed != false)
        ? Channel.fromPath(params.primer_bed, checkIfExists: true)
        : Channel.value(file("NO_FILE"))


    // ----------------------------------------------------------------
    // PATHWAY A: Raw reads → QC → map → iVar consensus + variants
    // ----------------------------------------------------------------

    if (params.raw_dir && params.raw_dir != false) {

        ch_raw = createRawChannel(params.raw_dir)

        FASTP_QC(ch_raw)

        MAP_BWA(
            FASTP_QC.out.trimmed,
            ch_reference
        )

        IVAR_CONSENSUS(
            MAP_BWA.out.bam,
            ch_reference
        )

        /*
         * Raw-read variant pathway.
         *
         * IVAR_CONSENSUS emits:
         *   - consensus FASTA
         *   - read-backed iVar variants.tsv
         *
         * The consensus files continue through the shared downstream
         * genome-QC / alignment / phylogeny pathway.
         *
         * The iVar variants are separately aggregated here so they are
         * retained as read-backed evidence and are not confused with
         * consensus-vs-reference SNPs.
         */

        ch_ivar_variant_files = IVAR_CONSENSUS.out.variants.collect()

        AGGREGATE_IVAR_VARIANTS(
            ch_ivar_variant_files
        )

        // Collect per-sample consensus FASTA files
        // (sample ID already preserved in FASTA headers)

        ch_cons_files = IVAR_CONSENSUS.out.consensus
            .map { sample, fasta -> fasta }
            .collect()

        MERGE_CONSENSUS(ch_cons_files)

        ch_merged = MERGE_CONSENSUS.out.merged
    }


    // ----------------------------------------------------------------
    // PATHWAY B: Pre-computed / public consensus sequences
    // ----------------------------------------------------------------

    else {

        def seq_path = params.sequences ?: params.curated_dir

        ch_seq_files = createSequenceChannel(seq_path)
            .collect()

        MERGE_CONSENSUS(ch_seq_files)

        ch_merged = MERGE_CONSENSUS.out.merged
    }


    // ----------------------------------------------------------------
    // Shared: genome QC → MSA → phylogeny
    // ----------------------------------------------------------------

    GENOME_QC(
        ch_merged,
        ch_reference
    )

    // Only PASS + WARN sequences proceed to alignment

    MAFFT_ALIGN(
        GENOME_QC.out.pass_fasta
    )

    IQTREE_PHYLO(
        MAFFT_ALIGN.out.alignment
    )


    /*
     * Consensus-vs-reference SNPs.
     *
     * This is intentionally separate from the read-backed iVar variant
     * aggregation above.
     *
     * source = consensus_vs_reference
     * depth / allele frequency / quality = not available from consensus
     */

    VARIANTS_FROM_CONSENSUS(
        GENOME_QC.out.pass_fasta,
        ch_reference
    )


    // ----------------------------------------------------------------
    // Metadata acquisition and validation
    //
    // Priority:
    //   1. User-supplied --metadata
    //   2. Automatic NCBI retrieval for public FASTA input
    //   3. No metadata for raw-read input without --metadata
    //
    // All metadata paths pass through VALIDATE_METADATA.
    // ----------------------------------------------------------------

    if (params.metadata && params.metadata != false) {

        /*
         * Explicit metadata always takes precedence.
         * This supports private datasets and custom metadata sources.
         */

        ch_meta = Channel.fromPath(
            params.metadata,
            checkIfExists: true
        )

        VALIDATE_METADATA(
            ch_meta,
            GENOME_QC.out.all_fasta
        )

        ch_meta_val = VALIDATE_METADATA.out.validated
        ch_meta_iss = VALIDATE_METADATA.out.issues
        ch_meta_for_fig = VALIDATE_METADATA.out.validated

    } else if (!params.raw_dir || params.raw_dir == false) {

        /*
         * Public / pre-computed FASTA pathway.
         *
         * GENOME_QC preserves the original FASTA identifiers in
         * all_sequences_with_status.fasta. Those identifiers are
         * used as NCBI accession/sample IDs for metadata retrieval.
         */

        FETCH_NCBI_METADATA(
            GENOME_QC.out.all_fasta
        )

        VALIDATE_METADATA(
            FETCH_NCBI_METADATA.out.metadata,
            GENOME_QC.out.all_fasta
        )

        ch_meta_val = VALIDATE_METADATA.out.validated
        ch_meta_iss = VALIDATE_METADATA.out.issues
        ch_meta_for_fig = VALIDATE_METADATA.out.validated

    } else {

        /*
         * Raw-read pathway without explicit metadata.
         *
         * FASTQ sample names are not assumed to be NCBI accessions,
         * so no automatic NCBI query is performed.
         */

        ch_meta_val = Channel.value([])
        ch_meta_iss = Channel.value([])
        ch_meta_for_fig = Channel.value([])
    }


    // ----------------------------------------------------------------
    // Figures
    //
    // Keep figures based on consensus-vs-reference variants.
    // Read-backed iVar variants are retained separately and are not
    // silently substituted into this existing analysis.
    // ----------------------------------------------------------------

    GENERATE_FIGURES(
        GENOME_QC.out.qc_table,
        IQTREE_PHYLO.out.tree,
        ch_meta_for_fig,
        VARIANTS_FROM_CONSENSUS.out.standardized,
        VARIANTS_FROM_CONSENSUS.out.frequency,
        VARIANTS_FROM_CONSENSUS.out.matrix
    )


    // ----------------------------------------------------------------
    // Provenance
    // ----------------------------------------------------------------

    input_mode_ch = Channel.value(
        (params.raw_dir && params.raw_dir != false)
            ? "raw_reads"
            : "consensus_fasta"
    )

    WRITE_PROVENANCE(
        ch_reference,
        GENOME_QC.out.qc_table,
        input_mode_ch,
        Channel.value("NA")
    )


    // ----------------------------------------------------------------
    // Report
    // ----------------------------------------------------------------

    GENERATE_REPORT(
        GENOME_QC.out.qc_table,
        IQTREE_PHYLO.out.tree,
        WRITE_PROVENANCE.out.provenance_tsv,
        GENERATE_FIGURES.out.manifest,
        VARIANTS_FROM_CONSENSUS.out.standardized,
        VARIANTS_FROM_CONSENSUS.out.qc,
        ch_meta_val,
        ch_meta_iss
    )
}


