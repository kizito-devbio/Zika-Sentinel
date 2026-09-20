# Zika-Sentinel v1.0

Nextflow DSL2 pipeline for **Zika virus (ZIKV)** genomic surveillance.

Refactored from [GBS-Sentinel](https://github.com/kizito-devbio/GBS-Sentinel) for viral (not bacterial) analysis.

## What it does

- **Consensus / public FASTA pathway**: merge → genome QC → MAFFT → IQ-TREE → consensus-vs-reference SNPs → optional metadata validation → figures → provenance → HTML report  
- **Raw-read pathway** (implemented): FASTP → BWA-MEM → iVar consensus/variants → same downstream stages  

Zika-Sentinel is **ZIKV-specific**. It does not implement bacterial MLST, Prokka, Panaroo, AMRFinderPlus, or GBS serotyping.

## Verified smoke test (5 public genomes)

Dataset: `test_data/zika_public.fasta` + `test_data/zika_public_metadata.tsv`  
Reference: `assets/reference/ZIKV_NC_012532.1.fasta` (NC_012532.1, 10794 bp)

Verified outside a full Docker process execution (see Limitations):

| Stage | Result |
|-------|--------|
| Genome QC | 5/5 PASS |
| Metadata validation | 5/5 matched; 0 issues |
| Consensus-vs-reference SNPs | 4717 SNP records; 1218 unique positions |
| MAFFT alignment | 5 sequences, length 10810 |
| IQ-TREE | Valid Newick; model TN+F+G4; tips = original sample IDs |
| Figures | 11 generated, 4 skipped with reasons |
| Provenance + HTML report | Generated from real tables |

**Real Newick (sample IDs preserved):**
```
(NC_012532.1:0.1590708917,KU501215.1:0.0009323112,((KU509998.1:0.0007517739,KX197192.1:0.0006527573)100:0.0008412084,MF574561.1:0.0025598138)74:0.0011397916);
```

## Requirements

- Nextflow ≥ 23 (tested with 26.04.6)
- Docker (recommended) **or** local installs of: fastp, bwa, samtools, iVar, MAFFT, IQ-TREE 2, Python 3 + Biopython/pandas/matplotlib

## Reproduce the consensus smoke test

```bash
# With Docker (when image is built)
docker build -t kizitodevbio/zika-sentinel:1.0.0 -f docker/Dockerfile .
nextflow run pipeline.nf -profile docker \
  --sequences test_data/zika_public.fasta \
  --metadata test_data/zika_public_metadata.tsv \
  --outdir results_smoke

# Local tools (no Docker)
nextflow run pipeline.nf -profile local \
  --sequences test_data/zika_public.fasta \
  --metadata test_data/zika_public_metadata.tsv \
  --outdir results_smoke \
  --max_cpus 2 --max_memory '4 GB'
```

## Input

Exactly one of:

- `--sequences` / `--curated_dir` — multi-FASTA or directory of FASTA  
- `--raw_dir` — paired-end FASTQ directory  

Optional: `--metadata` TSV/CSV with primary key `sample_id` matching sequence IDs.

Sample identifiers from the input are preserved end-to-end. The software name is never used as a sample label.

## Scientific limits (enforced)

- Dataset-level variant frequency is **not** population prevalence  
- Phylogenetic clustering is **not** confirmed transmission  
- No invented genotype/lineage system  
- No invented metadata  
- No functional/clinical mutation claims without validated annotation  
- Failed genomes are reported (PASS/WARN/FAIL), not silently dropped  

## Raw-read pathway status

Implemented (FASTP, BWA, iVar). **Not end-to-end validated** in this repository because no real ZIKV FASTQ test set is included. Do not claim raw-read validation until real FASTQ data are run.

## Repository layout

```
pipeline.nf
nextflow.config
docker/Dockerfile
modules/          # qc, mapping, consensus, merge, genome_qc, align, phylogeny,
                  # variants, metadata, visualization, reporting
bin/              # download_zika_sequences.py, generate_report.py
assets/reference/ # ZIKV_NC_012532.1.fasta
assets/schemas/   # metadata_schema.md
test_data/        # public genomes, metadata, stage verification artefacts
```

## Citation

See `CITATION.cff`.

## Licence

See `LICENSE`.
