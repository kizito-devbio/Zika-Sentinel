/*
 * Variant aggregation for Zika-Sentinel.
 *
 * Two sources (clearly distinguished):
 *   1. iVar variants.tsv from the raw-read pathway (read-backed)
 *   2. Consensus-vs-reference SNPs from pairwise alignment (public/consensus pathway)
 *
 * Sample identifiers are preserved exactly.
 * No invented variants. No functional/clinical annotation.
 * Dataset frequency is frequency within the analysed set, not population prevalence.
 */

process AGGREGATE_IVAR_VARIANTS {

    tag "aggregate_ivar"
    label 'low_compute'

    publishDir "${params.outdir}/variants/read_backed",
        mode: 'copy'

    publishDir "${params.outdir}/logs",
        mode: 'copy',
        pattern: "*_log.txt"

    input:
    path variant_files

    output:
    path "variants_standardized.tsv", emit: standardized
    path "variant_frequency.tsv", emit: frequency
    path "snp_matrix.tsv", emit: matrix
    path "variant_qc.tsv", emit: qc
    path "aggregate_ivar_log.txt", emit: log

    script:
    """
    #!/usr/bin/env bash
    set -euo pipefail
    LOG="aggregate_ivar_log.txt"

    {
        echo "========================================"
        echo "Stage:       Aggregate iVar variants"
        echo "Started:     \$(date +"%Y-%m-%d %H:%M:%S")"
        echo "Source:      iVar variants.tsv (read-backed)"
        echo "----------------------------------------"
    } > "\$LOG"

    python3 << 'PY'
import csv
from pathlib import Path
from collections import defaultdict, Counter

files = list(Path(".").glob("*.variants.tsv")) + list(Path(".").glob("*variants*.tsv"))
files = [f for f in files if f.name != "variants_standardized.tsv"]

print(f"Input files: {files}")

rows = []

for f in files:
    sid = f.name.replace(".variants.tsv", "").replace("_variants.tsv", "")

    with open(f, newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")

        for r in reader:
            pos = r.get("POS") or r.get("pos") or r.get("Position")
            ref = r.get("REF") or r.get("ref")
            alt = r.get("ALT") or r.get("alt")

            if not pos or not ref or not alt:
                continue

            try:
                pos = int(pos)
            except ValueError:
                continue

            depth = r.get("TOTAL_DP") or ""
            af = r.get("ALT_FREQ") or r.get("allele_frequency") or ""
            qual = r.get("ALT_QUAL") or r.get("quality") or ""

            rows.append({
                "sample_id": sid,
                "reference": r.get("REGION") or "",
                "position": pos,
                "reference_allele": ref,
                "alternate_allele": alt,
                "depth": depth,
                "allele_frequency": af,
                "quality": qual,
                "variant_type": "SNP" if (len(ref) == 1 and len(alt) == 1) else "OTHER",
                "source": "ivar"
            })

print(f"Total variant records: {len(rows)}")

fields = [
    "sample_id",
    "reference",
    "position",
    "reference_allele",
    "alternate_allele",
    "depth",
    "allele_frequency",
    "quality",
    "variant_type",
    "source"
]

with open("variants_standardized.tsv", "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=fields, delimiter="\t")
    w.writeheader()
    w.writerows(rows)

samples = sorted(set(r["sample_id"] for r in rows))
n_samples = len(samples) or 1

freq = Counter()

for r in rows:
    freq[
        (
            r["position"],
            r["reference_allele"],
            r["alternate_allele"]
        )
    ] += 1

with open("variant_frequency.tsv", "w", newline="") as fh:
    w = csv.writer(fh, delimiter="\t")

    w.writerow([
        "position",
        "reference_allele",
        "alternate_allele",
        "count_in_dataset",
        "n_samples_in_dataset",
        "frequency_in_dataset"
    ])

    for (pos, ra, aa), cnt in sorted(freq.items()):
        w.writerow([
            pos,
            ra,
            aa,
            cnt,
            n_samples,
            round(cnt / n_samples, 4)
        ])

pos_allele = defaultdict(dict)

for r in rows:
    pos_allele[r["sample_id"]][r["position"]] = r["alternate_allele"]

positions = sorted(set(r["position"] for r in rows))

with open("snp_matrix.tsv", "w", newline="") as fh:
    w = csv.writer(fh, delimiter="\t")

    w.writerow(["sample_id"] + [str(p) for p in positions])

    for sid in samples:
        w.writerow([
            sid
        ] + [
            pos_allele[sid].get(p, "0")
            for p in positions
        ])

with open("variant_qc.tsv", "w", newline="") as fh:
    w = csv.writer(fh, delimiter="\t")

    w.writerow(["metric", "value"])
    w.writerow(["n_samples", n_samples])
    w.writerow(["total_variant_records", len(rows)])
    w.writerow(["unique_positions", len(positions)])
    w.writerow(["source", "ivar"])

    for sid in samples:
        w.writerow([
            f"n_variants_{sid}",
            sum(
                1
                for r in rows
                if r["sample_id"] == sid
            )
        ])

with open("aggregate_ivar_log.txt", "a") as log:
    log.write(
        f"Records: {len(rows)}\\n"
        f"Unique positions: {len(positions)}\\n"
        f"Samples: {samples}\\n"
    )
PY

    echo "Completed: \$(date +"%Y-%m-%d %H:%M:%S")" >> "\$LOG"
    """
}


process VARIANTS_FROM_CONSENSUS {

    tag "variants_consensus"
    label 'high_compute'

    publishDir "${params.outdir}/variants/consensus_vs_reference",
        mode: 'copy'

    publishDir "${params.outdir}/logs",
        mode: 'copy',
        pattern: "*_log.txt"

    input:
    path consensus_fasta
    path reference

    output:
    path "variants_standardized.tsv", emit: standardized
    path "variant_frequency.tsv", emit: frequency
    path "snp_matrix.tsv", emit: matrix
    path "snp_allele_matrix.tsv", emit: allele_matrix
    path "variant_qc.tsv", emit: qc
    path "variants_consensus_log.txt", emit: log

    script:
    """
    #!/usr/bin/env bash
    set -euo pipefail
    LOG="variants_consensus_log.txt"

    {
        echo "========================================"
        echo "Stage:       Consensus-vs-reference SNP calling"
        echo "Started:     \$(date +"%Y-%m-%d %H:%M:%S")"
        echo "Method:      Pairwise global alignment to reference"
        echo "Note:        Not read-backed; depth/AF/quality are empty"
        echo "             Ambiguous bases (N) are skipped"
        echo "             Indels are not expanded in the SNP matrix"
        echo "----------------------------------------"
    } > "\$LOG"

    python3 << 'PY'
import csv
from collections import defaultdict, Counter
from Bio import SeqIO, Align

ref = SeqIO.read("${reference}", "fasta")
ref_seq = str(ref.seq).upper()
ref_id = ref.id

genomes = list(
    SeqIO.parse("${consensus_fasta}", "fasta")
)

aligner = Align.PairwiseAligner()
aligner.mode = "global"
aligner.match_score = 2
aligner.mismatch_score = -1
aligner.open_gap_score = -5
aligner.extend_gap_score = -0.5

all_rows = []
positions_by_sample = defaultdict(dict)

for g in genomes:

    sid = g.id
    seq = str(g.seq).upper()

    if sid == ref_id or seq == ref_seq:
        print(
            f"{sid}: reference or identical — 0 SNPs"
        )
        continue

    print(f"Aligning {sid} ...")

    aln = next(
        aligner.align(ref_seq, seq)
    )

    a_str, b_str = str(aln[0]), str(aln[1])

    ref_pos = 0
    snps = 0

    for a, b in zip(a_str, b_str):

        if a != '-':
            ref_pos += 1

        if (
            a in 'ACGT'
            and b in 'ACGT'
            and a != b
        ):

            snps += 1

            all_rows.append({
                "sample_id": sid,
                "reference": ref_id,
                "position": ref_pos,
                "reference_allele": a,
                "alternate_allele": b,
                "depth": "",
                "allele_frequency": "",
                "quality": "",
                "variant_type": "SNP",
                "source": "consensus_vs_reference"
            })

            positions_by_sample[sid][ref_pos] = b

    print(
        f"  {sid}: {snps} SNPs"
    )

fields = [
    "sample_id",
    "reference",
    "position",
    "reference_allele",
    "alternate_allele",
    "depth",
    "allele_frequency",
    "quality",
    "variant_type",
    "source"
]

with open(
    "variants_standardized.tsv",
    "w",
    newline=""
) as fh:

    w = csv.DictWriter(
        fh,
        fieldnames=fields,
        delimiter="\t"
    )

    w.writeheader()
    w.writerows(all_rows)

sample_ids = sorted(
    positions_by_sample.keys()
)

n_samples = len(sample_ids) or 1

unique_pos = sorted(
    set(
        r["position"]
        for r in all_rows
    )
)

freq = Counter(
    (
        r["position"],
        r["reference_allele"],
        r["alternate_allele"]
    )
    for r in all_rows
)

with open(
    "variant_frequency.tsv",
    "w",
    newline=""
) as fh:

    w = csv.writer(
        fh,
        delimiter="\t"
    )

    w.writerow([
        "position",
        "reference_allele",
        "alternate_allele",
        "count_in_dataset",
        "n_samples_compared",
        "frequency_in_dataset"
    ])

    for (
        pos,
        ra,
        aa
    ), cnt in sorted(freq.items()):

        w.writerow([
            pos,
            ra,
            aa,
            cnt,
            n_samples,
            round(
                cnt / n_samples,
                4
            )
        ])

with open(
    "snp_matrix.tsv",
    "w",
    newline=""
) as fh:

    w = csv.writer(
        fh,
        delimiter="\t"
    )

    w.writerow(
        ["sample_id"]
        + [
            str(p)
            for p in unique_pos
        ]
    )

    for sid in sample_ids:

        w.writerow([
            sid
        ] + [
            positions_by_sample[sid].get(
                p,
                "0"
            )
            for p in unique_pos
        ])

with open(
    "snp_allele_matrix.tsv",
    "w",
    newline=""
) as fh:

    w = csv.writer(
        fh,
        delimiter="\t"
    )

    w.writerow(
        ["sample_id"]
        + [
            str(p)
            for p in unique_pos
        ]
    )

    for sid in sample_ids:

        w.writerow([
            sid
        ] + [
            positions_by_sample[sid].get(
                p,
                "."
            )
            for p in unique_pos
        ])

with open(
    "variant_qc.tsv",
    "w",
    newline=""
) as fh:

    w = csv.writer(
        fh,
        delimiter="\t"
    )

    w.writerow([
        "metric",
        "value"
    ])

    w.writerow([
        "n_samples_compared",
        n_samples
    ])

    w.writerow([
        "total_snp_observations",
        len(all_rows)
    ])

    w.writerow([
        "unique_variant_positions",
        len(unique_pos)
    ])

    w.writerow([
        "source",
        "consensus_vs_reference_pairwise_alignment"
    ])

    w.writerow([
        "note",
        "Ambiguous bases skipped; indels not in SNP matrix; depth/AF empty"
    ])

    for sid in sample_ids:

        w.writerow([
            f"n_snps_{sid}",
            len(
                positions_by_sample[sid]
            )
        ])

print(
    f"Total SNPs: {len(all_rows)}; "
    f"unique positions: {len(unique_pos)}; "
    f"samples: {sample_ids}"
)

with open(
    "variants_consensus_log.txt",
    "a"
) as log:

    log.write(
        f"Total SNPs: {len(all_rows)}\\n"
        f"Unique positions: {len(unique_pos)}\\n"
        f"Samples: {sample_ids}\\n"
    )
PY

    echo "Completed: \$(date +"%Y-%m-%d %H:%M:%S")" >> "\$LOG"
    """
}

