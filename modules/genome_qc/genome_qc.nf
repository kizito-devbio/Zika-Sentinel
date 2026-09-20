/*
 * ZIKV genome quality control.
 *
 * Evaluates each consensus sequence against configurable thresholds:
 *   - sequence length relative to reference
 *   - fraction of ambiguous (N) bases
 *   - invalid characters
 *
 * Emits PASS / WARN / FAIL with recorded reason.
 * Failed sequences are NOT silently discarded; status is always written.
 * Sample identifiers are preserved exactly.
 */

process GENOME_QC {

    tag "genome_qc"
    label 'low_compute'

    publishDir "${params.outdir}/genome_qc",
        mode: 'copy'

    publishDir "${params.outdir}/logs",
        mode: 'copy',
        pattern: "genome_qc_log.txt"

    input:
    path consensus_fasta
    path reference

    output:
    path "genome_qc.tsv", emit: qc_table
    path "pass_sequences.fasta", emit: pass_fasta
    path "genome_qc_log.txt", emit: log
    path "all_sequences_with_status.fasta", emit: all_fasta

    script:
    """
    #!/usr/bin/env bash
    set -euo pipefail

    LOG="genome_qc_log.txt"
    START=\$(date +%s)

    {
        echo "========================================"
        echo "Stage:       ZIKV genome QC"
        echo "Started:     \$(date +"%Y-%m-%d %H:%M:%S")"
        echo "Reference:   ${reference}"
        echo "Min length fraction: ${params.qc_min_length_frac}"
        echo "Max N fraction:      ${params.qc_max_n_frac}"
        echo "Warn N fraction:     ${params.qc_warn_n_frac}"
        echo "----------------------------------------"
    } > "\$LOG"

    python3 << 'PY'
import sys
from pathlib import Path
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

ref_path = "${reference}"
cons_path = "${consensus_fasta}"
min_len_frac = float("${params.qc_min_length_frac}")
max_n_frac   = float("${params.qc_max_n_frac}")
warn_n_frac  = float("${params.qc_warn_n_frac}")

ref = SeqIO.read(ref_path, "fasta")
ref_len = len(ref.seq)
print(f"Reference length: {ref_len}", file=sys.stderr)

min_len = int(ref_len * min_len_frac)

rows = []
pass_recs = []
all_recs = []

seen_ids = set()
for rec in SeqIO.parse(cons_path, "fasta"):
    sid = rec.id
    if sid in seen_ids:
        status = "FAIL"
        reason = "duplicate_sample_id"
    else:
        seen_ids.add(sid)
        seq = str(rec.seq).upper()
        length = len(seq)
        # invalid characters (non-IUPAC)
        valid = set("ACGTNRYSWKMBDHV")
        invalid = sum(1 for c in seq if c not in valid)
        n_count = seq.count("N")
        n_frac = n_count / length if length > 0 else 1.0

        if length == 0:
            status, reason = "FAIL", "empty_sequence"
        elif invalid > 0:
            status, reason = "FAIL", f"invalid_characters:{invalid}"
        elif length < min_len:
            status, reason = "FAIL", f"too_short:{length}<{min_len}"
        elif n_frac > max_n_frac:
            status, reason = "FAIL", f"excess_N:{n_frac:.3f}>{max_n_frac}"
        elif n_frac > warn_n_frac:
            status, reason = "WARN", f"elevated_N:{n_frac:.3f}"
        else:
            status, reason = "PASS", "ok"

    rows.append({
        "sample_id": sid,
        "length": length if 'length' in dir() else 0,
        "ref_length": ref_len,
        "n_count": n_count if 'n_count' in dir() else 0,
        "n_fraction": round(n_frac, 4) if 'n_frac' in dir() else 1.0,
        "status": status,
        "reason": reason
    })

    # Preserve original identifier
    new_rec = SeqRecord(Seq(str(rec.seq)), id=sid, description=f"status={status} reason={reason}")
    all_recs.append(new_rec)
    if status in ("PASS", "WARN"):
        pass_recs.append(SeqRecord(Seq(str(rec.seq)), id=sid, description=""))

# Write QC table
with open("genome_qc.tsv", "w") as fh:
    fh.write("sample_id\\tlength\\tref_length\\tn_count\\tn_fraction\\tstatus\\treason\\n")
    for r in rows:
        fh.write(f"{r['sample_id']}\\t{r['length']}\\t{r['ref_length']}\\t{r['n_count']}\\t{r['n_fraction']}\\t{r['status']}\\t{r['reason']}\\n")

SeqIO.write(pass_recs, "pass_sequences.fasta", "fasta")
SeqIO.write(all_recs, "all_sequences_with_status.fasta", "fasta")

n_pass = sum(1 for r in rows if r["status"] == "PASS")
n_warn = sum(1 for r in rows if r["status"] == "WARN")
n_fail = sum(1 for r in rows if r["status"] == "FAIL")
print(f"PASS={n_pass} WARN={n_warn} FAIL={n_fail}", file=sys.stderr)

with open("genome_qc_log.txt", "a") as log:
    log.write(f"PASS={n_pass} WARN={n_warn} FAIL={n_fail}\\n")
    for r in rows:
        log.write(f"  {r['sample_id']}: {r['status']} ({r['reason']})\\n")
PY

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
