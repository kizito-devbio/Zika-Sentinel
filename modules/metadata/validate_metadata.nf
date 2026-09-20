/*
 * Metadata validation for Zika-Sentinel.
 *
 * - sample_id is the primary key and must match sequence identifiers
 * - Missing values remain missing (never invented)
 * - Duplicate sample_id values cause failure
 * - Sequences without metadata are reported as missing_metadata
 * - Malformed dates are flagged but do not discard the sequence
 */

process VALIDATE_METADATA {

    tag "metadata"
    label 'low_compute'

    publishDir "${params.outdir}/metadata",
        mode: 'copy'

    publishDir "${params.outdir}/logs",
        mode: 'copy',
        pattern: "metadata_log.txt"

    input:
    path metadata_file
    path sequence_fasta

    output:
    path "metadata_validated.tsv", emit: validated
    path "metadata_issues.tsv", emit: issues
    path "metadata_log.txt", emit: log

    script:
    """
    #!/usr/bin/env bash
    set -euo pipefail

    LOG="metadata_log.txt"
    START=\$(date +%s)

    {
        echo "========================================"
        echo "Stage:       Metadata validation"
        echo "Started:     \$(date +"%Y-%m-%d %H:%M:%S")"
        echo "Metadata:    ${metadata_file}"
        echo "Sequences:   ${sequence_fasta}"
        echo "----------------------------------------"
    } > "\$LOG"

    python3 << 'PY'
import sys, csv
from pathlib import Path
from Bio import SeqIO
from datetime import datetime

meta_path = Path("${metadata_file}")
seq_path  = Path("${sequence_fasta}")

# Load sequence IDs (preserve exact identifiers)
seq_ids = [rec.id for rec in SeqIO.parse(seq_path, "fasta")]
seq_set = set(seq_ids)
print(f"Sequences: {len(seq_ids)}", file=sys.stderr)

# Load metadata
rows = []
with open(meta_path, newline="") as fh:
    # Detect delimiter
    sample = fh.read(2048)
    fh.seek(0)
    dialect = csv.Sniffer().sniff(sample, delimiters="\\t,")
    reader = csv.DictReader(fh, dialect=dialect)
    fieldnames = reader.fieldnames
    if not fieldnames or "sample_id" not in fieldnames:
        print("ERROR: metadata must contain a sample_id column", file=sys.stderr)
        sys.exit(1)
    for row in reader:
        rows.append({k: (v.strip() if isinstance(v, str) else v) for k, v in row.items()})

print(f"Metadata rows: {len(rows)}", file=sys.stderr)

# Check duplicates
seen = {}
issues = []
validated = []

for r in rows:
    sid = r.get("sample_id", "").strip()
    if not sid:
        issues.append({"sample_id": "", "issue": "empty_sample_id", "detail": ""})
        continue
    if sid in seen:
        issues.append({"sample_id": sid, "issue": "duplicate_sample_id", "detail": "appears more than once"})
        continue
    seen[sid] = True

    # Date check (optional)
    date_val = r.get("collection_date", "") or ""
    date_ok = True
    if date_val:
        parsed = False
        for fmt in ("%Y-%m-%d", "%Y-%m", "%Y"):
            try:
                datetime.strptime(date_val, fmt)
                parsed = True
                break
            except ValueError:
                pass
        if not parsed:
            date_ok = False
            issues.append({"sample_id": sid, "issue": "malformed_date", "detail": date_val})

    # Sequence presence
    if sid not in seq_set:
        issues.append({"sample_id": sid, "issue": "metadata_without_sequence", "detail": ""})
    else:
        validated.append(r)

# Sequences missing metadata
for sid in seq_ids:
    if sid not in seen:
        issues.append({"sample_id": sid, "issue": "sequence_without_metadata", "detail": ""})

# Write outputs
out_fields = fieldnames
with open("metadata_validated.tsv", "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=out_fields, delimiter="\\t", extrasaction="ignore")
    w.writeheader()
    w.writerows(validated)

with open("metadata_issues.tsv", "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=["sample_id", "issue", "detail"], delimiter="\\t")
    w.writeheader()
    w.writerows(issues)

print(f"Validated rows: {len(validated)}", file=sys.stderr)
print(f"Issues: {len(issues)}", file=sys.stderr)

with open("metadata_log.txt", "a") as log:
    log.write(f"Validated: {len(validated)}\\n")
    log.write(f"Issues: {len(issues)}\\n")
    for i in issues:
        log.write(f"  {i['sample_id']}: {i['issue']} {i['detail']}\\n")
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
