/*
 * Provenance summary for Zika-Sentinel.
 */

process WRITE_PROVENANCE {

    tag "provenance"
    label 'low_compute'

    publishDir "${params.outdir}/provenance", mode: 'copy'
    publishDir "${params.outdir}/logs", mode: 'copy', pattern: "*_log.txt"

    input:
    path reference
    path qc_table
    val input_mode
    val sample_ids_str

    output:
    path "provenance.tsv", emit: provenance_tsv
    path "provenance.json", emit: provenance_json
    path "provenance_log.txt", emit: log

    script:
    """
    #!/usr/bin/env bash
    set -euo pipefail
    LOG="provenance_log.txt"
    echo "Stage: Provenance" > "\$LOG"
    echo "Started: \$(date -u +"%Y-%m-%dT%H:%M:%SZ")" >> "\$LOG"

    REF_MD5=\$(md5sum "${reference}" 2>/dev/null | awk '{print \$1}' || echo "unavailable")
    REF_LEN=\$(grep -v '^>' "${reference}" | tr -d '\\n' | wc -c)
    REF_ID=\$(head -1 "${reference}" | sed 's/^>//;s/ .*//')
    RUN_DATE=\$(date -u +"%Y-%m-%dT%H:%M:%SZ")

    capture_ver() {
        local name="\$1"; shift
        if command -v "\$1" >/dev/null 2>&1; then
            "\$@" 2>&1 | head -1 | tr '\\t' ' ' | head -c 200
        else
            echo "unavailable"
        fi
    }

    VER_FASTP=\$(capture_ver fastp fastp --version || true)
    VER_SAMTOOLS=\$(capture_ver samtools samtools --version || true)
    VER_IVAR=\$(capture_ver ivar ivar version || true)
    VER_MAFFT=\$(capture_ver mafft mafft --version || true)
    if command -v iqtree2 >/dev/null 2>&1; then
        VER_IQTREE=\$(iqtree2 --version 2>&1 | head -1)
    elif command -v iqtree >/dev/null 2>&1; then
        VER_IQTREE=\$(iqtree --version 2>&1 | head -1)
    else
        VER_IQTREE="unavailable"
    fi
    VER_PYTHON=\$(python3 --version 2>&1 || echo unavailable)
    VER_BIOPYTHON=\$(python3 -c 'import Bio; print(Bio.__version__)' 2>/dev/null || echo unavailable)
    VER_BWA=\$(bwa 2>&1 | head -1 || echo unavailable)
    VER_BCFTOOLS=\$(bcftools --version 2>&1 | head -1 || echo unavailable)
    NF_VER=\${NXF_VER:-unavailable}

    export REF_MD5 REF_LEN REF_ID RUN_DATE VER_FASTP VER_SAMTOOLS VER_IVAR VER_MAFFT VER_IQTREE VER_PYTHON VER_BIOPYTHON VER_BWA VER_BCFTOOLS NF_VER

    python3 - <<'PY'
import csv, json, os, sys
from pathlib import Path

qc_path = Path("${qc_table}")
sids = []
if qc_path.exists():
    import csv as csvmod
    with open(qc_path) as fh:
        for r in csvmod.DictReader(fh, delimiter="\\t"):
            if r.get("sample_id"):
                sids.append(r["sample_id"])

rows = [
    ("pipeline_name", "Zika-Sentinel"),
    ("pipeline_version", "${params.pipeline_version}"),
    ("run_date_utc", os.environ.get("RUN_DATE", "")),
    ("nextflow_version", os.environ.get("NF_VER", "unavailable")),
    ("container", "kizitodevbio/zika-sentinel:1.0.0"),
    ("input_mode", "${input_mode}"),
    ("reference_accession", os.environ.get("REF_ID", "")),
    ("reference_path", "${reference}"),
    ("reference_length_bp", os.environ.get("REF_LEN", "")),
    ("reference_md5", os.environ.get("REF_MD5", "")),
    ("sample_count", str(len(sids))),
    ("sample_ids", ";".join(sids)),
    ("param_min_depth", "${params.min_depth}"),
    ("param_min_freq", "${params.min_freq}"),
    ("param_min_qual", "${params.min_qual}"),
    ("param_qc_min_length_frac", "${params.qc_min_length_frac}"),
    ("param_qc_max_n_frac", "${params.qc_max_n_frac}"),
    ("param_qc_warn_n_frac", "${params.qc_warn_n_frac}"),
    ("param_iqtree_model", "${params.iqtree_model}"),
    ("param_iqtree_bootstrap", "${params.iqtree_bootstrap}"),
    ("tool_fastp", os.environ.get("VER_FASTP", "unavailable")),
    ("tool_bwa", os.environ.get("VER_BWA", "unavailable")),
    ("tool_samtools", os.environ.get("VER_SAMTOOLS", "unavailable")),
    ("tool_bcftools", os.environ.get("VER_BCFTOOLS", "unavailable")),
    ("tool_ivar", os.environ.get("VER_IVAR", "unavailable")),
    ("tool_mafft", os.environ.get("VER_MAFFT", "unavailable")),
    ("tool_iqtree", os.environ.get("VER_IQTREE", "unavailable")),
    ("tool_python", os.environ.get("VER_PYTHON", "unavailable")),
    ("tool_biopython", os.environ.get("VER_BIOPYTHON", "unavailable")),
]

with open("provenance.tsv", "w", newline="") as fh:
    w = csv.writer(fh, delimiter="\\t")
    w.writerow(["field", "value"])
    w.writerows(rows)
with open("provenance.json", "w") as fh:
    json.dump({k: v for k, v in rows}, fh, indent=2)
print("provenance written", len(sids), "samples", file=sys.stderr)
PY

    echo "Completed" >> "\$LOG"
    """
}
