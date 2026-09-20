#!/usr/bin/env python3
"""Generate integrated HTML report from Zika-Sentinel outputs."""
import argparse, csv, html as htmlmod, sys
from pathlib import Path
from datetime import datetime, timezone
from collections import Counter

def load_tsv(path):
    p = Path(path)
    if not p.exists() or p.stat().st_size == 0:
        return []
    if p.name in ("NO_FILE", "null", "NO_META", "NO_VAR", "NO_MANIFEST", "NO_PROV"):
        return []
    with open(p, newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))

def load_kv(path):
    rows = load_tsv(path)
    return {r.get("field", r.get("metric", "")): r.get("value", "") for r in rows if r}

def esc(s):
    return htmlmod.escape(str(s) if s is not None else "")

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--qc", required=True)
    ap.add_argument("--tree", required=True)
    ap.add_argument("--provenance", required=True)
    ap.add_argument("--manifest", required=True)
    ap.add_argument("--variants", required=True)
    ap.add_argument("--variant-qc", required=True)
    ap.add_argument("--metadata", required=True)
    ap.add_argument("--metadata-issues", required=True)
    ap.add_argument("--pipeline-version", default="1.0.0")
    args = ap.parse_args()

    qc_rows = load_tsv(args.qc)
    prov = load_kv(args.provenance)
    manifest = load_tsv(args.manifest)
    var_rows = load_tsv(args.variants)
    vqc = load_kv(args.variant_qc)
    meta_rows = load_tsv(args.metadata)
    meta_issues = load_tsv(args.metadata_issues)

    fields = ["sample_id", "length", "ref_length", "n_count", "n_fraction", "status", "reason"]
    with open("zika_summary.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields, extrasaction="ignore")
        w.writeheader()
        for r in qc_rows:
            w.writerow({k: r.get(k, "") for k in fields})

    n_pass = sum(1 for r in qc_rows if r.get("status") == "PASS")
    n_warn = sum(1 for r in qc_rows if r.get("status") == "WARN")
    n_fail = sum(1 for r in qc_rows if r.get("status") == "FAIL")
    sample_ids = [r.get("sample_id", "") for r in qc_rows]

    tree_path = Path(args.tree)
    has_tree = tree_path.exists() and tree_path.stat().st_size > 20
    tree_text = tree_path.read_text()[:200] if has_tree else ""
    tree_usable = has_tree and "insufficient" not in tree_text.lower() and tree_text.strip() not in ("()", ";")

    gen_figs = [m for m in manifest if m.get("status") == "generated"]
    skip_figs = [m for m in manifest if str(m.get("status", "")).startswith("skipped")]

    meta_fields_check = ["sample_id", "accession", "collection_date", "country", "location", "host"]
    meta_completeness = {}
    for col in meta_fields_check:
        present = sum(1 for r in meta_rows if (r.get(col) or "").strip())
        meta_completeness[col] = (present, max(0, len(meta_rows) - present))

    H = []
    H.append("<!DOCTYPE html><html><head><meta charset='utf-8'><title>Zika-Sentinel Report</title>")
    H.append("<style>")
    H.append("body{font-family:system-ui,sans-serif;margin:2rem;max-width:1100px;line-height:1.45;color:#222}")
    H.append("h1{color:#1a5276}h2{color:#2874a6;margin-top:1.8rem;border-bottom:1px solid #ddd;padding-bottom:0.3rem}")
    H.append("table{border-collapse:collapse;width:100%;margin:0.6rem 0;font-size:0.92rem}")
    H.append("th,td{border:1px solid #ccc;padding:0.35rem 0.55rem;text-align:left}th{background:#eaf2f8}")
    H.append(".pass{color:#196f3d;font-weight:600}.warn{color:#b7950b;font-weight:600}.fail{color:#922b21;font-weight:600}")
    H.append(".note{background:#fef9e7;border-left:4px solid #f4d03f;padding:0.6rem 1rem;margin:1rem 0}")
    H.append(".warnbox{background:#fdedec;border-left:4px solid #e74c3c;padding:0.6rem 1rem;margin:1rem 0}")
    H.append(".infobox{background:#eaf2f8;border-left:4px solid #5dade2;padding:0.6rem 1rem;margin:1rem 0}")
    H.append("code{background:#f4f6f7;padding:0.1rem 0.3rem}</style></head><body>")
    H.append("<h1>Zika-Sentinel — Genomic Surveillance Report</h1>")
    H.append(f"<p>Generated: {datetime.now(timezone.utc).strftime('%Y-%m-%d %H:%M UTC')}</p>")
    H.append(f"<p>Pipeline: Zika-Sentinel v{esc(prov.get('pipeline_version', args.pipeline_version))}</p>")
    H.append("<div class='warnbox'><b>Dataset scope:</b> Results are from the described run. "
             "A smoke-test dataset must not be interpreted as substantive Zika epidemiological findings.</div>")

    H.append("<h2>1. Run information</h2><table>")
    for label, key in [("Pipeline","pipeline_name"),("Version","pipeline_version"),("Run date (UTC)","run_date_utc"),
                       ("Nextflow version","nextflow_version"),("Container","container"),("Input mode","input_mode")]:
        H.append(f"<tr><th>{label}</th><td>{esc(prov.get(key, 'not recorded'))}</td></tr>")
    H.append("</table>")

    H.append("<h2>2. Reference genome</h2><table>")
    for label, key in [("Accession","reference_accession"),("Length (bp)","reference_length_bp"),
                       ("MD5","reference_md5"),("Path","reference_path")]:
        H.append(f"<tr><th>{label}</th><td>{esc(prov.get(key, 'not recorded'))}</td></tr>")
    H.append("</table>")

    H.append(f"<h2>3. Input samples</h2><p>Sample count: <b>{len(sample_ids)}</b></p><ul>")
    for sid in sample_ids:
        H.append(f"<li><code>{esc(sid)}</code></li>")
    H.append("</ul>")

    H.append("<h2>4. Metadata</h2>")
    if meta_rows:
        H.append(f"<p>Validated rows: <b>{len(meta_rows)}</b></p><table><tr><th>Field</th><th>Present</th><th>Missing</th></tr>")
        for col, (pres, miss) in meta_completeness.items():
            H.append(f"<tr><td>{esc(col)}</td><td>{pres}</td><td>{miss}</td></tr>")
        H.append("</table>")
        n_dates = sum(1 for r in meta_rows if (r.get("collection_date") or "").strip())
        if n_dates < len(meta_rows):
            H.append(f"<div class='note'>Temporal analysis limited: valid collection_date for {n_dates}/{len(meta_rows)} samples.</div>")
    else:
        H.append("<div class='note'>No metadata validated for this run.</div>")

    H.append(f"<h2>5. Genome QC</h2><p>PASS: <span class='pass'>{n_pass}</span> "
             f"WARN: <span class='warn'>{n_warn}</span> FAIL: <span class='fail'>{n_fail}</span></p>")
    H.append("<table><tr><th>Sample ID</th><th>Length</th><th>N fraction</th><th>Status</th><th>Reason</th></tr>")
    for r in qc_rows:
        st = r.get("status", "")
        cls = "pass" if st == "PASS" else ("warn" if st == "WARN" else "fail")
        H.append(f"<tr><td><code>{esc(r.get('sample_id',''))}</code></td><td>{esc(r.get('length',''))}</td>"
                 f"<td>{esc(r.get('n_fraction',''))}</td><td class='{cls}'>{esc(st)}</td><td>{esc(r.get('reason',''))}</td></tr>")
    H.append("</table>")

    H.append("<h2>6. Variant summary</h2>")
    H.append("<div class='infobox'><b>Definitions:</b> consensus_vs_reference SNPs are reference-relative differences "
             "(not read-backed). Dataset-level frequency is not population prevalence.</div>")
    if vqc:
        H.append("<table>")
        for k, v in vqc.items():
            if k:
                H.append(f"<tr><th>{esc(k)}</th><td>{esc(v)}</td></tr>")
        H.append("</table>")
    if var_rows:
        per = Counter(r.get("sample_id") for r in var_rows)
        H.append("<table><tr><th>Sample ID</th><th>N SNPs</th></tr>")
        for sid, n in sorted(per.items()):
            H.append(f"<tr><td><code>{esc(sid)}</code></td><td>{n}</td></tr>")
        H.append("</table>")
    else:
        H.append("<div class='note'>No variant table available.</div>")

    H.append("<h2>7. Phylogeny</h2>")
    if tree_usable:
        H.append("<p>Maximum-likelihood tree available under <code>phylogeny/</code>.</p>")
        H.append("<div class='note'>Phylogenetic relatedness does not establish transmission chains or direction.</div>")
    else:
        H.append("<div class='note'>No usable Newick tree for this run.</div>")

    H.append("<h2>8. Genotype / clade</h2>")
    H.append("<div class='note'>No validated ZIKV genotype/clade method is implemented.</div>")

    H.append(f"<h2>9. Figures</h2><p>Generated: <b>{len(gen_figs)}</b></p>")
    if gen_figs:
        H.append("<table><tr><th>ID</th><th>Title</th><th>File</th><th>Source data</th></tr>")
        for m in gen_figs:
            H.append(f"<tr><td>{esc(m.get('figure_id',''))}</td><td>{esc(m.get('title',''))}</td>"
                     f"<td><code>{esc(m.get('filename',''))}</code></td><td><code>{esc(m.get('source_data',''))}</code></td></tr>")
        H.append("</table>")
    if skip_figs:
        H.append("<p>Skipped:</p><ul>")
        for m in skip_figs:
            H.append(f"<li><code>{esc(m.get('figure_id',''))}</code>: {esc(m.get('status',''))}</li>")
        H.append("</ul>")

    H.append("<h2>10. Limitations</h2><ul>")
    H.append("<li>Smoke-test datasets are not substantive epidemiology.</li>")
    H.append("<li>Dataset-level variant frequency is not population prevalence.</li>")
    H.append("<li>Phylogenetic clustering is not transmission evidence.</li>")
    H.append("<li>No functional or clinical mutation annotation.</li>")
    H.append("<li>Genotype/clade assignment is not implemented.</li>")
    H.append("</ul>")

    H.append("<h2>11. Reproducibility</h2>")
    H.append("<p>See <code>provenance/provenance.tsv</code> and machine-readable outputs under the results directory.</p>")
    H.append("</body></html>")

    Path("zika_report.html").write_text("\n".join(H))
    print(f"Report written: samples={len(sample_ids)} figures={len(gen_figs)} variants={len(var_rows)}", file=sys.stderr)

if __name__ == "__main__":
    main()
