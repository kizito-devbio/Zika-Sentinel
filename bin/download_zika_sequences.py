#!/usr/bin/env python3
"""
Download Zika virus (ZIKV) complete genomes + basic metadata from NCBI Virus / Nucleotide.
Produces a multi-FASTA suitable for --sequences input to Zika-Sentinel.

Usage:
  python bin/download_zika_sequences.py --outdir test_data --max 20
  python bin/download_zika_sequences.py --accession-list accessions.txt --outdir data/
"""

import argparse
import sys
import time
from pathlib import Path
from urllib.request import urlopen
from urllib.error import URLError

try:
    from Bio import SeqIO
    from Bio.SeqRecord import SeqRecord
except ImportError:
    print("Biopython required: pip install biopython", file=sys.stderr)
    sys.exit(1)


def fetch_fasta(accession: str) -> SeqRecord | None:
    url = (
        "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
        f"?db=nuccore&id={accession}&rettype=fasta&retmode=text"
    )
    try:
        with urlopen(url, timeout=30) as resp:
            text = resp.read().decode("utf-8", errors="replace")
        from io import StringIO
        recs = list(SeqIO.parse(StringIO(text), "fasta"))
        if not recs:
            return None
        rec = recs[0]
        rec.id = accession
        rec.description = accession
        return rec
    except (URLError, Exception) as e:
        print(f"  WARNING: failed {accession}: {e}", file=sys.stderr)
        return None


def main():
    p = argparse.ArgumentParser(description="Download public ZIKV genomes for Zika-Sentinel")
    p.add_argument("--outdir", default="test_data", help="Output directory")
    p.add_argument("--max", type=int, default=10, help="Max sequences (when using default list)")
    p.add_argument("--accession-list", help="Text file with one accession per line")
    p.add_argument("--outfile", default="zika_public.fasta")
    args = p.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    # Curated set of complete / near-complete ZIKV genomes spanning Asian & African lineages
    default_ids = [
        "NC_012532.1",   # French Polynesia 2013 (Asian)
        "KU501215.1",    # Puerto Rico 2015
        "KU509998.1",    # Brazil 2015
        "KX197192.1",    # China
        "MF574561.1",    # Singapore
        "KY241758.1",    # Thailand
        "KX051563.1",    # Haiti
        "KX280026.1",    # Dominican Republic
        "AY632535.2",    # Uganda MR766 (African prototype)
        "KF268948.1",    # Central African Republic
        "KJ776791.2",    # French Polynesia
        "KU365777.1",    # Brazil
    ]

    if args.accession_list:
        ids = [ln.strip() for ln in open(args.accession_list) if ln.strip() and not ln.startswith("#")]
    else:
        ids = default_ids[: args.max]

    print(f"Fetching {len(ids)} ZIKV accessions ...", file=sys.stderr)
    records = []
    for i, acc in enumerate(ids, 1):
        print(f"  [{i}/{len(ids)}] {acc}", file=sys.stderr)
        rec = fetch_fasta(acc)
        if rec:
            records.append(rec)
        time.sleep(0.4)  # be polite to NCBI

    outpath = outdir / args.outfile
    SeqIO.write(records, outpath, "fasta")
    print(f"Wrote {len(records)} sequences → {outpath}", file=sys.stderr)


if __name__ == "__main__":
    main()
