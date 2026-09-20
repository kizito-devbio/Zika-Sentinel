#!/usr/bin/env python3

import argparse
import csv
import json
import sys
import time
from datetime import datetime
from io import StringIO
from pathlib import Path
from urllib.parse import urlencode
from urllib.request import Request, urlopen
from urllib.error import HTTPError, URLError

from Bio import SeqIO


EUTILS_BASE = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"

FIELDS = [
    "sample_id",
    "accession",
    "collection_date",
    "country",
    "location",
    "host",
    "sequencing_platform",
    "study_accession",
    "source",
]


def request_ncbi(endpoint, params):
    url = f"{EUTILS_BASE}/{endpoint}?{urlencode(params)}"

    request = Request(
        url,
        headers={
            "User-Agent": "Zika-Sentinel/1.1"
        },
    )

    for attempt in range(3):
        try:
            with urlopen(request, timeout=60) as response:
                return response.read().decode(
                    "utf-8",
                    errors="replace",
                )

        except (HTTPError, URLError, TimeoutError) as exc:
            if attempt == 2:
                raise RuntimeError(
                    f"NCBI request failed after 3 attempts: {exc}"
                ) from exc

            time.sleep(2 ** attempt)


def fasta_accessions(fasta):
    accessions = []

    for record in SeqIO.parse(fasta, "fasta"):
        accession = record.id.strip()

        if accession:
            accessions.append(accession)

    if len(accessions) != len(set(accessions)):
        duplicates = sorted(
            accession
            for accession in set(accessions)
            if accessions.count(accession) > 1
        )

        raise RuntimeError(
            "Duplicate FASTA identifiers detected: "
            + ", ".join(duplicates)
        )

    return accessions


def esearch(accessions):
    term = " OR ".join(
        f'"{accession}"[Accession]'
        for accession in accessions
    )

    text = request_ncbi(
        "esearch.fcgi",
        {
            "db": "nuccore",
            "term": term,
            "retmode": "json",
            "retmax": len(accessions),
        },
    )

    data = json.loads(text)

    return data["esearchresult"]["idlist"]


def efetch(uids):
    if not uids:
        return ""

    return request_ncbi(
        "efetch.fcgi",
        {
            "db": "nuccore",
            "id": ",".join(uids),
            "rettype": "gb",
            "retmode": "text",
        },
    )


def qualifier_values(record):
    values = {}

    for feature in record.features:
        for key, raw_values in feature.qualifiers.items():
            if key not in values:
                values[key] = []

            values[key].extend(
                str(value).strip()
                for value in raw_values
                if str(value).strip()
            )

    return values


def first(values, key):
    value = values.get(key, [])

    if not value:
        return ""

    return value[0]


def normalize_collection_date(value):
    """
    Convert common NCBI collection-date formats to the
    Zika-Sentinel metadata schema:

        YYYY-MM-DD
        YYYY-MM
        YYYY
    """

    if not value:
        return ""

    value = value.strip()

    formats = (
        "%d-%b-%Y",
        "%d-%B-%Y",
        "%b-%Y",
        "%B-%Y",
        "%Y-%m-%d",
        "%Y-%m",
        "%Y",
    )

    for fmt in formats:
        try:
            parsed = datetime.strptime(value, fmt)

            if fmt == "%Y":
                return parsed.strftime("%Y")

            if fmt in ("%b-%Y", "%B-%Y"):
                return parsed.strftime("%Y-%m")

            return parsed.strftime("%Y-%m-%d")

        except ValueError:
            continue

    # Never invent or alter an unrecognized date.
    # Leave it for the existing metadata validator to flag.
    return value


def parse_geo_location(value):
    """
    Normalize NCBI geographic metadata.

    Examples:

        French Polynesia
            -> country = French Polynesia
               location = ""

        Thailand:Bangkok
            -> country = Thailand
               location = Bangkok

        United States:Florida:Miami
            -> country = United States
               location = Florida:Miami
    """

    if not value:
        return "", ""

    value = value.strip()

    if ":" in value:
        country, location = value.split(":", 1)
        return country.strip(), location.strip()

    return value, ""


def extract_study_accession(record, qualifiers):
    """
    Retrieve BioProject information when it is exposed by
    the GenBank record.
    """

    # Feature-level db_xref values.
    for xref in qualifiers.get("db_xref", []):
        if xref.startswith("BioProject:"):
            return xref.split(":", 1)[1].strip()

    # Record-level dbxrefs.
    for xref in getattr(record, "dbxrefs", []):
        if xref.startswith("BioProject:"):
            return xref.split(":", 1)[1].strip()

    # NCBI records may expose the project in a structured
    # annotation value.
    for key in ("bioproject", "BioProject"):
        value = record.annotations.get(key)

        if value:
            if isinstance(value, (list, tuple)):
                return str(value[0]).strip()

            return str(value).strip()

    return ""


def extract_sequencing_platform(record, qualifiers):
    """
    Retrieve sequencing platform information when NCBI
    actually exposes it in the record.

    Missing platform information remains blank.
    """

    candidate_keys = (
        "sequencing_platform",
        "sequencing_technology",
        "platform",
        "technology",
    )

    for key in candidate_keys:
        value = first(qualifiers, key)

        if value:
            return value

    annotation_keys = (
        "sequencing_platform",
        "sequencing_technology",
        "platform",
        "technology",
    )

    for key in annotation_keys:
        value = record.annotations.get(key)

        if value:
            if isinstance(value, (list, tuple)):
                return str(value[0]).strip()

            return str(value).strip()

    return ""


def parse_records(genbank_text):
    records = {}

    for record in SeqIO.parse(
        StringIO(genbank_text),
        "genbank",
    ):
        accession = ""

        accessions = record.annotations.get(
            "accessions",
            [],
        )

        if accessions:
            accession = accessions[0]

        if not accession:
            accession = record.id

        accession = accession.strip()

        qualifiers = qualifier_values(record)

        # NCBI commonly uses geo_loc_name for geographic
        # information. Fall back to country when necessary.
        geo_value = first(
            qualifiers,
            "geo_loc_name",
        )

        if not geo_value:
            geo_value = first(
                qualifiers,
                "country",
            )

        country, location = parse_geo_location(
            geo_value
        )

        collection_date = normalize_collection_date(
            first(
                qualifiers,
                "collection_date",
            )
        )

        study_accession = extract_study_accession(
            record,
            qualifiers,
        )

        sequencing_platform = (
            extract_sequencing_platform(
                record,
                qualifiers,
            )
        )

        records[accession] = {
            "accession": accession,
            "collection_date": collection_date,
            "country": country,
            "location": location,
            "host": first(
                qualifiers,
                "host",
            ),
            "sequencing_platform": sequencing_platform,
            "study_accession": study_accession,
        }

    return records


def find_record(records, accession):
    if accession in records:
        return records[accession]

    versionless = accession.split(
        ".",
        1,
    )[0]

    for key, record in records.items():
        if key.split(".", 1)[0] == versionless:
            return record

    return None


def write_metadata(accessions, records, output):
    with open(
        output,
        "w",
        newline="",
        encoding="utf-8",
    ) as handle:

        writer = csv.DictWriter(
            handle,
            fieldnames=FIELDS,
            delimiter="\t",
        )

        writer.writeheader()

        for sample_id in accessions:
            record = find_record(
                records,
                sample_id,
            )

            if record is None:
                record = {
                    "accession": sample_id,
                    "collection_date": "",
                    "country": "",
                    "location": "",
                    "host": "",
                    "sequencing_platform": "",
                    "study_accession": "",
                }

            writer.writerow(
                {
                    # Preserve the FASTA identifier exactly.
                    "sample_id": sample_id,
                    "accession": sample_id,
                    "collection_date": record[
                        "collection_date"
                    ],
                    "country": record[
                        "country"
                    ],
                    "location": record[
                        "location"
                    ],
                    "host": record[
                        "host"
                    ],
                    "sequencing_platform": record[
                        "sequencing_platform"
                    ],
                    "study_accession": record[
                        "study_accession"
                    ],
                    "source": "NCBI Nucleotide",
                }
            )


def main():
    parser = argparse.ArgumentParser(
        description=(
            "Retrieve and normalize NCBI metadata for "
            "accessions contained in a FASTA file."
        )
    )

    parser.add_argument(
        "--fasta",
        required=True,
        help="FASTA containing NCBI accession identifiers.",
    )

    parser.add_argument(
        "--output",
        required=True,
        help="Output normalized metadata TSV.",
    )

    args = parser.parse_args()

    fasta = Path(args.fasta)

    if not fasta.exists():
        raise RuntimeError(
            f"FASTA does not exist: {fasta}"
        )

    accessions = fasta_accessions(fasta)

    if not accessions:
        raise RuntimeError(
            "No sequence identifiers were found in the FASTA."
        )

    print(
        f"FASTA sequences: {len(accessions)}",
        file=sys.stderr,
    )

    records = {}

    try:
        print(
            "Querying NCBI ESearch...",
            file=sys.stderr,
        )

        uids = esearch(accessions)

        print(
            f"NCBI records resolved: {len(uids)}",
            file=sys.stderr,
        )

        if len(uids) != len(accessions):
            print(
                "WARNING: not every FASTA identifier resolved in NCBI.",
                file=sys.stderr,
            )

        if uids:
            print(
                "Retrieving NCBI GenBank records...",
                file=sys.stderr,
            )

            genbank = efetch(uids)

            records = parse_records(genbank)

            print(
                f"NCBI records parsed: {len(records)}",
                file=sys.stderr,
            )
        else:
            print(
                "WARNING: NCBI returned no matching records. "
                "Continuing without NCBI metadata.",
                file=sys.stderr,
            )

    except RuntimeError as exc:
        print(
            f"WARNING: NCBI metadata retrieval failed: {exc}",
            file=sys.stderr,
        )
        print(
            "WARNING: continuing without NCBI metadata.",
            file=sys.stderr,
        )

    write_metadata(
        accessions,
        records,
        args.output,
    )

    print(
        f"Metadata written: {args.output}",
        file=sys.stderr,
    )


if __name__ == "__main__":
    main()
