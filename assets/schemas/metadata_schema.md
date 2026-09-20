# Zika-Sentinel metadata schema (v1)

## Primary key

`sample_id` — must match the biological sample identifier used as the FASTA header
(and as the sample name derived from FASTQ filenames in the raw-read pathway).

Zika-Sentinel never invents, replaces, or renames sample identifiers.

## Columns

| Column               | Required for analysis? | Description |
|----------------------|------------------------|-------------|
| sample_id            | **Yes** (always)       | Unique biological sample identifier |
| accession            | No                     | Public database accession if different from sample_id |
| collection_date      | Required for temporal figures | ISO-like date: YYYY-MM-DD, YYYY-MM, or YYYY |
| country              | Required for geographic figures | Country of collection |
| location             | No                     | Region / city / province |
| host                 | No                     | Host species |
| sequencing_platform  | No                     | e.g. Illumina |
| study_accession      | No                     | BioProject or study ID |
| source               | No                     | Publication or data source |

Only `sample_id` is required for the pipeline to run.
Temporal figures are generated only when valid `collection_date` values exist.
Geographic figures are generated only when non-empty `country` (or location) values exist.

## Validation rules

1. Empty `sample_id` → recorded as issue.
2. Duplicate `sample_id` → recorded as issue; not silently merged.
3. Metadata row with no matching sequence → `metadata_without_sequence`.
4. Sequence with no matching metadata row → `sequence_without_metadata`.
5. `collection_date` that cannot be parsed as YYYY / YYYY-MM / YYYY-MM-DD → `malformed_date`.
6. Missing optional fields remain empty. They are never filled from accession dates, upload dates, or any other source.

## Machine-readable outputs

- `metadata/metadata_validated.tsv` — rows that matched a sequence
- `metadata/metadata_issues.tsv` — all detected issues with sample_id, issue type, detail

## Example (test data)

See `test_data/zika_public_metadata.tsv`.
