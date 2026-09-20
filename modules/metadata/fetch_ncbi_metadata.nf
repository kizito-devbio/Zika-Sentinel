process FETCH_NCBI_METADATA {

    tag "NCBI metadata"

    label 'low_compute'

    input:
    path fasta

    output:
    path "metadata_ncbi.tsv", emit: metadata

    script:
    """
    python /opt/zika-sentinel/bin/fetch_ncbi_metadata.py \
        --fasta ${fasta} \
        --output metadata_ncbi.tsv
    """
}
