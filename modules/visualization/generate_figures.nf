/*
 * Data-driven visualization layer for Zika-Sentinel.
 *
 * Consumes standardized upstream outputs:
 *   - genome_qc.tsv
 *   - metadata (optional)
 *   - variants_standardized.tsv
 *   - variant_frequency.tsv
 *   - snp_matrix.tsv
 *   - tree Newick
 *
 * Generates only figures supported by available data.
 * No fabricated values. Sample IDs preserved exactly.
 * High-dimensional heatmaps are filtered for readability;
 * full matrices remain upstream.
 */

process GENERATE_FIGURES {

    tag "figures"
    label 'med_compute'

    publishDir "${params.outdir}/figures",
        mode: 'copy',
        pattern: "fig*.png"

    publishDir "${params.outdir}/figures/source_data",
        mode: 'copy',
        pattern: "fig*.csv"

    publishDir "${params.outdir}/figures",
        mode: 'copy',
        pattern: "figure_manifest.tsv"

    publishDir "${params.outdir}/logs",
        mode: 'copy',
        pattern: "figures_log.txt"

    input:
    path qc_table
    path tree
    path metadata_tsv
    path variants_tsv
    path variant_freq_tsv
    path snp_matrix_tsv

    output:
    path "fig*.png", optional: true, emit: figures
    path "fig*.csv", optional: true, emit: source_data
    path "figure_manifest.tsv", emit: manifest
    path "figures_log.txt", emit: log

    script:

    /*
     * Metadata is optional.
     *
     * When metadata is supplied, metadata_tsv contains the staged
     * metadata file.
     *
     * When metadata is absent, the workflow passes an empty list.
     * No fake NO_FILE file is created.
     */
    def has_meta = metadata_tsv != null && metadata_tsv.size() > 0
    def metadata_path = has_meta ? metadata_tsv[0].toString() : ''

    def has_var = variants_tsv != null
    def variants_path = has_var ? variants_tsv.toString() : ''

    def has_freq = variant_freq_tsv != null
    def frequency_path = has_freq ? variant_freq_tsv.toString() : ''

    def has_mat = snp_matrix_tsv != null
    def matrix_path = has_mat ? snp_matrix_tsv.toString() : ''

    """
    #!/usr/bin/env bash
    set -euo pipefail

    LOG="figures_log.txt"

    {
        echo "========================================"
        echo "Stage:       Visualization (Stage 4)"
        echo "Started:     \$(date +"%Y-%m-%d %H:%M:%S")"
        echo "QC:          ${qc_table}"

        if [ "${has_meta}" = "true" ]; then
            echo "Metadata:    ${metadata_path}"
        else
            echo "Metadata:    none supplied"
        fi

        if [ "${has_var}" = "true" ]; then
            echo "Variants:    ${variants_path}"
        else
            echo "Variants:    none supplied"
        fi

        if [ "${has_freq}" = "true" ]; then
            echo "Frequency:   ${frequency_path}"
        else
            echo "Frequency:   none supplied"
        fi

        if [ "${has_mat}" = "true" ]; then
            echo "SNP matrix:  ${matrix_path}"
        else
            echo "SNP matrix:  none supplied"
        fi

        echo "----------------------------------------"
    } > "\$LOG"

    python3 << 'PY'
import sys
import csv
from pathlib import Path
from collections import Counter
from datetime import datetime

import matplotlib
matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np


# ============================================================
# HELPERS
# ============================================================

manifest = []
skipped = []


def register(
    fig_id,
    title,
    question,
    filename,
    source,
    n_samples,
    meta_fields,
    filtering="none"
):
    manifest.append({
        "figure_id": fig_id,
        "title": title,
        "question": question,
        "filename": filename,
        "source_data": source,
        "n_samples": n_samples,
        "metadata_fields": meta_fields,
        "filtering": filtering,
        "status": "generated"
    })

    print(
        f"GENERATED {fig_id}: {filename}",
        file=sys.stderr
    )


def skip(fig_id, reason):
    skipped.append({
        "figure_id": fig_id,
        "reason": reason
    })

    manifest.append({
        "figure_id": fig_id,
        "title": "",
        "question": "",
        "filename": "",
        "source_data": "",
        "n_samples": 0,
        "metadata_fields": "",
        "filtering": "",
        "status": f"skipped: {reason}"
    })

    print(
        f"SKIPPED {fig_id}: {reason}",
        file=sys.stderr
    )


def load_tsv(path):
    # Load a TSV or CSV file if it exists and contains data.
    # The parser accepts either tab- or comma-delimited input.

    p = Path(path)

    if not p.exists() or p.stat().st_size == 0:
        return []

    with open(p, newline="") as fh:
        sample = fh.read(4096)
        fh.seek(0)

        try:
            dialect = csv.Sniffer().sniff(
                sample,
                delimiters="\\t,"
            )
        except csv.Error:
            dialect = csv.excel_tab

        return list(
            csv.DictReader(
                fh,
                dialect=dialect
            )
        )


# ============================================================
# LOAD INPUTS
# ============================================================

qc_rows = load_tsv("${qc_table}")

has_meta = ${has_meta ? 'True' : 'False'}
has_var = ${has_var ? 'True' : 'False'}
has_freq = ${has_freq ? 'True' : 'False'}
has_mat = ${has_mat ? 'True' : 'False'}

meta_rows = (
    load_tsv("${metadata_path}")
    if has_meta
    else []
)

var_rows = (
    load_tsv("${variants_path}")
    if has_var
    else []
)

freq_rows = (
    load_tsv("${frequency_path}")
    if has_freq
    else []
)


# Restrict metadata to QC sample IDs.
qc_ids = {
    r.get("sample_id")
    for r in qc_rows
    if r.get("sample_id")
}

meta_rows = [
    r
    for r in meta_rows
    if r.get("sample_id") in qc_ids
]

n_qc = len(qc_rows)


# ============================================================
# FIG01: SEQUENCE LENGTH
# ============================================================

if qc_rows:

    lengths = []
    labels = []

    for r in qc_rows:

        try:
            lengths.append(
                int(float(r.get("length", 0)))
            )

            labels.append(
                r.get("sample_id", "")
            )

        except (
            TypeError,
            ValueError
        ):
            pass

    if lengths:

        fig, ax = plt.subplots(
            figsize=(
                max(6, len(labels) * 0.6),
                4.5
            )
        )

        ax.bar(
            range(len(lengths)),
            lengths
        )

        ax.set_xticks(
            range(len(labels))
        )

        ax.set_xticklabels(
            labels,
            rotation=45,
            ha="right",
            fontsize=8
        )

        ax.set_ylabel(
            "Sequence length (bp)"
        )

        ax.set_xlabel(
            "Sample ID"
        )

        ax.set_title(
            "ZIKV genome length by sample"
        )

        fig.tight_layout()

        fig.savefig(
            "fig01_sequence_length.png",
            dpi=150
        )

        plt.close()

        with open(
            "fig01_sequence_length.csv",
            "w",
            newline=""
        ) as fh:

            w = csv.writer(fh)

            w.writerow([
                "sample_id",
                "length_bp"
            ])

            for sample, length in zip(
                labels,
                lengths
            ):
                w.writerow([
                    sample,
                    length
                ])

        register(
            "FIG01",
            "Sequence length by sample",
            "What is the length of each ZIKV consensus genome?",
            "fig01_sequence_length.png",
            "fig01_sequence_length.csv",
            n_qc,
            ""
        )

    else:

        skip(
            "FIG01",
            "no length values in QC table"
        )

else:

    skip(
        "FIG01",
        "no QC table"
    )


# ============================================================
# FIG02: N FRACTION
# ============================================================

if qc_rows:

    nfracs = []
    nlabels = []

    for r in qc_rows:

        try:

            nfracs.append(
                float(
                    r.get(
                        "n_fraction",
                        0
                    )
                )
            )

            nlabels.append(
                r.get(
                    "sample_id",
                    ""
                )
            )

        except (
            TypeError,
            ValueError
        ):
            pass

    if nfracs:

        fig, ax = plt.subplots(
            figsize=(
                max(6, len(nlabels) * 0.6),
                4.5
            )
        )

        ax.bar(
            range(len(nfracs)),
            nfracs
        )

        ax.set_xticks(
            range(len(nlabels))
        )

        ax.set_xticklabels(
            nlabels,
            rotation=45,
            ha="right",
            fontsize=8
        )

        ax.set_ylabel(
            "N fraction"
        )

        ax.set_xlabel(
            "Sample ID"
        )

        ax.set_title(
            "Ambiguous-base (N) fraction by sample"
        )

        ymax = max(nfracs)

        ax.set_ylim(
            0,
            max(
                0.05,
                ymax * 1.2
                if ymax > 0
                else 0.05
            )
        )

        fig.tight_layout()

        fig.savefig(
            "fig02_n_fraction.png",
            dpi=150
        )

        plt.close()

        with open(
            "fig02_n_fraction.csv",
            "w",
            newline=""
        ) as fh:

            w = csv.writer(fh)

            w.writerow([
                "sample_id",
                "n_fraction"
            ])

            for sample, value in zip(
                nlabels,
                nfracs
            ):

                w.writerow([
                    sample,
                    value
                ])

        register(
            "FIG02",
            "N fraction by sample",
            "What proportion of each genome is ambiguous (N)?",
            "fig02_n_fraction.png",
            "fig02_n_fraction.csv",
            n_qc,
            ""
        )

    else:

        skip(
            "FIG02",
            "no n_fraction values"
        )

else:

    skip(
        "FIG02",
        "no QC table"
    )


# ============================================================
# FIG03: QC STATUS
# ============================================================

if qc_rows:

    sc = Counter(
        r.get(
            "status",
            "UNKNOWN"
        )
        for r in qc_rows
    )

    order = [
        "PASS",
        "WARN",
        "FAIL"
    ]

    labels_s = (
        [
            s
            for s in order
            if s in sc
        ]
        +
        [
            s
            for s in sc
            if s not in order
        ]
    )

    vals = [
        sc[s]
        for s in labels_s
    ]

    fig, ax = plt.subplots(
        figsize=(5, 4)
    )

    ax.bar(
        labels_s,
        vals
    )

    ax.set_ylabel(
        "Number of genomes"
    )

    ax.set_title(
        "Genome QC status"
    )

    for i, value in enumerate(vals):

        ax.text(
            i,
            value + 0.05,
            str(value),
            ha="center"
        )

    fig.tight_layout()

    fig.savefig(
        "fig03_qc_status.png",
        dpi=150
    )

    plt.close()

    with open(
        "fig03_qc_status.csv",
        "w",
        newline=""
    ) as fh:

        w = csv.writer(fh)

        w.writerow([
            "status",
            "count"
        ])

        for status, value in zip(
            labels_s,
            vals
        ):

            w.writerow([
                status,
                value
            ])

    register(
        "FIG03",
        "QC status counts",
        "How many genomes passed, were warned, or failed QC?",
        "fig03_qc_status.png",
        "fig03_qc_status.csv",
        n_qc,
        ""
    )

else:

    skip(
        "FIG03",
        "no QC table"
    )


# ============================================================
# FIG04: SEQUENCES BY COUNTRY
# ============================================================

countries = [
    r.get(
        "country",
        ""
    ).strip()
    for r in meta_rows
    if r.get(
        "country",
        ""
    ).strip()
]

if countries:

    cc = Counter(countries)

    fig, ax = plt.subplots(
        figsize=(
            7,
            max(
                3.5,
                len(cc) * 0.4
            )
        )
    )

    names = list(
        cc.keys()
    )

    vals = [
        cc[n]
        for n in names
    ]

    ax.barh(
        names,
        vals
    )

    ax.set_xlabel(
        "Number of sequences"
    )

    ax.set_title(
        "ZIKV sequences by country (dataset)"
    )

    fig.tight_layout()

    fig.savefig(
        "fig04_sequences_by_country.png",
        dpi=150
    )

    plt.close()

    with open(
        "fig04_sequences_by_country.csv",
        "w",
        newline=""
    ) as fh:

        w = csv.writer(fh)

        w.writerow([
            "country",
            "count"
        ])

        for name, value in cc.most_common():

            w.writerow([
                name,
                value
            ])

    register(
        "FIG04",
        "Sequences by country",
        "How many sequences in this dataset were collected in each country?",
        "fig04_sequences_by_country.png",
        "fig04_sequences_by_country.csv",
        len(countries),
        "country"
    )

else:

    skip(
        "FIG04",
        "no country metadata"
    )


# ============================================================
# FIG05: SEQUENCES BY COLLECTION YEAR
# ============================================================

years = []

for r in meta_rows:

    date_value = (
        r.get(
            "collection_date",
            ""
        )
        or ""
    ).strip()

    if not date_value:
        continue

    for fmt in (
        "%Y-%m-%d",
        "%Y-%m",
        "%Y"
    ):

        try:

            years.append(
                datetime.strptime(
                    date_value,
                    fmt
                ).year
            )

            break

        except ValueError:

            pass


if years:

    yc = Counter(years)

    fig, ax = plt.subplots(
        figsize=(6, 4)
    )

    ys = sorted(
        yc.keys()
    )

    vals = [
        yc[y]
        for y in ys
    ]

    ax.bar(
        [
            str(y)
            for y in ys
        ],
        vals
    )

    ax.set_xlabel(
        "Collection year"
    )

    ax.set_ylabel(
        "Number of sequences"
    )

    ax.set_title(
        "ZIKV sequences by collection year (dataset)"
    )

    fig.tight_layout()

    fig.savefig(
        "fig05_sequences_by_year.png",
        dpi=150
    )

    plt.close()

    with open(
        "fig05_sequences_by_year.csv",
        "w",
        newline=""
    ) as fh:

        w = csv.writer(fh)

        w.writerow([
            "year",
            "count"
        ])

        for year in ys:

            w.writerow([
                year,
                yc[year]
            ])

    register(
        "FIG05",
        "Sequences by collection year",
        "How many sequences in this dataset were collected in each year?",
        "fig05_sequences_by_year.png",
        "fig05_sequences_by_year.csv",
        len(years),
        "collection_date"
    )

else:

    skip(
        "FIG05",
        "no valid collection_date values"
    )


# ============================================================
# FIG06 AND FIG07: VARIANTS
# ============================================================

if var_rows:

    # --------------------------------------------------------
    # FIG06: VARIANTS PER SAMPLE
    # --------------------------------------------------------

    per = Counter(
        r.get(
            "sample_id"
        )
        for r in var_rows
    )

    sids = sorted(
        per.keys()
    )

    vals = [
        per[s]
        for s in sids
    ]

    fig, ax = plt.subplots(
        figsize=(
            max(
                6,
                len(sids) * 0.7
            ),
            4.5
        )
    )

    ax.bar(
        range(len(sids)),
        vals
    )

    ax.set_xticks(
        range(len(sids))
    )

    ax.set_xticklabels(
        sids,
        rotation=45,
        ha="right",
        fontsize=8
    )

    ax.set_ylabel(
        "Number of SNPs vs reference"
    )

    ax.set_xlabel(
        "Sample ID"
    )

    ax.set_title(
        "SNPs per sample (vs reference)"
    )

    fig.tight_layout()

    fig.savefig(
        "fig06_variants_per_sample.png",
        dpi=150
    )

    plt.close()

    with open(
        "fig06_variants_per_sample.csv",
        "w",
        newline=""
    ) as fh:

        w = csv.writer(fh)

        w.writerow([
            "sample_id",
            "n_snps"
        ])

        for sample, value in zip(
            sids,
            vals
        ):

            w.writerow([
                sample,
                value
            ])

    register(
        "FIG06",
        "SNPs per sample",
        "How many SNPs does each sample have relative to the reference?",
        "fig06_variants_per_sample.png",
        "fig06_variants_per_sample.csv",
        len(sids),
        "",
        "source=variants_standardized"
    )


    # --------------------------------------------------------
    # FIG07: VARIANT POSITION DISTRIBUTION
    # --------------------------------------------------------

    positions = []

    for r in var_rows:

        try:

            positions.append(
                int(
                    r["position"]
                )
            )

        except (
            KeyError,
            ValueError,
            TypeError
        ):

            pass

    if positions:

        fig, ax = plt.subplots(
            figsize=(10, 3.5)
        )

        ax.hist(
            positions,
            bins=min(
                80,
                max(
                    10,
                    len(
                        set(
                            positions
                        )
                    ) // 5
                )
            ),
            edgecolor="none"
        )

        ax.set_xlabel(
            "Reference position (1-based)"
        )

        ax.set_ylabel(
            "SNP observations"
        )

        ax.set_title(
            "Distribution of SNP positions across the ZIKV reference"
        )

        fig.tight_layout()

        fig.savefig(
            "fig07_variant_position_distribution.png",
            dpi=150
        )

        plt.close()

        with open(
            "fig07_variant_position_distribution.csv",
            "w",
            newline=""
        ) as fh:

            w = csv.writer(fh)

            w.writerow([
                "position"
            ])

            for position in sorted(
                positions
            ):

                w.writerow([
                    position
                ])

        register(
            "FIG07",
            "SNP position distribution",
            "Where along the ZIKV reference do SNPs occur in this dataset?",
            "fig07_variant_position_distribution.png",
            "fig07_variant_position_distribution.csv",
            len(
                set(
                    r.get(
                        "sample_id"
                    )
                    for r in var_rows
                )
            ),
            ""
        )

    else:

        skip(
            "FIG07",
            "no parseable positions"
        )

else:

    skip(
        "FIG06",
        "no variant table"
    )

    skip(
        "FIG07",
        "no variant table"
    )


# ============================================================
# FIG08: MOST FREQUENT VARIANT POSITIONS
# ============================================================

TOP_N = 30

if freq_rows:

    scored = []

    for r in freq_rows:

        try:

            scored.append(
                (
                    int(
                        r.get(
                            "count_in_dataset"
                        )
                        or r.get(
                            "count"
                        )
                        or 0
                    ),

                    int(
                        r["position"]
                    ),

                    r.get(
                        "reference_allele",
                        ""
                    ),

                    r.get(
                        "alternate_allele",
                        ""
                    ),

                    float(
                        r.get(
                            "frequency_in_dataset"
                        )
                        or 0
                    )
                )
            )

        except (
            KeyError,
            ValueError,
            TypeError
        ):

            pass

    scored.sort(
        reverse=True
    )

    top = scored[
        :TOP_N
    ]

    if top:

        labels = [
            f"{position}:{ref}>{alt}"
            for (
                count,
                position,
                ref,
                alt,
                frequency
            )
            in top
        ]

        counts = [
            count
            for (
                count,
                position,
                ref,
                alt,
                frequency
            )
            in top
        ]

        fig, ax = plt.subplots(
            figsize=(
                10,
                max(
                    4,
                    len(top) * 0.25
                )
            )
        )

        ax.barh(
            range(
                len(top) - 1,
                -1,
                -1
            ),
            counts[::-1]
        )

        ax.set_yticks(
            range(len(top))
        )

        ax.set_yticklabels(
            labels[::-1],
            fontsize=7
        )

        ax.set_xlabel(
            "Count in dataset (not population prevalence)"
        )

        ax.set_title(
            f"Top {len(top)} most frequent SNP alleles in dataset"
        )

        fig.tight_layout()

        fig.savefig(
            "fig08_top_variant_frequency.png",
            dpi=150
        )

        plt.close()

        with open(
            "fig08_top_variant_frequency.csv",
            "w",
            newline=""
        ) as fh:

            w = csv.writer(fh)

            w.writerow([
                "position",
                "reference_allele",
                "alternate_allele",
                "count_in_dataset",
                "frequency_in_dataset"
            ])

            for (
                count,
                position,
                ref,
                alt,
                frequency
            ) in top:

                w.writerow([
                    position,
                    ref,
                    alt,
                    count,
                    frequency
                ])

        register(
            "FIG08",
            "Top frequent SNP alleles",
            "Which SNP alleles are most frequent within this dataset?",
            "fig08_top_variant_frequency.png",
            "fig08_top_variant_frequency.csv",
            n_qc,
            "",
            f"top_{TOP_N}_by_count_in_dataset"
        )

    else:

        skip(
            "FIG08",
            "frequency table empty after parse"
        )

else:

    skip(
        "FIG08",
        "no variant frequency table"
    )


# ============================================================
# FIG09: FILTERED SNP PRESENCE HEATMAP
# ============================================================

MAX_HEAT_POS = 40
MAX_HEAT_SAMPLES = 50

mat_rows = []

if (
    has_mat
    and Path(
        "${matrix_path}"
    ).exists()
):

    with open(
        "${matrix_path}",
        newline=""
    ) as fh:

        reader = csv.reader(
            fh,
            delimiter="\\t"
        )

        header = next(
            reader,
            []
        )

        pos_cols = header[1:]

        for row in reader:

            if row:

                mat_rows.append(
                    (
                        row[0],
                        row[1:]
                    )
                )

    if mat_rows and pos_cols:

        nonref = []

        for j, position in enumerate(
            pos_cols
        ):

            count = sum(
                1
                for _, values
                in mat_rows
                if (
                    j < len(values)
                    and values[j]
                    not in (
                        "0",
                        ".",
                        "",
                        "REF"
                    )
                )
            )

            nonref.append(
                (
                    count,
                    j,
                    position
                )
            )

        nonref.sort(
            reverse=True
        )

        selected = nonref[
            :MAX_HEAT_POS
        ]

        sel_idx = [
            j
            for _, j, _
            in selected
        ]

        sel_pos = [
            p
            for _, _, p
            in selected
        ]

        samples = [
            sid
            for sid, _
            in mat_rows[
                :MAX_HEAT_SAMPLES
            ]
        ]

        data = []

        for sid, values in mat_rows[
            :MAX_HEAT_SAMPLES
        ]:

            data.append([
                0
                if values[j]
                in (
                    "0",
                    ".",
                    ""
                )
                else 1
                for j in sel_idx
            ])

        data = np.array(
            data
        )

        fig, ax = plt.subplots(
            figsize=(
                max(
                    8,
                    len(sel_pos) * 0.22
                ),
                max(
                    4,
                    len(samples) * 0.35
                )
            )
        )

        im = ax.imshow(
            data,
            aspect="auto",
            cmap="Blues",
            interpolation="nearest"
        )

        ax.set_yticks(
            range(len(samples))
        )

        ax.set_yticklabels(
            samples,
            fontsize=7
        )

        ax.set_xticks(
            range(len(sel_pos))
        )

        ax.set_xticklabels(
            sel_pos,
            rotation=90,
            fontsize=6
        )

        ax.set_xlabel(
            "Reference position (filtered)"
        )

        ax.set_ylabel(
            "Sample ID"
        )

        ax.set_title(
            f"SNP presence heatmap "
            f"(top {len(sel_pos)} variable positions)"
        )

        fig.colorbar(
            im,
            ax=ax,
            fraction=0.02,
            label="0=ref/absent, 1=alt"
        )

        fig.tight_layout()

        fig.savefig(
            "fig09_snp_presence_heatmap.png",
            dpi=150
        )

        plt.close()

        with open(
            "fig09_snp_presence_heatmap.csv",
            "w",
            newline=""
        ) as fh:

            w = csv.writer(fh)

            w.writerow(
                ["sample_id"] + sel_pos
            )

            for sid, row in zip(
                samples,
                data
            ):

                w.writerow(
                    [sid] + list(row)
                )

        register(
            "FIG09",
            "SNP presence heatmap (filtered)",
            "Which of the most variable positions differ from the reference in each sample?",
            "fig09_snp_presence_heatmap.png",
            "fig09_snp_presence_heatmap.csv",
            len(samples),
            "",
            (
                f"top_{MAX_HEAT_POS}_positions_by_nonref_count; "
                f"max_{MAX_HEAT_SAMPLES}_samples; "
                "full matrix retained upstream"
            )
        )

    else:

        skip(
            "FIG09",
            "empty SNP matrix"
        )

else:

    skip(
        "FIG09",
        "no SNP matrix"
    )


# ============================================================
# FIG10: PAIRWISE SNP DISTANCE DISTRIBUTION
# ============================================================

if (
    has_mat
    and Path(
        "${matrix_path}"
    ).exists()
    and mat_rows
    and len(mat_rows) >= 2
):

    sids = [
        sid
        for sid, _
        in mat_rows
    ]

    binary = []

    for _, values in mat_rows:

        binary.append([
            0
            if value
            in (
                "0",
                ".",
                ""
            )
            else 1
            for value in values
        ])

    binary = np.array(
        binary
    )

    dists = []
    pairs = []

    for i in range(
        len(sids)
    ):

        for j in range(
            i + 1,
            len(sids)
        ):

            distance = int(
                np.sum(
                    binary[i]
                    != binary[j]
                )
            )

            dists.append(
                distance
            )

            pairs.append(
                (
                    sids[i],
                    sids[j],
                    distance
                )
            )

    if dists:

        fig, ax = plt.subplots(
            figsize=(6, 4)
        )

        ax.hist(
            dists,
            bins=min(
                20,
                max(
                    5,
                    len(
                        set(dists)
                    )
                )
            ),
            edgecolor="white"
        )

        ax.set_xlabel(
            "Pairwise SNP differences (binary presence)"
        )

        ax.set_ylabel(
            "Number of pairs"
        )

        ax.set_title(
            "Pairwise SNP-distance distribution (dataset)"
        )

        fig.tight_layout()

        fig.savefig(
            "fig10_pairwise_snp_distance.png",
            dpi=150
        )

        plt.close()

        with open(
            "fig10_pairwise_snp_distance.csv",
            "w",
            newline=""
        ) as fh:

            w = csv.writer(fh)

            w.writerow([
                "sample_a",
                "sample_b",
                "snp_differences"
            ])

            for (
                sample_a,
                sample_b,
                distance
            ) in pairs:

                w.writerow([
                    sample_a,
                    sample_b,
                    distance
                ])

        register(
            "FIG10",
            "Pairwise SNP-distance distribution",
            "How many SNP presence differences separate each pair of samples?",
            "fig10_pairwise_snp_distance.png",
            "fig10_pairwise_snp_distance.csv",
            len(sids),
            "",
            (
                "binary presence from full snp_matrix; "
                "not transmission evidence"
            )
        )

    else:

        skip(
            "FIG10",
            "could not compute distances"
        )

else:

    skip(
        "FIG10",
        "SNP matrix unavailable or <2 samples"
    )


# ============================================================
# FIG11: PHYLOGENETIC TREE
# ============================================================

tree_path = Path(
    "${tree}"
)

tree_ok = (
    tree_path.exists()
    and tree_path.stat().st_size > 20
)

tree_text = (
    tree_path.read_text()[:300]
    if tree_ok
    else ""
)

if (
    tree_ok
    and "insufficient"
    not in tree_text.lower()
    and tree_text.strip()
    not in (
        "()",
        ";"
    )
):

    # The validated Newick tree remains available from the
    # phylogeny stage. Graphical tree rendering is intentionally
    # not fabricated when a rendering dependency is unavailable.

    skip(
        "FIG11",
        "tree file present but graphical tree rendering not enabled in this environment (Newick available in phylogeny/)"
    )

else:

    skip(
        "FIG11",
        "no usable Newick tree"
    )


# ============================================================
# FIG12 AND FIG13: GENOTYPE / CLADE
# ============================================================

skip(
    "FIG12",
    "no validated ZIKV genotype/clade assignment implemented"
)

skip(
    "FIG13",
    "no validated ZIKV genotype/clade assignment implemented"
)


# ============================================================
# FIG14: RAW-READ DEPTH
# ============================================================

skip(
    "FIG14",
    "no raw-read depth tables in this consensus-pathway verification run"
)


# ============================================================
# FIG15: VARIANTS BY COUNTRY
# ============================================================

if var_rows and meta_rows:

    sid_to_country = {
        r["sample_id"]:
            r.get(
                "country",
                ""
            ).strip()
        for r in meta_rows
        if r.get(
            "country",
            ""
        ).strip()
    }

    if sid_to_country:

        cvar = Counter()

        for r in var_rows:

            country = sid_to_country.get(
                r.get(
                    "sample_id"
                ),
                ""
            )

            if country:

                cvar[country] += 1

        if cvar:

            fig, ax = plt.subplots(
                figsize=(
                    7,
                    max(
                        3.5,
                        len(cvar) * 0.4
                    )
                )
            )

            names = list(
                cvar.keys()
            )

            vals = [
                cvar[n]
                for n in names
            ]

            ax.barh(
                names,
                vals
            )

            ax.set_xlabel(
                "SNP observations in dataset"
            )

            ax.set_title(
                "SNP observations by country "
                "(dataset, not prevalence)"
            )

            fig.tight_layout()

            fig.savefig(
                "fig15_variants_by_country.png",
                dpi=150
            )

            plt.close()

            with open(
                "fig15_variants_by_country.csv",
                "w",
                newline=""
            ) as fh:

                w = csv.writer(fh)

                w.writerow([
                    "country",
                    "snp_observations"
                ])

                for name, value in cvar.most_common():

                    w.writerow([
                        name,
                        value
                    ])

            register(
                "FIG15",
                "SNP observations by country",
                "How many SNP observations in this dataset come from each country?",
                "fig15_variants_by_country.png",
                "fig15_variants_by_country.csv",
                len(sid_to_country),
                "country",
                "counts of SNP records, not population prevalence"
            )

        else:

            skip(
                "FIG15",
                "no country-linked variant records"
            )

    else:

        skip(
            "FIG15",
            "no country metadata linked to variant samples"
        )

else:

    if not var_rows:

        skip(
            "FIG15",
            "variants missing"
        )

    else:

        skip(
            "FIG15",
            "metadata missing"
        )


# ============================================================
# WRITE MANIFEST
# ============================================================

with open(
    "figure_manifest.tsv",
    "w",
    newline=""
) as fh:

    fields = [
        "figure_id",
        "title",
        "question",
        "filename",
        "source_data",
        "n_samples",
        "metadata_fields",
        "filtering",
        "status"
    ]

    w = csv.DictWriter(
        fh,
        fieldnames=fields,
        delimiter="\\t"
    )

    w.writeheader()
    w.writerows(
        manifest
    )


n_gen = sum(
    1
    for m in manifest
    if m["status"] == "generated"
)

print(
    f"Generated: {n_gen}; "
    f"Skipped entries: {len(skipped)}",
    file=sys.stderr
)


with open(
    "figures_log.txt",
    "a"
) as log:

    log.write(
        f"Generated: {n_gen}\\n"
    )

    for m in manifest:

        if m["status"] == "generated":

            log.write(
                f"  {m['figure_id']}: "
                f"{m['filename']}\\n"
            )

        else:

            log.write(
                f"  {m['figure_id']}: "
                f"{m['status']}\\n"
            )

PY

    echo "Completed: \$(date +"%Y-%m-%d %H:%M:%S")" >> "\$LOG"
    """
}
