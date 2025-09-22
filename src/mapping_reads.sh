#!/usr/bin/env bash

# Define base project path
PROJECT="/media/estandia/LaCie/projects/metabarcoding-gt/results"

# Define markers to process
MARKERS=("16s" "rbcl" "coi")

for marker in "${MARKERS[@]}"; do
    # Define paths
    PATH_ASV="${PROJECT}/${marker}/outputs/fasta_asv/"
    PATH_OTUS="${PROJECT}/${marker}/outputs/otus/"
    DB="${PROJECT}/${marker}/outputs/06_${marker}_reads_style.fasta"

    echo "=== Processing marker: $marker ==="

    # ---- Process ASV FASTA files ----
    for fasta in "${PATH_ASV}"/*.fasta; do
        [ -e "$fasta" ] || continue  # skip if no fasta files
        sample=$(basename "$fasta" .fasta)
        OUTPUT="${PROJECT}/${marker}/outputs/fasta_asv/mapping_reads/${sample}_${marker}_mapping_reads.tsv"

        echo "Running vsearch for ASV: $sample"
        vsearch --search_exact "$fasta" \
                --db "$DB" \
                --otutabout "$OUTPUT"
    done

    # ---- Process OTU FASTA files ----
    for fasta in "${PATH_OTUS}"/*.fasta; do
        [ -e "$fasta" ] || continue  # skip if no fasta files
        sample=$(basename "$fasta" .fasta)
        OUTPUT="${PROJECT}/${marker}/outputs/otus/mapping_reads/${sample}_${marker}_mapping_reads.tsv"

        echo "Running vsearch for OTU: $sample"
        vsearch --search_exact "$fasta" \
                --db "$DB" \
                --otutabout "$OUTPUT"
    done
done
