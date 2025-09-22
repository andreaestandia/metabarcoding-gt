#!/usr/bin/env bash

# Define base project path
PROJECT="/media/estandia/LaCie/projects/metabarcoding-gt/results"

# Define markers to process
MARKERS=("16s" "rbcl" "coi")

for marker in "${MARKERS[@]}"; do
    PATH_ASV="${PROJECT}/${marker}/outputs/fasta_asv/"
    PATH_OTUS="${PROJECT}/${marker}/outputs/otus/"

    echo "Processing marker: $marker"

    # Create output folder if missing
    mkdir -p "$PATH_OTUS"

    # Loop over all *_ASVs.fasta files for this marker
    for file in "${PATH_ASV}"*_ASVs.fasta; do
        [[ -e "$file" ]] || { echo "No ASV files found for $marker"; continue; }

        base=$(basename "$file" _ASVs.fasta)
        out="${PATH_OTUS}${base}_OTU.fasta"
        uc="${PATH_OTUS}${base}_OTU.uc"

        echo "  Clustering: $file -> $out + $uc"

        # Run vsearch and generate .uc file
        vsearch --cluster_size "$file" \
                --id 0.97 \
                --centroids "$out" \
                --uc "$uc" \
                --sizein \
                --relabel otu

        # Post-process UC file to replace numeric centroid labels with otu1, otu2...
        if [[ -s "$uc" ]]; then
            tmp_asvgroups="${uc%.uc}_asvgroups.temp"
            tmp_otus="${uc%.uc}_otus.temp"
            tmp_rest="${uc%.uc}_rest.temp"
            tmp_final="${uc%.uc}_final.temp"

            # Remove cluster summary lines
            grep -v "^C" "$uc" > "$tmp_asvgroups"

            # Replace centroid numbers with otu1, otu2, ...
            cut -f2 "$tmp_asvgroups" | while read -r l; do echo "otu$(($l + 1))"; done > "$tmp_otus"

            # Extract rest of table (columns 3-9)
            cut -f3-9 "$tmp_asvgroups" > "$tmp_rest"

            # Combine
            paste "$tmp_otus" "$tmp_rest" > "$tmp_final"

            # Replace original .uc file with processed version
            mv "$tmp_final" "$uc"

            # Cleanup
            rm "$tmp_asvgroups" "$tmp_otus" "$tmp_rest"
        else
            echo "  Warning: $uc was empty or not created!"
        fi
    done
done
