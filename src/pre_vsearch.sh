#!/usr/bin/env bash

# Define base project path
PROJECT="/media/estandia/LaCie/projects/metabarcoding-gt/results"

# Define markers to process
MARKERS=("16s" "rbcl" "coi")

# Loop over markers
for marker in "${MARKERS[@]}"; do
    PATH_ASV="${PROJECT}/${marker}/outputs/fasta_asv"
    echo "Processing marker: $marker in $PATH_ASV"

    # Loop over fasta files in the folder
    for file in "$PATH_ASV"/*.fasta; do
        [ -e "$file" ] || continue  # Skip if no fasta files
        echo "Processing $file ..."

        # Extract sample name (everything before last _ASV part in the first header)
        sample=$(grep -m1 "^>" "$file" | sed 's/^>//; s/_ASV[0-9]*//')

        # Create temporary file
        tmpfile=$(mktemp)

        # Rewrite headers with ASV numbering
        awk -v sample="$sample" '
        BEGIN { counter=1 }
        /^>/ {
            printf(">ASV%d;sample=%s\n", counter++, sample)
            next
        }
        { print }
        ' "$file" > "$tmpfile"

        # Overwrite original file
        mv "$tmpfile" "$file"
    done
done
