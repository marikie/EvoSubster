#!/bin/bash

set -e

# Function to generate TSV files
generate_tsv() {
    local input_file=$1
    local output2=$2
    local output3=$3
    local script=$4
    local complete_marker="${output2}.complete"

    if [ ! -s "$output2" ] || [ ! -s "$output3" ] || [ ! -e "$complete_marker" ]; then
        local tmp2="${output2}.tmp.$$"
        local tmp3="${output3}.tmp.$$"
        echo "time python3 $script $input_file -o2 ./$output2 -o3 ./$output3"
        if time python3 "$script" "$input_file" -o2 "$tmp2" -o3 "$tmp3" \
                && [ -s "$tmp2" ] && [ -s "$tmp3" ]; then
            mv "$tmp2" "$output2"
            mv "$tmp3" "$output3"
            : > "$complete_marker"
        else
            rm -f "$tmp2" "$tmp3"
            return 1
        fi
    else
        echo "$output2 and $output3 already exist"
    fi
}

# Get arguments
joinedFile=$1
org2tsv=$2
org3tsv=$3
scripts_analysis_path=$4

org2_dinuc_tsv=$5
org3_dinuc_tsv=$6

# Optional ncds arguments (only present when GFF file is available)
joinedFile_ncds=$7
org2tsv_ncds=$8
org3tsv_ncds=$9
org2_dinuc_tsv_ncds=${10}
org3_dinuc_tsv_ncds=${11}

# Single-nuc TSV files
echo "---making .tsv trinucleotide substitution files"
generate_tsv "$joinedFile" "$org2tsv" "$org3tsv" "$scripts_analysis_path/single_sbst_2TSVs.py"

# Dinuc TSV files
echo "---making .tsv dinucleotide substitution files"
generate_tsv "$joinedFile" "$org2_dinuc_tsv" "$org3_dinuc_tsv" "$scripts_analysis_path/disbst_2TSVs.py"

# Generate ncds TSV files if ncds arguments are provided
if [ -n "$joinedFile_ncds" ]; then
    # Single-nuc ncds TSV files
    echo "---making .tsv trinucleotide substitution files (ncds)"
    generate_tsv "$joinedFile_ncds" "$org2tsv_ncds" "$org3tsv_ncds" "$scripts_analysis_path/single_sbst_2TSVs.py"

    # Dinuc ncds TSV files
    echo "---making .tsv dinucleotide substitution files (ncds)"
    generate_tsv "$joinedFile_ncds" "$org2_dinuc_tsv_ncds" "$org3_dinuc_tsv_ncds" "$scripts_analysis_path/disbst_2TSVs.py"
fi
