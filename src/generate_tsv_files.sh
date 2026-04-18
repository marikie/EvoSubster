#!/bin/bash

# Function to generate TSV files
generate_tsv() {
    local input_file=$1
    local output2=$2
    local output3=$3
    local script=$4

    if [ ! -e "$output2" ] || [ ! -e "$output3" ]; then
        echo "time python $script $input_file -o2 ./$output2 -o3 ./$output3"
        time python "$script" "$input_file" -o2 "./$output2" -o3 "./$output3"
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

# Trinuc TSV files
echo "---making .tsv trinucleotide substitution files"
generate_tsv "$joinedFile" "$org2tsv" "$org3tsv" "$scripts_analysis_path/trisbst_2TSVs.py"

# Dinuc TSV files
echo "---making .tsv dinucleotide substitution files"
generate_tsv "$joinedFile" "$org2_dinuc_tsv" "$org3_dinuc_tsv" "$scripts_analysis_path/disbst_2TSVs.py"

# Generate ncds TSV files if ncds arguments are provided
if [ -n "$joinedFile_ncds" ]; then
    echo "---making .tsv trinucleotide substitution files (ncds)"
    generate_tsv "$joinedFile_ncds" "$org2tsv_ncds" "$org3tsv_ncds" "$scripts_analysis_path/trisbst_2TSVs.py"

    echo "---making .tsv dinucleotide substitution files (ncds)"
    generate_tsv "$joinedFile_ncds" "$org2_dinuc_tsv_ncds" "$org3_dinuc_tsv_ncds" "$scripts_analysis_path/disbst_2TSVs.py"
fi
