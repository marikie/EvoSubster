#!/bin/bash

# Get arguments
org2tsv=$1
org3tsv=$2
org2tsv_ncds=$3
org3tsv_ncds=$4

org2_dinuc_tsv=$5
org3_dinuc_tsv=$6
org2_dinuc_tsv_ncds=$7
org3_dinuc_tsv_ncds=$8

r_scripts_path=$9
figs_org2=${10}
figs_org3=${11}

org2_sn="${figs_org2}/singlenuc"
org3_sn="${figs_org3}/singlenuc"
org2_lr="${figs_org2}/singlenuc/log-ratio"
org3_lr="${figs_org3}/singlenuc/log-ratio"
org2_dn="${figs_org2}/dinuc"
org3_dn="${figs_org3}/dinuc"

# Trinucleotide substitutions graphs
echo "---making graphs of the trinucleotide substitutions"
time Rscript "$r_scripts_path/sbmut.R" "$org2tsv" "$org2_sn"
time Rscript "$r_scripts_path/sbmut.R" "$org3tsv" "$org3_sn"
time Rscript "$r_scripts_path/logRatioPlot.R" "$org2tsv" "$org2_lr"
time Rscript "$r_scripts_path/logRatioPlot.R" "$org3tsv" "$org3_lr"

# Trinucleotide substitutions graphs (ncds)
if [ -e "$org2tsv_ncds" ] && [ -e "$org3tsv_ncds" ]; then
    echo "---making graphs of the trinucleotide substitutions (ncds)"
    time Rscript "$r_scripts_path/sbmut.R" "$org2tsv_ncds" "$org2_sn"
    time Rscript "$r_scripts_path/sbmut.R" "$org3tsv_ncds" "$org3_sn"
    time Rscript "$r_scripts_path/logRatioPlot.R" "$org2tsv_ncds" "$org2_lr"
    time Rscript "$r_scripts_path/logRatioPlot.R" "$org3tsv_ncds" "$org3_lr"
fi

# Dinucleotide substitutions graphs
echo "---making graphs of the dinucleotide substitutions"
time Rscript "$r_scripts_path/dinucleotide-plot.R" "$org2_dinuc_tsv" "$org2_dn"
time Rscript "$r_scripts_path/dinucleotide-plot.R" "$org3_dinuc_tsv" "$org3_dn"

# Dinucleotide substitutions graphs (ncds)
if [ -e "$org2_dinuc_tsv_ncds" ] && [ -e "$org3_dinuc_tsv_ncds" ]; then
    echo "---making graphs of the dinucleotide substitutions (ncds)"
    time Rscript "$r_scripts_path/dinucleotide-plot.R" "$org2_dinuc_tsv_ncds" "$org2_dn"
    time Rscript "$r_scripts_path/dinucleotide-plot.R" "$org3_dinuc_tsv_ncds" "$org3_dn"
fi
