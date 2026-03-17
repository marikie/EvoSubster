#!/bin/bash

# Function to generate a graph; outputs PDF to out_dir and log to log_dir
generate_graph() {
    local input_file=$1
    local script=$2
    local out_dir=$3
    local log_dir=$4
    local base_name
    base_name=$(basename "${input_file%.*}")
    local output_log="${log_dir}/${base_name}.out"

    echo "time Rscript $script $input_file $out_dir > $output_log"
    time Rscript "$script" "$input_file" "$out_dir" > "$output_log"
}

# Get arguments
org2tsv=$1
org3tsv=$2
org2tsv_maflinked=$3
org3tsv_maflinked=$4
org2tsv_errprb=$5
org3tsv_errprb=$6
org2tsv_maflinked_errprb=$7
org3tsv_maflinked_errprb=$8

org2tsv_ncds=$9
org3tsv_ncds=${10}
org2tsv_maflinked_ncds=${11}
org3tsv_maflinked_ncds=${12}

org2_dinuc_tsv=${13}
org3_dinuc_tsv=${14}
org2_dinuc_tsv_maflinked=${15}
org3_dinuc_tsv_maflinked=${16}

org2_dinuc_tsv_ncds=${17}
org3_dinuc_tsv_ncds=${18}
org2_dinuc_tsv_maflinked_ncds=${19}
org3_dinuc_tsv_maflinked_ncds=${20}

r_scripts_path=${21}
figs_dir=${22:-figs}
misc_dir=${23:-statistics/misc}

figs_sn_wml="${figs_dir}/singlenuc/ratio/without_maflinked"
figs_sn_ml="${figs_dir}/singlenuc/ratio"
figs_lr_wml="${figs_dir}/singlenuc/log-ratio/without_maflinked"
figs_lr_ml="${figs_dir}/singlenuc/log-ratio"
figs_dn_wml="${figs_dir}/dinuc/without_maflinked"
figs_dn_ml="${figs_dir}/dinuc"


# Trinucleotide substitutions graphs
echo "---making graphs of the trinucleotide substitutions"
generate_graph "$org2tsv" "$r_scripts_path/sbmut.R" "$figs_sn_wml" "$misc_dir"
generate_graph "$org3tsv" "$r_scripts_path/sbmut.R" "$figs_sn_wml" "$misc_dir"
echo "time Rscript $r_scripts_path/logRatioPlot.R $org2tsv $figs_lr_wml"
time Rscript "$r_scripts_path/logRatioPlot.R" "$org2tsv" "$figs_lr_wml" \
    > "${misc_dir}/$(basename "${org2tsv%.*}").out"
echo "time Rscript $r_scripts_path/logRatioPlot.R $org3tsv $figs_lr_wml"
time Rscript "$r_scripts_path/logRatioPlot.R" "$org3tsv" "$figs_lr_wml" \
    > "${misc_dir}/$(basename "${org3tsv%.*}").out"

# Trinucleotide substitutions graphs (of ncds)
if [ -e "$org2tsv_ncds" ] && [ -e "$org3tsv_ncds" ]; then
    echo "---making graphs of the trinucleotide substitutions (of ncds)"
    generate_graph "$org2tsv_ncds" "$r_scripts_path/sbmut.R" "$figs_sn_wml" "$misc_dir"
    generate_graph "$org3tsv_ncds" "$r_scripts_path/sbmut.R" "$figs_sn_wml" "$misc_dir"
    echo "time Rscript $r_scripts_path/logRatioPlot.R $org2tsv_ncds $figs_lr_wml"
    time Rscript "$r_scripts_path/logRatioPlot.R" "$org2tsv_ncds" "$figs_lr_wml" \
        > "${misc_dir}/$(basename "${org2tsv_ncds%.*}").out"
    echo "time Rscript $r_scripts_path/logRatioPlot.R $org3tsv_ncds $figs_lr_wml"
    time Rscript "$r_scripts_path/logRatioPlot.R" "$org3tsv_ncds" "$figs_lr_wml" \
        > "${misc_dir}/$(basename "${org3tsv_ncds%.*}").out"
fi

# Trinucleotide substitutions graphs (with maf-linked)
echo "---making graphs of the trinucleotide substitutions (with maf-linked)"
generate_graph "$org2tsv_maflinked" "$r_scripts_path/sbmut.R" "$figs_sn_ml" "$misc_dir"
generate_graph "$org3tsv_maflinked" "$r_scripts_path/sbmut.R" "$figs_sn_ml" "$misc_dir"
echo "time Rscript $r_scripts_path/logRatioPlot.R $org2tsv_maflinked $figs_lr_ml"
time Rscript "$r_scripts_path/logRatioPlot.R" "$org2tsv_maflinked" "$figs_lr_ml" \
    > "${misc_dir}/$(basename "${org2tsv_maflinked%.*}").out"
echo "time Rscript $r_scripts_path/logRatioPlot.R $org3tsv_maflinked $figs_lr_ml"
time Rscript "$r_scripts_path/logRatioPlot.R" "$org3tsv_maflinked" "$figs_lr_ml" \
    > "${misc_dir}/$(basename "${org3tsv_maflinked%.*}").out"

# Trinucleotide substitutions graphs (of ncds and maf-linked)
if [ -e "$org2tsv_maflinked_ncds" ] && [ -e "$org3tsv_maflinked_ncds" ]; then
    echo "---making graphs of the trinucleotide substitutions (of maf-linked_ncds)"
    generate_graph "$org2tsv_maflinked_ncds" "$r_scripts_path/sbmut.R" "$figs_sn_ml" "$misc_dir"
    generate_graph "$org3tsv_maflinked_ncds" "$r_scripts_path/sbmut.R" "$figs_sn_ml" "$misc_dir"
    echo "time Rscript $r_scripts_path/logRatioPlot.R $org2tsv_maflinked_ncds $figs_lr_ml"
    time Rscript "$r_scripts_path/logRatioPlot.R" "$org2tsv_maflinked_ncds" "$figs_lr_ml" \
        > "${misc_dir}/$(basename "${org2tsv_maflinked_ncds%.*}").out"
    echo "time Rscript $r_scripts_path/logRatioPlot.R $org3tsv_maflinked_ncds $figs_lr_ml"
    time Rscript "$r_scripts_path/logRatioPlot.R" "$org3tsv_maflinked_ncds" "$figs_lr_ml" \
        > "${misc_dir}/$(basename "${org3tsv_maflinked_ncds%.*}").out"
fi

# Dinucleotide substitutions graphs
echo "---making graphs of the dinucleotide substitutions"
time Rscript "$r_scripts_path/dinucleotide-plot.R" "$org2_dinuc_tsv" "$figs_dn_wml" \
    > "${misc_dir}/$(basename "${org2_dinuc_tsv}").out"
time Rscript "$r_scripts_path/dinucleotide-plot.R" "$org3_dinuc_tsv" "$figs_dn_wml" \
    > "${misc_dir}/$(basename "${org3_dinuc_tsv}").out"

# Dinucleotide substitutions graphs (of ncds)
if [ -e "$org2_dinuc_tsv_ncds" ] && [ -e "$org3_dinuc_tsv_ncds" ]; then
    echo "---making graphs of the dinucleotide substitutions (of ncds)"
    time Rscript "$r_scripts_path/dinucleotide-plot.R" "$org2_dinuc_tsv_ncds" "$figs_dn_wml" \
        > "${misc_dir}/$(basename "${org2_dinuc_tsv_ncds}").out"
    time Rscript "$r_scripts_path/dinucleotide-plot.R" "$org3_dinuc_tsv_ncds" "$figs_dn_wml" \
        > "${misc_dir}/$(basename "${org3_dinuc_tsv_ncds}").out"
fi

# Dinucleotide substitutions graphs (with maf-linked)
echo "---making graphs of the dinucleotide substitutions (with maf-linked)"
time Rscript "$r_scripts_path/dinucleotide-plot.R" "$org2_dinuc_tsv_maflinked" "$figs_dn_ml" \
    > "${misc_dir}/$(basename "${org2_dinuc_tsv_maflinked}").out"
time Rscript "$r_scripts_path/dinucleotide-plot.R" "$org3_dinuc_tsv_maflinked" "$figs_dn_ml" \
    > "${misc_dir}/$(basename "${org3_dinuc_tsv_maflinked}").out"

# Dinucleotide substitutions graphs (of ncds and maf-linked)
if [ -e "$org2_dinuc_tsv_maflinked_ncds" ] && [ -e "$org3_dinuc_tsv_maflinked_ncds" ]; then
    echo "---making graphs of the dinucleotide substitutions (of ncds and maf-linked)"
    time Rscript "$r_scripts_path/dinucleotide-plot.R" "$org2_dinuc_tsv_maflinked_ncds" "$figs_dn_ml" \
        > "${misc_dir}/$(basename "${org2_dinuc_tsv_maflinked_ncds}").out"
    time Rscript "$r_scripts_path/dinucleotide-plot.R" "$org3_dinuc_tsv_maflinked_ncds" "$figs_dn_ml" \
        > "${misc_dir}/$(basename "${org3_dinuc_tsv_maflinked_ncds}").out"
fi

# Trinucleotide substitutions graphs (with error probability)
# echo "---making graphs of the trinucleotide substitutions (with error probability)"
# generate_graph "$org2tsv_errprb" "$r_scripts_path/sbmut.R" "$figs_sn_wml" "$misc_dir"
# generate_graph "$org3tsv_errprb" "$r_scripts_path/sbmut.R" "$figs_sn_wml" "$misc_dir"

# Trinucleotide substitutions graphs (with error probability and maf-linked)
# echo "---making graphs of the trinucleotide substitutions (with error probability and maf-linked)"
# generate_graph "$org2tsv_maflinked_errprb" "$r_scripts_path/sbmut.R" "$figs_sn_ml" "$misc_dir"
# generate_graph "$org3tsv_maflinked_errprb" "$r_scripts_path/sbmut.R" "$figs_sn_ml" "$misc_dir"
