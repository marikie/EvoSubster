#! /bin/bash

argNum=$#

if [ "$argNum" -ne 1 ]; then
    echo "Usage: $0 <Today's Date>"
    exit 1
fi

Date=$1

script=${HOME}/sbst/scripts/last/trisbst_3spc_fromDwl.sh
outDir=${HOME}/sbst/data/cnidaria
logDir=${HOME}/sbst/log

run_triple() {
    local log_suffix=$1
    shift
    local log_file="${logDir}/${Date}_${log_suffix}.log"
    if ! bash "${script}" --out-dir "${outDir}" "${Date}" "$@" &> "${log_file}"; then
        echo "Warning: run ${log_suffix} failed; check ${log_file}" >&2
    fi
}

# run_triple "chFus_chQui_chChe" GCA_009936425.2 GCA_014526335.2 GCA_011763395.2

# run_triple "caOrn_caXam_caAnd" GCA_964304725.1 GCA_964235115.1 GCA_018155075.1

# run_triple "heMag_stHad_stMer" GCA_011763375.2 GCA_049996035.1 GCA_011800005.2

# run_triple "heVer_phCru_phLol" GCA_044474045.1 GCA_009858155.3 GCA_044473265.1

# run_triple "acDig_acHyp_acMil" GCF_000222465.1 GCA_964291705.1 GCF_013753865.1

# run_triple "acDig_acHyp_acSpi" GCF_000222465.1 GCA_964291705.1 GCA_964261235.1

# run_triple "acMil_acNas_acMic" GCF_013753865.1 GCA_014634205.1 GCA_014634165.1

# run_triple "acAus_acTen_acSpa" GCA_964273435.1 GCA_014633955.1 GCA_031770025.1

# run_triple "moCap_moGri_moEff" GCA_949126865.1 GCA_043882275.1 GCA_014634505.1

# run_triple "duAx_deCr_tuCo" GCA_964258685.1 GCA_024195265.1 GCA_047759845.1

# run_triple "poRu_poAu_poCy" GCA_964035705.1 GCA_022179025.1 GCA_964035525.1

run_triple "poRu_poLo_poLu" GCA_964035705.1 GCA_942486035.1 GCF_958299795.1

# ------
# run_triple "poGr_poVer_poDa" GCA_964027065.2 GCF_036669915.1 GCF_003704095.1

# run_triple "stPi_poVer_poDa" GCF_002571385.2 GCF_036669915.1 GCF_003704095.1

# run_triple "ecHo_orFa_cySa" GCA_964199735.2 GCF_002042975.1 GCA_964194085.1

# run_triple "cySa_orFa_orFr" GCA_964194085.1 GCF_002042975.1 GCA_964199315.1

# run_triple "paAt_paDa_paSw" GCA_025448195.2 GCA_025403605.1 GCA_025403585.1

# run_triple "umNa_umCh_umPa" GCA_027574985.2 GCA_025448255.1 GCA_027575125.1

# run_triple "umKa_umCh_umPa" GCA_025403525.1 GCA_025448255.1 GCA_027575125.1

# run_triple "isTu_zoSo_zoPu" GCA_025447995.2 GCA_025331085.1 GCA_027575285.1

# run_triple "pyMi_pyCa_pyMu" GCA_042846405.1 GCA_965234985.1 GCA_027575235.1

# run_triple "pyCa_pyMu_pyGr" GCA_965234985.1 GCA_027575235.1 GCA_026546935.1

# run_triple "epIl_epRa_epPl" GCA_026546655.1 GCA_026413505.1 GCA_025388665.1

# run_triple "phMa_phSt_phSu" GCA_025400095.1 GCA_025400075.2 GCA_025434845.1

# run_triple "paPa_coRu_heIm" GCA_035771165.1 GCA_964265095.1 GCA_045791825.1
