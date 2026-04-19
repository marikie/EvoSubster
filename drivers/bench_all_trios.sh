#! /bin/bash
# Benchmark driver: run short-paper trios with /usr/bin/time -v
# Produces per-trio .time files in results/benchmark/ for later parsing.
#
# History:
#   * 2026-04-12 (bench_20260412): ran 19 fungi + 16 cnidaria trios (commented
#     out below for reference; their .time files are already in the benchmark
#     directory).
#   * 2026-04-16 (bench_20260416): added 5 Glomeromycetes fungi trios (active
#     below) that were missing from the original sweep.

argNum=$#

if [ "$argNum" -ne 1 ]; then
    echo "Usage: $0 <Date>"
    exit 1
fi

Date=$1

DRIVER_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
script="$DRIVER_DIR/../src/sbst_fromDwl.sh"
ROOT_DIR="$(cd "$DRIVER_DIR/.." && pwd)"
logDir=$ROOT_DIR/log
benchGenomeDir=$ROOT_DIR/bench_genomes
benchOutDir=$ROOT_DIR/bench_results
timeDir=$ROOT_DIR/bench_results/benchmark

rm -rf "$benchGenomeDir" "$benchOutDir"
mkdir -p "$benchGenomeDir" "$benchOutDir" "$timeDir"

TRIO_NUM=0
TOTAL=40

run_bench() {
    local label=$1
    local lineage=$2
    shift 2
    TRIO_NUM=$((TRIO_NUM + 1))
    echo "[${TRIO_NUM}/${TOTAL}] Starting: ${label} (${lineage})"
    /usr/bin/time -v -o "${timeDir}/${Date}_${label}.time" \
        bash "${script}" --genome-dir "${benchGenomeDir}" --out-dir "${benchOutDir}" "${Date}" "$@" \
        &> "${logDir}/${Date}_bench_${label}.log"
    local exit_code=$?
    echo "[${TRIO_NUM}/${TOTAL}] Done:     ${label} (exit ${exit_code})"
}

# ---- Fungi: Glomeromycetes (5 trios) ----

run_bench "denHet_gigRos_gigMar" fungi \
    GCA_910591775.1 GCA_003550325.1 GCA_009809945.1

run_bench "funGeo_funCal_funMos" fungi \
    GCA_946474995.1 GCA_910591825.1 GCA_910592005.1

run_bench "parBra_parOcc_parOcc" fungi \
    GCA_910592345.1 GCA_910592205.1 GCA_022605545.1

run_bench "rhiCla_rhiIrr_rhiPro" fungi \
    GCA_015698045.1 GCF_026210795.1 GCA_019425655.1

run_bench "rhiPro_rhiIrr_rhiIrr" fungi \
    GCA_019425655.1 GCF_026210795.1 GCA_020716745.1

# ---- Fungi (19 trios) ----

run_bench "agaBis_agaBit_agaSin" fungi \
    GCF_000300575.1 GCA_030246685.1 GCA_022315185.1

run_bench "armBor_armGal_armAlt" fungi \
    GCA_030435635.1 GCA_037576215.1 GCA_022818075.1

run_bench "bolBar_bolRet_bolNob" fungi \
    GCA_038088775.1 GCA_038093535.1 GCA_038088815.1

run_bench "bolNan_bolTom_bolMar" fungi \
    GCA_038093035.1 GCA_038092475.1 GCA_038093055.1

run_bench "bolRex_bolRet_bolEdu" fungi \
    GCA_038088795.1 GCA_018397855.1 GCA_015179015.1

run_bench "bolSem_bolTyl_bolPse" fungi \
    GCA_038090295.1 GCA_038093875.1 GCA_038092835.1

run_bench "bolVar_bolEdu_bolRex" fungi \
    GCA_038092395.1 GCA_015179015.1 GCA_038088795.1

run_bench "inoSue_inoTig_inoFlo" fungi \
    GCA_043168805.1 GCA_964248975.1 GCA_043167125.1

run_bench "lacAka_lacHen_lacPse" fungi \
    GCA_021524915.1 GCA_021525025.1 GCA_021525015.1

run_bench "lacAme_lacBic_lacTri" fungi \
    GCA_000827195.1 GCF_000143565.1 GCA_018417955.1

run_bench "lacSan_lacDel_lacHat" fungi \
    GCA_021527775.1 GCA_021525775.1 GCA_024734325.1

run_bench "lecGlu_lecDis_lecPro" fungi \
    GCA_038092055.1 GCA_038092035.1 GCA_038091695.1

run_bench "lecObs_lecIns_lecImi" fungi \
    GCA_038091815.1 GCA_038091875.1 GCA_038091955.1

run_bench "lenEdo_lenLat_lenNov" fungi \
    GCF_021015755.1 GCA_028011325.1 GCA_027921425.1

run_bench "pleOst_pleEry_pleTuo" fungi \
    GCF_014466165.1 GCA_029467805.1 GCA_036872985.1

run_bench "podMar_podPis_podRug" fungi \
    GCA_018524725.1 GCA_018524465.1 GCA_018524415.1

run_bench "podMin_podVer_podHum" fungi \
    GCA_016098005.1 GCA_000739165.1 GCA_025677895.1

run_bench "rusAbi_rusGri_rusLep" fungi \
    GCA_003313715.1 GCA_022884055.1 GCA_003316425.1

run_bench "strPac_strLuc_strSte" fungi \
    GCA_019915135.1 GCA_019915105.1 GCA_019915075.1

# ---- Cnidaria (16 trios) ----

run_bench "acrAus_acrTen_acrSpa" cnidaria \
    GCA_964273435.1 GCA_014633955.1 GCA_031770025.1

run_bench "acrDig_acrHya_acrMil" cnidaria \
    GCF_000222465.1 GCA_964291705.1 GCF_013753865.1

run_bench "acrDig_acrHya_acrSpi" cnidaria \
    GCF_000222465.1 GCA_964291705.1 GCA_964261235.1

run_bench "acrMil_acrNas_acrMic" cnidaria \
    GCF_013753865.1 GCA_014634205.1 GCA_014634165.1

run_bench "casOrn_casXam_casAnd" cnidaria \
    GCA_964304725.1 GCA_964235115.1 GCA_018155075.1

run_bench "cypSal_orbFav_orbFra" cnidaria \
    GCA_964194085.1 GCF_002042975.1 GCA_964199315.1

run_bench "dunAxi_denCri_tubCoc" cnidaria \
    GCA_964258685.1 GCA_024195265.1 GCA_047759845.1

run_bench "echHor_orbFav_cypSal" cnidaria \
    GCA_964199735.2 GCF_002042975.1 GCA_964194085.1

run_bench "hetMag_stiHad_stiMer" cnidaria \
    GCA_011763375.2 GCA_049996035.1 GCA_011800005.2

run_bench "monCap_monGri_monEff" cnidaria \
    GCA_949126865.1 GCA_043882275.1 GCA_014634505.1

run_bench "palCar_palMut_palGra" cnidaria \
    GCA_965234985.1 GCA_027575235.1 GCA_026546935.1

run_bench "palMiz_palCar_palMut" cnidaria \
    GCA_042846405.1 GCA_965234985.1 GCA_027575235.1

run_bench "pocGra_pocVer_pocDam" cnidaria \
    GCA_964027065.2 GCF_036669915.1 GCF_003704095.1

run_bench "porRus_porAus_porCyl" cnidaria \
    GCA_964035705.1 GCA_022179025.1 GCA_964035525.1

run_bench "porRus_porLob_porLut" cnidaria \
    GCA_964035705.1 GCA_942486035.1 GCF_958299795.1

run_bench "styPis_pocVer_pocDam" cnidaria \
    GCF_002571385.2 GCF_036669915.1 GCF_003704095.1

echo "All ${TOTAL} trios complete. Time files in: ${timeDir}"

# ---- Post-processing: generate bench_table.tsv and bench_figure.pdf ----

echo "Running post-processing..."

python "${ROOT_DIR}/src/report/parse_bench_times.py" "${timeDir}"

python "${ROOT_DIR}/src/report/collect_run_summary.py" "${benchOutDir}" \
    --output "${timeDir}/benchmark_summary.json"

python "${ROOT_DIR}/src/report/collect_bench_metadata.py" "${benchOutDir}" "${benchGenomeDir}"

Rscript "${ROOT_DIR}/src/visualize/bench_figure.R" \
    "${timeDir}/bench_summary.tsv" \
    "${timeDir}/bench_metadata.tsv" \
    "${timeDir}"

echo "Done. bench_table.tsv and bench_figure.pdf in: ${timeDir}"
