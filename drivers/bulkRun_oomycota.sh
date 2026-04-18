#! /bin/bash

argNum=$#

if [ $argNum -ne 1 ]; then
    echo "Usage: $0 <Today's Date>"
    exit 1
fi

Date=$1

DRIVER_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
script="$DRIVER_DIR/../src/sbst_fromDwl.sh"
ROOT_DIR="$(cd "$DRIVER_DIR/.." && pwd)"
outDir=$ROOT_DIR/results/oomycota
logDir=$ROOT_DIR/log

bash ${script} --out-dir ${outDir} ${Date} GCA_023338075.1 GCA_023338115.1 GCA_023334395.1 &> ${logDir}/${Date}_phyLit_phyBor_phySin.log

bash ${script} --out-dir ${outDir} ${Date} GCA_023338175.1 GCA_023338125.1 GCA_023338075.1 &> ${logDir}/${Date}_phyDel_phyAic_phyLit.log

bash ${script} --out-dir ${outDir} ${Date} GCA_023338055.1 GCA_023338135.1 GCA_023338025.1 &> ${logDir}/${Date}_phyMir_phyKen_phyOed.log

bash ${script} --out-dir ${outDir} ${Date} GCA_023335835.1 GCA_023336395.1 GCA_023334555.1 &> ${logDir}/${Date}_eloSen_eloDim_eloUnd.log

bash ${script} --out-dir ${outDir} ${Date} GCA_023333865.1 GCA_023334105.1 GCA_023334365.1 &> ${logDir}/${Date}_gloPer_gloNod_gloAca.log

bash ${script} --out-dir ${outDir} ${Date} GCA_023334105.1 GCA_023334075.1 GCA_023334005.1 &> ${logDir}/${Date}_gloNod_gloNun_gloOth.log

bash ${script} --out-dir ${outDir} ${Date} GCA_023334305.1 GCA_023334245.1 GCA_023336215.1 &> ${logDir}/${Date}_gloCan_gloCed_gloIwa.log

bash ${script} --out-dir ${outDir} ${Date} GCA_031305395.1 GCA_018806915.1 GCA_032432875.1 &> ${logDir}/${Date}_phyCit_phyCol_phyMea.log

bash ${script} --out-dir ${outDir} ${Date} GCA_024679075.1 GCA_042082425.1 GCA_030463285.1 &> ${logDir}/${Date}_phyUli_phyFle_phyEur.log

bash ${script} --out-dir ${outDir} ${Date} GCA_024679115.1 GCA_024679175.1 GCA_016169955.1 &> ${logDir}/${Date}_phyPis_phyCaj_phyVig.log

bash ${script} --out-dir ${outDir} ${Date} GCA_024679275.1 GCF_000149755.1 GCA_024679225.1 &> ${logDir}/${Date}_phyNie_phySoj_phyPis.log

bash ${script} --out-dir ${outDir} ${Date} GCA_016880985.1 GCA_016864655.1 GCA_042082605.1 &> ${logDir}/${Date}_phyIda_phyCac_phyHed.log

