#! /bin/bash

argNum=$#

if [ $argNum -ne 1 ]; then
    echo "Usage: $0 <Today's Date>"
    exit 1
fi

Date=$1

script=${HOME}/sbst/scripts/last/trisbst_3spc_fromDwl.sh
outDir=${HOME}/sbst/data/arthropoda
logDir=${HOME}/sbst/log

bash ${script} --out-dir ${outDir} ${Date} GCF_000001215.4 GCF_016746395.2 GCF_004382195.2 &> ${logDir}/${Date}_droMel_droSim_droSec.log

bash ${script} --out-dir ${outDir} ${Date} GCF_000469605.1 GCF_003254395.2 GCF_029169275.1 &> ${logDir}/${Date}_apDor_apMel_apCer.log

# bash ${script} --out-dir ${outDir} ${Date} GCF_016920715.1 GCF_943734735.2 GCF_017562075.2 &> ${logDir}/${Date}_anoAra_anoGam_anoMer.log

bash ${script} --out-dir ${outDir} ${Date} GCA_002197625.1 GCF_030269925.1 GCF_003987935.1 &> ${logDir}/${Date}_bomHut_bomMor_bomMan.log

# bash ${script} --out-dir ${outDir} ${Date} GCA_965638035.1 GCA_046254815.1 GCA_038098605.1 &> ${logDir}/${Date}_gryBim_gryAss_gryLon.log

# bash ${script} --out-dir ${outDir} ${Date} GCF_015345945.1 GCF_031307605.1 GCA_939628115.1 &> ${logDir}/${Date}_triMad_triCas_triFre.log

# bash ${script} --out-dir ${outDir} ${Date} GCA_032361865.1 GCA_032273845.1 GCA_032361485.1 &> ${logDir}/${Date}_panPas_panLon_panCyg.log

# bash ${script} --out-dir ${outDir} ${Date} GCA_032361705.1 GCA_043589495.1 GCF_036320965.1 &> ${logDir}/${Date}_panVer_panHom_panOrn.log

# bash ${script} --out-dir ${outDir} ${Date} GCF_017591435.1 GCA_965230865.2 GCA_035046885.1 &> ${logDir}/${Date}_porTri_porSeg_porPel.log

# bash ${script} --out-dir ${outDir} ${Date} GCA_034695945.1 GCA_023279145.1 GCA_034696005.1 &> ${logDir}/${Date}_gloHex_gloMae_gloMar.log

# bash ${script} --out-dir ${outDir} ${Date} GCA_000239455.1 GCA_034683725.1 GCA_949358305.1 &> ${logDir}/${Date}_strMar_strCra_strAcu.log

