#! /bin/bash

argNum=$#

if [ $argNum -ne 1 ]; then
    echo "Usage: $0 <Today's Date>"
    exit 1
fi

Date=$1

DRIVER_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
script="$DRIVER_DIR/../src/trisbst_3spc_fromDwl.sh"
ROOT_DIR="$(cd "$DRIVER_DIR/.." && pwd)"
outDir=$ROOT_DIR/results/actinopteri
logDir=$ROOT_DIR/log

bash ${script} --out-dir ${outDir} ${Date}  GCF_007364275.1 GCA_015108585.1 GCA_013435755.1 &> ${logDir}/${Date}_Arccen1_Ampzal2_Ampcit3.log

