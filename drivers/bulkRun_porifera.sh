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
outDir=$ROOT_DIR/results/porifera
logDir=$ROOT_DIR/log

bash ${script} --out-dir ${outDir} ${Date} GCA_965643675.1 GCA_949841015.1 GCA_964659575.1 &> ${logDir}/${Date}_apCau_apAer_apCav.log

