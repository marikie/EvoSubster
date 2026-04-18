#!/bin/bash

argNum=3
if [ $# -ne $argNum ]; then
	echo "You need $argNum arguments" 1>&2
	echo "You'll get a joined .maf file from two .maf files" 1>&2
	echo "- path to the org1-org2 .maf file" 1>&2 # $1
	echo "- path to the org1-org3 .maf file" 1>&2 # $2
	echo "- path to the output joined .maf file" 1>&2 # $3
	exit 1
fi

org1_org2=$1
org1_org3=$2
org1_org2_base=$(basename "$org1_org2" .maf)
org1_org2_dir=$(dirname "$org1_org2")
org1_org2_sorted="${org1_org2_dir}/${org1_org2_base}_sorted.maf"
org1_org3_base=$(basename "$org1_org3" .maf)
org1_org3_dir=$(dirname "$org1_org3")
org1_org3_sorted="${org1_org3_dir}/${org1_org3_base}_sorted.maf"
joinedFile=$3

echo "---maf-joining the two .maf files"
if [ ! -e $org1_org2_sorted ]; then
	echo "---sorting $org1_org2"
	echo "time maf-sort $org1_org2 >$org1_org2_sorted"
	time maf-sort $org1_org2 >$org1_org2_sorted
else
	echo "$org1_org2_sorted already exists"
fi
# cleanup: remove unsorted MAF once the sorted file is available
if [ -s "$org1_org2_sorted" ] && [ -e "$org1_org2" ]; then
	echo "cleanup: rm -f $org1_org2"
	rm -f "$org1_org2"
fi
if [ ! -e $org1_org3_sorted ]; then
	echo "---sorting $org1_org3"
	echo "time maf-sort $org1_org3 >$org1_org3_sorted"
	time maf-sort $org1_org3 >$org1_org3_sorted
else
	echo "$org1_org3_sorted already exists"
fi
if [ -s "$org1_org3_sorted" ] && [ -e "$org1_org3" ]; then
	echo "cleanup: rm -f $org1_org3"
	rm -f "$org1_org3"
fi
if [ ! -e $joinedFile ]; then
	echo "---joining $org1_org2_sorted and $org1_org3_sorted"
	echo "time maf-join $org1_org2_sorted $org1_org3_sorted >$joinedFile"
	time maf-join $org1_org2_sorted $org1_org3_sorted >$joinedFile
else
	echo "$joinedFile already exists"
fi