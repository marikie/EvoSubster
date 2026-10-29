#!/bin/bash

if [ $# -lt 5 ] || [ $# -gt 6 ]; then
    echo "Usage: $0 <Today's Date> <org2FASTA> <org3FASTA> <org2ShortName> <org3ShortName> [out_dir]"
    exit 1
fi
# Get arguments
DATE=$1
org2FASTA=$2
org3FASTA=$3
org2ShortName=$4
org3ShortName=$5
out_dir="${6:-.}"

dbBasename="${org2ShortName}db_${DATE}"
dbDir="${out_dir}/${dbBasename}"
trainFile="${out_dir}/${org2ShortName}2${org3ShortName}_${DATE}.train"

# Validate settings before reusing a DB or profile.
source "$(cd "$(dirname "${BASH_SOURCE[0]}")/../lib" && pwd)/last_options.sh"
last_options_load || exit 1
last_options_guard "$out_dir" training || exit 1

# lastdb
echo "---lastdb"
if [ ! -d "$dbDir" ]; then
	echo "making lastdb"
	dbTmp="${dbDir}.tmp.$$"
	rm -rf "$dbTmp"
	mkdir -p "$dbTmp"
	if (cd "$dbTmp" && echo "time lastdb ${LASTDB_TRAIN_ARGS[*]} $dbBasename $org2FASTA" && \
		time lastdb "${LASTDB_TRAIN_ARGS[@]}" "$dbBasename" "$org2FASTA"); then
		mv "$dbTmp" "$dbDir"
	else
		rm -rf "$dbTmp"
		exit 1
	fi
else
	echo "$dbDir already exists"
fi
# -P4: makes it faster by using 4 threads (This has no effect on the results.)
# -c: Soft-mask lowercase letters.  This means that, when we compare
#     these sequences to some other sequences using lastal, lowercase
#     letters will be excluded from initial matches.  This will apply
#     to lowercase letters in both sets of sequences.
# -uRY4: selects a seeding scheme that reduces the run time and memory use, but also reduces sensitivity.

# last-train
echo "--last-train"
if [ ! -e "$trainFile" ]; then
	echo "time last-train ${LAST_TRAIN_ARGS[*]} $dbDir/$dbBasename $org3FASTA >$trainFile"
	trainTmp="${trainFile}.tmp.$$"
	if time last-train "${LAST_TRAIN_ARGS[@]}" "$dbDir/$dbBasename" "$org3FASTA" >"$trainTmp" \
			&& [ -s "$trainTmp" ]; then
		mv "$trainTmp" "$trainFile"
	else
		rm -f "$trainTmp"
		exit 1
	fi
else
	echo "$trainFile already exists"
fi
# --revsym: Force the substitution scores to have reverse-complement
#           symmetry, e.g. score(A→G) = score(T→C).  This is often
#           appropriate, if neither strand is "special".
# -C COUNT: Before extending gapped alignments, discard any gapless alignment whose query range lies in COUNT other gapless alignments with higher score-per-length. This aims to reduce run time. -C2 may reduce run time with little effect on accuracy.
