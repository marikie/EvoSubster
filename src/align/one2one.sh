#!/bin/bash

# Get arguments
DATE=$1
org1FASTA=$2
org2FASTA=$3
dbDir=$4
trainFile=$5
o2omaf=$6
dbBasename=$(basename "$dbDir")

threadNum=${THREAD_NUM:-8}
# lastdb
echo "---lastdb"
if [ ! -d "$dbDir" ]; then
	echo "making lastdb"
	mkdir -p "$dbDir"
	(cd "$dbDir" && echo "time lastdb -P${threadNum} -c -uRY4 $dbBasename $org1FASTA" && \
		time lastdb -P${threadNum} -c -uRY4 "$dbBasename" "$org1FASTA") || exit 1
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
	echo "time last-train -P${threadNum} --revsym -C2 $dbDir/$dbBasename $org2FASTA >$trainFile"
	time last-train -P${threadNum} --revsym -C2 "$dbDir/$dbBasename" "$org2FASTA" >"$trainFile"
else
	echo "$trainFile already exists"
fi
# --revsym: Force the substitution scores to have reverse-complement
#           symmetry, e.g. score(A→G) = score(T→C).  This is often
#           appropriate, if neither strand is "special".
# -C COUNT: Before extending gapped alignments, discard any gapless alignment whose query range lies in COUNT other gapless alignments with higher score-per-length. This aims to reduce run time. -C2 may reduce run time with little effect on accuracy.

# lastal | last-split -r (one-to-one alignments), then maf-linked
# maf-linked requires a file argument (doesn't accept stdin), so write to
# a temporary MAF first, then filter and delete the temporary.
echo "---lastal | last-split -r, then maf-linked"
o2omaf_base=$(basename "$o2omaf" .maf)
o2omaf_dir=$(dirname "$o2omaf")
o2omaf_sorted="${o2omaf_dir}/${o2omaf_base}_sorted.maf"
if [ ! -e "$o2omaf" ] && [ ! -e "$o2omaf_sorted" ]; then
	o2omaf_tmp="${o2omaf}.raw"
	echo "time lastal -P${threadNum} -H1 -C2 --split-f=MAF+ -p $trainFile $dbDir/$dbBasename $org2FASTA | last-split -r >$o2omaf_tmp"
	time lastal -P"${threadNum}" -H1 -C2 --split-f=MAF+ -p "$trainFile" "$dbDir/$dbBasename" "$org2FASTA" \
		| last-split -r >"$o2omaf_tmp"
	echo "time maf-linked $o2omaf_tmp >$o2omaf"
	time maf-linked "$o2omaf_tmp" >"$o2omaf"
	if [ -s "$o2omaf" ]; then
		rm -f "$o2omaf_tmp"
	fi
elif [ -e "$o2omaf_sorted" ]; then
	echo "$o2omaf_sorted already exists"
else
	echo "$o2omaf already exists"
fi
# -H EXPECT: report alignments that are expected by chance at most EXPECT times, in all the sequences. This option requires reading the queries twice (to get their lengths before finding alignments), so it doesn't allow piped-in queries.
# -j4 and --split-f=MAF+: lastal can optionally write "p" lines, indicating the probability that each base is misaligned due to wrong gap placement. last-split, on the other hand, writes "p" lines indicating the probability that each base is aligned to the wrong genomic locus. You can combine both sources of error (roughly) by taking the maximum of the two error probabilities for each base.
# -r: reverse the roles of the two sequences in each alignment: use the 1st(top) sequence as the query and the 2nd(bottom) sequence as the reference.
# maf-linked: reads pair-wise sequence alignments in MAF format, and omits isolated alignments. It keeps groups of alignments that are nearby in both sequences. It may be useful for genome-to-genome alignments: It removes alignments between non-homologous insertions of homologous transposons.
