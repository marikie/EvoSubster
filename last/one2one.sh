#!/bin/bash

# Get arguments
DATE=$1
org1FASTA=$2
org2FASTA=$3
dbName=$4
trainFile=$5
m2omaf=$6
o2omaf=$7
o2omaf_maflinked=$8

threadNum=${THREAD_NUM:-8}
# lastdb
echo "---lastdb"
if [ ! -d $dbName ]; then
	echo "making lastdb"
	mkdir $dbName
	cd $dbName
	echo "time lastdb -P${threadNum} -c -uRY4 $dbName $org1FASTA"
	time lastdb -P${threadNum} -c -uRY4 $dbName $org1FASTA
	cd ..
else
	echo "$dbName already exists"
fi
# -P4: makes it faster by using 4 threads (This has no effect on the results.)
# -c: Soft-mask lowercase letters.  This means that, when we compare
#     these sequences to some other sequences using lastal, lowercase
#     letters will be excluded from initial matches.  This will apply
#     to lowercase letters in both sets of sequences.
# -uRY4: selects a seeding scheme that reduces the run time and memory use, but also reduces sensitivity.

# last-train
echo "--last-train"
if [ ! -e $trainFile ]; then
	echo "time last-train -P${threadNum} --revsym -C2 $dbName/$dbName $org2FASTA >$trainFile"
	time last-train -P${threadNum} --revsym -C2 $dbName/$dbName $org2FASTA >$trainFile
else
	echo "$trainFile already exists"
fi
# --revsym: Force the substitution scores to have reverse-complement
#           symmetry, e.g. score(A→G) = score(T→C).  This is often
#           appropriate, if neither strand is "special".
# -C COUNT: Before extending gapped alignments, discard any gapless alignment whose query range lies in COUNT other gapless alignments with higher score-per-length. This aims to reduce run time. -C2 may reduce run time with little effect on accuracy.

# lastal
echo "---lastal"
if [ ! -e $m2omaf ]; then
	echo "time lastal -P${threadNum} -H1 -C2 --split-f=MAF+ -p $trainFile $dbName/$dbName $org2FASTA >$m2omaf"
	time lastal -P${threadNum} -H1 -C2 --split-f=MAF+ -p $trainFile $dbName/$dbName $org2FASTA >$m2omaf
else
	echo "$m2omaf already exists"
fi
#(-j4: show the confidence of each alignment column)
# -H EXPECT: report alignments that are expected by chance at most EXPECT times, in all the sequences. This option requires reading the queries twice (to get their lengths before finding alignments), so it doesn't allow piped-in queries.
# -j4 and --split-f=MAF+: lastal can optionally write "p" lines, indicating the probability that each base is misaligned due to wrong gap placement. last-split, on the other hand, writes "p" lines indicating the probability that each base is aligned to the wrong genomic locus. You can combine both sources of error (roughly) by taking the maximum of the two error probabilities for each base.

# last-split (without maf-linked)
echo "---last-split (without maf-linked)" 
if [ ! -e $o2omaf ]; then
	echo "time last-split -r $m2omaf >$o2omaf"
	time last-split -r $m2omaf >$o2omaf
else
	echo "$o2omaf already exists"
fi
# -r: reverse the roles of the two sequences in each alignment: use the 1st(top) sequence as the query and the 2nd(bottom) sequence as the reference.

# last-split (with maf-linked)
echo "---last-split (with maf-linked)"
if [ ! -e $o2omaf_maflinked ]; then
	echo "maf-linked $o2omaf >$o2omaf_maflinked"
	time maf-linked $o2omaf >$o2omaf_maflinked
else
	echo "$o2omaf_maflinked already exists"
fi
# maf-linked: maf-linked reads pair-wise sequence alignments in MAF format, and omits isolated alignments. It keeps groups of alignments that are nearby in both sequences. It may be useful for genome-to-genome alignments: It removes alignments between non-homologous insertions of homologous transposons

# Cleanup: remove many-to-one MAF after successful downstream generation.
if [ -s "$o2omaf" ] && [ -s "$o2omaf_maflinked" ]; then
	if [ -e "$m2omaf" ]; then
		echo "cleanup: rm -f $m2omaf"
		rm -f "$m2omaf"
	fi
fi
