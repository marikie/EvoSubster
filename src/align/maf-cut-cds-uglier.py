#! /usr/bin/env python3

"""Get parts of alignments (in MAF format) where the top sequence is
outside CDS regions.

The output is MAF format, except it doesn't bother to line up the columns.
"""

import argparse
import bisect
import collections
import gzip
import itertools
import sys

def openFile(fileName):
    if fileName == "-":
        return sys.stdin
    if fileName.endswith(".gz"):
        return gzip.open(fileName, "rt")
    return open(fileName)

def shrunkRange(beg, end):  # try to shrink the range, but ensure length >= 0
    return beg + 1, max(end - 1, beg + 1)

def mergedRanges(ranges):  # yield ranges with overlapping/touching ones merged
    beg = end = 0
    for b, e in sorted(ranges):
        if b > end:
            if end > 0:
                yield beg, end
            beg = b
        end = max(end, e)
    yield beg, end

def readAnnotations(typeOfThing, lines):
    for line in lines:
        f = line.split("\t")
        if f and f[0][0] != "#" and f[2] == typeOfThing:
            seqName = f[0]
            beg = int(f[3]) - 1  # convert to in-between coordinates
            end = int(f[4])
            if beg >= end:
                raise RuntimeError("bad GFF line: begin > end")
            yield seqName, beg, end

def readMafAlignments(lines):
    rows = []
    for line in lines:
        f = line.split()
        if not f:  # blank line: separates alignments
            if rows: yield rows
            rows = []
        elif f[0] == "s":
            s, seqName, start, span, strand, seqLen, seq = f
            rows.append((seqName, int(start), strand, seqLen, seq))
    if rows: yield rows

def printMaf(seqNames, begs, ends, strands, seqLengths, seqs, alnBeg, alnEnd):
    print("a")
    z = zip(seqNames, begs, ends, strands, seqLengths, seqs)
    for seqName, beg, end, strand, seqLen, seq in z:
        print("s", seqName, beg, end - beg, strand, seqLen, seq[alnBeg:alnEnd])
    print()

def cutAlignment(alignmentRows, rangesPerSequence):
    seqNames, starts, strands, seqLengths, seqs = zip(*alignmentRows)
    if strands[0] == "-":
        raise RuntimeError("sorry, can't handle '-' strand"
                           "for the top sequence in an alignment")
    ranges = rangesPerSequence[seqNames[0]]
    n = len(ranges)
    j = bisect.bisect(ranges, (starts[0], starts[0]))
    coords = list(starts)
    alnBeg = 0
    alignmentColumns = zip(*seqs)
    for i, col in enumerate(alignmentColumns):
        if j < n and coords[0] == ranges[j][0] and col[0] != "-":
            if i > 0:
                printMaf(seqNames, starts, coords, strands, seqLengths, seqs,
                         alnBeg, i)
            j += 1
        if j > 0 and coords[0] == ranges[j-1][1] and col[0] != "-":
            starts = coords.copy()
            alnBeg = i
        for k, x in enumerate(col):
            if x != "-":
                coords[k] += 1
    if j < 1 or coords[0] > ranges[j-1][1]:
        printMaf(seqNames, starts, coords, strands, seqLengths, seqs,
                 alnBeg, i+1)

def main(args):
    rangesPerSequence = collections.defaultdict(list)
    for seqName, beg, end in readAnnotations(args.type, openFile(args.gff)):
        rangesPerSequence[seqName].append((beg, end))
    for k, v in rangesPerSequence.items():
        v = mergedRanges(v)
        v = itertools.starmap(shrunkRange, v)
        rangesPerSequence[k] = list(v)

    gffSeqNames = frozenset(rangesPerSequence)
    totalAlns = 0
    matchedAlns = 0
    sampleMafSeqName = None
    for a in readMafAlignments(openFile(args.maf)):
        totalAlns += 1
        topSeqName = a[0][0]
        if sampleMafSeqName is None:
            sampleMafSeqName = topSeqName
        if topSeqName in gffSeqNames:
            matchedAlns += 1
        cutAlignment(a, rangesPerSequence)

    # If GFF had annotations but no MAF alignment referenced a matching
    # sequence, the reference FASTA and GFF almost certainly describe different
    # assemblies (e.g. cds_from_genomic.fna was used instead of the genomic
    # FASTA). Fail loudly -- otherwise the "ncds" MAF is identical to the input.
    if gffSeqNames and totalAlns > 0 and matchedAlns == 0:
        sampleGffSeqName = next(iter(gffSeqNames))
        sys.stderr.write(
            "Error: No MAF top-sequence matched any GFF sequence;"
            " CDS filtering is a no-op.\n"
            "  Sample MAF top-seq name: {!r}\n"
            "  Sample GFF seq name:     {!r}\n"
            "  Common cause: the reference alignment used a non-genomic FASTA\n"
            "  (e.g. cds_from_genomic.fna) or a GFF from a different assembly.\n"
            .format(sampleMafSeqName, sampleGffSeqName))
        sys.exit(1)

if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=
                                 argparse.ArgumentDefaultsHelpFormatter)
    ap.add_argument("-t", "--type", default="CDS", help="type of annotation")
    ap.add_argument("gff", help="GFF-format file of annotations")
    ap.add_argument("maf", help="MAF-format file of alignments")
    args = ap.parse_args()
    main(args)
