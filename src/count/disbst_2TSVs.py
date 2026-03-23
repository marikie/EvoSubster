#! /usr/bin/env python3

# Written by: Martin C. Frith

"""
Input:
    - a 3-genome joined alignment .maf file
      (the top sequence should be the outgroup)
    - outputFilePath1
    - outputFilePath2
Output:
    - two tsv files with the following columns:
        - mutation type
        - mutNum
        - totalRootNum
"""

import argparse
import collections
import csv
import itertools
import os
from Util import getJoinedAlignmentObj


#############
# functions
#############


revDict = {"A": "T", "T": "A", "C": "G", "G": "C"}

# We don't want to count any original doublet separately from its
# reverse complement.  So use the reverse-complement of these (arbitrary?) ones
reversedOriginalDoublets = "AA GT AG CA GG GA"


def initialize_mut_dict():
    """
    key: mutType (e.g. ACGA which means AC -> GA)
    value: 0
    """
    mutDict = {}

    for x, y, b, c in itertools.product("ACGT", repeat=4):
        originalDoublet = x + y
        mutatedDoublet = b + c
        # require both bases to be different:
        if x == b or y == c:
            continue
        # arbitrarily(?) exclude cases whose reverse-complement is included:
        if (
            originalDoublet in reversedOriginalDoublets
            or originalDoublet == "AT"
            and mutatedDoublet in "GG TC TG"
            or originalDoublet == "CG"
            and mutatedDoublet in "AA AC GA"
            or originalDoublet == "GC"
            and mutatedDoublet in "CT TG TT"
            or originalDoublet == "TA"
            and mutatedDoublet in "AC AG CC"
        ):
            continue
        mutType = originalDoublet + mutatedDoublet
        mutDict[mutType] = 0
    return mutDict


def write_output_file(outputFilePath, mutDict, totDict):
    with open(outputFilePath, "w") as tsvfile:
        writer = csv.writer(tsvfile, delimiter="\t", lineterminator="\n")
        writer.writerow(["mutType", "mutNum", "totalRootNum"])
        for mutType in sorted(mutDict):
            originalDoublet = mutType[0:2]
            writer.writerow(
                [
                    originalDoublet + ">" + mutType[2:4],
                    mutDict[mutType],
                    totDict[originalDoublet],
                ]
            )


def mutDictFromCounts(counts):
    mutDict = initialize_mut_dict()
    for key, count in counts.items():
        if key not in mutDict:
            x, y, b, c = key
            key = revDict[y] + revDict[x] + revDict[c] + revDict[b]
        mutDict[key] += count
    return mutDict


def totDictFromCounts(counts):
    totDict = collections.Counter()
    for doublet, count in counts.items():
        if doublet in reversedOriginalDoublets:
            x, y = doublet
            doublet = revDict[y] + revDict[x]
        totDict[doublet] += count
    return totDict


###################
# main procedures
###################
def main(alnFileHandle, outputFilePath2, outputFilePath3):
    originalDoubletCounts = collections.Counter()
    mutCounts2 = collections.Counter()
    mutCounts3 = collections.Counter()

    for aln in getJoinedAlignmentObj(alnFileHandle):
        gSeq1 = aln.gSeq1.upper()
        gSeq2 = aln.gSeq2.upper()
        gSeq3 = aln.gSeq3.upper()
        assert (
            len(gSeq1) == len(gSeq2)
            and len(gSeq1) == len(gSeq3)
            and len(gSeq2) == len(gSeq3)
        ), "gSeq1, gSeq2, and gSeq3 should have the same length"
        for i in range(len(gSeq1) - 3):
            # Maybe we could ignore edge bases, but this code checks
            # they're the same, which might avoid unreliable alignment(?)
            w, x, y, z = gSeq1[i : i + 4]
            a, b, c, d = gSeq2[i : i + 4]
            e, f, g, h = gSeq3[i : i + 4]
            # if both edge bases same
            # and outgroup middle bases not unique
            # and there are no gaps or non-ACGT symbols:
            if (
                w == a
                and w == e
                and z == d
                and z == h
                and (x + y == b + c or x + y == f + g)
                and w in "ACGT"
                and z in "ACGT"
                and b in "ACGT"
                and c in "ACGT"
                and f in "ACGT"
                and g in "ACGT"
            ):
                originalDoublet = x + y
                originalDoubletCounts[originalDoublet] += 1
                if b != x and c != y:
                    mutCounts2[originalDoublet + b + c] += 1
                if f != x and g != y:
                    mutCounts3[originalDoublet + f + g] += 1

    alnFileHandle.close()

    mutDict2 = mutDictFromCounts(mutCounts2)
    mutDict3 = mutDictFromCounts(mutCounts3)
    totDict = totDictFromCounts(originalDoubletCounts)

    write_output_file(outputFilePath2, mutDict2, totDict)
    write_output_file(outputFilePath3, mutDict3, totDict)


def get_default_output_file_names(joinedAlnFile):
    # file name without extension
    filename = os.path.splitext(os.path.basename(joinedAlnFile))[0]
    # Get the path before the filename at the end
    path_before_filename = os.path.dirname(joinedAlnFile)

    filename_parts = filename.split("_")
    if len(filename_parts) < 3:
        raise ValueError("The file name should be org1_org2_org3_*.maf")
    org2 = filename_parts[1]
    org3 = filename_parts[2]
    rest = "_".join(filename_parts[3:])
    if rest == "":
        outFile2 = f"{org2}_dinuc.tsv"
        outFile3 = f"{org3}_dinuc.tsv"
    else:
        outFile2 = f"{org2}_{rest}_dinuc.tsv"
        outFile3 = f"{org3}_{rest}_dinuc.tsv"
    outputFilePath2 = os.path.join(path_before_filename, outFile2)
    outputFilePath3 = os.path.join(path_before_filename, outFile3)
    return outputFilePath2, outputFilePath3


if __name__ == "__main__":
    ###################
    # parse arguments
    ###################
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "joinedAlignmentFile",
        help="a 3-genome joined alignment .maf file (the top sequence should be the outgroup), the file name should be org1_org2_org3_*.maf",
    )
    parser.add_argument(
        "-o2",
        "--outputFilePath2",
        help="output file path for organism2",
    )
    parser.add_argument(
        "-o3",
        "--outputFilePath3",
        help="output file path for organism3",
    )
    args = parser.parse_args()
    joinedAlnFile = args.joinedAlignmentFile
    outputFilePath2 = (
        args.outputFilePath2 or get_default_output_file_names(joinedAlnFile)[0]
    )
    outputFilePath3 = (
        args.outputFilePath3 or get_default_output_file_names(joinedAlnFile)[1]
    )

    # print(f"joinedAlnFile: {joinedAlnFile}")
    # print(f"outputFilePath2: {outputFilePath2}")
    # print(f"outputFilePath3: {outputFilePath3}")
    # exit(1)

    ###################
    # main
    ###################
    alnFileHandle = open(joinedAlnFile)
    main(alnFileHandle, outputFilePath2, outputFilePath3)
