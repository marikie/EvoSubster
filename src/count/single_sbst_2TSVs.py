#! /usr/bin/env python3

# Modified by: Martin C. Frith

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
import os
from Util import getJoinedAlignmentObj


#############
# functions
#############


revDict = {"A": "T", "T": "A", "C": "G", "G": "C"}


def initialize_mut_dict():
    """
    key: mutType (e.g. ACGA which means ACG -> AAG)
    value: 0
    """
    mutDict = {}
    letters = "ACGT"
    midLetters = "CT"
    for i in letters:
        for j in midLetters:
            for k in letters:
                for l in letters:
                    if l != j:
                        mutType = i + j + k + l
                        mutDict[mutType] = 0
    return mutDict


def write_output_file(outputFilePath, mutDict, totDict):
    with open(outputFilePath, "w") as tsvfile:
        writer = csv.writer(tsvfile, delimiter="\t", lineterminator="\n")
        writer.writerow(["mutType", "mutNum", "totalRootNum"])
        mutTypeList = sorted(
            list(mutDict.keys()), key=lambda x: (x[1], x[3], x[0], x[2])
        )
        for mutType in mutTypeList:
            originalTriplet = mutType[0:3]
            writer.writerow(
                [
                    mutType[0] + "[" + mutType[1] + ">" + mutType[3] + "]" + mutType[2],
                    mutDict[mutType],
                    totDict[originalTriplet],
                ]
            )


def mutDictFromCounts(counts):
    mutDict = initialize_mut_dict()
    for key, count in counts.items():
        x, y, z, b = key
        if y == "A" or y == "G":
            key = revDict[z] + revDict[y] + revDict[x] + revDict[b]
        mutDict[key] += count
    return mutDict


def totDictFromCounts(counts):
    totDict = collections.Counter()
    for triplet, count in counts.items():
        x, y, z = triplet
        if y == "A" or y == "G":
            triplet = revDict[z] + revDict[y] + revDict[x]
        totDict[triplet] += count
    return totDict


###################
# main procedures
###################
def main(alnFileHandle, outputFilePath2, outputFilePath3):
    originalTripletCounts = collections.Counter()
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
        for i in range(len(gSeq1) - 2):
            x, y, z = gSeq1[i : i + 3]
            a, b, c = gSeq2[i : i + 3]
            d, e, f = gSeq3[i : i + 3]
            # if both edge bases same, and outgroup middle base not unique:
            if x == a and x == d and z == c and z == f and (y == b or y == e):
                # skip cases with gaps or any other non-ACGT symbols:
                if x in "ACGT" and z in "ACGT" and b in "ACGT" and e in "ACGT":
                    originalTriplet = x + y + z
                    originalTripletCounts[originalTriplet] += 1
                    if b != y:
                        mutCounts2[originalTriplet + b] += 1
                    if e != y:
                        mutCounts3[originalTriplet + e] += 1

    alnFileHandle.close()

    mutDict2 = mutDictFromCounts(mutCounts2)
    mutDict3 = mutDictFromCounts(mutCounts3)
    totDict = totDictFromCounts(originalTripletCounts)

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
        outFile2 = f"{org2}.tsv"
        outFile3 = f"{org3}.tsv"
    else:
        outFile2 = f"{org2}_{rest}.tsv"
        outFile3 = f"{org3}_{rest}.tsv"
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
