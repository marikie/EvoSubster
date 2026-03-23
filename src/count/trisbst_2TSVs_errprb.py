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
import csv
import os
from Util import getJoinedAlignmentObj

oriDict = {
    "ACA": ["ACAA", "ACAG", "ACAT"],
    "ACC": ["ACCA", "ACCG", "ACCT"],
    "ACG": ["ACGA", "ACGG", "ACGT"],
    "ACT": ["ACTA", "ACTG", "ACTT"],
    "CCA": ["CCAA", "CCAG", "CCAT"],
    "CCC": ["CCCA", "CCCG", "CCCT"],
    "CCG": ["CCGA", "CCGG", "CCGT"],
    "CCT": ["CCTA", "CCTG", "CCTT"],
    "GCA": ["GCAA", "GCAG", "GCAT"],
    "GCC": ["GCCA", "GCCG", "GCCT"],
    "GCG": ["GCGA", "GCGG", "GCGT"],
    "GCT": ["GCTA", "GCTG", "GCTT"],
    "TCA": ["TCAA", "TCAG", "TCAT"],
    "TCC": ["TCCA", "TCCG", "TCCT"],
    "TCG": ["TCGA", "TCGG", "TCGT"],
    "TCT": ["TCTA", "TCTG", "TCTT"],
    "ATA": ["ATAA", "ATAC", "ATAG"],
    "ATC": ["ATCA", "ATCC", "ATCG"],
    "ATG": ["ATGA", "ATGC", "ATGG"],
    "ATT": ["ATTA", "ATTC", "ATTG"],
    "CTA": ["CTAA", "CTAC", "CTAG"],
    "CTC": ["CTCA", "CTCC", "CTCG"],
    "CTG": ["CTGA", "CTGC", "CTGG"],
    "CTT": ["CTTA", "CTTC", "CTTG"],
    "GTA": ["GTAA", "GTAC", "GTAG"],
    "GTC": ["GTCA", "GTCC", "GTCG"],
    "GTG": ["GTGA", "GTGC", "GTGG"],
    "GTT": ["GTTA", "GTTC", "GTTG"],
    "TTA": ["TTAA", "TTAC", "TTAG"],
    "TTC": ["TTCA", "TTCC", "TTCG"],
    "TTG": ["TTGA", "TTGC", "TTGG"],
    "TTT": ["TTTA", "TTTC", "TTTG"],
}

revDict = {"A": "T", "T": "A", "C": "G", "G": "C"}


#############
# functions
#############
def bothEdgeBasesSame(g1Tri, g2Tri, g3Tri):
    return (g1Tri[0] == g2Tri[0] and g2Tri[0] == g3Tri[0]) and (
        g1Tri[2] == g2Tri[2] and g2Tri[2] == g3Tri[2]
    )


def lowErrProb(pLines, i):
    asciiThreshold = 53  # ord(char) >= 53 (char >= 5) -> error probability =< 0.01
    if all(
        ord(char) >= asciiThreshold for pLine in pLines for char in pLine[i : i + 3]
    ):
        return True
    else:
        return False


def rev(triNuc):
    return revDict[triNuc[2]] + revDict[triNuc[1]] + revDict[triNuc[0]]


def mutType(gTri, ori, mut):
    if ori in set(["C", "T"]):
        return gTri[0] + ori + gTri[2] + mut
    else:
        return revDict[gTri[2]] + revDict[ori] + revDict[gTri[0]] + revDict[mut]


# ACTG -> ACT
# ACGG -> ACG
def ori(mutType):
    return mutType[0:3]


def add2totalNum(mutDict, triNuc):
    if triNuc[1] == "A" or triNuc[1] == "G":
        oriNuc = rev(triNuc)
    else:
        oriNuc = triNuc

    for mutType in oriDict[oriNuc]:
        mutDict[mutType]["totalRootNum"] += 1


def add2MutDict(g1Tri, g2Tri, g3Tri, mutDict2, mutDict3):
    """
    If the middle bases are all different to each other,
    do nothing.
    If the mutation happened on genome2,
    add the mutation count to mutDict2 (N[majority > minority]N),
    add the total count to mutDict2 (N(majority)N),
    and add the total count to mutDict3 (N(majority)N).
    If the mutation happened on genome3,
    add the mutation count to mutDict3 (N[majority > minority]N),
    add the total count to mutDict3 (N(majority)N),
    and add the total count to mutDict2 (N(majority)N).
    If the mutation happend on genome1,
    do nothing.
    """
    # print("add2MutDict called")
    # print("g1Tri: ", g1Tri)
    # print("g2Tri: ", g2Tri)
    # print("g3Tri: ", g3Tri)

    middleList = [g1Tri[1], g2Tri[1], g3Tri[1]]
    majority = ""
    minority = ""
    # print("middleList: ", middleList)
    # print("set(middleList): ", set(middleList))

    if len(set(middleList)) == 1:
        # only original count
        add2totalNum(mutDict2, g1Tri)
        add2totalNum(mutDict3, g1Tri)
    # if there is a mutation on org1, org2, or org3
    elif len(set(middleList)) == 2:
        for base in set(middleList):
            if middleList.count(base) == 2:
                majority = base
            elif middleList.count(base) == 1:
                minority = base
            else:
                raise (Exception)
        # print("majority: ", majority, "minority: ", minority)
        # if the mutation is on org1
        # ambiguous: minority > majority or majority > minority
        if g1Tri[1] == minority:
            # print("g1Tri[1] == minority")
            pass
        # if the mutation is on org2
        # majority > minority on mutDict2
        elif g2Tri[1] == minority:
            # print("g2Tri[1] == minority")
            # mutation count
            # print(
            #     "mutType(g2Tri, majority, minority): ",
            #     mutType(g2Tri, majority, minority),
            # )
            mutDict2[mutType(g2Tri, majority, minority)]["mutNum"] += 1
            # total count (original is g1Tri) to mutDict2 and mutDict3
            add2totalNum(mutDict2, g1Tri)
            add2totalNum(mutDict3, g1Tri)
        # if the mutation is on org3
        # majority > minority on mutDict3
        elif g3Tri[1] == minority:
            # print("g3Tri[1] == minority")
            # mutation count
            mutDict3[mutType(g3Tri, majority, minority)]["mutNum"] += 1
            # total count (original is g1Tri) to mutDict2 and mutDict3
            add2totalNum(mutDict2, g1Tri)
            add2totalNum(mutDict3, g1Tri)


def initialize_mut_dict():
    """
    key: mutType (e.g. ACGA which means ACG -> AAG)
    value: {
        "mutNum": 0,
        "totalRootNum": 0
    }
    """
    mutDict = {}
    letters = ["A", "C", "G", "T"]
    midLetters = ["C", "T"]
    cSubs = ["A", "G", "T"]
    tSubs = ["A", "C", "G"]
    for i in letters:
        for j in midLetters:
            if j == "C":
                for k in cSubs:
                    for l in letters:
                        mtype = i + j + l + k
                        mutDict[mtype] = {"mutNum": 0, "totalRootNum": 0}
            if j == "T":
                for k in tSubs:
                    for l in letters:
                        mtype = i + j + l + k
                        mutDict[mtype] = {"mutNum": 0, "totalRootNum": 0}
    return mutDict


def write_output_file(outputFilePath, mutDict):
    with open(outputFilePath, "w") as tsvfile:
        writer = csv.writer(tsvfile, delimiter="\t", lineterminator="\n")
        writer.writerow(["mutType", "mutNum", "totalRootNum"])
        mutTypeList = sorted(
            list(mutDict.keys()), key=lambda x: (x[1], x[3], x[0], x[2])
        )
        for mutType in mutTypeList:
            writer.writerow(
                [
                    mutType[0] + "[" + mutType[1] + ">" + mutType[3] + "]" + mutType[2],
                    mutDict[mutType]["mutNum"],
                    mutDict[mutType]["totalRootNum"],
                ]
            )


###################
# main procedures
###################
def main(alnFileHandle, outputFilePath2, outputFilePath3):
    mutDict2 = initialize_mut_dict()
    mutDict3 = initialize_mut_dict()
    # print("mutDict2: ", len(mutDict2.keys()))

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
            g1Tri = gSeq1[i : i + 3]
            g2Tri = gSeq2[i : i + 3]
            g3Tri = gSeq3[i : i + 3]
            # if indels are included, go next
            # if not all (i in "ACGT" for i in gTri)
            if not (
                all(list(map(lambda b: b in set(["A", "C", "G", "T"]), g1Tri)))
                and all(list(map(lambda b: b in set(["A", "C", "G", "T"]), g2Tri)))
                and all(list(map(lambda b: b in set(["A", "C", "G", "T"]), g3Tri)))
            ):
                continue
            # if there are no indels
            # and the edge bases are the same,
            # and the column error probability is less than 0.01,
            # count original trinucs as total and count mutations if there is a mutation
            if bothEdgeBasesSame(g1Tri, g2Tri, g3Tri) and lowErrProb(aln.pLines, i):
                add2MutDict(g1Tri, g2Tri, g3Tri, mutDict2, mutDict3)

    alnFileHandle.close()

    # write to outputFilePath2
    write_output_file(outputFilePath2, mutDict2)
    # write to outputFilePath3
    write_output_file(outputFilePath3, mutDict3)


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

    ###################
    # main
    ###################
    alnFileHandle = open(joinedAlnFile)
    main(alnFileHandle, outputFilePath2, outputFilePath3)
