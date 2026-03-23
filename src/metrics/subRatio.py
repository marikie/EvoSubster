"""
Input:
    - a 3-genome joined alignment .maf file
      (the top sequence should be the outgroup)
Output:
    - the percentage of single-base substitutions (without considering neighboring bases) for the second and third species respectively
"""

import argparse
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "count"))
from Util import getJoinedAlignmentObj


def count(alnFile):
    alnFileHandle = open(alnFile)
    subB = 0
    nonSubB = 0
    subC = 0
    nonSubC = 0
    for aln in getJoinedAlignmentObj(alnFileHandle):
        gSeqA = aln.gSeq1.upper()
        gSeqB = aln.gSeq2.upper()
        gSeqC = aln.gSeq3.upper()
        assert (
            len(gSeqA) == len(gSeqB)
            and len(gSeqA) == len(gSeqC)
            and len(gSeqB) == len(gSeqC)
        ), "gSeqA, gSeqB, and gSeqC should have the same length"
        for i in range(len(gSeqA)):
            gA_base = gSeqA[i]
            gB_base = gSeqB[i]
            gC_base = gSeqC[i]
            # ignore indels
            if gA_base == "-" or gB_base == "-" or gC_base == "-":
                continue
            # count substitutions
            if gA_base == gB_base and gB_base == gC_base:
                nonSubB += 1
                nonSubC += 1
            elif gA_base != gB_base and gB_base == gC_base:
                continue
            elif gA_base != gB_base and gA_base == gC_base:
                subB += 1
                nonSubC += 1
            elif gA_base == gB_base and gB_base != gC_base:
                nonSubB += 1
                subC += 1
            elif gA_base != gB_base and gB_base != gC_base:
                continue
            else:
                error_message = f"Unexpected base combination at position {i}: gSeqA: {gA_base}, gSeqB: {gB_base}, gSeqC: {gC_base}"
                raise ValueError(error_message)
    return subB, nonSubB, subC, nonSubC


def printPercentage(subB, nonSubB, subC, nonSubC):
    # Calculate and print the substitution percentage for species B
    if subB + nonSubB > 0:
        subPercentageB = (subB / (subB + nonSubB)) * 100
        print(f"Substitution percentage for species B: {subPercentageB:.2f}%")
    else:
        error_message = (
            "No valid bases for species B to calculate substitution percentage."
        )
        raise ValueError(error_message)
    # Calculate and print the substitution percentage for species C
    if subC + nonSubC > 0:
        subPercentageC = (subC / (subC + nonSubC)) * 100
        print(f"Substitution percentage for species C: {subPercentageC:.2f}%")
    else:
        error_message = (
            "No valid bases for species C to calculate substitution percentage."
        )
        raise ValueError(error_message)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "joinedAlignmentFile",
        help="a 3-genome joined alignment .maf file (the top sequence should be the outgroup)",
    )
    args = parser.parse_args()
    joinedAlnFile = args.joinedAlignmentFile
    subB, nonSubB, subC, nonSubC = count(joinedAlnFile)
    printPercentage(subB, nonSubB, subC, nonSubC)
