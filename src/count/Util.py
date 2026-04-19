from Alignment import Alignment
from Alignment import JoinedAlignment
from itertools import groupby


def resetCoordToItsStrand(aln, readStrand):
    """
    reset coordinates of the read to
    + strand coordinates if its rStrand == '+',
    - strand coordinates if its rStrand == '-'
    """
    rStart = None
    rEnd = None
    if readStrand == "+":
        if aln.rStrand == "+":
            pass
        # aln.rStrand == '-'
        else:
            rStart, rEnd = convert2CoorOnOppositeStrand(aln)
    # readStrand == '-'
    else:
        if aln.rStrand == "+":
            rStart, rEnd = convert2CoorOnOppositeStrand(aln)
        else:
            pass
    return rStart, rEnd


def filterComments(lines):
    """
    from Alignments.py in JSA. written by Anish Shrestha.
    ignores lines in a maf file begining with #
    """
    for line in lines:
        if line.startswith("#"):
            pass
        else:
            yield line


def getMAFBlock(mafLines):
    """
    from Alignments.py in JSA. written by Anish Shrestha.
    takes in raw lines from a mafFile, filters comment lines,
    return a list containing lines of a maf block
    i.e. ['a ...', 's ...', 's ...']
    """
    for k, groupedLines in groupby(filterComments(mafLines), str.isspace):
        if not k:
            mafBlock = [line for line in groupedLines]
            # print(mafBlock)
            yield mafBlock


def getAln(mafLines):
    for k, groupedLines in groupby(filterComments(mafLines), str.isspace):
        if not k:
            mafBlock = [line for line in groupedLines]
            # print(mafBlock)
            yield Alignment.fromMAFEntry(mafBlock)


def getJoinedAlignmentObj(alnFileHandle):
    """
    takes an alignmentFile and yield a JoinedAlignment object
    """
    for mafEntry in getMAFBlock(alnFileHandle):
        yield JoinedAlignment.fromMAFEntry(mafEntry)


def getAlignmentObjsOneByOne(alnFileHandle):
    """
    takes an alignmentFile and yield Alignment objects one by one
    """
    for mafEntry in getMAFBlock(alnFileHandle):
        yield Alignment.fromMAFEntry(mafEntry)


def getAlignmentObjs(alignmentFile):
    """
    takes an alignmentFile and returns a list of Alignment objects
    """
    alnFileHandle = open(alignmentFile)
    alnObjList = []
    for mafEntry in getMAFBlock(alnFileHandle):
        alnObjList.append(Alignment.fromMAFEntry(mafEntry))
    return alnObjList


def getMultiMAFEntries_all(alignmentFile):
    """
    from Alignments.py in JSA. written by Anish Shrestha.
    takes in mafEntries in list of strings (line) form.
    yield all Alignment objs.
    """
    alnFileHandle = open(alignmentFile, "r")
    # alignments must be in MAF format.
    for rid, mafEntries in groupby(
        getMAFBlock(alnFileHandle), lambda x: x[2].split()[1]
    ):
        mafEntries = list(mafEntries)
        alnObjectList = []
        for mafEntry in mafEntries:
            alnObjectList.append(Alignment.fromMAFEntry(mafEntry))
        yield alnObjectList


def getMultiMAFEntries(alignmentFile):
    """
    ONLY SPLICED READS
    from Alignments.py in JSA. written by Anish Shrestha.
    takes in mafEntries in list of strings (line) form.
    If a query has multiple entries:
    constructs Alignment objects, and yields them.
    """
    alnFileHandle = open(alignmentFile, "r")
    # alignments must be in MAF format.
    for rid, mafEntries in groupby(
        getMAFBlock(alnFileHandle), lambda x: x[2].split()[1]
    ):
        mafEntries = list(mafEntries)
        alnObjectList = []
        if len(mafEntries) > 1:
            for mafEntry in mafEntries:
                alnObjectList.append(Alignment.fromMAFEntry(mafEntry))
            print(type(alnObjectList))
            yield alnObjectList


def convert2CoorOnOppositeStrand(alnObj):
    """
    Convert start and end coordinates of a read
    to start and end coordinates on the reverse strand of read's coordinates
    """
    start = alnObj.rLength - alnObj.rEnd
    end = alnObj.rLength - alnObj.rStart
    return (start, end)


def plusStrandStart(aln):
    if aln.rStrand == "+":
        return aln.rStart
    else:
        return convert2CoorOnOppositeStrand(aln)[0]


def plusStrandEnd(aln):
    if aln.rStrand == "+":
        return aln.rEnd
    else:
        return convert2CoorOnOppositeStrand(aln)[1]


def setToPlusCoord(aln):
    if aln.rStrand == "-":
        # convert coord to + strand coord
        alnStart, alnEnd = convert2CoorOnOppositeStrand(aln)
    else:
        # do nothing
        alnStart, alnEnd = (aln.rStart, aln.rEnd)
    return alnStart, alnEnd


def getIntronCoord(readStrand, aln1, aln2):
    """
    Get which strand of read is aligned (readStrand) and its two alignments.
    Return intron start coordinates and intron end coodinates.
    The coordinates include strand info too.
    # aln1 is 5' side alignment (having donor),
      aln2 is 3' side alignment (having acceptor).
    """
    # assuming aln.gStrand == '+'
    try:
        if not (aln1.gStrand == "+" and aln2.gStrand == "+"):
            raise Exception
    except Exception:
        print('< gStrand = "-" >')
        print(aln1)
        print(aln2)
    # assuming the order of [aln1, aln2]:
    # aln1: 5' side (having donor),
    # aln2: 3' side (having acceptor)
    if (readStrand == "+" and aln1.rStrand == "+") or (
        readStrand == "-" and aln1.rStrand == "-"
    ):
        intronStart = (aln1.gChr, aln1.gEnd, "+")
    else:
        # ((readStrand == '+' and aln1.rStrand == '-') or
        #   (readStrand == '-' and aln1.rStrand == '+'))
        intronStart = (aln1.gChr, aln1.gStart, "-")

    if (readStrand == "+" and aln2.rStrand == "+") or (
        readStrand == "-" and aln2.rStrand == "-"
    ):
        intronEnd = (aln2.gChr, aln2.gStart, "+")
    else:
        # ((readStrand == '+' and aln1.rStrand == '-') or
        #   (readStrand == '-' and aln1.rStrand == '+'))
        intronEnd = (aln2.gChr, aln2.gEnd, "-")

    return (intronStart, intronEnd)


def toSTR(intronCoords):
    """
    convert intronCoords into string
    intronCoords == ((sChr(str), sPos(int), sStrand(str)),
                     (eChr(str), ePos(int), eStrand(str))
    -> sChr_sPos_sStrand_eChr_ePos_eStrand
    """
    startCoord = "_".join(
        [intronCoords[0][0], str(intronCoords[0][1]), intronCoords[0][2]]
    )
    endCoord = "_".join(
        [intronCoords[1][0], str(intronCoords[1][1]), intronCoords[1][2]]
    )
    intronCoords_str = "_".join([startCoord, endCoord])
    return intronCoords_str


def isExactSplit(readStrand, aln1, aln2):
    """
    Return True if the read is exactly spliced
    and there is no unaligned bases between two aligned regions.
    Otherwise return False.

    transcript order is aln1 -> aln2
    AlignmentObj: in-between coord
    """
    # calculate with + strand coord
    if readStrand == "+":
        # calculate aln1End
        if aln1.rStrand == "+":
            aln1End = aln1.rEnd
        else:
            aln1End = aln1.rLength - aln1.rStart

        # calculate aln2Start
        if aln2.rStrand == "+":
            aln2Start = aln2.rStart
        else:
            aln2Start = aln2.rLength - aln2.rEnd

    # calculate with - strand coord
    else:
        # calculate aln1End
        if aln1.rStrand == "-":
            aln1End = aln1.rEnd
        else:
            aln1End = aln1.rLength - aln1.rStart

        # calculate aln2Start
        if aln2.rStrand == "-":
            aln2Start = aln2.rStart
        else:
            aln2Start = aln2.rLength - aln2.rEnd

    diff = aln2Start - aln1End

    if diff == 0:
        return True
    else:
        return False


def getSplitAlignmentData(alignmentFile):
    """
    Input: an alignment .maf file
    Output:
        return only 'Exact Splits' and
        alignments with don-acc pairs
        return a list of ((readStrand, aln1, aln2),
                          (intronLeft, intronRight),
                          (intronStart, intronEnd))
        readStrand: '+'/'-'
        aln: Alignment()
        intron: (aln.gChr, aln.gStart/aln.gEnd, '+'/'-')
    """
    print("--- Reading alignmentFile")
    alignments_list = []
    for alignments in getMultiMAFEntries(alignmentFile):
        # prerequisite:
        # alignments are already sorted
        # according to + strand's coordinates

        readStrand = None
        # get the order of alignments
        # if the first alignment has donor and doesn't have acceptor
        # or the last alignment doesn't have donor and has acceptor
        if (alignments[0].don and not alignments[0].acc) or (
            not alignments[-1].don and alignments[-1].acc
        ):
            # set readStrand to '+'
            readStrand = "+"
        # if the last alignment has donor and doesn't have acceptor
        # or the first alignment doesn't have donor and has acceptor
        elif (alignments[-1].don and not alignments[-1].acc) or (
            not alignments[0].don and alignments[0].acc
        ):
            # set readStrand to '-'
            readStrand = "-"
            # reverse the alignments list
            alignments.reverse()
        else:
            # go to next readID
            # (do NOT append to alignments_list)
            continue

        for aln1, aln2 in zip(alignments, alignments[1:]):
            # if two separate alignments are continuous on the reaad
            # (checking only "Exact Splits")
            # do NOT append alignments with inexact splits
            if isExactSplit(readStrand, aln1, aln2):
                intronStart, intronEnd = getIntronCoord(readStrand, aln1, aln2)
                # print('intronStart: ', intronStart[0],
                #      str(intronStart[1]), intronStart[2])
                # print('intronEnd: ', intronEnd[0],
                #      str(intronEnd[1]), intronEnd[2])
                intronLeft = min(
                    [intronStart, intronEnd], key=lambda c: (c[0], c[1], c[2])
                )
                intronRight = max(
                    [intronStart, intronEnd], key=lambda c: (c[0], c[1], c[2])
                )
                # print('intronLeft: ', intronLeft[0],
                #      str(intronLeft[1]), intronLeft[2])
                # print('intronRight: ', intronRight[0],
                #      str(intronRight[1]), intronRight[2])
                # alignments_list: list of tuple
                alignments_list.append(
                    (
                        (readStrand, aln1, aln2),
                        (intronLeft, intronRight),
                        (intronStart, intronEnd),
                    )
                )

    # sort alignments_list according to
    # intronLeft and intronRight chromosom name and coord
    print("--- Sorting alignments")
    alignments_list.sort(
        key=lambda a_c: ((a_c[1][0][0], a_c[1][0][1]), (a_c[1][1][0], a_c[1][1][1]))
    )
    # print sorted alignments_list
    # for alignments in alignments_list:
    #    intronLeft = alignments[1][0]
    #    intronRight = alignments[1][1]
    #    print('intronLeft:', intronLeft[0],
    #          str(intronLeft[1]), intronLeft[2])
    #    print('intronRight:', intronRight[0],
    #          str(intronRight[1]), intronRight[2])

    return alignments_list


def start(aln):
    return (aln.gChr, aln.gStart)


def end(aln):
    return (aln.gChr, aln.gEnd)


def s1(a1, b1, a2, b2):
    return (a1 - a2) * (b1 - b2)


def s2(a1, b1, a2, b2):
    return (a1 - b2) * (b1 - a2)


def noOverlap(chr1, a1, b1, chr2, a2, b2):
    """
    compare reference coordinates
    check if they don't overlap at all
    """
    if chr1 == chr2:
        if s1(a1, b1, a2, b2) > 0 and s2(a1, b1, a2, b2) > 0:
            return True
        else:
            return False
    else:
        return True


def meetAtPoint(chr1, a1, b1, chr2, a2, b2):
    """
    compare reference coordinates
    check if they meet at a point
    """
    if chr1 == chr2:
        if s1(a1, b1, a2, b2) > 0 and s2(a1, b1, a2, b2) == 0:
            return True
        else:
            return False
    else:
        return False


def exactMatch(chr1, a1, b1, chr2, a2, b2):
    """
    compare reference coordinates
    check if they match exactly
    """
    if chr1 == chr2:
        if s1(a1, b1, a2, b2) == 0 and s2(a1, b1, a2, b2) < 0:
            return True
        else:
            return False
    else:
        return False


def oneIncludesTheOther(chr1, a1, b1, chr2, a2, b2):
    """
    compare reference coordinates
    check if one includes the other
    """
    if chr1 == chr2:
        if s1(a1, b1, a2, b2) <= 0 and s2(a1, b1, a2, b2) < 0:
            return True
        else:
            return False
    else:
        return False


def overlap(chr1, a1, b1, chr2, a2, b2):
    """
    compare reference coordinates
    check if they overlap
    """
    if chr1 == chr2:
        if s1(a1, b1, a2, b2) > 0 and s2(a1, b1, a2, b2) <= 0:
            return True
        else:
            return False
    else:
        return False


def bothEdgeBasesSame(g1Tri, g2Tri, g3Tri):
    """
    Return True if both edge bases are the same.
    """
    return (g1Tri[0] == g2Tri[0] and g2Tri[0] == g3Tri[0]) and (
        g1Tri[2] == g2Tri[2] and g2Tri[2] == g3Tri[2]
    )


def rev(triNuc):
    revDict = {"A": "T", "T": "A", "C": "G", "G": "C"}
    return revDict[triNuc[2]] + revDict[triNuc[1]] + revDict[triNuc[0]]
