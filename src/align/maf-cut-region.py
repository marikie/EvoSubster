#! /usr/bin/env python3

"""Split a 3-way (or pairwise) MAF alignment by genomic regions of the TOP
sequence (the outgroup), keeping only the columns whose top-sequence
coordinate falls INSIDE -- or OUTSIDE -- a set of regions given as BED or GFF.

This generalises ``maf-cut-cds-uglier.py`` (which only removes CDS, i.e.
``--keep outside --format gff --type CDS``) to arbitrary functional regions
such as transposable-element (TE) annotations, so the downstream
trinucleotide / tetranucleotide substitution counters can be run separately on
each region class (e.g. TE vs non-TE).

Coordinates
-----------
Regions are interpreted on the TOP sequence's forward strand, in 0-based
half-open ("in-between") coordinates -- the same convention used internally by
the counting code.
  * BED  : columns are (chrom, start, end) already 0-based half-open.
  * GFF  : columns 1/4/5 are (seqid, start, end) 1-based inclusive; converted
           to 0-based half-open as ``beg = start-1, end = end``. A feature-type
           filter (``--type``) selects rows by column 3.

The output is MAF format (it does not bother to line up the columns), written
to stdout, ready to pipe into ``src/count/single_sbst_2TSVs.py`` etc.

Design note (why per-column classification): each alignment column is labelled
inside/outside by the current top-sequence coordinate, and maximal runs with
the requested membership are emitted as sub-alignments. Trinucleotide windows
that straddle a region boundary are split across two emitted blocks and are
therefore dropped from both classes -- exactly as the CDS cutter drops windows
at CDS boundaries -- so ``inside`` and ``outside`` are a clean partition of the
countable windows minus a small, well-defined set of boundary windows.
"""

import argparse
import bisect
import collections
import sys

try:
    from src.align.maf_io import iter_maf_blocks, open_text
except ModuleNotFoundError:
    from maf_io import iter_maf_blocks, open_text


def open_file(file_name):
    return open_text(file_name)


def read_gff_regions(lines, want_type):
    """Yield (seqName, beg, end) from GFF rows of the requested feature type."""
    for line in lines:
        if not line or line[0] == "#":
            continue
        f = line.rstrip("\n").split("\t")
        if len(f) < 5:
            continue
        if want_type is not None and f[2] != want_type:
            continue
        seq_name = f[0]
        beg = int(f[3]) - 1  # 1-based inclusive -> 0-based half-open
        end = int(f[4])
        if beg >= end:
            raise RuntimeError("bad GFF line (begin >= end): " + line.rstrip())
        yield seq_name, beg, end


def read_bed_regions(lines):
    """Yield (seqName, beg, end) from BED rows (already 0-based half-open)."""
    for line in lines:
        if not line or line[0] == "#":
            continue
        if line.startswith(("track", "browser")):
            continue
        f = line.rstrip("\n").split("\t")
        if len(f) < 3:
            continue
        seq_name = f[0]
        beg = int(f[1])
        end = int(f[2])
        if beg >= end:
            raise RuntimeError("bad BED line (begin >= end): " + line.rstrip())
        yield seq_name, beg, end


def merged_ranges(ranges):
    """Yield non-overlapping (beg, end) with overlapping/touching ones merged."""
    it = sorted(ranges)
    if not it:
        return
    cur_beg, cur_end = it[0]
    for b, e in it[1:]:
        if b <= cur_end:  # overlapping or touching
            cur_end = max(cur_end, e)
        else:
            yield cur_beg, cur_end
            cur_beg, cur_end = b, e
    yield cur_beg, cur_end


def padded(beg, end, pad):
    return max(0, beg - pad), end + pad


def shrunk_range(beg, end):
    """Shrink a range by 1 bp on each side (length kept >= 0). Reproduces the
    boundary policy of the legacy maf-cut-cds-uglier.py when --shrink is set."""
    return beg + 1, max(end - 1, beg + 1)


def read_maf_alignments(lines):
    """Yield a MAF block as a list of (seqName, start, strand, seqLen, seq)."""
    for block in iter_maf_blocks(lines):
        yield [
            (row.seqid, row.start0, row.strand, str(row.seqlen), row.sequence)
            for row in block
        ]


def print_maf(seq_names, begs, ends, strands, seq_lengths, seqs, aln_beg, aln_end,
              out=sys.stdout):
    out.write("a\n")
    for name, beg, end, strand, seq_len, seq in zip(
        seq_names, begs, ends, strands, seq_lengths, seqs
    ):
        out.write(
            "s {} {} {} {} {} {}\n".format(
                name, beg, end - beg, strand, seq_len, seq[aln_beg:aln_end]
            )
        )
    out.write("\n")


def cut_alignment(rows, ranges_by_seq, keep, out=sys.stdout):
    """Emit sub-alignments whose TOP-sequence coordinate has the wanted
    membership. ``keep`` is "inside" or "outside"."""
    seq_names, starts, strands, seq_lengths, seqs = zip(*rows)
    if strands[0] == "-":
        raise RuntimeError(
            "can't handle '-' strand for the top sequence in an alignment"
        )
    ranges = ranges_by_seq.get(seq_names[0], [])
    begs = [r[0] for r in ranges]  # for bisect; ranges are sorted, non-overlapping
    want_inside = keep == "inside"
    n_col = len(seqs[0])

    def is_inside(coord):
        k = bisect.bisect_right(begs, coord) - 1
        return k >= 0 and ranges[k][0] <= coord < ranges[k][1]

    coords = list(starts)
    run_start_col = 0
    run_start_coords = list(coords)
    cur_member = is_inside(coords[0]) if ranges else False

    for i in range(n_col):
        top_base = seqs[0][i]
        member_here = is_inside(coords[0]) if (ranges and top_base != "-") else cur_member
        if member_here != cur_member:
            if cur_member == want_inside and i > run_start_col:
                print_maf(seq_names, run_start_coords, coords, strands,
                          seq_lengths, seqs, run_start_col, i, out)
            run_start_col = i
            run_start_coords = list(coords)
            cur_member = member_here
        for k in range(len(seqs)):
            if seqs[k][i] != "-":
                coords[k] += 1

    if cur_member == want_inside and n_col > run_start_col:
        print_maf(seq_names, run_start_coords, coords, strands,
                  seq_lengths, seqs, run_start_col, n_col, out)


def load_regions(args):
    fmt = args.format
    if fmt == "auto":
        low = args.regions.lower()
        if low.endswith((".bed", ".bed.gz")):
            fmt = "bed"
        elif low.endswith((".gff", ".gff3", ".gtf", ".gff.gz", ".gff3.gz")):
            fmt = "gff"
        else:
            fmt = "gff" if args.type else "bed"

    with open_file(args.regions) as fh:
        if fmt == "bed":
            raw = list(read_bed_regions(fh))
        else:
            raw = list(read_gff_regions(fh, args.type))

    ranges_by_seq = collections.defaultdict(list)
    for seq_name, beg, end in raw:
        if args.pad:
            beg, end = padded(beg, end, args.pad)
        ranges_by_seq[seq_name].append((beg, end))
    for k in list(ranges_by_seq):
        ranges = list(merged_ranges(ranges_by_seq[k]))
        if getattr(args, "shrink", False):
            ranges = [shrunk_range(b, e) for b, e in ranges]
        ranges_by_seq[k] = ranges
    return ranges_by_seq, fmt


def main(args):
    ranges_by_seq, fmt = load_regions(args)
    region_seq_names = frozenset(ranges_by_seq)

    total_alns = 0
    matched_alns = 0
    sample_maf_seq_name = None
    maf_fh = open_file(args.maf)
    try:
        for rows in read_maf_alignments(maf_fh):
            total_alns += 1
            top = rows[0][0]
            if sample_maf_seq_name is None:
                sample_maf_seq_name = top
            if top in region_seq_names:
                matched_alns += 1
            cut_alignment(rows, ranges_by_seq, args.keep)
    finally:
        if maf_fh is not sys.stdin:
            maf_fh.close()

    # Guard against a silent no-op: regions present but no MAF top-sequence name
    # matched any region sequence -> almost certainly mismatched assemblies.
    if region_seq_names and total_alns > 0 and matched_alns == 0:
        sample_region = next(iter(region_seq_names))
        sys.stderr.write(
            "Error: no MAF top-sequence name matched any region sequence; "
            "region filtering would be a no-op.\n"
            "  Sample MAF top-seq name: {!r}\n"
            "  Sample region seq name:  {!r}\n"
            "  Common cause: region file and reference FASTA are different "
            "assemblies, or chromosome naming differs (e.g. 'chr1' vs '1').\n"
            .format(sample_maf_seq_name, sample_region)
        )
        sys.exit(1)


if __name__ == "__main__":
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    ap.add_argument("regions", help="BED or GFF file of regions on the top sequence")
    ap.add_argument("maf", help="MAF-format alignment file ('-' for stdin)")
    ap.add_argument("--keep", choices=["inside", "outside"], default="inside",
                    help="keep columns inside the regions (e.g. TE) or outside "
                         "(e.g. non-TE)")
    ap.add_argument("--format", choices=["auto", "bed", "gff"], default="auto",
                    help="region file format (auto: by extension)")
    ap.add_argument("-t", "--type", default=None,
                    help="for GFF, keep only this feature type (column 3), "
                         "e.g. CDS")
    ap.add_argument("--pad", type=int, default=0,
                    help="extend each region by this many bp on both sides")
    ap.add_argument("--shrink", action="store_true",
                    help="shrink each region by 1 bp on both sides; reproduces "
                         "the boundary policy of maf-cut-cds-uglier.py")
    args = ap.parse_args()
    main(args)
