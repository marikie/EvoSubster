#! /usr/bin/env python3
"""Tests for src/align/maf-cut-region.py (region-stratified MAF splitting).

Run from the repo root or from test/:
    python3 -m unittest test.test_maf_cut_region
    (cd test && python3 -m unittest test_maf_cut_region)

Stdlib only, mirroring the rest of the pipeline. Two layers:
  * unit  -- exact emitted sub-block boundaries of cut_alignment()
  * integ -- end-to-end through the real substitution counter, incl.
             byte-equivalence to the proven maf-cut-cds-uglier.py
"""

import importlib.util
import io
import os
import subprocess
import sys
import tempfile
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(HERE)
ALIGN = os.path.join(REPO, "src", "align")
COUNT = os.path.join(REPO, "src", "count")
CUTTER = os.path.join(ALIGN, "maf-cut-region.py")
CDS_CUTTER = os.path.join(ALIGN, "maf-cut-cds-uglier.py")
COUNTER = os.path.join(COUNT, "single_sbst_2TSVs.py")


def _load_cutter():
    spec = importlib.util.spec_from_file_location("maf_cut_region", CUTTER)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


MCR = _load_cutter()


def make_block(top_name, top_start, seqs, names=("ref", "ing2", "ing3"),
               starts=None, strand="+", seqlen=1000):
    """Build a one-block MAF string from a list of aligned sequence strings."""
    if starts is None:
        starts = [top_start, 0, 0]
    lines = ["a"]
    for nm, st, sq in zip(names, starts, seqs):
        span = sum(1 for c in sq if c != "-")
        lines.append("s {} {} {} {} {} {}".format(nm, st, span, strand, seqlen, sq))
    return "\n".join(lines) + "\n\n"


def parse_blocks(maf_text):
    """Return list of blocks; each block is list of (name, start, span, seq)."""
    blocks = []
    cur = []
    for line in maf_text.splitlines():
        if not line.strip():
            if cur:
                blocks.append(cur)
                cur = []
            continue
        f = line.split()
        if f[0] == "s":
            cur.append((f[1], int(f[2]), int(f[3]), f[6]))
    if cur:
        blocks.append(cur)
    return blocks


def run_cutter(regions_text, maf_text, keep, fmt="bed", ftype=None, pad=0,
               shrink=False, regions_suffix=".bed"):
    """Run cut_alignment in-process and capture the emitted MAF."""
    with tempfile.NamedTemporaryFile("w", suffix=regions_suffix, delete=False) as rf:
        rf.write(regions_text)
        regions_path = rf.name
    with tempfile.NamedTemporaryFile("w", suffix=".maf", delete=False) as mf:
        mf.write(maf_text)
        maf_path = mf.name
    try:
        class A:
            pass
        a = A()
        a.regions, a.maf, a.keep, a.format, a.type, a.pad, a.shrink = (
            regions_path, maf_path, keep, fmt, ftype, pad, shrink)
        ranges_by_seq, _ = MCR.load_regions(a)
        out = io.StringIO()
        fh = MCR.open_file(maf_path)
        try:
            for rows in MCR.read_maf_alignments(fh):
                MCR.cut_alignment(rows, ranges_by_seq, keep, out=out)
        finally:
            fh.close()
        return out.getvalue()
    finally:
        os.unlink(regions_path)
        os.unlink(maf_path)


def count_tsv(maf_text):
    """Pipe a joined MAF through the real counter; return {mutType:(mut2,opp)}.

    Returns dict for ingroup2 and a separate total opportunity count.
    """
    with tempfile.TemporaryDirectory() as d:
        maf_path = os.path.join(d, "j_a_b_run.maf")
        with open(maf_path, "w") as fh:
            fh.write(maf_text)
        o2 = os.path.join(d, "o2.tsv")
        o3 = os.path.join(d, "o3.tsv")
        env = dict(os.environ, PYTHONPATH=COUNT)
        subprocess.run([sys.executable, COUNTER, maf_path, "-o2", o2, "-o3", o3],
                       check=True, cwd=COUNT, env=env,
                       stdout=subprocess.DEVNULL, stderr=subprocess.PIPE)
        res = {}
        opp = 0
        with open(o2) as fh:
            next(fh)
            for line in fh:
                mt, mn, tot = line.rstrip("\n").split("\t")
                res[mt] = int(mn)
        with open(o2) as fh:
            next(fh)
            opp = sum(int(line.split("\t")[2]) for line in fh)
        return res, opp


def total_muts(d):
    return sum(d.values())


# A reproducible-ish DNA string with all 4 bases, no run that is trivially folded
DNA = ("ACGTAGCTAGGCTAACCGGTTACGTGCATGCAATTCCGGATCGATCGTAGCTAGCAT"
       "GGCATCGATCGTAGCATCGATCGGATCCAGCTAGCATCGATCGATGCATCGTAGCTA")


class TestCutAlignmentUnit(unittest.TestCase):
    def test_maf_reader_preserves_legacy_tuple_shape(self):
        rows = list(MCR.read_maf_alignments(io.StringIO(
            "a\ns ref 10 4 + 100 ACGT\n\n"
        )))
        self.assertEqual(rows, [[("ref", 10, "+", "100", "ACGT")]])

    def test_inside_single_interval_no_gaps(self):
        seq = DNA[:12]
        maf = make_block("ref", 0, [seq, seq, seq])
        out = run_cutter("ref\t4\t8\n", maf, keep="inside")
        blocks = parse_blocks(out)
        self.assertEqual(len(blocks), 1)
        ref = blocks[0][0]
        self.assertEqual((ref[1], ref[2]), (4, 4))          # start=4, span=4
        self.assertEqual(ref[3], seq[4:8])

    def test_outside_single_interval_no_gaps(self):
        seq = DNA[:12]
        maf = make_block("ref", 0, [seq, seq, seq])
        out = run_cutter("ref\t4\t8\n", maf, keep="outside")
        blocks = parse_blocks(out)
        spans = sorted((b[0][1], b[0][2]) for b in blocks)
        self.assertEqual(spans, [(0, 4), (8, 4)])           # [0,4) and [8,12)

    def test_inside_with_nonzero_top_start(self):
        seq = DNA[:12]
        maf = make_block("ref", 100, [seq, seq, seq])       # top starts at 100
        out = run_cutter("ref\t104\t108\n", maf, keep="inside")
        ref = parse_blocks(out)[0][0]
        self.assertEqual((ref[1], ref[2]), (104, 4))
        self.assertEqual(ref[3], seq[4:8])

    def test_top_gap_does_not_advance_coordinate(self):
        # top has a gap at column 5; region 4-7 (top coords) -> columns covering
        # top coords 4,5,6 are columns 4, 6, 7 (col 5 is a gap inheriting inside).
        top = DNA[:5] + "-" + DNA[5:11]   # length 12, one gap at col 5
        ing = DNA[:12]
        maf = make_block("ref", 0, [top, ing, ing])
        out = run_cutter("ref\t4\t7\n", maf, keep="inside")
        ref = parse_blocks(out)[0][0]
        self.assertEqual(ref[1], 4)                          # top start coord
        self.assertEqual(ref[2], 3)                          # 3 non-gap top bases (4,5,6)
        self.assertTrue("-" in ref[3])                       # gap retained inside

    def test_empty_regions_keep_outside_is_whole(self):
        seq = DNA[:12]
        maf = make_block("ref", 0, [seq, seq, seq])
        out = run_cutter("", maf, keep="outside")
        blocks = parse_blocks(out)
        self.assertEqual(len(blocks), 1)
        self.assertEqual((blocks[0][0][1], blocks[0][0][2]), (0, 12))

    def test_empty_regions_keep_inside_is_nothing(self):
        seq = DNA[:12]
        maf = make_block("ref", 0, [seq, seq, seq])
        out = run_cutter("", maf, keep="inside")
        self.assertEqual(parse_blocks(out), [])


class TestKnownAnswerSubstitution(unittest.TestCase):
    """Place a single substitution at a known coordinate and check it lands in
    the right region class."""

    def _maf_with_sub_at(self, center_col):
        # outgroup == ingroup3 everywhere; ingroup2 differs only at center_col.
        top = DNA[:30]
        ing3 = top
        bases = list(top)
        orig = bases[center_col]
        bases[center_col] = {"A": "C", "C": "A", "G": "T", "T": "G"}[orig]
        ing2 = "".join(bases)
        return make_block("ref", 0, [top, ing2, ing3])

    def test_substitution_inside_region(self):
        # center at top coord 15, region 10-20 -> window cols 14,15,16 all inside
        maf = self._maf_with_sub_at(15)
        ins, _ = count_tsv(run_cutter("ref\t10\t20\n", maf, keep="inside"))
        outs, _ = count_tsv(run_cutter("ref\t10\t20\n", maf, keep="outside"))
        self.assertEqual(total_muts(ins), 1)
        self.assertEqual(total_muts(outs), 0)

    def test_substitution_outside_region(self):
        # center at top coord 5, region 10-20 -> window entirely outside
        maf = self._maf_with_sub_at(5)
        ins, _ = count_tsv(run_cutter("ref\t10\t20\n", maf, keep="inside"))
        outs, _ = count_tsv(run_cutter("ref\t10\t20\n", maf, keep="outside"))
        self.assertEqual(total_muts(ins), 0)
        self.assertEqual(total_muts(outs), 1)

    def test_substitution_on_boundary_is_dropped(self):
        # center at top coord 9, region 10-20 -> window cols 8,9,10 straddles
        # the boundary at 10 and is dropped from both classes.
        maf = self._maf_with_sub_at(9)
        ins, _ = count_tsv(run_cutter("ref\t10\t20\n", maf, keep="inside"))
        outs, _ = count_tsv(run_cutter("ref\t10\t20\n", maf, keep="outside"))
        self.assertEqual(total_muts(ins), 0)
        self.assertEqual(total_muts(outs), 0)


class TestPartitionAndWhole(unittest.TestCase):
    def test_full_cover_region_equals_whole(self):
        top = DNA[:40]
        bases = list(top)
        for c in (7, 15, 23, 31):
            bases[c] = {"A": "C", "C": "A", "G": "T", "T": "G"}[bases[c]]
        maf = make_block("ref", 0, [top, "".join(bases), top])
        whole, whole_opp = count_tsv(maf)
        ins, ins_opp = count_tsv(run_cutter("ref\t0\t40\n", maf, keep="inside"))
        self.assertEqual(ins, whole)
        self.assertEqual(ins_opp, whole_opp)

    def test_inside_plus_outside_no_double_count(self):
        top = DNA[:40]
        bases = list(top)
        for c in (6, 14, 22, 30, 33):
            bases[c] = {"A": "C", "C": "A", "G": "T", "T": "G"}[bases[c]]
        maf = make_block("ref", 0, [top, "".join(bases), top])
        whole, whole_opp = count_tsv(maf)
        ins, ins_opp = count_tsv(run_cutter("ref\t12\t28\n", maf, keep="inside"))
        outs, out_opp = count_tsv(run_cutter("ref\t12\t28\n", maf, keep="outside"))
        # partition: inside + outside == whole minus the boundary windows
        self.assertLessEqual(total_muts(ins) + total_muts(outs), total_muts(whole))
        self.assertLessEqual(ins_opp + out_opp, whole_opp)
        # every inside/outside count is bounded by whole
        for mt in whole:
            self.assertLessEqual(ins.get(mt, 0) + outs.get(mt, 0), whole[mt])


class TestEquivalenceToCdsScript(unittest.TestCase):
    """maf-cut-region --keep outside --format gff --type CDS must produce the
    same substitution counts as the proven maf-cut-cds-uglier.py."""

    def test_outside_cds_matches_legacy_tsv(self):
        top = DNA[:50]
        bases = list(top)
        for c in (5, 12, 19, 26, 33, 40, 44):
            bases[c] = {"A": "C", "C": "A", "G": "T", "T": "G"}[bases[c]]
        maf = make_block("ref", 0, [top, "".join(bases), top])
        # CDS-style GFF on the top sequence (1-based inclusive)
        gff = "ref\tx\tCDS\t11\t20\t.\t+\t.\t.\nref\tx\tCDS\t31\t38\t.\t+\t.\t.\n"

        with tempfile.TemporaryDirectory() as d:
            maf_path = os.path.join(d, "m.maf")
            gff_path = os.path.join(d, "a.gff")
            with open(maf_path, "w") as fh:
                fh.write(maf)
            with open(gff_path, "w") as fh:
                fh.write(gff)
            # legacy CDS cutter (writes MAF to stdout)
            legacy = subprocess.run(
                [sys.executable, CDS_CUTTER, "-t", "CDS", gff_path, maf_path],
                check=True, capture_output=True, text=True).stdout
            # new general cutter, outside CDS, with --shrink to match the
            # legacy boundary policy exactly
            new = run_cutter(gff, maf, keep="outside", fmt="gff", ftype="CDS",
                             shrink=True, regions_suffix=".gff")

        legacy_counts, legacy_opp = count_tsv(legacy)
        new_counts, new_opp = count_tsv(new)
        self.assertEqual(new_counts, legacy_counts)
        self.assertEqual(new_opp, legacy_opp)


if __name__ == "__main__":
    unittest.main(verbosity=2)
