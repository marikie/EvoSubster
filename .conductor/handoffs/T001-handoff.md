# Handoff: T001 - Shared strict MAF I/O implemented

task_id: T001
status: partial
claim_boundary_ok: true

## Summary

Implemented the T003-facing MAF interface in `src/align/maf_io.py` using only
the Python 3.8+ standard library. The module exposes immutable
`MafRow(seqid, start0, span, strand, seqlen, sequence)`, transparent plain/gzip
`open_text(path)`, streaming `iter_maf_blocks(lines)`, `forward_coord(row,
ungapped_offset)`, and `aligned_forward_coords(row)`. Parsed rows retain all
fields needed to verify declared spans and aligned lengths. Sequence rows are
rejected when their field count, integer fields, coordinate bounds, strand, or
declared ungapped span are invalid. A block is rejected when its aligned row
lengths differ. Blank lines, a new `a` header, and EOF delimit non-empty blocks.

Forward coordinates use the source-relative MAF start directly for `+` rows
and `seqlen - 1 - (start0 + ungapped_offset)` for `-` rows. Aligned coordinate
lists preserve alignment columns with `None` at gaps.

Refactored `src/align/maf-cut-region.py` to reuse `open_text` and
`iter_maf_blocks`. Its existing `open_file` and `read_maf_alignments` entry
points remain available, and the latter preserves the legacy tuple shape,
including string-valued `seqLen`. The existing cutter suite remains green and
the CLI can still run directly from the repository root.

TDD order was preserved. The new MAF tests were added and observed RED before
the production module existed. A separate regression test was then observed
RED when the first wrapper version changed legacy `seqLen` from string to
integer; the wrapper was corrected before the focused suite was rerun.

The required full discovery command remains red for non-T001 work: two active
T003 tests require the not-yet-implemented `analyze_trio`, `test_chi2_context`
imports a missing module, and three legacy tests cannot find relative fixture
paths. No classifier, CDS interval, output schema, or data file was modified.

## Changed files

- `src/align/maf_io.py`
- `src/align/maf-cut-region.py`
- `test/test_maf_io.py`
- `test/test_maf_cut_region.py`

## Machine checks

- TDD MAF-module RED observation: PASS
- TDD legacy-tuple regression RED observation: PASS
- `python3 -m unittest test.test_maf_io test.test_maf_cut_region`: PASS
- `python3 -m unittest discover -s test -p 'test_*.py'`: FAIL
- `python3 -m py_compile src/align/maf_io.py src/align/maf-cut-region.py`: PASS
- `python3 src/align/maf-cut-region.py --help`: PASS
- owned-file whitespace checks: PASS

## Unresolved risk / escalation

- Full discovery has six errors outside T001 ownership: two concurrent T003
  RED tests and four legacy test-infrastructure errors. T001 focused tests pass.
- No public-interface judgment remains for T001.

## Next action

Have the conductor inspect the T001 diff and rerun full discovery after T003
finishes; separately resolve or waive the four legacy baseline errors.

delegation_effect: difficulty=simple; agents=0; useful=bounded shared parser and
coordinate API with compatibility evidence; misses=full-suite baseline remained
red; next_mix=one conductor reviewer plus the T003 and legacy-test owners
