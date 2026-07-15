# Task card: T001 - Shared MAF I/O and coordinate helpers

## Goal
Create an importable stdlib MAF module with strict parsing and correct forward
coordinates for both strands. Refactor the existing region cutter to reuse it
without changing cutter behavior.

## Authority / inputs
- Shared context: `.conductor/context-packet.md`
- Own: `src/align/maf_io.py`, `test/test_maf_io.py`
- May modify only as needed for reuse: `src/align/maf-cut-region.py`,
  `test/test_maf_cut_region.py`

## CLAIM BOUNDARY (hard)
- Allowed: define `MafRow`, transparent gzip/plain opening, block iteration,
  row validation, `forward_coord`, and aligned coordinate iteration.
- MUST NOT change classifier, CDS interval modules, output schemas, or data files.
- Preserve current `maf-cut-region.py` output byte-for-byte for existing tests.
- No git commit. A worker proposes; the conductor promotes.

## Steps
1. Add failing tests for plus/minus coordinates, gaps, gzip input, malformed rows,
   span mismatch, unequal aligned lengths, and block boundaries.
2. Run the focused tests and record the expected RED result.
3. Implement the minimum shared module and refactor the cutter to import it.
4. Run focused and full test suites.
5. Write `.conductor/handoffs/T001-handoff.md`.

## Machine checks
- `python3 -m unittest test.test_maf_io test.test_maf_cut_region`
- `python3 -m unittest discover -s test -p 'test_*.py'`
- `python3 -m py_compile src/align/maf_io.py src/align/maf-cut-region.py`

## Stop conditions
- Stop on any behavior regression or judgment question.
