# Task card: T002 - CDS annotation segmentation

## Goal
Implement deterministic CDS GFF parsing and disjoint strand-state segments so
same-strand isoform overlap is not double-counted and opposite-strand overlap is
explicitly ambiguous.

## Authority / inputs
- Shared context: `.conductor/context-packet.md`
- Own: `src/analysis/cds_intervals.py`, `test/test_cds_intervals.py`

## CLAIM BOUNDARY (hard)
- Allowed: GFF CDS parsing, attribute parsing, segment construction, interval
  querying, and per-seqid metadata needed by the classifier.
- MUST NOT modify MAF modules, classifier, CLI, schemas, docs, or data files.
- Invalid coordinates/strand values must fail explicitly; no silent coercion.
- No git commit. A worker proposes; the conductor promotes.

## Steps
1. Add failing tests for half-open conversion, touching/non-touching boundaries,
   same-strand overlap, duplicate Parent records, opposite-strand ambiguity, and
   absent seqids.
2. Run focused tests and record the expected RED result.
3. Implement a compact disjoint segment index with binary-search entry lookup.
4. Run focused and full test suites.
5. Write `.conductor/handoffs/T002-handoff.md`.

## Machine checks
- `python3 -m unittest test.test_cds_intervals`
- `python3 -m unittest discover -s test -p 'test_*.py'`
- `python3 -m py_compile src/analysis/cds_intervals.py`

## Stop conditions
- Stop on any unclear annotation semantics or failed machine check.
