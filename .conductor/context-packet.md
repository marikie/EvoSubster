# Context packet - CDS annotation and strand consistency QC

## What this project/phase is
Implement a reproducible, streaming QC that compares independent CDS annotations
for the outgroup and both ingroups in three-way MAF alignments. The implementation
lives on `feature/protein-coding-region`; large inputs and raw results remain in the
external co-author package.

## Binding constraints
- Python 3.8+ standard library only; code and comments are English.
- Use TDD: add a focused failing test and observe the expected failure before production code.
- The three supported trios are `lenEdo`, `puffer`, and `styPis`.
- The four classes are `nonCodingIncluded`, `codingStrandConsistent`,
  `strandContradictions`, and `insufficientEvidence`.
- The primary minimum support is 10 non-gap bases in each ingroup; sensitivity
  thresholds are 1, 10, and 30.
- Same-strand overlapping CDS annotations are a genomic union. Opposite-strand
  overlap is ambiguous and never coerced to a boolean strand match.
- This QC establishes coding/strand consistency only. It does not prove orthology,
  synteny, or absence of paralogy.
- Raw output must stream; never retain all classified runs in memory.
- Do not revert or overwrite unrelated changes. Multiple workers are active.

## Machine checks every worker must leave green
- `python3 -m unittest discover -s test -p 'test_*.py'`
- `python3 -m py_compile <changed Python files>`

## Hard stops
- No git commit; the conductor reviews and commits.
- Escalate scientific or public-interface judgment questions.
- On machine-check failure, stop and report the exact failing check.
