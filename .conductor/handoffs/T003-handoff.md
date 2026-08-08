# Handoff: T003 - Adversarial cleanup review and resolution verification

task_id: T003
Status: done
claim_boundary_ok: true

## Summary

Reviewed the complete uncommitted implementation and test diff against
`.conductor/context-packet.md` and
`docs/superpowers/specs/2026-08-08-unused-trio-species-cleanup-design.md`.
The focused suites are green, and adversarial probes confirmed that traversal
is rejected before deletion, a parent symlink escaping the configured root is
rejected, a symlink inside a recursively removed owned tree is not followed,
and a retained child under an unused recursive target causes a fatal preflight
error while the selected file survives. Existing focused tests also establish
selected/unlisted-data retention, direct-target symlink unlinking, opt-out
forwarding, and downstream blocking after selector failure.

Concrete findings:

- **[P1] Zero-trio tree mode reports failure after successful cleanup.**
  `src/select/trio_selection.R:969` returns a zero-column `data.frame()`, and
  `src/select/trio_selection.R:1000` serializes it as a one-byte newline. The
  wrapper's size test at `src/sbst_fromDwl.sh:375-379` therefore treats the file
  as nonempty, then `src/sbst_fromDwl.sh:381-394` fails because no accession
  header exists. The adversarial probe measured one byte and one empty line.
  Thus zero-trio owned artifacts are deleted, but the documented wrapper path
  exits nonzero instead of taking its explicit “selected no trios” branch.
  Add an integration regression for a real headerless zero-row selector result
  and either emit a typed zero-row trio table or make the wrapper recognize it.

- **[P1] Cleanup preflight failures leave no current failure audit and can
  preserve a stale success audit.** `src/select/trio_selection.R:293-309`
  creates only the audit directory before path preparation; conflict and type
  checks at `src/select/trio_selection.R:337-373` can stop before any write.
  Only an unlink failure at `src/select/trio_selection.R:386-389` or overall
  success at `src/select/trio_selection.R:395` writes the audit. A seeded prior
  audit remained byte-for-byte `stale-success` after an unsafe `../outside.txt`
  manifest failed, although the outside sentinel correctly survived. This is
  misleading for a destructive default. Initialize/atomically replace the
  audit for each attempt and record unsafe/conflict/type-validation failures
  before propagating the fatal error.

- **[P2] Current-run LAST database directories have no ownership record.**
  `src/align/last_train.sh:15-29` creates `<acc>db_<date>` and its index files,
  while `src/select/trio_selection.R:843-851` records only the `.train` file.
  Zero-trio cleanup therefore leaves current-run database artifacts for unused
  accessions. This conflicts with the overall unused-data/“train caches” claim,
  although the narrow design paragraph names only a train file. The conductor
  should either extend ownership and safe cleanup to newly created lastdb trees
  or narrow README/design claims explicitly.

## Resolution verification

The zero-trio fix at `src/sbst_fromDwl.sh:375-383` now treats the selector's
one-byte blank table as no trios, exits successfully, and starts no downstream
pipeline; both the exact file-shape counterexample and the new wrapper
regression pass. LAST database ownership at
`src/select/trio_selection.R:870-901` is now keyed to the reference accession
and only recorded when the database directory was absent before training. An
exact zero-trio cleanup probe deleted an owned database tree while preserving
an unlisted pre-existing database tree. The new focused ownership and cleanup
assertions pass.

The original unused-row stale-audit counterexample is also fixed: an unsafe
path atomically replaces the prior report with a current `failed` row and the
outside sentinel survives. A follow-up adversarial case initially found that a
selected unsafe artifact remained labeled `retained`. The final fix at
`src/select/trio_selection.R:338-438` tracks the active validation/deletion row
and marks that cause row failed with the fatal reason. The exact selected-row
counterexample now reports `failed`, includes `Unsafe cleanup path`, preserves
the outside sentinel, and the new regression at
`test/test_trio_selection_helpers.R:343-365` passes.

release_blocker_status: CLEARED - no remaining blocker found within T003 scope.

delegation_effect: difficulty=hard; agents=1; useful=confirmed all three
original fixes, found the selected-row audit hole, and verified its final fix;
misses=live NCBI/LAST was not exercised; next_mix=one focused adversarial
reviewer remains useful for destructive-default changes.

## Changed files

- `.conductor/handoffs/T003-handoff.md`

## Machine checks

- `bash test/test_dwl_organism.sh`: PASS
- `Rscript test/test_trio_selection_helpers.R`: PASS
- `bash test/test_sbst_from_dwl_cache.sh`: PASS
- `bash -n src/dwl_organism.sh src/sbst_fromDwl.sh`: PASS
- `git diff --check`: PASS
- Adversarial traversal, root, symlink, retained-overlap probes: PASS
- Adversarial zero-trio wrapper-file probe: PASS
- Adversarial unused-row cleanup-failure audit probe: PASS
- Adversarial LAST database ownership/cleanup probe: PASS
- Adversarial retained-row cleanup-failure audit probe: PASS

## Unresolved risk / escalation

- No unresolved release blocker was found within the assigned scope. Live
  NCBI/LAST execution remains outside this focused fixture-based verification.

## Next action

The conductor should independently confirm the final diff and machine-check
state, then make the release decision and commit if its gate is satisfied.
