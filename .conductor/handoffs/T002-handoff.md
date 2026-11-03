# Handoff: T002 - Selector cleanup and wrapper propagation

task_id: T002
Status: done
claim_boundary_ok: true

## Summary

Implemented current-run artifact cleanup in `src/select/trio_selection.R` and
the tree-mode opt-out propagation in `src/sbst_fromDwl.sh`. The selector now
requests the fixed downloader manifest contract only when cleanup is enabled,
validates the exact `accession`, `artifact_type`, and `relative_path` columns,
and aggregates manifest rows by accession without changing the downloader's
six-field stdout handling. `make_fetchers()` separately records a train cache
only after the current run creates or replaces it and
`parse_train_identity()` confirms a complete valid result. Valid pre-existing
cache hits are never added to the owned-artifact table.

After `selected_trios.tsv` is written, the selector derives the selected set
from `out_acc`, `in1_acc`, and `in2_acc`. Current-run genome artifacts are
retained only for selected accessions, and current-run train files are retained
only when both accessions are selected. A zero-row trio result therefore
deletes all current-run-owned artifacts. Unlisted pre-existing paths are never
cleanup targets. The selector writes `cleanup_unused_species.tsv` with
accessions, artifact type, normalized path, action, and reason.

Cleanup rejects empty, absolute, dotted, traversal, root-resolving, and
parent-symlink escape paths before deleting anything. Targets are constructed
under normalized genome or train-cache roots; target types are verified without
following symlinks. Files and symlinks are removed before deepest-first
directories, while only `directory_tree` records permit recursive removal.
Conflicts between retained artifacts and recursive targets are fatal. Any
deletion status failure or surviving target writes a failed audit entry and
raises an error, so the existing wrapper selector-failure branch exits before
any downstream `sbst.sh` call.

Added `--keep-unused-species-data` to the selector and wrapper help and
forwarding path. It disables manifest requests, train ownership tracking, and
post-selection cleanup. `README.md` now documents the default cleanup,
preservation rules, audit location, opt-out, and failure behavior. The fetched
Stage 0 baseline regression was also fixed in `test/test_stage0_gate.R` by
explicitly loading `dplyr` before sourcing the selector.

TDD evidence was observed before production edits: the helper test failed
because the ownership/cleanup API and help option did not exist; the wrapper
test failed because the opt-out was rejected as unknown; and the Stage 0
baseline failed because `%>%` was unavailable. Focused tests now cover selected
retention, unused deletion, pre-existing retention, generated cache cleanup,
zero trios, traversal rejection, symlink non-following, incomplete deletion,
opt-out forwarding, and selector-failure propagation.

delegation_effect: difficulty=hard; agents=1; useful=isolated selector and
wrapper implementation with destructive-path tests; misses=no known misses
within T002 scope; next_mix=one adversarial integration reviewer for downloader
manifest consumption and cleanup overlap cases.

## Changed files

- `src/select/trio_selection.R`
- `src/sbst_fromDwl.sh`
- `README.md`
- `test/test_trio_selection_helpers.R`
- `test/test_sbst_from_dwl_cache.sh`
- `test/test_stage0_gate.R`
- `.conductor/handoffs/T002-handoff.md`

## Machine checks

- TDD helper and wrapper RED observation: PASS
- Stage 0 deferred-dplyr RED observation: PASS
- `bash test/test_dwl_organism.sh`: PASS
- `Rscript test/test_trio_selection_helpers.R`: PASS
- `Rscript test/test_stage0_gate.R`: PASS
- `Rscript test/test_outgroup_selection.R`: PASS
- `bash test/test_sbst_from_dwl_cache.sh`: PASS
- `bash -n src/dwl_organism.sh src/sbst_fromDwl.sh`: PASS
- `git diff --check`: PASS

## Unresolved risk / escalation

- No judgment escalation is required. Live NCBI/LAST execution was not run;
  focused process-boundary fixtures and the downloader manifest suite cover the
  implemented contracts.

## Next action

The conductor should have T003 adversarially review the combined T001/T002 diff,
then rerun the load-bearing checks and make the release decision. Do not commit
T002 independently from the downloader manifest implementation.
