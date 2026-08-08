# Handoff: T007 - Integrate and verify the Stage 0 refactor

task_id: T007
status: done
claim_boundary_ok: true

## Summary

The Stage 0 assembly-selection refactor is complete on `feature/stage0-reference-quality-ranking`, relative to baseline `96787c57a0e18baf69f4e48f59be2e7316afffbd`. The implementation changes only internal structure and characterization tests. It does not change the CLI, external metadata inputs, output schemas, ranking fields, selection policy, validation rules, or reason-code vocabulary.

`src/select/assembly_selection.R` now keeps `rank_assembly_candidates()` focused on orchestration: candidate preparation, per-species shortlisting, baseline selection, final-rank assignment, and result projection. Ranking-state initialization moved to `stage0_initialize_ranking_state()`. The eukaryote audit/optional-override state transitions moved behind `stage0_process_eukaryote_qc()`, with focused helpers for QC preference, distinct assembly-equivalence classes, audit reasons, blocker calculation, and strict override application. The five-tier reference-first baseline functions and the prokaryote ranking path were not modified.

`test/test_stage0_quality_selection.R` adds two characterization cases recommended by the design review. One fixes the behavior when two alternative assemblies both satisfy strict override criteria: the strongest BUSCO candidate is selected and the other is recorded as `not_top_override_candidate`. The other fixes the exact stable order of simultaneous blocker codes: `not_higher_single_copy;higher_duplicated;higher_fragmented`.

The independent T006 adversarial reviewer reported no Critical, Important, or Minor findings. In addition to the named regression suites, that reviewer compared the complete list returned by the baseline and refactored modules for 1,800 generated combinations spanning eukaryote/prokaryote profiles, one-to-five candidates, three shortlist sizes, override on/off, references, annotations, hard exclusions, complete BUSCO evidence, and paired accessions. Every comparison was identical. The conductor then reran the complete relevant test set on the final working-tree state.

The final diff was also inspected to confirm that validation and ordering functions outside the extracted region remained unchanged.

Unrelated untracked paths `analysis/`, `stage0-processing-flow.html`, and `stage1-matching-algorithm.html` remain untouched and are outside this handoff. No remote state was changed.

delegation_effect: difficulty=medium; agents=2; useful=design review identified two missing characterization cases and adversarial review supplied 1800 exact baseline-current comparisons; misses=no production-data replay was performed; next_mix=one design worker and one differential reviewer remains appropriate for mutation-heavy ranking refactors.

## Changed files

- `src/select/assembly_selection.R`
- `test/test_stage0_quality_selection.R`
- `.conductor/context-refactor-stage0.md`
- `.conductor/task-cards/T005-stage0-refactor-design.md`
- `.conductor/task-cards/T006-stage0-refactor-review.md`
- `.conductor/task-cards/T007-stage0-refactor-integration.md`
- `.conductor/handoffs/T005-handoff.md`
- `.conductor/handoffs/T006-handoff.md`
- `.conductor/handoffs/T007-handoff.md`

## Machine checks

- `Rscript test/test_stage0_quality_selection.R` PASS
- `Rscript test/test_stage0_gate.R` PASS
- `Rscript test/test_trio_selection_helpers.R` PASS
- `python3 -m unittest test.test_fetch_metadata_taxon` PASS (7 tests)
- `python3 -m unittest test.test_trio_filter` PASS (25 tests)
- `bash test/test_sbst_from_dwl_cache.sh` PASS
- `bash -n src/sbst_fromDwl.sh` PASS
- `git diff --check` PASS
- T005 and T006 handoff lint PASS
- Baseline/current randomized differential comparison PASS (1,800 cases)

## Next action

Commit this bounded refactor. Push only when explicitly requested.
