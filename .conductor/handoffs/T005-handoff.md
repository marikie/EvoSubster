# Handoff: T005 - Isolate eukaryote QC state transitions from Stage 0 orchestration

task_id: T005
status: done
claim_boundary_ok: true

## Summary

I inspected the baseline at `96787c57a0e18baf69f4e48f59be2e7316afffbd`, focusing on `src/select/assembly_selection.R` and the four Stage 0 regression suites named in the context packet. The highest-value maintainability problem is `rank_assembly_candidates()` at lines 485-719. It currently combines five responsibilities: audit-state initialization (497-511), per-species shortlisting (513-529), baseline ranking/selection (531-540), eukaryote QC review and optional override (542-703), and result projection (709-718). The eukaryote block is especially risky because the same full candidate table is mutated through `qc_preferred`, review-reason, blocker, selection, basis, and override fields while global row indexes and paired-accession equivalence must remain aligned.

The smallest useful refactor is extraction, not a ranking rewrite. First extract `stage0_initialize_ranking_state(candidates)` from lines 497-511; it must append the same columns in the same order and construct `assembly_equivalence_key` last. Then extract `stage0_apply_eukaryote_qc(candidates, shortlist_index, baseline_index, baseline_order_index, allow_qc_override)` from lines 542-703. It should return `list(candidates = candidates, selected_index = selected_index)` and own every QC-specific state transition, while the existing species loop continues to own shortlisting, baseline ranks, and final-rank assignment. Within that helper, extract one pure `stage0_override_blockers(candidates, baseline_index, alternative_index, complete_index)` from lines 645-677. It must return blocker codes in the current order: baseline completeness, alternative completeness, run comparability, Single-copy, Duplicated, Fragmented, Missing. Existing `stage0_add_review_reason()` and `stage0_qc_order()` should remain the single writers/order authority.

Do not combine this with cleanup of `stage0_phase1_order()` lines 433-453 or `stage0_final_order()` lines 455-483, even though some profile branches appear unreachable through the current coordinator. Removing them is a separate behavior/API decision and adds no value to the targeted extraction. Likewise, leave validation/error text in lines 68-401 untouched.

Regression risks are exact audit-column order; accidentally assigning `qc_preferred` to a paired GCA/GCF alias rather than the baseline representative; reviewing equivalent paired accessions; changing semicolon-separated reason/blocker order; allowing one-sided or incomplete BUSCO evidence to override; changing `baseline_rank` when an override receives `final_rank == 1`; and disturbing metadata tie-breaking in `stage0_qc_order()`. Existing tests pin most of these: `test/test_stage0_quality_selection.R:130-208`, `256-410`, and `480-739`; gate behavior is pinned by `test/test_stage0_gate.R:47-68`; integration and exact review schema are pinned by `test/test_trio_selection_helpers.R:114-278`. Before or with implementation, add focused characterization for multiple simultaneously valid override candidates (`not_top_override_candidate`) and a multi-code blocker string, because those transitions are not directly asserted today.

delegation_effect: difficulty=medium; agents=1; useful=isolated the mutation-heavy QC boundary and two missing characterization cases; misses=did not inspect a completed implementation diff; next_mix=one implementer plus one independent reviewer.

## Changed files

- `.conductor/handoffs/T005-handoff.md`

## Machine checks

- `Rscript test/test_stage0_quality_selection.R` PASS
- `Rscript test/test_stage0_gate.R` PASS
- `Rscript test/test_trio_selection_helpers.R` PASS
- `bash test/test_sbst_from_dwl_cache.sh` PASS

## Unresolved risk / escalation

- The recommendation is based on the unchanged baseline; it does not establish that a future implementation preserves behavior.

## Next action

Implement the three bounded extractions, add the two characterization cases, then independently compare the final diff and rerun all four checks plus `git diff --check` before making any behavior-preservation claim.
