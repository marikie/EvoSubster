# Handoff: T006 - Adversarial review of the Stage 0 refactor

task_id: T006
status: done
claim_boundary_ok: true

## Verdict

APPROVE. I found no Critical, Important, or Minor correctness findings in the refactor diff relative to `96787c57a0e18baf69f4e48f59be2e7316afffbd`. This is an independent review verdict, not the conductor's final behavior-preservation claim.

## Summary

The production diff is confined to `src/select/assembly_selection.R` and extracts the mutation-heavy portion of `rank_assembly_candidates()` into cohesive helpers. The reference-first five-tier ranking is unchanged because `stage0_baseline_order()` and all of its ordering fields remain byte-for-byte untouched. The species loop still computes the shortlist and metadata baseline before invoking the eukaryote-only QC helper. Prokaryote candidates continue to bypass that helper and retain their ANI, CheckM, Merqury, assembly-quality ordering.

Default QC remains audit-only: `selected` is initialized to the baseline row, and the extracted override path is entered only when `isTRUE(allow_qc_override)` and a distinct assembly-equivalence class exists. The blocker predicates and their stable order are preserved: baseline completeness, alternative completeness, run comparability, Single-copy, Duplicated, Fragmented, and Missing. Paired GCA/GCF rows are still collapsed through the unchanged equivalence-key construction, and QC preference is assigned to the baseline-ranked representative.

Audit column initialization remains in the original order, and the final `audit`, `shortlist`, `best`, and `review` projections and ordering code are unchanged. Validation before ranking, including exact versioned-accession matching, duplicate rejection, BUSCO genome-mode requirements, lineage/version comparability, component validation, and Single-copy derivation, is outside the refactored region and unchanged. The two added characterization cases correctly pin multiple valid override selection plus `not_top_override_candidate`, and the stable ordering of multiple blocker codes.

I also loaded the baseline and current modules into separate R environments and compared the complete returned lists for 300 generated metadata/QC cases under `shortlist_k` values 1, 2, and `Inf`, with override both disabled and enabled. All 1,800 comparisons were `identical()`. Those cases varied eukaryote/prokaryote profiles, one-to-five candidates, reference/annotation status, hard exclusions, assembly metadata, complete BUSCO evidence, and paired accessions.

## Findings

- Critical: none.
- Important: none.
- Minor: none.

## Changed files

- `.conductor/handoffs/T006-handoff.md`

## Machine checks

- `Rscript test/test_stage0_quality_selection.R` PASS
- `Rscript test/test_stage0_gate.R` PASS
- `Rscript test/test_trio_selection_helpers.R` PASS
- `bash test/test_sbst_from_dwl_cache.sh` PASS
- `git diff --check 96787c57a0e18baf69f4e48f59be2e7316afffbd` PASS
- Baseline/current randomized differential comparison, 1,800 cases PASS

## Claim boundary

The evidence supports this review verdict, but final verification remains owned by the conductor and requires its fresh test rerun against the final working-tree state.

delegation_effect: difficulty=medium; agents=1; useful=confirmed exact output equivalence across targeted tests and 1,800 differential cases; misses=no independent production-data replay; next_mix=main conductor only for final diff inspection and machine-check rerun.

## Next action

The conductor should inspect the final diff and working-tree scope, rerun the complete relevant checks on the final state, and only then promote the behavior-preservation claim.
