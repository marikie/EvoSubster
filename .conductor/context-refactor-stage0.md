# Context packet: Stage 0 refactor

- Repository: `/Users/marikonakagawa/biohazard/EvoSubster`
- Branch: `feature/stage0-reference-quality-ranking`
- Baseline commit: `96787c57a0e18baf69f4e48f59be2e7316afffbd`
- Goal: improve the internal structure and maintainability of Stage 0 assembly selection without changing its documented behavior, CLI, output schemas, or ranking results.
- Primary implementation: `src/select/assembly_selection.R`
- Integration boundary: `src/select/trio_selection.R`, `src/sbst_fromDwl.sh`
- Primary tests: `test/test_stage0_quality_selection.R`, `test/test_stage0_gate.R`, `test/test_trio_selection_helpers.R`, `test/test_sbst_from_dwl_cache.sh`
- Preserve unrelated untracked files: `analysis/`, `stage0-processing-flow.html`, `stage1-matching-algorithm.html`.
- Code and comments must be English. User-facing explanation is Japanese.
- Do not push or alter remote state.
- Claim boundary: this phase may claim behavior preservation only after the conductor reruns the complete relevant test set and reviews the final diff.
