# Context packet - Eupercaria external BUSCO experiment

## What this project/phase is
Evaluate real assembly candidates represented by `/Users/marikonakagawa/biohazard/TheTreeOfLife/Eupercaria_tree.nwk`, create an external BUSCO genome-mode TSV compatible with EvoSubster Stage 0, and determine whether BUSCO evidence changes any Reference-first selection.

## Binding constraints
- Work in `/Users/marikonakagawa/biohazard/EvoSubster` on `feature/stage0-reference-quality-ranking`.
- Preserve all pre-existing user changes and untracked files, especially `analysis/`, `stage0-processing-flow.html`, and `stage1-matching-algorithm.html`.
- Use exact versioned GCA/GCF accessions and record all BUSCO software, lineage-dataset, and command details needed for reproducibility.
- BUSCO comparisons are valid only for genome FASTA inputs run with the same BUSCO version, exact lineage dataset/version, and materially identical settings.
- A worker may produce evidence and draft interpretations; only the conductor may accept the final claim that a Reference is or is not overturned.
- Do not commit or push. The conductor owns git history.

## Calibration
- The final claim is release-facing for the Stage 0 design, so workers leave it as `draft` and the conductor independently verifies it.
- An honest negative result is acceptable; do not search selectively for a positive override.

## Machine checks every worker must leave green
- `git status --short --branch`
- Confirm every reported artifact path exists.

## Hard stops
- No git commit.
- Do not modify application source files.
- Escalate judgment questions rather than silently changing the comparison rule.
