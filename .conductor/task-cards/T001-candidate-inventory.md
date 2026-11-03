# Task card: T001 - Inventory the real Eupercaria Stage 0 candidates

## Goal
Determine how the supplied Newick tree maps to the current Stage 0 pipeline and identify the existing or reproducible candidate/selection artifacts, especially species with multiple viable assemblies and a selected Reference.

## Authority / inputs
- Context: `.conductor/context-packet.md`
- Read: `src/select/*`, `test/*stage0*`, `README.md`, `analysis/`, the supplied Newick file, and existing generated outputs under the repo.

## CLAIM BOUNDARY (hard)
- Allowed: report exact commands, paths, species/accessions, and current ranking behavior supported by inspected artifacts.
- MUST NOT: run large downloads, edit application source, claim that BUSCO changes a selection, commit, or alter user files.
- A worker proposes findings; the conductor promotes them.

## Steps
1. Read the context packet and current Stage 0 code/tests.
2. Locate prior Eupercaria Stage 0 outputs or derive the exact command to regenerate them.
3. Count tree tips/species, identify assembly candidate tables and baseline winners, and recommend a bounded but scientifically fair BUSCO comparison set.
4. Write a compact handoff with artifact paths and pass/fail checks.

## Machine checks
- `git status --short --branch`
- `test -f /Users/marikonakagawa/biohazard/TheTreeOfLife/Eupercaria_tree.nwk`

## Stop conditions
- If required outputs cannot be located, report the missing artifact and exact regeneration command.
- Escalate any scientific judgment about sampling scope.

## Handoff
Write `.conductor/handoffs/T001-handoff.md` (400-800 words, no raw logs).
