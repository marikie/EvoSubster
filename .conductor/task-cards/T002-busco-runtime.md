# Task card: T002 - Audit BUSCO runtime and reproducible comparison conditions

## Goal
Audit the local environment for BUSCO and related download/runtime dependencies, identify the appropriate currently available lineage dataset for Eupercaria fish assemblies, and estimate the practical execution scope without modifying source code.

## Authority / inputs
- Context: `.conductor/context-packet.md`
- Read: local executables/environments/caches, Stage 0 external-QC schema/tests, official BUSCO/NCBI metadata when needed.

## CLAIM BOUNDARY (hard)
- Allowed: perform read-only environment checks and small metadata/network probes; report exact versions, cache paths, commands, and realistic resource constraints.
- MUST NOT: install packages system-wide, launch a large all-species BUSCO batch, edit application source, claim an override result, commit, or alter user files.
- A worker proposes the execution recipe; the conductor decides and runs the experiment.

## Steps
1. Read the context packet and the external-QC schema in code/tests.
2. Check for BUSCO, conda/mamba, Docker, Augustus/MetaEuk/Miniprot, and existing lineage datasets.
3. Determine a valid same-lineage genome-mode recipe for Eupercaria and note dataset-version compatibility.
4. Estimate runtime/storage for the candidate set suggested by available metadata.
5. Write a compact handoff with machine-check lines only.

## Machine checks
- `git status --short --branch`
- Confirm every reported executable/version with a command.

## Stop conditions
- Do not install or run a large batch; report the least invasive viable route.
- Escalate if no reproducible local runtime exists.

## Handoff
Write `.conductor/handoffs/T002-handoff.md` (400-800 words, no raw logs).
