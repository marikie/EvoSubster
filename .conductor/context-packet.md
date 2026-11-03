# Context packet - unused trio species cleanup

## What this project/phase is

Implement automatic cleanup on the latest `feature/newick-input` branch in the
isolated worktree. Trio selection downloads candidate genomes and generates
LAST cache files; only artifacts created by this run for accessions absent from
all selected trios may be deleted.

## Binding constraints

- Do not touch the primary worktree or its unrelated untracked files.
- Preserve all pre-existing genome and cache paths.
- Preserve selected accession data and all selection audit TSVs.
- Downloader stdout remains six pipe-separated fields.
- Manifest header is `accession\tartifact_type\trelative_path`.
- Manifest paths are relative to `--out-dir`; allowed types are
  `directory_tree`, `directory`, `file`, and `symlink`.
- Any unsafe path or deletion failure is fatal and blocks downstream analysis.
- Code and comments are English. Workers do not commit or revert unrelated work.

## Machine checks every worker must leave green

- `bash -n src/dwl_organism.sh src/sbst_fromDwl.sh`
- `git diff --check`
- The focused tests named in the assigned task card.

## Hard stops

- No git commit; the conductor commits.
- Do not change the manifest contract without escalation.
- On a machine-check failure, stop and report it.
- Escalate judgment questions with `ESCALATE: JUDGMENT_REQUIRED`.
