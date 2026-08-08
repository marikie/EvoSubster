# Unused Trio Species Cleanup Design

## Goal

After successful tree-mode trio selection, automatically remove data created by
the current selection run for accessions that are not present in any selected
trio. Preserve all data that existed before the run, all data required by a
selected trio, and the selection audit tables.

## Ownership model

`dwl_organism.sh` gains an optional `--artifact-manifest PATH` interface. Its
existing six-field stdout contract remains unchanged. On a successful call, the
manifest contains a tab-separated header and only paths created by that call:

```text
accession	artifact_type	relative_path
```

Paths are relative to the requested `--out-dir`. If the organism directory did
not exist before the call, the manifest contains one `directory_tree` row for
that directory. If it already existed, the manifest contains only newly created
`file`, `directory`, or `symlink` entries. Existing or overwritten paths are not
owned by the run and are never cleanup targets.

`trio_selection.R` aggregates these manifests by accession. It also records a
LAST train file only when the current run had to create or replace it and the
result passes the existing complete-train validator. A valid cache hit is not
owned by the run.

## Cleanup decision

After `selected_trios.tsv` is written, the selector builds the selected
accession set from `out_acc`, `in1_acc`, and `in2_acc`.

- A newly created genome artifact is retained when its accession is selected;
  otherwise it is deleted.
- A newly generated train cache is retained when both accessions are selected;
  otherwise it is deleted.
- With zero selected trios, every artifact owned by the run is deleted.
- `selected_assemblies.tsv`, `selected_trios.tsv`, train-selection metadata, and
  all pre-existing genome/cache data are retained.

The selector writes `cleanup_unused_species.tsv` with accession(s), artifact
type, path, action, and reason. `--keep-unused-species-data` disables tracking
and cleanup for debugging; cleanup remains the default in both direct selector
and wrapper tree mode.

## Safety and failure behavior

Manifest paths must be non-empty relative paths without absolute prefixes,
empty components, or `..`. The cleanup implementation constructs targets only
under the normalized genome or cache root, rejects the root itself, and checks
the normalized parent before deleting. Symlinks are unlinked without following
their targets. Files and symlinks are removed before directories; newly owned
directory trees may be removed recursively only after the same root checks.

Any unsafe path, deletion failure, or incomplete cleanup is fatal. The selector
returns non-zero and `sbst_fromDwl.sh` does not start downstream trio analysis.
Selection failure before `selected_trios.tsv` is finalized is outside this
cleanup contract. Concurrent processes must not mutate the same genome
directory during selection.

