# Handoff: T001 - Downloader artifact manifest

task_id: T001
Status: done
claim_boundary_ok: true

## Summary

Implemented the downloader-owned artifact manifest contract in
`src/dwl_organism.sh` and its focused regression coverage in
`test/test_dwl_organism.sh`. The optional `--artifact-manifest PATH` argument
now writes a tab-separated manifest with the required
`accession`, `artifact_type`, and `relative_path` header. Manifest publication
is deferred until the downloader has completed successfully. The output is
built in a temporary file adjacent to the requested destination and finalized
with `mv`, so a completed manifest replaces its destination atomically. The
normal downloader result remains exactly six pipe-separated stdout fields.

Ownership is based on the organism directory state before the downloader moves
the fetched summary into its final location. If that directory was absent, the
manifest contains exactly one `directory_tree` row for the sanitized organism
directory relative to `--out-dir`; descendants are intentionally represented by
that single owned tree. If the directory already existed, the implementation
records a sorted before snapshot, records a sorted successful after snapshot,
and emits only paths absent from the before snapshot. Each surviving new entry
is classified as `symlink`, `directory`, or `file`, in that order, without
following symlinks during classification. Existing names remain unowned even
when the downloader overwrites their contents. Temporary snapshot and
incomplete manifest files are removed by the exit trap.

The test first exercised the fetched retry implementation on the host Bash
3.2 and reproduced the existing empty-array failure before production changes.
The retry state now uses one optional scalar instead of expanding an empty Bash
array under `set -u`; the initial request still has no HTTP-version override,
and retryable failures still switch subsequent requests to HTTP/1.1. After that
baseline became green, manifest assertions were added and observed failing on
the previously unknown option. The fresh-directory case asserts the exact
two-line manifest and verifies that stdout still has six fields. The
pre-existing-directory case seeds both a marker and an accession summary,
then verifies that only the newly downloaded FASTA, retained archive, and
taxonomy JSON are listed. It also verifies that the marker and overwritten
summary are absent from the manifest and that the marker remains on disk.

No selector, wrapper, README, or unrelated test files were edited by T001. No
commit was created. Other workers' concurrent changes visible in the shared
worktree were left untouched.

delegation_effect: difficulty=medium; agents=1; useful=isolated downloader implementation and Bash 3.2 verification; misses=none observed within the assigned claim boundary; next_mix=one integration reviewer should verify selector consumption and destructive cleanup ordering.

## Changed files

- `src/dwl_organism.sh`
- `test/test_dwl_organism.sh`
- `.conductor/handoffs/T001-handoff.md`

## Machine checks

- Bash 3.2 empty-array regression RED observation: PASS
- Manifest option RED observation: PASS
- `bash test/test_dwl_organism.sh`: PASS
- `bash -n src/dwl_organism.sh src/sbst_fromDwl.sh`: PASS
- `git diff --check -- src/dwl_organism.sh test/test_dwl_organism.sh`: PASS
- `git diff --check`: PASS
- `check-handoff.py .conductor/handoffs/T001-handoff.md`: PASS

## Unresolved risk / escalation

- The focused downloader tests cover fresh-tree ownership and new files in an
  existing directory; downstream deletion safety and cleanup ordering belong to
  T002/T003.
- No judgment escalation is required.

## Next action

The conductor should review the scoped diff, integrate it with T002's manifest
consumer, and run the combined selector/wrapper tests before the release gate.
