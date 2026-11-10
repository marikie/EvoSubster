# Pipeline manifest v2

New analyses use `metadata/metadata_manifest.json` with `schema_version: 2`.
Existing schema-2 results remain usable in place, including the earlier hashed
directory layout. Older results without schema-2 manifests must be rerun to use
pipeline reports or caches. There is no migration or legacy filename fallback.

`organisms` records exactly `org1` (outgroup), `org2`, and `org3` (ingroups).
Each input has a content SHA256, an optional versioned assembly accession, and
an explicit GFF hash or null. Display names do not determine input identity.
Metadata and taxonomy copies have relative paths and hashes.

`trio_key` hashes the ordered input identities; `run_key` additionally hashes
effective LAST settings, threads, configuration, and computation/tool versions.
New runs live under
`<output>/<name1>_<ACC1>__<name2>_<ACC2>__<name3>_<ACC3>/<label>-<run_key[:12]>/`.
Slots are ordered outgroup/org1, ingroup/org2, ingroup/org3. Names come from existing
`raw_organism_name`, then `ncbi_full_name`, then the local FASTA parent directory;
no metadata lookup is added. Invalid ASCII filename character blocks become `_`,
repeated underscores collapse, and leading/trailing `_`, `.`, and `-` are removed.
Empty names become `orgN`. Versioned accessions remain complete; missing accessions
use `local-<fasta_sha256[:12]>`. The trio component is bounded to the filesystem's
`NAME_MAX` or 255 bytes, whichever is smaller. After reserving complete tokens and
separators, names share equal length budgets and only overlong names are shortened.
Tokens that cannot fit cause an error. The run label remains intact; overflow
causes an error. A valid occupied manifest with another full key extends the
suffix by four characters, up to all 64; exhausted collisions cause an error.

Before allocating, a bounded scan checks only
`<output>/*/*/metadata/metadata_manifest.json`, including incomplete and old-layout
runs. The same full `run_key` and exact stored label/date reuse that physical run,
even after display or input-path renames. Symlink aliases are deduplicated; several
distinct physical matches cause an error listing candidates. Only physical runs
at the output root's trio/run depth are eligible writer resumes; aliases targeting
external or deeper directories are ignored. New destinations are checked before
creating directories or locks. Unknown or broken
occupied candidates are preserved and rejected. Empty or `.run.lock`-only
directories may be initialized. An output-root allocation lock covers selection
through manifest preparation; the per-run lock covers the entire analysis.
Preparation receives the selected explicit path and checks identity and label
before recovery cleanup. A readable trio parent may group different full
`trio_key` values when contents differ despite matching accessions. Discovery still
selects the latest completed run per full `trio_key`. Aggregate link naming and
legacy/mixed input handling are unchanged.

Shared training profiles
live under `<cache>/v2/<training-key>/`, keyed by ordered FASTA content, accessions,
settings, threads, and LAST version. Profiles and DBs have completion hashes.
Assembly downloads use `assemblies/<accession>/versions/` with an atomically
published `genome` pointer and a verified `assembly_bundle.json`.

`artifacts` records `kind`, `slot`, run-relative `path`, and `sha256`.
Training artifacts additionally record final identity and reference/query slots.
The main kinds are `train`, `gc`, `coding`, `sbst_ratio`, `single`, `dinuc`,
`single_ncds`, `dinuc_ncds`, and the corresponding `pdf_*` plot kinds.
Reports and PCA resolve those records rather than parsing names or directories.

Statuses progress through `pending`, `identity_complete`, `counts_complete`, and
`complete`. Training-only execution is not a complete analysis. A failed stage
cannot publish full success. Recovery removes unverified generated files;
verified earlier stages remain reusable. The run lock prevents concurrent writers.
The complete plot set is staged privately, validated, then registered with hashes.
A plot/font failure is a pipeline failure. The automatic frozen regression verifies
all eight count TSVs and all required plots before declaring success.

## Review repair coverage

| Finding | Repair | Regression coverage |
| --- | --- | --- |
| Short-name cache collisions | Content-addressed manifest and directed caches | `test_run_manifest.py`, `test_train_cache.sh` |
| Assembly/GFF mixing | Assembly bundles, explicit annotation absence | `test_dwl_organism.sh`, `test_run_manifest.py` |
| TSV self-link deletion | Full preflight, collision refusal | `test_review_repairs.py` |
| Hidden plot failures | Required fonts, private plots, completion inventory | `test_pipeline_review_failures.py`, full frozen regression |
| DBS null probabilities | Folding multiplicities, categorical GOF | `test_pipeline_review_failures.py` |
| Premature outgroup stop | Evaluate following candidates | `test_outgroup_selection.R` |
| Ambiguous subRatio columns | All three bases must be A/C/G/T | `test_subRatio.py`, pipeline review reproduction |
| Inconsistent identity | Complete final LAST block, shared parser | `test_trio_selection_helpers.R`, report/PCA manifest tests |
| Ignored training failures | Check exits and complete profiles | `test_pipeline_review_failures.py`, `test_last_train_atomic.sh` |

## Validation record (2026-11-03)

- Full Python discovery passed **236 tests** at the recorded broad-suite snapshot
  (`tmp/codex/review-repairs/python-final.log`). Later additions and fixes were
  checked separately: 10 manifest audit tests, 54 filter/manifest/report tests,
  11 final report tests, and the executable PCA main regression. These overlapping
  targeted counts are not additional coverage totals or part of the 236 snapshot.
- All nine existing standalone R suites and seven shell suites passed on the final
  rerun; the added manifest plot-label R test also passed. Tests cover renamed
  inputs, optional local taxonomy, explicit annotations, failure propagation,
  assembly/version isolation, and cache integrity.
- Frozen full runs `pipeline-zjrlsmja` and `pipeline-0klt9jhd` passed, with all four
  input hashes/FASTA ordering verified and all eight expected count TSVs matching
  byte for byte. Both used eight threads and explicit `--lastdb-options '-uRY4'`.
  Expected TSVs and frozen inputs were not changed.
- The latest full run completed in 469 seconds, within the 600-second limit.
  Its 20 PDFs passed `pdfinfo` page checks. A subsequent full resume regenerated
  plots through the new label resolver while reusing all 11 training/count files
  with unchanged modification times; completed manifest/hash validation and all
  20 PDF page checks passed again.
- Actual PCA on the latest completed run succeeded for both log-ratio and observed
  variants, filtered and unfiltered: four CSVs, two samples each, finite PC scores.
  The main regression exercises manifest discovery and zero-variance channels.
  Missing scientific taxonomy for local-only inputs remains unknown; two-sample
  plotting can emit nonfatal metric/line warnings.

Full regression logs/results are retained under `tmp/codex-checks/`; focused test
and resume logs are retained under `tmp/codex/review-repairs/`.
