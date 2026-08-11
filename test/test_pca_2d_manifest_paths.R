#!/usr/bin/env Rscript

# Regression tests for metadata manifests copied with stale absolute paths.

this_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
root <- dirname(dirname(normalizePath(this_file)))
source(file.path(root, "src", "visualize", "pca_2d.R"))

fail <- 0L
check <- function(cond, msg) {
  if (isTRUE(cond)) {
    cat("ok:", msg, "\n")
  } else {
    cat("FAIL:", msg, "\n")
    fail <<- fail + 1L
  }
}

make_entry <- function(slot, role, accession, short_name, metadata_json) {
  list(
    slot = slot,
    role = role,
    accession = accession,
    short_name = short_name,
    metadata_json = metadata_json
  )
}

write_manifest <- function(manifest_dir, metadata_dir, org2_json, org3_json) {
  dir.create(manifest_dir, recursive = TRUE, showWarnings = FALSE)
  manifest_path <- file.path(manifest_dir, "metadata_manifest.json")
  manifest <- list(
    metadata_dir = metadata_dir,
    organisms = list(
      make_entry("org1", "outgroup", "GCA_000000001.1", "Out1", file.path(metadata_dir, "Out1_GCA_000000001.1.json")),
      make_entry("org2", "ingroup", "GCA_000000002.1", "In2", org2_json),
      make_entry("org3", "ingroup", "GCA_000000003.1", "In3", org3_json)
    )
  )
  jsonlite::write_json(manifest, manifest_path, auto_unbox = TRUE, pretty = TRUE)
  manifest_path
}

write_fixture_file <- function(path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  writeLines("{}", path)
}

write_identity <- function(path, value) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  writeLines(sprintf("# substitution percent identity: %s%%", value), path)
}

run_tests <- function() {
  test_root <- tempfile("pca-manifest-paths-")
  dir.create(test_root)
  on.exit(unlink(test_root, recursive = TRUE), add = TRUE)

  relocated_dir <- file.path(test_root, "relocated", "metadata")
  old_dir <- file.path(test_root, "old", "metadata")
  local_org2 <- file.path(relocated_dir, "In2_GCA_000000002.1.json")
  local_org3 <- file.path(relocated_dir, "In3_GCA_000000003.1.json")
  old_org2 <- file.path(old_dir, basename(local_org2))
  old_org3 <- file.path(old_dir, basename(local_org3))
  for (path in c(local_org2, local_org3, old_org2, old_org3)) write_fixture_file(path)
  relocated_manifest <- write_manifest(relocated_dir, old_dir, old_org2, old_org3)
  relocated <- load_manifest_lookup(relocated_manifest)
  check(
    identical(relocated$lookup[["GCA_000000002.1::In2"]]$metadata_path, normalizePath(local_org2)),
    "relocated metadata beside the manifest takes precedence over the recorded path"
  )
  check(
    identical(relocated$lookup[["GCA_000000003.1::In3"]]$taxonomy_path,
              file.path(normalizePath(relocated_dir), "taxonomy_In3_GCA_000000003.1.json")),
    "taxonomy path follows the relocated metadata directory"
  )
  check(
    length(relocated$lookup) == 2 && length(relocated$slot_map) == 3 &&
      identical(relocated$slot_map$org1$role, "outgroup"),
    "outgroup remains in the slot map but is excluded from the ingroup lookup"
  )

  fallback_manifest_dir <- file.path(test_root, "fallback-manifest")
  fallback_dir <- file.path(test_root, "fallback-metadata")
  fallback_org2 <- file.path(fallback_dir, "In2_GCA_000000002.1.json")
  fallback_org3 <- file.path(fallback_dir, "In3_GCA_000000003.1.json")
  write_fixture_file(fallback_org2)
  write_fixture_file(fallback_org3)
  fallback_manifest <- write_manifest(fallback_manifest_dir, fallback_dir, fallback_org2, fallback_org3)
  fallback <- load_manifest_lookup(fallback_manifest)
  check(
    identical(fallback$lookup[["GCA_000000002.1::In2"]]$metadata_path, normalizePath(fallback_org2)),
    "valid recorded metadata path remains a fallback"
  )
  check(
    identical(fallback$lookup[["GCA_000000002.1::In2"]]$taxonomy_path,
              file.path(normalizePath(fallback_dir), "taxonomy_In2_GCA_000000002.1.json")),
    "taxonomy path follows fallback metadata when no relocated copy exists"
  )

  relative_dir <- file.path(test_root, "relative", "metadata")
  relative_org2 <- file.path(relative_dir, "nested", "In2_GCA_000000002.1.json")
  relative_org3 <- file.path(relative_dir, "nested", "In3_GCA_000000003.1.json")
  write_fixture_file(relative_org2)
  write_fixture_file(relative_org3)
  relative_manifest <- write_manifest(
    relative_dir,
    ".",
    file.path("nested", basename(relative_org2)),
    file.path("nested", basename(relative_org3))
  )
  relative <- load_manifest_lookup(relative_manifest)
  check(
    identical(relative$lookup[["GCA_000000002.1::In2"]]$metadata_path, normalizePath(relative_org2)),
    "relative metadata path remains resolved from the manifest directory"
  )

  missing_manifest_dir <- file.path(test_root, "missing", "metadata")
  missing_dir <- file.path(test_root, "missing-old", "metadata")
  missing_manifest <- write_manifest(
    missing_manifest_dir,
    missing_dir,
    file.path(missing_dir, "In2_GCA_000000002.1.json"),
    file.path(missing_dir, "In3_GCA_000000003.1.json")
  )
  error_message <- tryCatch(
    {
      load_manifest_lookup(missing_manifest)
      NULL
    },
    error = function(e) e$message
  )
  check(
    !is.null(error_message) && grepl("checked relocated path", error_message, fixed = TRUE) &&
      grepl("recorded path", error_message, fixed = TRUE),
    "missing metadata reports both relocated and recorded candidates"
  )

  empty_manifest <- write_manifest(
    file.path(test_root, "empty-metadata-dir"),
    "",
    "In2_GCA_000000002.1.json",
    "In3_GCA_000000003.1.json"
  )
  empty_error <- tryCatch(
    {
      load_manifest_lookup(empty_manifest)
      NULL
    },
    error = function(e) e$message
  )
  check(
    !is.null(empty_error) && grepl("'metadata_dir' missing or empty", empty_error, fixed = TRUE),
    "empty metadata_dir fails with a clear manifest error"
  )

  artifact_run <- file.path(test_root, "artifacts", "20260810")
  train_dir <- file.path(artifact_run, "intermediateFiles")
  write_identity(file.path(artifact_run, "Out12In2_20260801.train"), 71)
  write_identity(file.path(train_dir, "Out12In2_20260809.train"), 89)
  write_identity(file.path(train_dir, "Out12In2_20260810.train"), 91)
  write_identity(file.path(train_dir, "Out12In3_20260810.train"), 92)
  write_identity(file.path(train_dir, "In22In3_20260810.train"), 97)
  identity_result <- collect_identity_strings(artifact_run, relocated$slot_map)
  check(
    length(identity_result$issues) == 0 &&
      identical(identity_result$identity, list(idt_12 = "89%", idt_13 = "92%", idt_23 = "97%")),
    "identity discovery prefers intermediateFiles and selects its first sorted match"
  )

  legacy_train_run <- file.path(test_root, "legacy-train", "20260810")
  write_identity(file.path(legacy_train_run, "Out12In2_20260810.train"), 81)
  write_identity(file.path(legacy_train_run, "Out12In3_20260810.train"), 82)
  write_identity(file.path(legacy_train_run, "In22In3_20260810.train"), 87)
  legacy_identity <- collect_identity_strings(legacy_train_run, relocated$slot_map)
  check(
    length(legacy_identity$issues) == 0 &&
      identical(legacy_identity$identity, list(idt_12 = "81%", idt_13 = "82%", idt_23 = "87%")),
    "identity discovery falls back to the legacy run-root layout"
  )

  ingroup2_tsv <- file.path(artifact_run, "statistics", "In2", "singlenuc", "GCA_000000002.1_In2_20260810_ncds.tsv")
  ingroup3_tsv <- file.path(artifact_run, "statistics", "In3", "singlenuc", "GCA_000000003.1_In3_20260810_ncds.tsv")
  dinuc_tsv <- file.path(artifact_run, "statistics", "In2", "dinuc", "GCA_000000002.1_In2_20260810_dinuc_ncds.tsv")
  misplaced_dinuc_tsv <- file.path(artifact_run, "statistics", "In2", "singlenuc", "GCA_000000002.1_In2_20260810_dinuc_ncds.tsv")
  misplaced_outgroup_tsv <- file.path(artifact_run, "statistics", "In2", "singlenuc", "GCA_000000001.1_Out1_20260810_ncds.tsv")
  legacy_tsv <- file.path(artifact_run, "GCA_000000002.1_In2_20260810_ncds.tsv")
  for (path in c(ingroup2_tsv, ingroup3_tsv, dinuc_tsv, misplaced_dinuc_tsv, misplaced_outgroup_tsv, legacy_tsv)) write_fixture_file(path)
  nested_tsvs <- collect_substitution_tsvs(artifact_run, "*_ncds.tsv", relocated$slot_map)
  check(
    identical(nested_tsvs, sort(c(ingroup2_tsv, ingroup3_tsv))),
    "TSV discovery uses ingroup singlenuc directories and excludes dinuc and legacy duplicates"
  )

  legacy_run <- file.path(test_root, "legacy-artifacts", "20260810")
  legacy_only_tsv <- file.path(legacy_run, "GCA_000000002.1_In2_20260810_ncds.tsv")
  legacy_dinuc_tsv <- file.path(legacy_run, "GCA_000000002.1_In2_20260810_dinuc_ncds.tsv")
  legacy_outgroup_tsv <- file.path(legacy_run, "GCA_000000001.1_Out1_20260810_ncds.tsv")
  for (path in c(legacy_only_tsv, legacy_dinuc_tsv, legacy_outgroup_tsv)) write_fixture_file(path)
  check(
    identical(collect_substitution_tsvs(legacy_run, "*_ncds.tsv", relocated$slot_map), legacy_only_tsv),
    "legacy TSV fallback keeps only ingroup singlenuc files"
  )

  mixed_run <- file.path(test_root, "mixed-artifacts", "20260810")
  mixed_nested_tsv <- file.path(mixed_run, "statistics", "In2", "singlenuc", "GCA_000000002.1_In2_20260810_ncds.tsv")
  mixed_legacy_tsv <- file.path(mixed_run, "GCA_000000003.1_In3_20260810_ncds.tsv")
  write_fixture_file(mixed_nested_tsv)
  write_fixture_file(mixed_legacy_tsv)
  check(
    identical(
      collect_substitution_tsvs(mixed_run, "*_ncds.tsv", relocated$slot_map),
      sort(c(mixed_nested_tsv, mixed_legacy_tsv))
    ),
    "TSV discovery falls back independently for each ingroup in a mixed layout"
  )

  if (fail > 0L) quit(status = 1)
  cat("All PCA manifest path tests passed.\n")
}

run_tests()