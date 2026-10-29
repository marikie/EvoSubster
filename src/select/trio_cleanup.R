# Helpers loaded by trio_selection.R into its calling environment.

empty_owned_artifacts <- function() {
  data.frame(
    accession_1 = character(),
    accession_2 = character(),
    storage = character(),
    artifact_type = character(),
    relative_path = character(),
    stringsAsFactors = FALSE
  )
}

read_artifact_manifest <- function(path, requested_accession) {
  if (!file.exists(path)) {
    stop("dwl_organism.sh did not create its artifact manifest for ", requested_accession,
         call. = FALSE)
  }
  manifest <- read_tsv(path)
  expected <- c("accession", "artifact_type", "relative_path")
  if (!identical(names(manifest), expected)) {
    stop("Invalid downloader artifact manifest header: ", path, call. = FALSE)
  }
  if (!nrow(manifest)) return(empty_owned_artifacts())

  allowed_types <- c("directory_tree", "directory", "file", "symlink")
  invalid <- is.na(manifest$accession) | !nzchar(manifest$accession) |
    manifest$accession != requested_accession |
    is.na(manifest$artifact_type) | !manifest$artifact_type %in% allowed_types |
    is.na(manifest$relative_path) | !nzchar(manifest$relative_path)
  if (any(invalid)) {
    stop("Invalid row in downloader artifact manifest: ", path, call. = FALSE)
  }

  data.frame(
    accession_1 = manifest$accession,
    accession_2 = "",
    storage = "genome",
    artifact_type = manifest$artifact_type,
    relative_path = manifest$relative_path,
    stringsAsFactors = FALSE
  )
}

path_exists_without_following <- function(path) {
  link_target <- Sys.readlink(path)
  file.exists(path) || dir.exists(path) || (!is.na(link_target) && nzchar(link_target))
}

prepare_cleanup_target <- function(relative_path, root) {
  if (is.na(relative_path) || !nzchar(relative_path) ||
      startsWith(relative_path, "/") ||
      grepl("^[A-Za-z]:[/\\\\]", relative_path)) {
    stop("Unsafe cleanup path: ", relative_path, call. = FALSE)
  }
  components <- strsplit(relative_path, "/", fixed = TRUE)[[1]]
  if (any(!nzchar(components)) || any(components %in% c(".", ".."))) {
    stop("Unsafe cleanup path: ", relative_path, call. = FALSE)
  }

  normalized_root <- normalizePath(root, mustWork = TRUE)
  normalized_parent <- normalizePath(
    file.path(normalized_root, dirname(relative_path)), mustWork = TRUE
  )
  root_prefix <- paste0(normalized_root, .Platform$file.sep)
  if (normalized_parent != normalized_root && !startsWith(normalized_parent, root_prefix)) {
    stop("Cleanup path escapes configured root: ", relative_path, call. = FALSE)
  }

  target <- file.path(normalized_parent, basename(relative_path))
  if (identical(target, normalized_root)) {
    stop("Refusing to delete a cleanup root", call. = FALSE)
  }
  target
}

build_cleanup_audit <- function(artifacts, selected) {
  required_columns <- c(
    "accession_1", "accession_2", "storage", "artifact_type", "relative_path"
  )
  if (!all(required_columns %in% names(artifacts))) {
    stop("Owned artifact table is missing required columns", call. = FALSE)
  }

  selected_accessions <- character()
  if (nrow(selected)) {
    selected_columns <- c("out_acc", "in1_acc", "in2_acc")
    if (!all(selected_columns %in% names(selected))) {
      stop("Selected trio table is missing accession columns", call. = FALSE)
    }
    selected_accessions <- unique(unname(unlist(selected[selected_columns], use.names = FALSE)))
    selected_accessions <- selected_accessions[!is.na(selected_accessions) &
                                                 nzchar(selected_accessions)]
  }

  audit <- data.frame(
    accessions = character(), artifact_type = character(), path = character(),
    action = character(), reason = character(), stringsAsFactors = FALSE
  )
  if (!nrow(artifacts)) {
    return(list(audit = audit, retained = logical()))
  }

  accessions <- vapply(seq_len(nrow(artifacts)), function(i) {
    values <- c(artifacts$accession_1[i], artifacts$accession_2[i])
    paste(values[!is.na(values) & nzchar(values)], collapse = ",")
  }, character(1))
  pair_artifact <- !is.na(artifacts$accession_2) & nzchar(artifacts$accession_2)
  retained <- artifacts$accession_1 %in% selected_accessions &
    (artifacts$storage == "genome" | !pair_artifact |
       artifacts$accession_2 %in% selected_accessions)
  selected_reason <- ifelse(
    artifacts$storage == "genome" | !pair_artifact,
    "selected_accession", "selected_accession_pair"
  )
  unused_reason <- ifelse(
    artifacts$storage == "genome" | !pair_artifact,
    "unused_accession", "unused_accession_pair"
  )
  audit <- data.frame(
    accessions = accessions,
    artifact_type = artifacts$artifact_type,
    path = artifacts$relative_path,
    action = ifelse(retained, "retained", "pending"),
    reason = ifelse(retained, selected_reason, unused_reason),
    stringsAsFactors = FALSE
  )
  list(audit = audit, retained = retained)
}

cleanup_owned_artifacts <- function(artifacts, selected, genome_root, train_cache_root,
                                    audit_path) {
  initial <- build_cleanup_audit(artifacts, selected)
  audit <- initial$audit
  retained <- initial$retained
  dir.create(dirname(audit_path), recursive = TRUE, showWarnings = FALSE)
  if (!nrow(artifacts)) {
    write_tsv_atomic(audit, audit_path)
    return(invisible(audit))
  }
  # Replace a previous run's report before any destructive preflight begins.
  write_tsv_atomic(audit, audit_path)

  active_index <- NA_integer_
  # These phases share audit progress so failures retain the exact active row.
  preflight_targets <- function() {
    allowed_types <- c("directory_tree", "directory", "file", "symlink")
    invalid_train_type <- artifacts$storage == "train_cache" &
      !artifacts$artifact_type %in% c("file", "directory_tree")
    if (any(!artifacts$storage %in% c("genome", "train_cache")) ||
        any(!artifacts$artifact_type %in% allowed_types) || any(invalid_train_type)) {
      stop("Owned artifact table contains an invalid storage or artifact type", call. = FALSE)
    }

    roots <- c(genome = genome_root, train_cache = train_cache_root)
    targets <- character(nrow(artifacts))
    for (i in seq_len(nrow(artifacts))) {
      active_index <<- i
      targets[i] <- prepare_cleanup_target(
        artifacts$relative_path[i], roots[[artifacts$storage[i]]]
      )
    }
    active_index <<- NA_integer_
    audit$path <<- targets
    write_tsv_atomic(audit, audit_path)

    delete_indices <- which(!retained)
    if (length(delete_indices)) {
      duplicate_targets <- targets[duplicated(targets) | duplicated(targets, fromLast = TRUE)]
      if (length(duplicate_targets)) {
        duplicate_targets <- unique(duplicate_targets)
        conflict <- any(vapply(duplicate_targets, function(target) {
          states <- retained[targets == target]
          any(states) && any(!states)
        }, logical(1)))
        if (conflict) {
          stop("Cleanup target is both retained and scheduled for deletion", call. = FALSE)
        }
        stop("Owned artifact table contains a duplicate cleanup target", call. = FALSE)
      }

      recursive_targets <- targets[
        !retained & artifacts$artifact_type == "directory_tree"
      ]
      retained_targets <- targets[retained]
      for (target in recursive_targets) {
        prefix <- paste0(target, .Platform$file.sep)
        if (any(retained_targets == target | startsWith(retained_targets, prefix))) {
          stop("Recursive cleanup target contains a retained artifact: ", target,
               call. = FALSE)
        }
      }

      # Validate every target before deleting any of them.
      for (i in delete_indices) {
        active_index <<- i
        target <- targets[i]
        link_target <- Sys.readlink(target)
        is_link <- !is.na(link_target) && nzchar(link_target)
        valid_type <- switch(
          artifacts$artifact_type[i],
          symlink = is_link,
          file = !is_link && file.exists(target) && !dir.exists(target),
          directory = !is_link && dir.exists(target),
          directory_tree = !is_link && dir.exists(target),
          FALSE
        )
        if (!valid_type) {
          stop("Cleanup target is missing or has an unexpected type: ", target,
               call. = FALSE)
        }
      }
      active_index <<- NA_integer_
    }
    delete_indices
  }
  delete_targets <- function(delete_indices, targets) {
    depth <- lengths(strsplit(artifacts$relative_path, "/", fixed = TRUE))
    type_order <- match(
      artifacts$artifact_type,
      c("file", "symlink", "directory", "directory_tree")
    )
    delete_indices <- delete_indices[order(type_order[delete_indices],
                                           -depth[delete_indices])]
    for (i in delete_indices) {
      active_index <<- i
      recursive <- artifacts$artifact_type[i] == "directory_tree"
      status <- unlink(targets[i], recursive = recursive, force = FALSE)
      if (status != 0 || path_exists_without_following(targets[i])) {
        stop("Failed to delete current-run artifact: ", targets[i], call. = FALSE)
      }
      audit$action[i] <<- "deleted"
    }
    active_index <<- NA_integer_
  }
  tryCatch({
    delete_indices <- preflight_targets()
    delete_targets(delete_indices, audit$path)

    write_tsv_atomic(audit, audit_path)
    invisible(audit)
  }, error = function(error) {
    failed_rows <- audit$action == "pending"
    if (!is.na(active_index)) failed_rows[active_index] <- TRUE
    if (!any(failed_rows)) failed_rows <- audit$action != "deleted"
    audit$action[failed_rows] <- "failed"
    audit$reason[failed_rows] <- paste0(
      audit$reason[failed_rows], "; ", conditionMessage(error)
    )
    write_tsv_atomic(audit, audit_path)
    stop(error)
  })
}
