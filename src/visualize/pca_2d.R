#!/usr/bin/env Rscript
.pca_source_file <- local({
  sourced <- Filter(Negate(is.null), lapply(sys.frames(), function(x) x$ofile))
  if (length(sourced)) tail(sourced, 1)[[1]] else sub("^--file=", "", grep("^--file=", commandArgs(), value=TRUE)[1])
})
.pca_manifest_helper <- normalizePath(file.path(dirname(.pca_source_file), "..", "lib", "pca_input_adapter.py"))

# Keep this entry point sourceable; all modules share the caller's environment.
.pca_module_dir <- file.path(dirname(normalizePath(.pca_source_file)), "pca")
source(file.path(dirname(.pca_module_dir), "transparent_background.R"), local = TRUE)
for (.pca_module in c("config.R", "input.R", "taxonomy.R", "identity.R",
                      "compute.R", "plot.R", "output.R", "pipeline.R")) {
  source(file.path(.pca_module_dir, .pca_module), local = TRUE)
}
rm(.pca_module)


if (sys.nframe() == 0L) {
  tryCatch(
    main(),
    error = function(e) {
      message(sprintf("Error: %s", e$message))
      quit(status = 1)
    }
  )
}
