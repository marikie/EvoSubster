#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 7) {
    stop("Expected: template, input JSON, output, format, preview dir, images dir, pdftools flag")
}

rmd_path <- args[[1]]
input_json <- args[[2]]
output_file <- args[[3]]
output_format <- args[[4]]
preview_dir <- args[[5]]
images_dir <- args[[6]]
pdftools_available <- identical(args[[7]], "TRUE")

adjust_docx_image_widths <- function(docx_path, target_cm = 16.5) {
    if (!requireNamespace("xml2", quietly = TRUE)) {
        message("Skipping width normalization; 'xml2' package not available.")
        return(invisible(FALSE))
    }

    target_emu <- as.integer(round(target_cm / 2.54 * 914400))
    work_dir <- tempfile("docx_fix_")
    dir.create(work_dir, recursive = TRUE, showWarnings = FALSE)
    on.exit(unlink(work_dir, recursive = TRUE, force = TRUE), add = TRUE)

    utils::unzip(docx_path, exdir = work_dir)
    doc_xml_path <- file.path(work_dir, "word", "document.xml")
    if (!file.exists(doc_xml_path)) {
        return(invisible(FALSE))
    }

    document <- xml2::read_xml(doc_xml_path)
    namespaces <- xml2::xml_ns(document)
    extent_nodes <- xml2::xml_find_all(document, ".//wp:extent", namespaces)
    if (length(extent_nodes) == 0) {
        return(invisible(FALSE))
    }

    changed <- FALSE
    for (node in extent_nodes) {
        parent_inline <- xml2::xml_parent(node)
        picture <- xml2::xml_find_first(
            parent_inline,
            ".//pic:pic/pic:nvPicPr/pic:cNvPr",
            namespaces
        )
        description <- if (inherits(picture, "xml_missing")) {
            NA_character_
        } else {
            xml2::xml_attr(picture, "descr")
        }
        if (is.na(description) || !startsWith(description, "images/")) {
            next
        }

        old_width <- suppressWarnings(as.numeric(xml2::xml_attr(node, "cx")))
        old_height <- suppressWarnings(as.numeric(xml2::xml_attr(node, "cy")))
        if (!is.finite(old_width) || old_width <= 0) {
            next
        }

        scale <- target_emu / old_width
        xml2::xml_set_attr(node, "cx", as.character(target_emu))
        if (is.finite(old_height) && old_height > 0) {
            xml2::xml_set_attr(node, "cy", as.character(round(old_height * scale)))
        }
        changed <- TRUE
    }

    if (!changed) {
        return(invisible(FALSE))
    }

    zip_binary <- Sys.which("zip")
    if (!nzchar(zip_binary)) {
        stop("Cannot normalize DOCX image widths because 'zip' is not available.")
    }

    output_path <- normalizePath(docx_path, mustWork = TRUE)
    staged_docx <- tempfile(
        pattern = paste0(".", basename(output_path), "."),
        tmpdir = dirname(output_path),
        fileext = ".docx"
    )
    on.exit(unlink(staged_docx, force = TRUE), add = TRUE)

    xml2::write_xml(document, doc_xml_path)
    old_working_dir <- getwd()
    on.exit(setwd(old_working_dir), add = TRUE)
    setwd(work_dir)
    files_to_zip <- list.files(".", recursive = TRUE, all.files = TRUE)
    zip_output <- system2(
        zip_binary,
        c("-r9Xq", shQuote(staged_docx), shQuote(files_to_zip)),
        stdout = TRUE,
        stderr = TRUE
    )
    zip_status <- attr(zip_output, "status")
    staged_size <- if (file.exists(staged_docx)) file.info(staged_docx)$size else 0
    if ((!is.null(zip_status) && zip_status != 0) || staged_size <= 0) {
        stop("Failed to repack DOCX: ", paste(zip_output, collapse = "\n"))
    }
    if (!file.rename(staged_docx, output_path)) {
        stop("Failed to replace DOCX with normalized output: ", output_path)
    }
    invisible(TRUE)
}

for (path in c(rmd_path, input_json)) {
    if (!file.exists(path)) {
        stop(sprintf("Input file not found: %s", path))
    }
}
for (path in c(preview_dir, images_dir)) {
    dir.create(path, recursive = TRUE, showWarnings = FALSE)
}

if (!requireNamespace("rmarkdown", quietly = TRUE)) {
    stop("Missing required R package: rmarkdown")
}

rmarkdown::render(
    input = rmd_path,
    params = list(
        input_json = input_json,
        preview_dir = preview_dir,
        images_dir = images_dir,
        pdftools_available = pdftools_available,
        page_break_between_datasets = TRUE
    ),
    output_file = output_file,
    output_format = output_format
)

if (identical(output_format, "word_document")) {
    adjust_docx_image_widths(output_file)
}
