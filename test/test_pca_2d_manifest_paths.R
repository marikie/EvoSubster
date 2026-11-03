#!/usr/bin/env Rscript
this_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value=TRUE)[1])
root <- dirname(dirname(normalizePath(this_file)))
source(file.path(root, "src", "visualize", "pca_2d.R"))
work <- tempfile("pca-v2-", tmpdir=file.path(root,"tmp","codex"))
dir.create(work)
on.exit <- NULL
code <- paste0("import sys; from pathlib import Path; sys.path.insert(0,sys.argv[1]); ",
               "from test_collect_run_summary import create_run; ",
               "print(create_run(Path(sys.argv[2])/'renamable'))")
run <- system2("python3", c("-c",shQuote(code),shQuote(file.path(root,"test")),shQuote(work)),stdout=TRUE)
stopifnot(is.null(attr(run,"status")))
run <- tail(run,1)
bundle <- load_manifest_lookup(file.path(run,"metadata","metadata_manifest.json"))
stopifnot(length(bundle$slot_map)==3, length(bundle$lookup)==2)
found <- collect_substitution_tsvs(run,"*.tsv",bundle$slot_map)
stopifnot(identical(sort(basename(found)),c("org2.tsv","org3.tsv")))
identities <- collect_identity_strings(run,bundle$slot_map)
stopifnot(identical(identities$identity,list(idt_12="90 %",idt_13="91 %",idt_23="97 %")))
renamed <- file.path(work,"unrelated-name")
stopifnot(file.rename(run,renamed))
stopifnot(identical(discover_dataset_dirs(work), normalizePath(renamed)))
relocated <- load_manifest_lookup(file.path(renamed,"metadata","metadata_manifest.json"))
stopifnot(all(file.exists(vapply(relocated$lookup,function(x)x$metadata_path,character(1)))))
stopifnot(identical(unname(build_classification_palette(c("Apicomplexa","Arthropoda"),"phylum")),c("#E69F00","#0072B2")))
payload <- jsonlite::fromJSON(file.path(renamed,"metadata","metadata_manifest.json"),simplifyVector=FALSE)
for (i in seq_along(payload$organisms)) {
  payload$organisms[[i]]$accession <- NULL
  payload$organisms[[i]]$taxonomy_json <- NULL
  payload$organisms[[i]]$taxonomy_json_sha256 <- NULL
}
jsonlite::write_json(payload,file.path(renamed,"metadata","metadata_manifest.json"),auto_unbox=TRUE,null="null")
local_bundle <- load_manifest_lookup(file.path(renamed,"metadata","metadata_manifest.json"))
stopifnot(length(local_bundle$lookup)==2, all(is.na(vapply(local_bundle$lookup,function(x)x$taxonomy_path,character(1)))))
unlink(work,recursive=TRUE)
cat("All PCA v2 manifest tests passed.\n")
