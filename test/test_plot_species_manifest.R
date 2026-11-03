#!/usr/bin/env Rscript
this <- sub("^--file=","",grep("^--file=",commandArgs(),value=TRUE)[1])
root <- dirname(dirname(normalizePath(this)))
source(file.path(root,"src/visualize/plot_species_label.R"))
work <- tempfile("plot-label-v2-",tmpdir=file.path(root,"tmp/codex"))
dir.create(work)
dir.create(file.path(work,"metadata"))
dir.create(file.path(work,"arbitrary-parent"))
input <- file.path(work,"arbitrary-parent","unrelated-name.tsv")
writeLines("counts",input)
manifest <- list(schema_version=2,organisms=list(
 list(slot="org2",raw_organism_name="Genus_species"),
 list(slot="org3",raw_organism_name="Other_species")),
 artifacts=list(list(path="arbitrary-parent/unrelated-name.tsv",slot="org2",kind="single")))
path <- file.path(work,"metadata/metadata_manifest.json")
jsonlite::write_json(manifest,path,auto_unbox=TRUE)
stopifnot(identical(plot_species_name(input),"Genus species"))
renamed <- file.path(work,"even-more-arbitrary.tsv")
stopifnot(file.rename(input,renamed))
manifest$artifacts[[1]]$path <- "even-more-arbitrary.tsv"
jsonlite::write_json(manifest,path,auto_unbox=TRUE)
stopifnot(identical(plot_species_name(renamed),"Genus species"))
manifest$organisms[[1]]$raw_organism_name <- NULL
jsonlite::write_json(manifest,path,auto_unbox=TRUE)
stopifnot(identical(suppressWarnings(plot_species_name(renamed)),""))
stopifnot(identical(plot_species_name(renamed,"Explicit species"),"Explicit species"))
unlink(work,recursive=TRUE)
cat("All manifest plot label tests passed.\n")
