# Helpers loaded by trio_selection.R into its calling environment.

defaults <- list(
  tree = NULL,
  idt_threshold = 80,
  ingroup_idt_threshold = NULL,
  outgroup_idt_threshold = NULL,
  min_idt_gap_between_out_in = 0,
  max_outgroup_tries = 5,
  ingroup_pairing = "matching",
  min_rel_contig_n50 = 0,
  stage0_top_k = 3L,
  assembly_qc = NULL,
  allow_qc_override = FALSE,
  genome_dir = "./genomes",
  out_dir = "./results/trio_selection",
  train_cache_dir = NULL,
  date = format(Sys.Date(), "%Y%m%d"),
  threads = 8,
  dry_run = FALSE,
  keep_unused_species_data = FALSE
)

usage <- function() {
  cat(
    "Usage: Rscript src/select/trio_selection.R --tree <file.nwk> [options]",
    "",
    "Options:",
    "  --tree FILE             Newick tree; leaf labels are accession-free complete taxon names",
    "                          with spaces encoded as underscores (for example,",
    "                          Chaunax_sp._Z400). Legacy accession-suffixed trees must first",
    "                          be converted with src/select/strip_newick_accessions.py.",
    "  --idt-threshold N       Common identity threshold (default: 80); individual values override it.",
    "  --ingroup-idt-threshold N  Ingroup identity must strictly exceed this (default: 80).",
    "  --outgroup-idt-threshold N  Both outgroup-ingroup identities must strictly exceed this (default: 80).",
    "  --min-idt-gap-between-out-in N  Both ingroup minus outgroup-ingroup identity gaps must",
    "                          strictly exceed N percentage points (default: 0). Exempt when",
    "                          ingroups share a genus distinct from the outgroup.",
    "  --max-outgroup-tries N  Give up on an ingroup pair after training this many outgroup",
    "                          candidates without a thesis-rule pass (default: 5). Caps the",
    "                          per-pair last-train cost.",
    "  --ingroup-pairing MODE  Which ingroup couples to consider (default: matching).",
    "                            matching  greedy closest-first disjoint matching -- each",
    "                                      species used at most once, ~n/2 independent",
    "                                      sister-pair couples (phylogenetic-pairs design).",
    "                            all       every tip pair (exhaustive; C(n,2) couples).",
    "  --min-rel-contig-n50 F  Drop a species unless some current assembly has relative contig",
    "                          N50 (contig_n50 / total_ungapped_length) >= F. Size-normalized,",
    "                          so one threshold works across lineages -- absolute N50 is",
    "                          meaningless across taxa (a perfect ~12 Mb yeast reference has a",
    "                          0.92 Mb contig N50). Default 0 = off (keep every species);",
    "                          e.g. 0.005 = 0.5%.",
    "  --stage0-top-k N        Keep N candidates per species in the metadata shortlist for QC audit and optional explicit strict override (default: 3).",
    "                          The eligible NCBI reference",
    "                          genome, when present, is always anchored in this shortlist.",
    "  --assembly-qc FILE      Optional TSV with accession plus comparable external QC:",
    "                          BUSCO genome-mode fields (qc_busco_*) and optional Merqury",
    "                          fields. NCBI annotation BUSCO is recorded but is not treated",
    "                          as assembly-level BUSCO genome-mode evidence.",
    "  --allow-qc-override    Explicitly allow a strict BUSCO-based override of the metadata",
    "                          baseline. Off by default; external QC is otherwise audit-only.",
    "  --genome-dir DIR        Genome storage (default: ./genomes).",
    "  --out-dir DIR           Output directory (default: ./results/trio_selection).",
    "  --train-cache-dir DIR   Directory for cached last-train files (default: <out-dir>/train_cache).",
    "  --date YYYYMMDD         Run date, used in cache filenames (default: today).",
    "  --threads N             Threads for lastdb/last-train (default: 8).",
    "  --dry-run               Enumerate ingroup pairs + candidate-outgroup order only;",
    "                          download nothing and run no last-train.",
    "  --keep-unused-species-data",
    "                          Keep genome and train artifacts created for candidates that are",
    "                          absent from every selected trio (default: remove them).",
    "",
    sep = "\n"
  )
}

parse_args <- function(argv) {
  opts <- defaults
  specified <- character()
  i <- 1
  while (i <= length(argv)) {
    key <- argv[i]
    specified <- c(specified, gsub("-", "_", sub("^--", "", key)))
    take <- function() {
      if (i + 1 > length(argv)) stop("Missing value for ", key, call. = FALSE)
      argv[i + 1]
    }
    switch(key,
      "--tree"               = { opts$tree <- take(); i <- i + 2 },
      "--idt-threshold"      = { opts$idt_threshold <- as.numeric(take()); i <- i + 2 },
      "--ingroup-idt-threshold" = { opts$ingroup_idt_threshold <- as.numeric(take()); i <- i + 2 },
      "--outgroup-idt-threshold" = { opts$outgroup_idt_threshold <- as.numeric(take()); i <- i + 2 },
      "--min-idt-gap-between-out-in" = { opts$min_idt_gap_between_out_in <- as.numeric(take()); i <- i + 2 },
      "--max-outgroup-tries" = { opts$max_outgroup_tries <- as.integer(take()); i <- i + 2 },
      "--ingroup-pairing"    = { opts$ingroup_pairing <- take(); i <- i + 2 },
      "--min-rel-contig-n50" = { opts$min_rel_contig_n50 <- as.numeric(take()); i <- i + 2 },
      "--stage0-top-k"       = { opts$stage0_top_k <- as.integer(take()); i <- i + 2 },
      "--assembly-qc"        = { opts$assembly_qc <- take(); i <- i + 2 },
      "--allow-qc-override"  = { opts$allow_qc_override <- TRUE; i <- i + 1 },
      "--genome-dir"         = { opts$genome_dir <- take(); i <- i + 2 },
      "--out-dir"            = { opts$out_dir <- take(); i <- i + 2 },
      "--train-cache-dir"    = { opts$train_cache_dir <- take(); i <- i + 2 },
      "--date"               = { opts$date <- take(); i <- i + 2 },
      "--threads"            = { opts$threads <- as.integer(take()); i <- i + 2 },
      "--dry-run"            = { opts$dry_run <- TRUE; i <- i + 1 },
      "--keep-unused-species-data" = { opts$keep_unused_species_data <- TRUE; i <- i + 1 },
      "-h"                   = { usage(); quit(status = 0) },
      "--help"               = { usage(); quit(status = 0) },
      stop("Unknown option: ", key, call. = FALSE)
    )
  }
  if (is.null(opts$tree)) {
    usage()
    stop("--tree is required.", call. = FALSE)
  }
  if (!opts$ingroup_pairing %in% c("matching", "all")) {
    stop("--ingroup-pairing must be 'matching' or 'all' (got: ", opts$ingroup_pairing, ")",
         call. = FALSE)
  }
  if (is.na(opts$min_rel_contig_n50) || opts$min_rel_contig_n50 < 0 ||
      opts$min_rel_contig_n50 > 1) {
    stop("--min-rel-contig-n50 must be a fraction in [0, 1] (got: ", opts$min_rel_contig_n50,
         "); it is contig_n50 / total_ungapped_length, e.g. 0.005 for 0.5%.", call. = FALSE)
  }
  if (is.na(opts$stage0_top_k) || opts$stage0_top_k < 1) {
    stop("--stage0-top-k must be a positive integer (got: ", opts$stage0_top_k, ")",
         call. = FALSE)
  }
  opts <- normalize_identity_options(opts)
  attr(opts, "specified_options") <- unique(specified)
  opts
}

normalize_identity_options <- function(opts) {
  # Resolve only after parsing: individual settings win regardless of CLI order.
  if (is.null(opts$ingroup_idt_threshold)) opts$ingroup_idt_threshold <- opts$idt_threshold
  if (is.null(opts$outgroup_idt_threshold)) opts$outgroup_idt_threshold <- opts$idt_threshold
  for (key in c("idt_threshold", "ingroup_idt_threshold", "outgroup_idt_threshold",
                "min_idt_gap_between_out_in")) {
    value <- opts[[key]]
    if (length(value) != 1L || !is.finite(value) || value < 0 || value > 100) {
      stop("--", gsub("_", "-", key), " must be a finite number in [0, 100].", call. = FALSE)
    }
  }
  opts
}
