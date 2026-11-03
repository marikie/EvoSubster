# Required plot fonts. No substitute family is accepted.
load_required_fonts <- function() {
  cache <- Sys.getenv("EVOSUBSTER_FONT_CACHE", unset = "")
  if (!nzchar(cache)) {
    cache <- file.path(path.expand("~"), ".cache", "evo-subster", "fonts")
  }
  dir.create(cache, recursive = TRUE, showWarnings = FALSE)
  specs <- list(mn = list(name = "Courier Prime", regular = "700", bold = "700"),
                os = list(name = "Roboto", regular = "regular", bold = "700"))
  database <- NULL
  for (family in names(specs)) {
    spec <- specs[[family]]
    paths <- list()
    for (face in c("regular", "bold")) {
      target <- file.path(cache, paste0(gsub(" ", "-", spec$name), "-", spec[[face]], ".ttf"))
      if (!file.exists(target) || file.info(target)$size == 0) {
        if (is.null(database)) database <- sysfonts:::google_font_db(TRUE, curl::new_handle())
        records <- database[[2]]
        matched <- Filter(function(x) identical(x$family, spec$name), records)
        if (length(matched) != 1) stop("Required Google font not found: ", spec$name)
        url <- matched[[1]]$files[[spec[[face]]]]
        if (is.null(url)) stop("Required font face unavailable: ", spec$name, " ", face)
        url <- sub("^http:", "https:", url)
        temporary <- tempfile(pattern = ".font-", tmpdir = cache)
        tryCatch({
          download.file(url, temporary, mode = "wb", quiet = TRUE, method = "libcurl")
          if (!file.exists(temporary) || file.info(temporary)$size == 0) stop("Empty font download")
          # Loading verifies that the downloaded file is a readable font.
          sysfonts::font_add(paste0("check-", family, "-", face), regular = temporary)
          if (!file.rename(temporary, target)) stop("Cannot publish required font cache")
        }, error = function(e) stop("Required font unavailable: ", spec$name, " ", face, ": ", e$message),
        finally = unlink(temporary))
      }
      paths[[face]] <- target
    }
    sysfonts::font_add(family, regular = paths$regular, bold = paths$bold)
  }
  showtext::showtext_auto()
}
