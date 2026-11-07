shape_channels <- function() {
  unlist(lapply(c("C>A", "C>G", "C>T", "T>A", "T>C", "T>G"), function(s) {
    unlist(lapply(c("A", "C", "G", "T"), function(a) paste0(a, "[", s, "]", c("A", "C", "G", "T"))))
  }), use.names = FALSE)
}
shape_vector <- function(df) {
  if (!all(c("mutType", "mutNum", "totalRootNum") %in% names(df))) stop("missing required columns")
  expected <- shape_channels()
  bad <- c(setdiff(expected, df$mutType), setdiff(df$mutType, expected), df$mutType[duplicated(df$mutType)])
  if (length(bad)) stop("invalid/missing/duplicate channels: ", paste(unique(bad), collapse = ","))
  df <- df[match(expected, df$mutType), ]
  count <- suppressWarnings(as.numeric(df$mutNum))
  total <- suppressWarnings(as.numeric(df$totalRootNum))
  rate <- count / total
  invalid <- !is.finite(count) | count <= 0 | !is.finite(total) | total <= 0 | !is.finite(rate) | rate <= 0
  if (any(invalid)) stop("zero/missing/invalid frequency: ", paste(expected[invalid], collapse = ","))
  setNames(log2(rate / mean(rate)), expected)
}
shape_cluster <- function(x, method = "ward", clusters = "auto", max_clusters = 10L) {
  if (nrow(x) < 2) stop("At least two valid samples are required")
  d <- if (method == "ward") {
    dist(x)
  } else {
    if (any(apply(x, 1, sd) == 0)) stop("Correlation undefined for constant spectrum")
    z <- 1 - cor(t(x))
    z[] <- pmax(0, pmin(2, z))
    diag(z) <- 0
    as.dist(z)
  }
  tree <- hclust(d, method = if (method == "ward") "ward.D2" else "average")
  scores <- data.frame(k = integer(), silhouette = numeric())
  if (clusters == "auto") {
    k <- 1L
    if (nrow(x) > 2 && any(d > 0)) {
      for (candidate in 2:min(max_clusters, nrow(x) - 1L)) {
        s <- mean(cluster::silhouette(cutree(tree, candidate), d)[, "sil_width"])
        scores <- rbind(scores, data.frame(k = candidate, silhouette = s))
      }
      k <- scores$k[which.max(scores$silhouette)]
    }
  } else {
    k <- suppressWarnings(as.integer(clusters))
    if (is.na(k) || as.character(k) != clusters || k < 1 || k > nrow(x)) stop("Invalid --clusters")
  }
  membership <- cutree(tree, k)
  dm <- as.matrix(d)
  representatives <- vapply(sort(unique(membership)), function(g) {
    ids <- sort(names(membership)[membership == g])
    ids[which.min(rowSums(dm[ids, ids, drop = FALSE]))]
  }, character(1))
  list(distance = d, tree = tree, membership = membership, scores = scores, representatives = representatives, k = k)
}
shape_palette <- function(labels, level) {
  labels <- sort(unique(as.character(labels)))
  colors <- setNames(grDevices::hcl.colors(length(labels), "Dark 3"), labels)
  known <- c(Apicomplexa = "#E69F00", Arthropoda = "#0072B2", Ascomycota = "#009E73", Basidiomycota = "#F0E442", Chordata = "#D55E00", Cnidaria = "#CC79A7", Mucoromycota = "#56B4E9", Oomycota = "#000000", Unknown = "#999999")
  if (level == "phylum") colors[intersect(labels, names(known))] <- known[intersect(labels, names(known))]
  if ("Unknown" %in% labels) colors["Unknown"] <- "#999999"
  colors
}
shape_umap <- function(x, method, seed) {
  if (nrow(x) < 4) {
    return(NULL)
  }
  if (method == "correlation") {
    x <- x - rowMeans(x)
    x <- x / sqrt(rowSums(x * x))
  }
  set.seed(seed)
  uwot::umap(x,
    n_neighbors = min(15L, nrow(x) - 1L), min_dist = 0.1, n_components = 2,
    metric = "euclidean", n_threads = 1, n_sgd_threads = 1, init = "random", verbose = FALSE
  )
}
