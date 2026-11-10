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
# CIELCh hue is undefined near the neutral axis; exclude chroma < 5.
shape_color_coordinates <- function(colors) {
  lab <- grDevices::convertColor(t(grDevices::col2rgb(colors)) / 255, from = "sRGB", to = "Lab")
  chroma <- sqrt(rowSums(lab[, 2:3, drop = FALSE]^2))
  hue <- (atan2(lab[, 3], lab[, 2]) * 180 / pi) %% 360
  hue[chroma < 5] <- NA_real_
  list(lab = lab, hue = hue, chroma = chroma)
}
shape_annotation_edges <- function(annotations, keys) {
  edges <- matrix(integer(), ncol = 2)
  if (ncol(annotations) < 2) return(edges)
  for (j in seq_len(ncol(annotations) - 1L)) {
    a <- match(as.character(annotations[[j]]), keys$label[keys$level == names(annotations)[j]])
    b <- match(as.character(annotations[[j + 1L]]), keys$label[keys$level == names(annotations)[j + 1L]])
    edges <- rbind(edges, cbind(which(keys$level == names(annotations)[j])[a], which(keys$level == names(annotations)[j + 1L])[b]))
  }
  edges <- unique(t(apply(edges, 1, sort)))
  edges[order(edges[, 1], edges[, 2]), , drop = FALSE]
}
shape_color_score <- function(selected, edges, distance, hue_bad) {
  pairs <- if (length(edges)) cbind(selected[edges[, 1]], selected[edges[, 2]]) else matrix(integer(), ncol = 2)
  c(-sum(hue_bad[pairs]), if (nrow(pairs)) min(distance[pairs]) else Inf,
    if (length(selected) > 1) min(distance[selected, selected][upper.tri(diag(length(selected)))]) else Inf)
}
shape_score_better <- function(a, b) {
  different <- which(a != b)
  length(different) > 0 && a[different[1]] > b[different[1]]
}
shape_optimize_colors <- function(keys, edges, candidates, fixed, reserved) {
  pool <- unique(c(keys$color, candidates))
  coord <- shape_color_coordinates(pool)
  distance <- as.matrix(dist(coord$lab))
  delta <- abs(outer(coord$hue, coord$hue, "-"))
  hue_bad <- pmin(delta, 360 - delta) < 60
  hue_bad[is.na(hue_bad)] <- FALSE
  selected <- match(keys$color, pool)
  available <- which(!pool %in% reserved)
  movable <- which(!fixed)
  score <- shape_color_score(selected, edges, distance, hue_bad)
  repeat {
    changed <- FALSE
    for (i in movable) {
      others <- selected[-i]
      incident <- which(edges[, 1] == i | edges[, 2] == i)
      neighbors <- if (length(incident)) ifelse(edges[incident, 1] == i, edges[incident, 2], edges[incident, 1]) else integer()
      remaining <- if (length(incident)) edges[-incident, , drop = FALSE] else edges
      base <- shape_color_score(others, matrix(integer(), ncol = 2), distance, hue_bad)[3]
      base_pairs <- if (nrow(remaining)) cbind(selected[remaining[, 1]], selected[remaining[, 2]]) else matrix(integer(), ncol = 2)
      choices <- setdiff(available, others)
      bad <- rep(-sum(hue_bad[base_pairs]), length(choices))
      adjacent <- rep(if (nrow(base_pairs)) min(distance[base_pairs]) else Inf, length(choices))
      if (length(neighbors)) {
        bad <- bad - rowSums(hue_bad[choices, selected[neighbors], drop = FALSE])
        adjacent <- pmin(adjacent, apply(distance[choices, selected[neighbors], drop = FALSE], 1, min))
      }
      global <- pmin(base, apply(distance[choices, others, drop = FALSE], 1, min))
      best <- order(-bad, -adjacent, -global, choices)[1]
      proposed <- c(bad[best], adjacent[best], global[best])
      if (shape_score_better(proposed, score)) {
        selected[i] <- choices[best]
        score <- proposed
        changed <- TRUE
      }
    }
    if (length(movable) > 1) {
      for (a in seq_len(length(movable) - 1L)) for (b in (a + 1L):length(movable)) {
        trial <- selected
        pair <- movable[c(a, b)]
        trial[pair] <- rev(trial[pair])
        proposed <- shape_color_score(trial, edges, distance, hue_bad)
        if (shape_score_better(proposed, score)) {
          selected <- trial
          score <- proposed
          changed <- TRUE
        }
      }
    }
    if (!changed) break
  }
  pool[selected]
}
# Allocate once across annotation levels; RGB uniqueness includes Unknown labels.
shape_annotation_colors <- function(annotations) {
  keys <- do.call(rbind, lapply(sort(names(annotations)), function(level) {
    data.frame(level = level, label = sort(unique(as.character(annotations[[level]]))), stringsAsFactors = FALSE)
  }))
  rownames(keys) <- NULL
  known <- shape_palette(c("Apicomplexa", "Arthropoda", "Ascomycota", "Basidiomycota", "Chordata", "Cnidaria", "Mucoromycota", "Oomycota", "Unknown"), "phylum")
  keys$color <- NA_character_
  fixed <- keys$level == "phylum" & keys$label %in% names(known)
  keys$color[fixed] <- known[keys$label[fixed]]
  grid <- expand.grid(h = seq(0, 355, 5), c = c(40, 60, 80), l = c(35, 50, 65, 80))
  candidates <- sort(unique(na.omit(grDevices::hcl(grid$h, grid$c, grid$l, fixup = FALSE))))
  candidates <- setdiff(candidates, unname(known))
  remaining <- candidates
  reserved <- shape_color_coordinates(unname(known))$lab
  for (i in which(!fixed)) {
    if (!length(remaining)) stop("Annotation color candidates exhausted")
    candidate_lab <- shape_color_coordinates(remaining)$lab
    distances <- vapply(seq_len(nrow(candidate_lab)), function(j) min(rowSums(sweep(reserved, 2, candidate_lab[j, ])^2)), numeric(1))
    chosen <- which.max(distances)
    keys$color[i] <- remaining[chosen]
    reserved <- rbind(reserved, candidate_lab[chosen, ])
    remaining <- remaining[-chosen]
  }
  keys$color <- shape_optimize_colors(keys, shape_annotation_edges(annotations, keys), candidates, fixed, unname(known))
  keys$identifier <- character(nrow(keys))
  prefixes <- c(phylum = "P", class = "C", order = "O", cluster = "K")
  for (level in unique(keys$level)) {
    rows <- which(keys$level == level)
    prefix <- if (level %in% names(prefixes)) prefixes[[level]] else level
    keys$identifier[rows] <- paste0(prefix, if (level == "cluster") keys$label[rows] else seq_along(rows))
  }
  stopifnot(!anyDuplicated(keys$color), !anyDuplicated(keys$identifier))
  keys
}
shape_annotation_runs <- function(labels, index = seq_along(labels)) {
  runs <- rle(as.character(labels[index]))
  end <- cumsum(runs$lengths)
  start <- end - runs$lengths + 1L
  data.frame(label = runs$values, start = start, end = end, center = (start + end) / 2,
    y = 1 - ((start + end) / 2 - 0.5) / length(index))
}
# Pure geometry: text widths are measured at 8 pt by the active device.
shape_cluster_bracket_layout <- function(labels, index, body_height_mm, text_width_mm) {
  runs <- shape_annotation_runs(labels, index)
  fraction <- (runs$end - runs$start + 1) / length(index)
  runs$text <- paste("Cluster", runs$label)
  runs$top <- 1 - (runs$start - 1) / length(index)
  runs$bottom <- 1 - runs$end / length(index)
  available <- pmax(0, body_height_mm * fraction - 2)
  runs$fontsize <- pmax(6, pmin(8, 8 * available / text_width_mm))
  runs$text_height_mm <- text_width_mm * runs$fontsize / 8
  runs$fits <- runs$text_height_mm + 2 <= body_height_mm * fraction + 1e-8
  list(runs = runs, required_body_mm = max((text_width_mm * 6 / 8 + 2) / fraction))
}
shape_contrast_text <- function(colors) {
  rgb <- t(grDevices::col2rgb(colors)) / 255
  linear <- ifelse(rgb <= 0.04045, rgb / 12.92, ((rgb + 0.055) / 1.055)^2.4)
  luminance <- as.vector(linear %*% c(0.2126, 0.7152, 0.0722))
  ifelse((luminance + 0.05) / 0.05 >= 1.05 / (luminance + 0.05), "black", "white")
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
