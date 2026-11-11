# Build feature vectors and choose the k-means cluster count.

compute_observed_vector <- function(df) {
  if (!all(c("mutType", "mutNum", "totalRootNum") %in% names(df))) {
    stop("TSV file is missing required columns (mutType, mutNum, totalRootNum)")
  }
  totals <- as.numeric(df$totalRootNum)
  mut_counts <- as.numeric(df$mutNum)
  observed <- ifelse(totals != 0, mut_counts / totals, 0)
  observed
}

compute_logratio_vector <- function(df, observed = NULL) {
  if (is.null(observed)) observed <- compute_observed_vector(df)
  mean_observed <- ifelse(length(observed) > 0, mean(observed), 0)
  adjusted_mean <- max(mean_observed, EPS)
  ratios <- (observed + EPS) / adjusted_mean
  log2(ratios)
}

determine_clusters <- function(data, max_clusters, seed) {
  sample_count <- nrow(data)
  if (sample_count < 2) stop("At least two samples are required for clustering.")
  if (sample_count == 2) {
    return(list(best_k = 2, inertia = numeric(), silhouette = numeric()))
  }
  adjusted_max <- min(max_clusters, sample_count - 1)
  if (adjusted_max < 2) adjusted_max <- 2
  set.seed(seed)
  inertia <- numeric()
  silhouette_scores <- numeric()
  for (k in 2:adjusted_max) {
    model <- kmeans(data, centers = k, nstart = 10)
    inertia <- c(inertia, model$tot.withinss)
    silhouette_scores <- c(silhouette_scores, mean(cluster::silhouette(model$cluster, stats::dist(data))[, "sil_width"]))
  }
  best_k <- which.max(silhouette_scores) + 2 - 1
  list(best_k = best_k, inertia = inertia, silhouette = silhouette_scores)
}

compute_pca_projection <- function(vectors, feature_order, indices) {
  subset_vectors <- vectors[indices]
  matrix_data <- do.call(rbind, subset_vectors)
  feature_sd <- apply(matrix_data,2,sd)
  feature_sd[!is.finite(feature_sd) | feature_sd == 0] <- 1
  scaled_data <- scale(matrix_data,scale=feature_sd)
  pca <- prcomp(scaled_data, center = FALSE, scale. = FALSE)
  pca_data <- pca$x[, 1:2, drop = FALSE]
  if (is.null(colnames(pca$rotation))) {
    colnames(pca$rotation) <- feature_order
  }
  list(pca = pca, data = pca_data)
}
