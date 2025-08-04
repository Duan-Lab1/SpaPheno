#' @title Compute Local Cell Type Composition in the Spatial Neighborhood
#'
#' @description
#' Computes the local cell type composition for each spot or cell based on spatial
#' proximity. Supports both single-cell resolution using nearest neighbors and
#' spot-level resolution using distance thresholding with cell type proportions.
#'
#' @details
#' This function calculates the distribution of cell types around each spatial unit
#' using one of the following modes:
#' \itemize{
#'   \item \code{"single_cell"}: Uses discrete cell-type labels and K-nearest neighbors.
#'   \item \code{"spot"}: Uses continuous deconvolution proportions and radius threshold.
#' }
#'
#' The returned matrix represents normalized cell-type compositions in each local neighborhood.
#'
#' @author Bin Duan
#'
#' @param sample_information_coordinate A data.frame or matrix of spatial coordinates (x, y).
#'   Row names must correspond to cell or spot IDs.
#' @param resolution One of \code{"single_cell"} or \code{"spot"}.
#' @param sample_information_cellType Named vector of cell types for single-cell resolution.
#'   Names must match rownames of \code{sample_information_coordinate}.
#' @param sample_information_decon A matrix of cell type proportions per spot (for \code{"spot"} mode).
#'   Row names must match those of \code{sample_information_coordinate}.
#' @param k Integer, number of nearest neighbors (used only in \code{"single_cell"} mode). Default is 50.
#' @param r Numeric, distance threshold for neighborhood (used only in \code{"spot"} mode). Default is 4.
#'
#' @return A matrix of normalized cell type proportions for each spatial unit (rows) across cell types (columns).
#'
#' @importFrom FNN get.knn
#'
#' @examples
#' \dontrun{
#' data("osmFISH_metadata_cellType")
#' data("osmFISH_bulk_decon")
#' data("osmFISH_bulk_pheno")
#' data("osmFISH_bulk_coordinate")
#'
#' Cell_type_neighborhood(
#'   sample_information_coordinate = osmFISH_bulk_coordinate,
#'   resolution = "single_cell",
#'   sample_information_cellType = osmFISH_metadata_cellType
#' )
#'
#' }
#'
#' @export
#'
Cell_type_neighborhood <- function(
    sample_information_coordinate,
    resolution = c("single_cell", "spot"),
    sample_information_cellType = NULL,
    sample_information_decon = NULL,
    k = 50,
    r = 4) {

  # Match the resolution argument
  resolution <- match.arg(resolution)

  # Check for required dependency
  if (!requireNamespace("FNN", quietly = TRUE)) {
    stop("Package 'FNN' is required. Please install it using install.packages('FNN').")
  }

  # Convert coordinates to matrix and ensure rownames are consistent
  coords <- as.matrix(sample_information_coordinate)
  rownames(coords) <- rownames(sample_information_coordinate)

  # Use SpatialKNN function to get nearest neighbor indices and distances
  knn_result <- SpatialKNN(sample_information_coordinate = coords, k = max(k, 1))
  knn_sample <- knn_result$knn_sample  # matrix of neighbor names (Self + k nearest)
  knn_value <- knn_result$knn_value    # matrix of distances (0 for self in first column)

  # Determine neighbors for each spatial unit
  knn_neighbors <- lapply(seq_len(nrow(knn_sample)), function(i) {
    if (resolution == "single_cell") {
      knn_sample[i, 1:(k + 1)]  # include Self and top-k neighbors
    } else {
      knn_sample[i, which(knn_value[i, ] <= r)]  # neighbors within distance r
    }
  })
  names(knn_neighbors) <- rownames(knn_sample)

  # Compute local composition for single-cell mode
  if (resolution == "single_cell") {
    if (is.null(sample_information_cellType)) {
      stop("Please provide 'sample_information_cellType' for single_cell mode.")
    }

    # Identify unique cell types
    cell_types <- sort(unique(sample_information_cellType))

    # Initialize result matrix
    dist_mat <- matrix(0, nrow = length(knn_neighbors), ncol = length(cell_types))
    rownames(dist_mat) <- names(knn_neighbors)
    colnames(dist_mat) <- cell_types

    # Count cell types in neighbors
    for (i in seq_along(knn_neighbors)) {
      neighbor_types <- sample_information_cellType[knn_neighbors[[i]]]
      ct_tab <- table(factor(neighbor_types, levels = cell_types))
      dist_mat[i, ] <- ct_tab / sum(ct_tab)  # normalize by total count
    }

  } else if (resolution == "spot") {
    if (is.null(sample_information_decon)) {
      stop("Please provide 'sample_information_decon' for spot mode.")
    }

    # Convert deconvolution matrix to standard form
    decon_mat <- as.matrix(sample_information_decon)
    cell_types <- colnames(decon_mat)

    # Initialize result matrix
    dist_mat <- matrix(0, nrow = length(knn_neighbors), ncol = length(cell_types))
    rownames(dist_mat) <- names(knn_neighbors)
    colnames(dist_mat) <- cell_types

    # Aggregate average cell type proportions from neighbors
    for (i in seq_along(knn_neighbors)) {
      neighbor_ids <- knn_neighbors[[i]]
      neighbor_mat <- decon_mat[neighbor_ids, , drop = FALSE]
      dist_mat[i, ] <- if (nrow(neighbor_mat) == 1) neighbor_mat else colMeans(neighbor_mat)
    }

  } else {
    stop("Invalid resolution type.")
  }

  # Normalize each row to sum to 1
  dist_mat <- t(apply(dist_mat, 1, function(x) x / sum(x)))
  return(dist_mat)
}
