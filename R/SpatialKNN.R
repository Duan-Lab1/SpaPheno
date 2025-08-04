#' @title Compute K-Nearest Neighbors for Spatial Coordinates
#'
#' @description
#' Given a matrix of spatial coordinates, this function computes the k-nearest neighbors (KNN)
#' for each point using Euclidean distance. It returns the neighbor indices and distances.
#'
#' @details
#' - The input coordinate matrix must have rownames corresponding to spot or cell IDs.
#' - The output includes a neighbor matrix with the first column as the self point.
#' - Distances to self are set to 0 by definition.
#'
#' @author Bin Duan
#'
#' @param sample_information_coordinate A data.frame or matrix. Rows are spatial units (e.g., cells or spots), columns are coordinates (e.g., x and y).
#' @param k Integer. Number of nearest neighbors to compute (default = 50).
#'
#' @return A list with:
#' \item{knn_value}{A matrix of distances to the k nearest neighbors (first column is self, with distance 0).}
#' \item{knn_sample}{A matrix of neighbor names (cell/spot IDs), with self in the first column.}
#'
#' @importFrom FNN get.knn
#'
#' @examples
#' \dontrun{
#' coords <- data.frame(x = runif(10), y = runif(10))
#' rownames(coords) <- paste0("cell", 1:10)
#' result <- SpatialKNN(coords, k = 3)
#' result$knn_sample
#' result$knn_value
#' }
#'
#' @export
#'
SpatialKNN <- function(sample_information_coordinate, k = 50) {
  if (!requireNamespace("FNN", quietly = TRUE)) {
    stop("Package 'FNN' is required. Please install it using install.packages('FNN').")
  }

  coord_mat <- as.matrix(sample_information_coordinate)
  rownames(coord_mat) <- rownames(sample_information_coordinate)

  knn_result <- FNN::get.knn(coord_mat, k = k)

  knn_sample <- apply(knn_result$nn.index, 1, function(idx) rownames(coord_mat)[idx])
  knn_sample <- t(knn_sample)
  rownames(knn_sample) <- rownames(coord_mat)
  knn_sample <- cbind(Self = rownames(coord_mat), knn_sample)

  knn_value <- cbind(Self = rep(0, nrow(coord_mat)), knn_result$nn.dist)

  return(list(
    knn_value = knn_value,
    knn_sample = knn_sample
  ))
}

