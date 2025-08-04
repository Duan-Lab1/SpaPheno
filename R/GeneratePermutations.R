#' @title Generate Permutation Scores for Model Significance Testing
#'
#' @description
#' Generates phenotype prediction scores based on permuted spatial coordinates,
#' enabling permutation-based significance testing of spatial features in glmnet models.
#'
#' @details
#' This function permutes the row names of the coordinate matrix to shuffle spatial positions,
#' recalculates the local cell type neighborhood features via \code{Cell_type_neighborhood()},
#' and predicts phenotype scores using the input glmnet model.
#' The process is repeated \code{n} times to generate a null distribution of prediction scores.
#'
#' The function relies on:
#' \itemize{
#'   \item \code{sample_information_decon}: global matrix of cell type proportions (for spot mode).
#'   \item \code{sample_information_cellType}: global named vector of cell types (for single-cell mode).
#' }
#'
#' @author Bin Duan
#'
#' @param model A fitted glmnet model object, typically from \code{cv.glmnet()}, containing \code{model} and \code{lambda}.
#' @param resolution Character; one of \code{"spot"} or \code{"single_cell"}.
#' @param coordinate A matrix or data.frame of spatial coordinates. Row names must match spot or cell IDs.
#' @param family Character; glmnet family used for prediction (\code{"cox"}, \code{"binomial"}, or \code{"gaussian"}).
#' @param r Numeric; radius for neighborhood definition (used for \code{"spot"}). Default is 2.
#' @param k Integer; number of neighbors (used for \code{"single_cell"}). Default is 50.
#' @param n Integer; number of permutations to run. Default is 50.
#' @param seed Integer; random seed for reproducibility. Default is 2025.
#' @param sample_information_decon A list of spot-level cell type proportions (only used if \code{resolution = "spot"}). Each entry is a data.frame with rows as spots and columns as cell types.
#' @param sample_information_cellType A list of single-cell annotations (only used if \code{resolution = "single_cell"}). Each entry is a character vector of cell types.
#'
#' @return A numeric vector of predicted phenotype scores from each permutation replicate.
#'
#' @importFrom stats predict
#'
#' @examples
#' \dontrun{
#'
#' data("osmFISH_bulk_decon")
#' data("osmFISH_bulk_pheno")
#' data("osmFISH_bulk_coordinate")
#'
#' model <- BuildPhenoModelAutoAlpha(
#'   expr = osmFISH_bulk_decon,
#'   pheno = osmFISH_bulk_pheno,
#'   family = "binomial")
#'
#' data("osmFISH_metadata_cellType")
#' perm_scores <- GeneratePermutations(
#'   model = model,
#'   resolution = "single_cell",
#'   coordinate = osmFISH_bulk_coordinate,
#'   family = "binomial",
#'   r = 2,
#'   n = 100,
#'   sample_information_cellType = osmFISH_metadata_cellType
#' )
#' hist(perm_scores, breaks = 20, main = "Permutation Null Distribution")
#' }
#'
#' @export
#'
GeneratePermutations <- function(
    model,
    resolution = c("spot", "single_cell"),
    coordinate = NULL,
    family = c("cox", "binomial", "gaussian"),
    r = 2,
    k = 50,
    n = 50,
    seed = 2025,
    sample_information_decon = NULL,
    sample_information_cellType = NULL) {

  set.seed(seed)
  resolution <- match.arg(resolution)
  family <- match.arg(family)

  # Validate required data
  if (resolution == "spot" && is.null(sample_information_decon)) {
    stop("Please provide 'sample_information_decon' for spot resolution.")
  }
  if (resolution == "single_cell" && is.null(sample_information_cellType)) {
    stop("Please provide 'sample_information_cellType' for single_cell resolution.")
  }

  # Permutation loop
  perm_scores <- replicate(n, {
    # Shuffle spatial positions
    perm_coord <- coordinate
    rownames(perm_coord) <- sample(rownames(coordinate))

    # Recalculate neighborhood feature matrix
    feature_matrix <- Cell_type_neighborhood(
      sample_information_coordinate = perm_coord,
      resolution = resolution,
      sample_information_decon = sample_information_decon,
      sample_information_cellType = sample_information_cellType,
      k = k,
      r = r
    )

    # Predict phenotype scores using permuted features
    predict(
      object = model$model,
      family = family,
      newx = as.matrix(feature_matrix),
      s = model$lambda
    )
  })

  return(as.numeric(perm_scores))
}
