#' @title Map Bulk Phenotype to Spatial Transcriptomics via Interpretable Modeling
#'
#' @description
#' Perform interpretable phenotype projection from bulk transcriptomic phenotypes
#' to spatial transcriptomics data using elastic net modeling and SHAP-based spatial mapping.
#'
#' @author Bin Duan
#'
#' @param bulk_decon A matrix or data.frame with bulk samples (rows) and cell type proportions (columns), typically from deconvolution.
#' @param bulk_pheno A named vector of bulk phenotypes. Can be binary labels (e.g., 0/1), survival time, or continuous traits. Names must match \code{rownames(bulk_decon)}.
#' @param family Modeling family. One of \code{"cox"} (survival), \code{"binomial"} (classification), or \code{"gaussian"} (regression).
#' @param coord A named list of spatial coordinate matrices for each spatial sample (data.frame with x, y).
#' @param resolution Spatial resolution. One of \code{"spot"} (e.g., 10x Visium) or \code{"single_cell"} (e.g., MERFISH, CosMx).
#' @param sample_information_decon A list of spot-level cell type proportions (only used if \code{resolution = "spot"}). Each entry is a data.frame with rows as spots and columns as cell types.
#' @param sample_information_cellType A list of single-cell annotations (only used if \code{resolution = "single_cell"}). Each entry is a character vector of cell types.
#' @param k Number of spatial neighbors to consider for KNN-based feature construction. Default is 50.
#' @param r Number of spatial shells/layers for context expansion. Default is 2.
#' @param p SHAP value threshold for phenotype labeling and significance testing. Default is 0.
#' @param size Point size used in visualization plots. Default is 1.
#' @param n_perm Number of permutations for significance estimation. Default is 10.
#' @param seed Random seed for reproducibility. Default is 2025.
#' @param nfolds Number of cross-validation folds for elastic net modeling. Default is 10.
#'
#' @return A list containing:
#' \item{pred_score}{A data.frame of prediction scores and assigned phenotype labels.}
#' \item{cell_type_distribution}{Feature matrix used for prediction.}
#' \item{model}{Trained elastic net model object.}
#'
#' @examples
#' \dontrun{
#'
#' data("osmFISH_metadata_cellType")
#' data("osmFISH_bulk_decon")
#' data("osmFISH_bulk_pheno")
#' data("osmFISH_bulk_coordinate")
#'
#' result <- SpatialPhenoMap(
#'   bulk_decon = osmFISH_bulk_decon,
#'   bulk_pheno = osmFISH_bulk_pheno,
#'   family = "binomial",
#'   coord = osmFISH_bulk_coordinate,
#'   resolution = "single_cell",
#'   sample_information_cellType = osmFISH_metadata_cellType
#' )
#' }
#'
#' @export
#'
SpatialPhenoMap <- function(
    bulk_decon,
    bulk_pheno,
    family = c("cox", "binomial", "gaussian"),
    coord,
    resolution = "spot",
    sample_information_decon = NULL,
    sample_information_cellType = NULL,
    k = 50,
    r = 2,
    p = 0,
    size = 1,
    n_perm = 10,
    seed = 2025,
    nfolds = 10) {

  stopifnot(resolution %in% c("spot", "single_cell"))
  stopifnot(!is.null(bulk_decon), !is.null(bulk_pheno), !is.null(coord))

  message("Building bulk phenotype model (with automatic alpha)...")
  family <- family[1]
  bulk_decon <- t(apply(bulk_decon, 1, function(x) (x / sum(x))))
  colnames(bulk_decon) <- make.names(colnames(bulk_decon))

  model <- BuildPhenoModelAutoAlpha(
    expr = bulk_decon,
    pheno = bulk_pheno,
    family = family,
    seed = seed,
    nfolds = nfolds
  )

  message("Constructing spatial cell type distribution...")
  if (resolution == "spot") {
    if (is.null(sample_information_decon)) stop("Missing sample_information_decon.")
    colnames(sample_information_decon) <- make.names(colnames(sample_information_decon))
    sample_information_decon <- sample_information_decon[, colnames(bulk_decon)]
    feature_matrix <- Cell_type_neighborhood(
      sample_information_coordinate = coord,
      resolution = "spot",
      sample_information_decon = sample_information_decon,
      k = k,
      r = r
    )
  } else {
    if (is.null(sample_information_cellType)) stop("Missing sample_information_cellType.")
    sample_name <- names(sample_information_cellType)
    sample_information_cellType <- make.names(sample_information_cellType)
    names(sample_information_cellType) <- sample_name
    feature_matrix <- Cell_type_neighborhood(
      sample_information_coordinate = coord,
      resolution = "single_cell",
      sample_information_cellType = sample_information_cellType,
      k = k,
      r = r
    )
  }

  message("Predicting phenotype scores...")
  pred_score <- stats::predict(
    object = model$model,
    family = family,
    s = model$lambda,
    newx = as.matrix(feature_matrix)
  )

  message("Generating permutation scores...")
  perm_scores <- GeneratePermutations(
    model = model,
    family = family,
    coordinate = coord,
    resolution = resolution,
    r = r,
    k = k,
    n = n_perm,
    seed = seed,
    sample_information_decon = sample_information_decon,
    sample_information_cellType = sample_information_cellType
  )

  message("Drawing plots...")
  FeatureImportanceCoef(model)
  ThresholdPlot(pred_score, perm_scores, p = p)
  PhenotypePlot(coord, pred_score, perm_scores, size = size, p = p)

  # Assign labels
  high_thre <- quantile(perm_scores, 1 - p)
  low_thre <- quantile(perm_scores, p)
  pred_score <- as.data.frame(pred_score)
  colnames(pred_score) <- "pred_score"
  pred_score$label <- ifelse(
    pred_score$pred_score > high_thre, "phenotype+",
    ifelse(pred_score$pred_score < low_thre, "phenotype-", "background")
  )

  print(PhenoAbundancePlot(feature_matrix, pred_score))

  return(list(
    pred_score = pred_score,
    cell_type_distribution = feature_matrix,
    model = model
  ))
}
