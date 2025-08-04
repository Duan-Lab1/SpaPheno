#' @title Compute SHAP Values for Spatial Transcriptomics Data
#'
#' @description
#' This function computes SHAP (SHapley Additive exPlanations) values for spatial
#' transcriptomics data using a predictive model trained on bulk RNA-seq data.
#' It returns a long-format data frame with one row per sample-feature pair.
#'
#' @details
#' \itemize{
#'   \item Normalizes \code{X_bulk} by row-wise sum to ensure comparability across samples.
#'   \item Replaces \code{NA} and \code{NaN} in \code{X_spatial} with 0 to ensure model stability.
#'   \item Ensures feature alignment and column order consistency between bulk and spatial datasets.
#'   \item Computes SHAP values using the \code{iml::Shapley} class for each spatial unit.
#' }
#'
#' @author Bin Duan
#'
#' @param model A fitted glmnet model object (e.g., from \code{cv.glmnet}) with named elements \code{$model} and \code{$lambda}.
#' @param X_bulk A numeric matrix or data.frame of bulk RNA-seq gene expression (samples × features).
#' @param y_bulk A numeric vector of phenotype labels corresponding to \code{X_bulk}.
#' @param X_spatial A numeric matrix or data.frame of spatial transcriptomics expression (samples × features).
#'
#' @return A data.frame of SHAP values for each sample-feature pair with columns:
#' \itemize{
#'   \item \code{feature}: Feature (gene) name.
#'   \item \code{phi}: The SHAP value of that feature.
#'   \item \code{sample}: The spatial sample identifier.
#'   \item \code{feature.value}: Original feature expression value.
#' }
#'
#' @importFrom iml Predictor Shapley
#' @importFrom stats predict
#' @importFrom utils head
#'
#' @examples
#' \dontrun{
#'
#' data("osmFISH_metadata_cellType")
#' data("osmFISH_bulk_decon")
#' data("osmFISH_bulk_pheno")
#' data("osmFISH_bulk_coordinate")
#'
#' PhenoResult <- SpatialPhenoMap(
#'   bulk_decon = osmFISH_bulk_decon,
#'   bulk_pheno = osmFISH_bulk_pheno,
#'   family = "binomial",
#'   coord = osmFISH_bulk_coordinate,
#'   resolution = "single_cell",
#'   sample_information_cellType = osmFISH_metadata_cellType
#' )
#'
#' pred_result <- PhenoResult$pred_score
#' phenoPlus <- row.names(pred_result[pred_result$label %in% "phenotype+", ])
#' model <- PhenoResult$model
#' X <- as.data.frame(PhenoResult$cell_type_distribution[phenoPlus, ])
#'
#' # This step took a very long time
#' shap_test_plus <- compute_shap_spatial(
#'   model = model,
#'   X_bulk = as.data.frame(osmFISH_bulk_decon),
#'   y_bulk = osmFISH_bulk_pheno,
#'   X_spatial = X)
#' }
#'
#' @export
#'
compute_shap_spatial <- function(model, X_bulk, y_bulk, X_spatial) {

  # Replace NA and NaN values in spatial expression with zeros
  X_spatial[is.na(X_spatial)] <- 0
  X_spatial[is.nan(as.matrix(X_spatial))] <- 0
  X_spatial <- as.data.frame(X_spatial)

  # Ensure syntactic column names
  colnames(X_bulk) <- make.names(colnames(X_bulk))
  colnames(X_spatial) <- make.names(colnames(X_spatial))

  # Align spatial features with bulk features
  common_genes <- intersect(colnames(X_bulk), colnames(X_spatial))
  missing_genes <- setdiff(colnames(X_bulk), common_genes)

  # Add missing genes to spatial data with value 0
  for (gene in missing_genes) {
    X_spatial[[gene]] <- 0
  }

  # Reorder columns to match bulk
  X_spatial <- X_spatial[, colnames(X_bulk), drop = FALSE]

  # Normalize bulk data row-wise
  X_bulk <- t(apply(X_bulk, 1, function(x) x / sum(x)))

  # Create Predictor object from iml package
  predictor <- iml::Predictor$new(
    model = model$model,
    data = as.data.frame(X_bulk),
    y = y_bulk,
    predict.function = function(m, newdata) {
      stats::predict(m, newx = as.matrix(newdata), s = model$lambda, type = "response") |> 
        as.numeric()
    }
  )

  # Compute SHAP values for each spatial sample
  shap_list <- lapply(1:nrow(X_spatial), function(i) {
    sample_i <- rownames(X_spatial)[i]
    x_i <- X_spatial[i, , drop = FALSE]

    shap_obj <- iml::Shapley$new(predictor, x.interest = x_i)
    shap_df <- shap_obj$results

    # Add sample name and feature value
    shap_df$sample <- sample_i
    shap_df$feature.value <- as.numeric(x_i[1, shap_df$feature])

    return(shap_df)
  })

  # Combine all SHAP values into a single data frame
  result <- do.call(rbind, shap_list)
  return(result)
}
