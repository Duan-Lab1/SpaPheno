#' @title Residual-based SHAP Analysis and Visualization in Spatial Transcriptomics
#'
#' @description
#' Performs residual-based SHAP analysis for a selected feature, identifying
#' outliers based on SHAP residuals and visualizing them in both spatial and dependence plots.
#'
#' @details
#' This function regresses SHAP values on feature values to obtain residuals, then uses
#' Z-score normalization to identify spatial units with unusually high or low SHAP effects
#' that are not explained by the raw expression. The results are visualized both spatially
#' and as enhanced SHAP dependence plots.
#'
#' Residual outliers are defined as:
#' \itemize{
#'   \item Z ≥ 2: High residual
#'   \item Z ≤ -2: Low residual
#'   \item Otherwise: Normal
#' }
#'
#' @author Bin Duan
#'
#' @param shap_df A data.frame containing SHAP results with at least the columns:
#'   \code{feature}, \code{phi}, \code{feature.value}, and \code{sample}.
#' @param feature_name Character; the feature (e.g., gene) to analyze.
#' @param coordinate_df A data.frame of spatial coordinates with rownames matching \code{sample} IDs;
#'   must contain columns named \code{X} and \code{Y}.
#' @param size Numeric; point size for plotting. Default is 1.
#' @param title Character; title for the spatial residual plot. Default is \code{"SHAP Residual Distribution"}.
#'
#' @return A list containing:
#' \itemize{
#'   \item \code{spatial_plot}: A \code{ggplot2} object showing spatial distribution of SHAP residuals.
#'   \item \code{dependence_plot}: A SHAP dependence plot with residual outlier groups.
#'   \item \code{residual_table}: A data.frame with residuals, Z-scores, and group annotations.
#' }
#'
#' @importFrom ggplot2 ggplot aes geom_point geom_smooth scale_color_manual labs coord_fixed
#' @importFrom ggplot2 theme_minimal theme element_text ggtitle xlab ylab
#' @importFrom dplyr filter mutate case_when
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
#'
#' resi_result <- SpaPheno_SHAP_residual_analysis(
#' shap_df = shap_test_plus,
#' feature_name = "Perivascular.Macrophages",
#' coordinate_df = test_coordinate, size = 0.8
#' )
#' resi_hot <- resi_result$residual_table
#' head(resi_hot[order(abs(resi_hot$phi_resid_z), decreasing = T), ], 5)
#' SpaPheno_SHAP_waterfall_plot(shap_test_plus, "cell_5593", top_n = 10)
#' resi_result$dependence_plot
#' resi_result$spatial_plot
#'
#' }
#'
#' @export
#'
SpaPheno_SHAP_residual_analysis <- function(
    shap_df,
    feature_name,
    coordinate_df,
    size = 1,
    title = "SHAP Residual Distribution") {

  # 1. Filter SHAP results for the specified feature
  df_feature <- dplyr::filter(shap_df, feature == feature_name)

  # 2. Fit linear model: SHAP ~ feature value
  model <- stats::lm(phi ~ feature.value, data = df_feature)
  df_feature$phi_residual <- stats::residuals(model)

  # 3. Compute Z-scores for residuals and classify into groups
  df_feature <- dplyr::mutate(df_feature,
                              phi_resid_z = scale(phi_residual),
                              resid_group = dplyr::case_when(
                                phi_resid_z >= 2 ~ "High residual",
                                phi_resid_z <= -2 ~ "Low residual",
                                TRUE ~ "Normal"
                              )
  )

  # 4. Assign residual groups to spatial coordinates
  coordinate_df$resid_group <- "Background"
  coordinate_df[df_feature$sample, "resid_group"] <- df_feature$resid_group

  spatial_plot <- ggplot2::ggplot(coordinate_df, ggplot2::aes(X, Y, color = resid_group)) +
    ggplot2::geom_point(size = size) +
    ggplot2::scale_color_manual(values = c(
      "High residual" = "#E41A1C",
      "Low residual" = "#377EB8",
      "Normal" = "#4D4D4D",
      "Background" = "#D9D9D9"
    )) +
    ggplot2::coord_fixed() +
    ggplot2::theme_minimal() +
    ggplot2::labs(title = title)

  # 5. Dependence plot with group coloring
  df_feature$color_group <- ifelse(
    df_feature$resid_group == "High residual", "High",
    ifelse(df_feature$resid_group == "Low residual", "Low", "Normal")
  )

  dependence_plot <- ggplot2::ggplot(df_feature, ggplot2::aes(x = feature.value, y = phi)) +
    ggplot2::geom_point(ggplot2::aes(color = color_group), size = size, alpha = 0.7) +
    ggplot2::geom_smooth(method = "loess", se = TRUE, color = "black") +
    ggplot2::scale_color_manual(values = c(
      "High" = "#E41A1C",
      "Low" = "#377EB8",
      "Normal" = "#4D4D4D"
    )) +
    ggplot2::theme_minimal(base_size = 14) +
    ggplot2::xlab("Feature Value") +
    ggplot2::ylab("SHAP Value") +
    ggplot2::ggtitle(paste0("SHAP Dependence Plot - ", feature_name)) +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"))

  # 6. Return plots and data
  return(list(
    spatial_plot = spatial_plot,
    dependence_plot = dependence_plot,
    residual_table = df_feature
  ))
}
