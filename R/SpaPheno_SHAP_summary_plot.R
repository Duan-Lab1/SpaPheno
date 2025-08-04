#' @title SHAP Summary Plot
#'
#' @description
#' Generate a SHAP summary plot to visualize the distribution, direction, and magnitude of SHAP values
#' for the top contributing features in a spatial transcriptomics model.
#'
#' @details
#' The function selects the top `n` features ranked by mean absolute SHAP value, and displays
#' a scatter plot of their SHAP contributions. Each point is colored by the corresponding
#' feature value, allowing joint interpretation of importance and directionality.
#'
#' @author Bin Duan
#'
#' @param shap_df A data.frame containing SHAP values. Must include at least the columns:
#'   \code{feature}, \code{phi} (SHAP value), and \code{feature.value}.
#' @param top_n Integer; number of top features to include in the plot. Default is 10.
#'
#' @return A \code{ggplot2} object displaying a SHAP summary plot for top features.
#'
#' @importFrom ggplot2 ggplot aes geom_jitter labs theme_minimal
#' @importFrom dplyr group_by summarise arrange slice_head pull filter %>%
#' @importFrom viridis scale_color_viridis
#'
#' @examples
#' \dontrun{
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
#' SpaPheno_SHAP_summary_plot(shap_test_plus, top_n = 31)
#'
#' }
#'
#' @export
#'
SpaPheno_SHAP_summary_plot <- function(shap_df, top_n = 10) {
  # Select top N features based on mean absolute SHAP values
  top_features <- dplyr::group_by(shap_df, feature) %>%
    dplyr::summarise(mean_abs_shap = mean(abs(phi)), .groups = "drop") %>%
    dplyr::arrange(desc(mean_abs_shap)) %>%
    dplyr::slice_head(n = top_n) %>%
    dplyr::pull(feature)

  plot_df <- dplyr::filter(shap_df, feature %in% top_features)

  ggplot2::ggplot(plot_df, ggplot2::aes(
    x = phi,
    y = reorder(feature, phi, FUN = median),
    color = feature.value
  )) +
    ggplot2::geom_jitter(height = 0.2, size = 1.2, alpha = 0.6) +
    viridis::scale_color_viridis(option = "plasma", direction = -1) +
    ggplot2::labs(
      x = "SHAP value",
      y = "Feature",
      color = "Feature value",
      title = "SHAP Summary Plot"
    ) +
    ggplot2::theme_minimal()
}
