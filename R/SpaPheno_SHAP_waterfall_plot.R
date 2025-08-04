#' @title SHAP Waterfall Plot for a Single Sample
#'
#' @description
#' Generate a SHAP waterfall plot to visualize the most influential features
#' contributing to a model's prediction for a given sample. Features are ordered
#' by the magnitude of their SHAP values, and color-coded by direction (positive/negative).
#'
#' @details
#' - Filters SHAP values for the selected \code{sample_id}.
#' - Selects top \code{top_n} features based on absolute SHAP value.
#' - Plots bars horizontally using \code{ggplot2}, colored by contribution sign.
#'
#' @author Bin Duan
#'
#' @param shap_df A data.frame of SHAP results. Must include:
#'   \code{sample}, \code{feature}, and \code{phi} columns.
#' @param sample_id A string. The ID of the sample to visualize (must match \code{sample} in \code{shap_df}).
#' @param top_n Integer; number of top features to include. Default is 10.
#'
#' @return A \code{ggplot2} object displaying the SHAP waterfall plot.
#'
#' @importFrom dplyr filter arrange slice_head mutate
#' @importFrom ggplot2 ggplot aes geom_bar coord_flip scale_fill_manual labs theme_minimal
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
#' SpaPheno_SHAP_waterfall_plot(shap_test_plus, sample_id = "cell_3477", top_n = 10)
#'
#' }
#'
#' @export
#'
SpaPheno_SHAP_waterfall_plot <- function(shap_df, sample_id, top_n = 10) {

  df <- shap_df %>%
    dplyr::filter(sample == sample_id) %>%
    dplyr::arrange(desc(abs(phi))) %>%
    dplyr::slice_head(n = top_n) %>%
    dplyr::mutate(
      direction = ifelse(phi >= 0, "Positive", "Negative"),
      feature = reorder(feature, phi)
    )

  p <- ggplot2::ggplot(df, ggplot2::aes(x = feature, y = phi, fill = direction)) +
    ggplot2::geom_bar(stat = "identity") +
    ggplot2::coord_flip() +
    ggplot2::scale_fill_manual(values = c("Positive" = "#E41A1C", "Negative" = "#377EB8")) +
    ggplot2::labs(
      title = paste("SHAP Waterfall for", sample_id),
      y = "SHAP value", x = "Top Features"
    ) +
    ggplot2::theme_minimal()

  return(p)
}
