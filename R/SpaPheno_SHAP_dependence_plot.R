#' @title SHAP Dependence Plot for Selected Features
#'
#' @description
#' Visualizes the relationship between SHAP values and corresponding feature values
#' across samples, stratified by selected features.
#'
#' @details
#' This function generates a scatter plot of SHAP values versus feature values,
#' overlaid with a LOESS smoothed curve. Each feature specified in \code{feature_names}
#' is plotted in a separate facet panel. Useful for interpreting the directionality
#' and interaction effects of SHAP-based model explanations.
#'
#' @author Bin Duan
#'
#' @param shap_df A data.frame with SHAP output. Must contain columns:
#'   \itemize{
#'     \item \code{feature}: Name of the feature.
#'     \item \code{feature.value}: The actual value of the feature in a sample.
#'     \item \code{phi}: The SHAP value associated with that feature/sample.
#'   }
#' @param feature_names Character vector of feature names to include in the plot.
#' @param title Optional character string for the plot title. Default is \code{"SHAP Dependence Plot"}.
#'
#' @return A \code{ggplot2} object showing SHAP dependence plots for selected features.
#'
#' @importFrom ggplot2 ggplot aes geom_point geom_smooth facet_wrap
#' @importFrom ggplot2 theme_minimal xlab ylab ggtitle theme element_text
#' @importFrom dplyr filter
#'
#' @examples
#' \dontrun{
#' shap_df <- data.frame(
#'   feature = rep(c("GeneA", "GeneB"), each = 50),
#'   feature.value = rnorm(100),
#'   phi = rnorm(100)
#' )
#' SpaPheno_SHAP_dependence_plot(shap_df, c("GeneA", "GeneB"), title = "Example SHAP Plot")
#' }
#'
#' @export
#'
SpaPheno_SHAP_dependence_plot <- function(shap_df, feature_names, title = NULL) {

  # Filter for selected features
  plot_df <- dplyr::filter(shap_df, feature %in% feature_names)

  # Set default title
  if (is.null(title)) {
    title <- "SHAP Dependence Plot"
  }

  # Create ggplot
  p <- ggplot2::ggplot(plot_df, ggplot2::aes(x = feature.value, y = phi)) +
    ggplot2::geom_point(alpha = 0.5, size = 1) +
    ggplot2::geom_smooth(method = "loess", se = TRUE, color = "firebrick") +
    ggplot2::facet_wrap(~ feature, scales = "free_x") +
    ggplot2::theme_minimal(base_size = 14) +
    ggplot2::xlab("Feature Value") +
    ggplot2::ylab("SHAP Value") +
    ggplot2::ggtitle(title) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5, face = "bold")
    )

  return(p)
}
