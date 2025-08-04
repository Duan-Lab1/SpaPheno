#' @title Visualize Average Cell Type Abundance by Phenotype Label
#'
#' @description
#' This function computes the average abundance of each cell type stratified by phenotype labels
#' and visualizes the result using a line plot. Labels are based on prediction scores or user-provided values.
#'
#' @details
#' This plot is useful for interpreting how cell type composition differs across predicted phenotypic groups.
#' Input data should contain cell type proportions and prediction scores with associated phenotype labels.
#'
#' @author Bin Duan
#'
#' @param cell_type_distribution A data.frame or matrix of cell type proportions per spot (rows = spots, columns = cell types).
#' @param pred_score A data.frame with prediction scores. Must contain at least column \code{"pred_score"} and optional \code{"label"}.
#'                   Row names should match those of \code{cell_type_distribution}.
#'
#' @return A \code{ggplot2} object showing average cell type abundance across phenotype groups.
#'
#' @importFrom dplyr inner_join group_by summarise
#' @importFrom tidyr pivot_longer
#' @importFrom ggplot2 ggplot aes geom_point geom_line scale_color_manual labs theme_bw theme element_text
#'
#' @examples
#' \dontrun{
#' cell_type_distribution <- matrix(runif(300), nrow = 100,
#'   dimnames = list(paste0("spot", 1:100), paste0("cell", 1:3)))
#' pred_score <- data.frame(
#'   pred_score = rnorm(100),
#'   label = sample(c("phenotype+", "phenotype-", "background"), 100, replace = TRUE),
#'   row.names = paste0("spot", 1:100)
#' )
#' PhenoAbundancePlot(cell_type_distribution, pred_score)
#' }
#'
#' @export
#'
PhenoAbundancePlot <- function(cell_type_distribution, pred_score) {

  # Convert inputs to data.frames
  cell_type_distribution <- as.data.frame(cell_type_distribution)
  pred_score <- as.data.frame(pred_score)

  # Add spot IDs for merging
  pred_score$spot_id <- rownames(pred_score)
  cell_type_distribution$spot_id <- rownames(cell_type_distribution)

  # Merge by spot ID
  merged_data <- dplyr::inner_join(pred_score, cell_type_distribution, by = "spot_id")

  # Identify cell type columns by exclusion
  cell_type_cols <- setdiff(colnames(merged_data), c("spot_id", "pred_score", "label"))

  # Reshape to long format for plotting
  long_data <- tidyr::pivot_longer(
    data = merged_data,
    cols = dplyr::all_of(cell_type_cols),
    names_to = "cell_type",
    values_to = "abundance"
  )

  # Compute mean abundance per phenotype and cell type
  mean_abundance <- long_data |>
    dplyr::group_by(label, cell_type) |>
    dplyr::summarise(
      mean_abundance = mean(abundance, na.rm = TRUE),
      .groups = "drop"
    )

  # Plot with ggplot2
  p <- ggplot2::ggplot(mean_abundance, ggplot2::aes(
    x = cell_type,
    y = mean_abundance,
    group = label,
    color = label
  )) +
    ggplot2::geom_point(size = 3) +
    ggplot2::geom_line(linewidth = 1) +
    ggplot2::scale_color_manual(
      values = c("phenotype+" = "#F8766D", "phenotype-" = "lightblue", "background" = "lightgray")
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)
    ) +
    ggplot2::labs(
      title = "Average Cell Type Abundance by Phenotype Label",
      x = "Cell Type",
      y = "Average Abundance"
    )

  return(p)
}
