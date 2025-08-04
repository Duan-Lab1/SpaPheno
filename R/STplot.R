#' @title Spatial Transcriptomics Scatter Plot with Labels
#'
#' @description
#' This function creates a scatter plot of spatial coordinates colored by labels.
#' It is useful for visualizing spatial transcriptomics or spatial data with sample annotations.
#'
#' @author Bin Duan
#'
#' @param coordinate A data frame containing spatial coordinates. Must have columns named 'X' and 'Y'.
#'                   Row names should correspond to spots or cells.
#' @param label A named vector or factor of labels corresponding to the rows of `coordinate`.
#'              The function will reorder the labels to match the coordinate row order.
#' @param size Numeric. The size of the points in the plot. Default is 1.
#'
#' @return A ggplot2 scatter plot object showing spatial distribution colored by label.
#'
#' @import ggplot2
#' @importFrom viridis scale_color_viridis
#'
#' @examples
#' \dontrun{
#' coord_df <- data.frame(X = runif(100), Y = runif(100))
#' rownames(coord_df) <- paste0("spot", 1:100)
#' labels <- sample(c("TypeA", "TypeB"), 100, replace = TRUE)
#' names(labels) <- rownames(coord_df)
#' STplot(coord_df, labels)
#' }
#'
#' @export
#'
STplot <- function(coordinate, label, size = 1) {

  # Reorder labels to match coordinate row names
  label <- label[row.names(coordinate)]

  # Generate scatter plot with viridis color scale
  pl <- ggplot(coordinate, aes(x = X, y = Y, color = label)) +
    geom_point(size = size) +
    viridis::scale_color_viridis(option = "D", direction = 1, discrete = TRUE) +
    theme_minimal()

  return(pl)
}
