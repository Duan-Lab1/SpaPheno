#' @title Plot Spatial Phenotype Significance Map with Permutation Thresholds
#'
#' @description
#' Visualizes phenotype prediction scores in spatial coordinates and highlights spatial units
#' with statistically significant scores based on permutation-derived null distributions.
#'
#' @details
#' Points with scores above the upper threshold (top \code{1 - p}) are shown in shades of red,
#' and those below the lower threshold (bottom \code{p}) in shades of blue. Non-significant points
#' are shown in gray. A separate side legend visualizes the score gradient and labels for phenotype interpretation.
#'
#' @author Bin Duan
#'
#' @param coordinate A data.frame or matrix with two columns (x, y). Row names should match prediction order.
#' @param prediction A numeric vector of phenotype prediction scores for each spatial unit.
#' @param permutation A numeric vector from the permutation null distribution.
#' @param size Numeric; point size in the spatial plot. Default is 1.
#' @param p Numeric; two-sided significance level (between 0 and 1). Default is 0.01.
#'
#' @return A combined \code{ggplot2} visualization with spatial points and a legend,
#' rendered using \code{gridExtra::grid.arrange}.
#'
#' @importFrom ggplot2 ggplot aes geom_point scale_color_identity theme_minimal
#' @importFrom ggplot2 geom_raster scale_fill_identity annotate theme element_blank
#' @importFrom gridExtra grid.arrange
#'
#' @examples
#' \dontrun{
#' set.seed(123)
#' coords <- data.frame(x = runif(100), y = runif(100))
#' pred <- rnorm(100)
#' perm <- rnorm(1000)
#' PhenotypePlot(coords, pred, perm)
#' }
#'
#' @export
#'
PhenotypePlot <- function(
    coordinate,
    prediction,
    permutation,
    size = 1,
    p = 0.01) {

  # Compute significance thresholds from permutation distribution
  upper_threshold <- quantile(permutation, 1 - p)
  lower_threshold <- quantile(permutation, p)

  # Color assignment function
  assign_color <- function(values, upper, lower) {
    colors <- rep("lightgray", length(values))

    if (any(values > upper)) {
      red_palette <- colorRampPalette(c("pink", "red"))(100)
      colors[values > upper] <- red_palette[cut(
        values[values > upper],
        breaks = 100,
        labels = FALSE
      )]
    }

    if (any(values < lower)) {
      blue_palette <- colorRampPalette(c("blue", "lightblue"))(100)
      colors[values < lower] <- blue_palette[cut(
        values[values < lower],
        breaks = 100,
        labels = FALSE
      )]
    }

    return(colors)
  }

  # Assign colors to prediction values
  prediction_color <- assign_color(prediction, upper_threshold, lower_threshold)

  # Prepare coordinate data
  coordinate <- as.data.frame(coordinate)
  colnames(coordinate) <- c("x", "y")

  # Spatial scatter plot
  scatter_plot <- ggplot2::ggplot(coordinate, ggplot2::aes(x = x, y = y)) +
    ggplot2::geom_point(ggplot2::aes(color = prediction_color), size = size, show.legend = FALSE) +
    ggplot2::scale_color_identity() +
    ggplot2::theme_minimal()

  # Legend (vertical gradient of scores)
  legend_data <- data.frame(
    value = seq(min(prediction), max(prediction), length.out = 300),
    x = 0.5,
    y = seq(0, 5, length.out = 300)
  )

  # Assign fill color for legend based on score range
  legend_data$fill_color <- sapply(legend_data$value, function(val) {
    if (val > upper_threshold) {
      colorRampPalette(c("pink", "red"))(100)[
        cut(val, breaks = seq(upper_threshold, max(prediction), length.out = 101), labels = FALSE)
      ]
    } else if (val < lower_threshold) {
      colorRampPalette(c("blue", "lightblue"))(100)[
        cut(val, breaks = seq(min(prediction), lower_threshold, length.out = 101), labels = FALSE)
      ]
    } else {
      "lightgray"
    }
  })

  # Legend plot
  legend_plot <- ggplot2::ggplot(legend_data, ggplot2::aes(x = x, y = value, fill = fill_color)) +
    ggplot2::geom_raster() +
    ggplot2::scale_fill_identity() +
    ggplot2::annotate("text", x = 1, y = min(prediction), label = round(min(prediction), 2), hjust = 0) +
    ggplot2::annotate("text", x = 1, y = lower_threshold, label = paste0("Low (", round(lower_threshold, 2), ")"), hjust = 0, color = "blue") +
    ggplot2::annotate("text", x = 1, y = upper_threshold, label = paste0("High (", round(upper_threshold, 2), ")"), hjust = 0, color = "red") +
    ggplot2::annotate("text", x = 1, y = max(prediction), label = round(max(prediction), 2), hjust = 0) +
    ggplot2::annotate("text", x = 1.1, y = (upper_threshold + max(prediction)) / 2, label = "Phenotype+", color = "red", fontface = "bold", size = 4, hjust = 0) +
    ggplot2::annotate("text", x = 1.1, y = (upper_threshold + lower_threshold) / 2, label = "Background", color = "gray40", fontface = "bold", size = 4, hjust = 0) +
    ggplot2::annotate("text", x = 1.1, y = (min(prediction) + lower_threshold) / 2, label = "Phenotype-", color = "blue", fontface = "bold", size = 4, hjust = 0) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      axis.title = ggplot2::element_blank(),
      axis.text = ggplot2::element_blank(),
      axis.ticks = ggplot2::element_blank(),
      legend.position = "none",
      plot.title = ggplot2::element_text(size = 14, face = "bold", hjust = 0.5)
    ) +
    ggplot2::labs(title = "Prediction")

  # Combine scatter and legend plots
  gridExtra::grid.arrange(
    scatter_plot + ggplot2::theme(plot.margin = unit(c(1, 1, 1, 1), "lines")),
    legend_plot + ggplot2::theme(plot.margin = unit(c(1, 1, 1, 1), "lines")),
    ncol = 2,
    widths = c(3, 1)
  )
}
