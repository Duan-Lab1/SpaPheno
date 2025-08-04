#' @title Plot Feature Coefficients from a glmnet Model
#'
#' @description
#' Visualizes feature coefficients from a penalized regression model fitted by \code{glmnet}.
#' Positive and negative coefficients are color-coded for better interpretation.
#'
#' @details
#' The function extracts non-zero coefficients at the optimal lambda from a glmnet model object.
#' It removes the intercept term and displays a horizontal bar plot using ggplot2.
#' Coefficients are sorted to enhance readability.
#'
#' @author Bin Duan
#'
#' @param model A fitted \code{glmnet} model object that includes components \code{model} and \code{lambda},
#'   such as from the output of \code{BuildPhenoModelAutoAlpha()} or \code{cv.glmnet()}.
#'
#' @return A \code{ggplot2} object showing feature coefficients with direction-based color coding.
#'
#' @importFrom ggplot2 ggplot aes geom_col labs scale_fill_manual geom_vline theme_minimal
#' @importFrom ggplot2 theme element_text element_blank margin
#' @importFrom dplyr filter mutate
#'
#' @examples
#' \dontrun{
#' library(glmnet)
#' x <- matrix(rnorm(100 * 20), nrow = 100, ncol = 20)
#' y <- sample(c(0, 1), 100, replace = TRUE)
#' fit <- cv.glmnet(x, y, family = "binomial")
#' FeatureImportanceCoef(model = model)
#' }
#'
#' @export
#'
FeatureImportanceCoef <- function(model) {

  # Extract coefficients at selected lambda
  model_coef <- as.matrix(
    stats::coef(
      object = model$model,
      s = model$lambda
      )
    )

  # Assign column name and remove intercept
  colnames(model_coef) <- "Coef"
  model_coef <- model_coef[!(row.names(model_coef) %in% "(Intercept)"), , drop = FALSE]

  # Build data frame with coefficient values
  data <- data.frame(
    feature = rownames(model_coef),
    coefficient = model_coef[, 1],
    stringsAsFactors = FALSE
  )

  # Clean and format for plotting
  data <- dplyr::filter(data, !is.na(coefficient))
  data <- dplyr::mutate(data,
                        direction = ifelse(coefficient > 0, "Positive", "Negative"),
                        feature = factor(feature, levels = feature[order(coefficient)])
  )

  # Generate horizontal bar plot
  p <- ggplot2::ggplot(data, ggplot2::aes(x = coefficient, y = feature, fill = direction)) +
    ggplot2::geom_col(width = 0.8) +
    ggplot2::scale_fill_manual(values = c(Negative = "#377EB8", Positive = "#E41A1C")) +
    ggplot2::geom_vline(xintercept = 0, color = "gray50", linewidth = 0.5) +
    ggplot2::labs(
      x = "Coefficient Value",
      y = NULL,
      title = "Feature Coefficients Visualization",
      fill = "Direction"
    ) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      axis.text.y = ggplot2::element_text(size = 10, margin = ggplot2::margin(r = 15)),
      plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"),
      legend.position = "top",
      panel.grid.major.y = ggplot2::element_blank()
    )

  return(p)
}
