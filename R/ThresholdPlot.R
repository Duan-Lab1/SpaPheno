#' @title Plot Prediction vs. Permutation Distributions with Significance Thresholds
#'
#' @description
#' Visualizes the density distribution of prediction scores in comparison to a null distribution
#' (e.g., from permutation), and highlights statistically significant regions based on a two-tailed
#' threshold. This is useful for visually identifying predictions that deviate meaningfully
#' from chance expectation.
#'
#' @author Bin Duan
#'
#' @param prediction A numeric vector of predicted scores.
#' @param permutation A numeric vector representing the null distribution (e.g., scores from permutations).
#' @param p Numeric, the significance level for defining extreme regions (two-tailed). Default is 0.01.
#'
#' @return No return value. This function generates a base R plot showing the density curves
#' of the prediction and permutation distributions, with shaded areas representing values
#' beyond the significance thresholds.
#'
#' @examples
#' set.seed(123)
#' pred <- rnorm(1000, mean = 0.5)
#' perm <- rnorm(1000)
#' ThresholdPlot(prediction = pred, permutation = perm, p = 0.01)
#'
#' @export
#'
ThresholdPlot <- function(
    prediction,
    permutation,
    p = 0.01) {
  x <- prediction
  y <- permutation

  high_thre <- quantile(y, 1 - p)
  low_thre <- quantile(y, p)

  density_x <- density(x)

  plot(density_x, col = NA, xlab = "Value", main = "Density Plot of Prediction and Permutation", lwd = 2)

  # Shaded significant regions
  polygon(c(density_x$x[density_x$x < low_thre], low_thre),
          c(density_x$y[density_x$x < low_thre], 0),
          col = "gray", border = NA
  )

  polygon(c(density_x$x[density_x$x > high_thre], high_thre),
          c(density_x$y[density_x$x > high_thre], 0),
          col = "gray", border = NA
  )

  # Non-significant region (white background)
  polygon(density_x$x[density_x$x >= low_thre & density_x$x <= high_thre],
          density_x$y[density_x$x >= low_thre & density_x$x <= high_thre],
          col = "white", border = NA
  )

  # Overlay lines
  lines(density(y), col = "red", lwd = 2)
  lines(density(x), col = "blue", lwd = 2)

  # Significance thresholds
  abline(v = c(low_thre, high_thre), col = "red", lwd = 2, lty = 2)

  # Legend
  legend("topright",
         legend = c("Prediction (Significant)", "Prediction", "Permutation"),
         fill = c("gray", "blue", "red"),
         border = NA
  )
}
