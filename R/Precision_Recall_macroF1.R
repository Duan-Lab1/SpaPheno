#' @title Compute Macro-Averaged Precision, Recall, and F1 Score for Multiclass Classification
#'
#' @description
#' Computes macro-averaged precision, recall, and F1 score using a confusion matrix
#' between true labels and predicted labels. Suitable for multi-class classification tasks.
#'
#' @details
#' Macro-averaging computes the mean of the metric (precision, recall, F1) over all classes,
#' treating all classes equally regardless of their support. Missing values due to undefined metrics
#' are excluded from the averaging (e.g., divisions by zero).
#'
#' @author Bin Duan
#'
#' @param actual A vector of true class labels.
#' @param predicted A vector of predicted class labels (must have same length as \code{actual}).
#'
#' @return A named list containing:
#' \itemize{
#'   \item \code{Precision}: Macro-averaged precision across all classes.
#'   \item \code{Recall}: Macro-averaged recall across all classes.
#'   \item \code{f1_score}: Macro-averaged F1 score.
#' }
#'
#' @importFrom caret confusionMatrix
#'
#' @examples
#' \dontrun{
#' actual <- c("A", "B", "A", "C", "B", "C")
#' predicted <- c("A", "B", "C", "C", "B", "A")
#' Precision_Recall_macroF1(actual, predicted)
#' }
#'
#' @export
#'
Precision_Recall_macroF1 <- function(
    actual,
    predicted) {

  # Compute confusion matrix using caret
  cm <- caret::confusionMatrix(
    data = factor(predicted),
    reference = factor(actual)
  )

  # Extract per-class precision (Positive Predictive Value) and recall (Sensitivity)
  precision <- cm$byClass[, "Pos Pred Value"]
  recall <- cm$byClass[, "Sensitivity"]

  # Compute per-class F1 scores
  f1_scores <- 2 * (precision * recall) / (precision + recall)

  # Compute macro-averaged metrics
  macro_precision <- mean(precision, na.rm = TRUE)
  macro_recall <- mean(recall, na.rm = TRUE)
  macro_f1 <- mean(f1_scores, na.rm = TRUE)

  # Return result as named list
  return(list(
    Precision = macro_precision,
    Recall = macro_recall,
    f1_score = macro_f1
  ))
}
