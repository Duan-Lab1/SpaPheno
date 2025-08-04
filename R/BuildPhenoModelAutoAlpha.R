#' @title Build Elastic Net Model with Automatic Alpha Selection
#'
#' @description
#' This function builds an Elastic Net model using \code{glmnet}, with automatic
#' alpha selection via cross-validation. It supports generalized linear models
#' for survival analysis (Cox), binary classification, and regression tasks.
#'
#' @details
#' The function performs K-fold cross-validation across a range of alpha values
#' to identify the optimal alpha that minimizes the cross-validation error.
#' Once the best alpha is selected, the model is re-trained using that alpha.
#'
#' Applicable families:
#' \itemize{
#'   \item \code{"cox"}: Cox proportional hazards model (survival)
#'   \item \code{"binomial"}: Logistic regression (binary classification)
#'   \item \code{"gaussian"}: Linear regression (continuous outcomes)
#' }
#'
#' @author Bin Duan
#'
#' @param expr A numeric matrix or data.frame of gene expression or cell-type proportion data (samples × features).
#' @param pheno A response vector of phenotypes: survival object, binary labels, or continuous outcomes.
#' @param family Character string indicating the type of outcome. Must be one of \code{"cox"}, \code{"binomial"}, or \code{"gaussian"}.
#' @param nfolds Number of folds to use in cross-validation. Default is 10.
#' @param seed Random seed for reproducibility. Default is 2025.
#'
#' @return A list with the following components:
#' \describe{
#'   \item{\code{model}}{The fitted \code{cv.glmnet} object with best alpha.}
#'   \item{\code{lambda}}{The lambda.min value corresponding to the best alpha.}
#'   \item{\code{alpha}}{The optimal alpha value selected.}
#' }
#'
#' @importFrom glmnet cv.glmnet
#'
#' @examples
#' \dontrun{
#' # Binary phenotype example
#' data("osmFISH_bulk_decon")
#' data("osmFISH_bulk_pheno")
#' model <- BuildPhenoModelAutoAlpha(
#'   expr = osmFISH_bulk_decon,
#'   pheno = osmFISH_bulk_pheno,
#'   family = "binomial")
#'
#' # Survival example
#' library(survival)
#' surv_obj <- Surv(time = c(5, 10, 15), event = c(1, 0, 1))
#' model <- BuildPhenoModelAutoAlpha(
#'   expr = matrix(rnorm(30), nrow = 3),
#'   pheno = surv_obj,
#'   family = "cox")
#' }
#'
#' @export
#'
BuildPhenoModelAutoAlpha <- function(
    expr,
    pheno,
    family = c("cox", "binomial", "gaussian"),
    nfolds = 10,
    seed = 2025) {

  family <- match.arg(family)
  alpha_list <- c(
    0.005, 0.01, 0.05, 0.1,
    0.2, 0.3, 0.4, 0.5, 0.6,
    0.7, 0.8, 0.9)

  best_cvm <- Inf
  best_alpha <- NULL
  best_fit <- NULL
  set.seed(seed)

  for (a in alpha_list) {
    fit <- cv.glmnet(x = as.matrix(expr), y = pheno, family = family,
                     alpha = a, nfolds = nfolds)
    this_cvm <- fit$cvm[fit$lambda == fit$lambda.min]
    if (this_cvm < best_cvm) {
      best_cvm <- this_cvm
      best_alpha <- a
      best_fit <- fit
    }
  }

  message(sprintf("Optimal alpha selected: %.3f, with lambda.min: %.5f",
                  best_alpha, best_fit$lambda.min))

  return(list(
    model = best_fit,
    lambda = best_fit$lambda.min,
    alpha = best_alpha
  ))
}
