#' @title Generate Simulated Pseudo-Bulk Data for Two Phenotypes
#'
#' @description
#' Simulates pseudo-bulk data from single-cell or spatial transcriptomics input
#' by aggregating data across regions/phenotypes, with optional perturbation.
#'
#' @details
#' The function supports two modes:
#' \itemize{
#'   \item \code{"expression"}: Uses numeric expression matrix (e.g., genes × spots).
#'   \item \code{"proportion"}: Uses categorical input (e.g., cell type labels) and converts to proportions.
#' }
#' For each phenotype group, the function randomly perturbs the data and averages columns to simulate
#' multiple pseudo-bulk samples. This is useful for benchmarking or downstream analysis.
#'
#' @author Bin Duan
#'
#' @param input_data A matrix or data.frame. For \code{mode = "expression"}, it should be a numeric matrix with rows as genes and columns as spots. For \code{"proportion"}, it can be a vector of cell types.
#' @param region_labels A named vector with phenotype labels. Names must match column names (spots) of \code{input_data}.
#' @param phenotypes A character vector of two phenotype labels to compare, e.g., \code{c("A", "B")}.
#' @param perturbation_percent Numeric; percentage (between 0 and 1) of random noise to add. Default is 0.1.
#' @param num_samples Integer; number of pseudo-bulk samples to generate for each phenotype. Default is 50.
#' @param mode Character; either \code{"expression"} or \code{"proportion"}. Default is \code{"proportion"}.
#' @param seed Integer; random seed for reproducibility. Default is 123.
#'
#' @return A list with two data.frames:
#' \itemize{
#'   \item First element: pseudo-bulk samples for phenotype I
#'   \item Second element: pseudo-bulk samples for phenotype II
#' }
#'
#' @examples
#' \dontrun{
#' # Simulate from gene expression matrix
#' gene_mat <- matrix(runif(1000), nrow = 100, ncol = 10)
#' colnames(gene_mat) <- paste0("Spot", 1:10)
#' labels <- setNames(rep(c("A", "B"), each = 5), colnames(gene_mat))
#' result <- generate_simulated_bulk_data(
#'   input_data = gene_mat,
#'   region_labels = labels,
#'   phenotypes = c("A", "B"),
#'   mode = "expression"
#' )
#'
#' data("osmFISH_metadata_cellType")
#' data("osmFISH_metadata_region")
#' data("osmFISH_phenotype_simu")
#' pseudo_bulk_simi <- generate_simulated_bulk_data(
#'   input_data = osmFISH_metadata_cellType,
#'   region_labels = osmFISH_metadata_region,
#'   phenotypes = osmFISH_phenotype_simu,
#'   perturbation_percent = 0.1,
#'   num_samples = 50,
#'   mode = "proportion")
#'
#' head(pseudo_bulk_simi[[1]])
#' }
#'
#' @export
#'
generate_simulated_bulk_data <- function(
    input_data,
    region_labels,
    phenotypes,
    perturbation_percent = 0.1,
    num_samples = 50,
    mode = c("proportion", "expression"),
    seed = 123) {

  mode <- match.arg(mode)

  phenotype_I <- phenotypes[1]
  phenotype_II <- phenotypes[2]

  # Select column names for each phenotype group
  idx_I <- names(region_labels)[!(region_labels %in% phenotype_II)]
  idx_II <- names(region_labels)[!(region_labels %in% phenotype_I)]

  # If input is a vector and in proportion mode, compute table-based frequency matrix
  if (mode == "proportion" && is.null(dim(input_data))) {
    data_I_p <- prop.table(table(input_data[idx_I]))
    data_II_p <- prop.table(table(input_data[idx_II]))
    data_I <- matrix(data_I_p, ncol = 1)
    rownames(data_I) <- names(data_I_p)
    data_II <- matrix(data_II_p, ncol = 1)
    rownames(data_II) <- names(data_II_p)
  } else {
    # Subset expression matrix columns based on phenotype labels
    data_I <- input_data[, idx_I, drop = FALSE]
    data_II <- input_data[, idx_II, drop = FALSE]
  }

  # Function to generate simulated pseudo-bulk samples
  generate_bulk <- function(mat, phenotype_label) {
    n_row <- nrow(mat)
    pseudo_bulk <- matrix(nrow = n_row, ncol = num_samples)

    for (i in seq_len(num_samples)) {
      # Generate perturbation matrix
      perturb <- 1 + runif(n_row * ncol(mat), -perturbation_percent, perturbation_percent)
      perturbed <- mat * matrix(perturb, nrow = n_row, ncol = ncol(mat))
      perturbed[perturbed < 0] <- 0  # Clamp negative values to zero
      pseudo_bulk[, i] <- rowMeans(perturbed)
    }

    colnames(pseudo_bulk) <- paste0("PseudoBulk_", phenotype_label, "_", seq_len(num_samples))
    rownames(pseudo_bulk) <- rownames(mat)
    return(as.data.frame(pseudo_bulk))
  }

  # Set seed for reproducibility
  set.seed(seed)

  # Generate pseudo-bulk samples for both phenotypes
  pseudo_bulk_I <- generate_bulk(as.matrix(data_I), "PhenotypeI")
  pseudo_bulk_II <- generate_bulk(as.matrix(data_II), "PhenotypeII")

  return(list(pseudo_bulk_I, pseudo_bulk_II))
}
