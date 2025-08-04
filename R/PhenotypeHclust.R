#' @title Phenotype-based Hierarchical Clustering and Visualization of Cell Type Distributions
#'
#' @description
#' This function performs phenotype-stratified hierarchical clustering of cell type proportions
#' using Jensen-Shannon Divergence (JSD) as the distance metric. It then visualizes average
#' cell type distributions across identified clusters/domains.
#'
#' @details
#' The function filters spots/cells by a specified phenotype label, computes pairwise JSD distances,
#' applies hierarchical clustering, and visualizes the mean cell type profiles of each domain.
#' Requires the \code{SpaDo} package for clustering and plotting domains. devtools::install_github("bm2-lab/SpaDo").
#'
#' @author Bin Duan
#'
#' @param PhenotypeMap_result A list with at least:
#'   \itemize{
#'     \item \code{cell_type_distribution}: Matrix or data.frame of cell type proportions (rows = spots).
#'     \item \code{pred_score}: Data.frame with phenotype labels. Row names must match the distribution matrix.
#'   }
#' @param phenotype Character string indicating phenotype label to subset. Default is \code{"phenotype-"}.
#' @param coordinate A data.frame of spatial coordinates for each spot/cell. Used for spatial plotting.
#' @param k Integer; number of clusters/domains to define. Default is 2.
#' @param size Numeric; point size used in \code{SpaDo::DomainPlot}. Default is 2.
#'
#' @return A \code{ggplot2} object displaying mean cell type proportions per identified domain.
#'
#' @importFrom SpaDo DistributionDistance DomainHclust DomainPlot
#' @importFrom reshape2 melt
#' @importFrom ggplot2 ggplot aes geom_line geom_point theme element_text
#'
#' @examples
#' \dontrun{
#' res <- list(
#'   cell_type_distribution = cell_type_df,
#'   pred_score = phenotype_df
#' )
#' PhenotypeHclust(res, phenotype = "phenotype+", coordinate = coord_df, k = 3)
#' }
#'
#' @export
#'
PhenotypeHclust <- function(
    PhenotypeMap_result,
    phenotype = c("phenotype-", "phenotype+"),
    coordinate,
    k = 2,
    size = 2) {

  phenotype <- match.arg(phenotype)

  # Extract input components
  cell_type_distribution <- PhenotypeMap_result$cell_type_distribution
  label_df <- PhenotypeMap_result$pred_score

  # Filter spots/cells by selected phenotype label
  spot_ids <- rownames(label_df[label_df$label == phenotype, ])
  filtered_distribution <- cell_type_distribution[spot_ids, , drop = FALSE]

  # Compute Jensen-Shannon Divergence distance matrix
  distance_matrix <- SpaDo::DistributionDistance(
    filtered_distribution,
    distance = "JSD",
    no_cores = 1
  )

  # Perform hierarchical clustering (fixed number of clusters)
  domain_hclust <- SpaDo::DomainHclust(
    distribution_distance = distance_matrix,
    autoselection = FALSE,
    domain_num = k
  )

  # Plot spatial domains
  SpaDo::DomainPlot(
    domain_hclust = domain_hclust,
    distribution_distance = distance_matrix,
    sample_information_coordinate = coordinate,
    size = size
  )

  # Assign clusters to cells/spots
  cluster_assignments <- domain_hclust[[1]][, ncol(domain_hclust[[1]])]
  names(cluster_assignments) <- rownames(domain_hclust[[1]])

  # Compute mean cell type proportions per cluster
  cell_means <- lapply(unique(cluster_assignments), function(clust_id) {
    sample_ids <- names(cluster_assignments[cluster_assignments == clust_id])
    cluster_mat <- filtered_distribution[sample_ids, , drop = FALSE]
    colMeans(cluster_mat, na.rm = TRUE)
  })

  # Combine into one data.frame
  df <- as.data.frame(do.call(cbind, cell_means))
  colnames(df) <- paste0("Domain_", unique(cluster_assignments))
  df$Cell_type <- rownames(df)

  # Reshape for plotting
  df_long <- reshape2::melt(df, id.vars = "Cell_type")
  colnames(df_long) <- c("Cell_type", "Domain", "Proportion")

  # Generate plot
  p <- ggplot2::ggplot(df_long, ggplot2::aes(
    x = Cell_type,
    y = Proportion,
    group = Domain,
    color = Domain
  )) +
    ggplot2::geom_line() +
    ggplot2::geom_point(size = 3) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, vjust = 0.5, hjust = 0.5),
      aspect.ratio = 1
    )

  return(p)
}
