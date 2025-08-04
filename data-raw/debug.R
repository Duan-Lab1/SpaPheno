# build package
roxygen2::roxygenise(getwd())
devtools::check(document = FALSE)
devtools::build(binary = FALSE, manual = TRUE, quiet = FALSE)

# build pkgdown
pkgdown::build_site() # Run to build the website

library(usethis)
library(SpaPheno)

load("../../data/raw/Simulation_osmFISH.RData")
region_all <- names(table(sample_information_region))
phenotype_simu <- region_all[c(3, 4)]
sample_information_region_choose <- sample_information_region
sample_information_region_choose[!sample_information_region_choose %in% phenotype_simu] <- NA

custom_colors <- c("red", "lightgray", "blue")
names(custom_colors) <- c(phenotype_simu[2], "Background", phenotype_simu[1])
sample_information_region_choose[!sample_information_region_choose %in% phenotype_simu] <- "Background"
Ground_truth <- factor(
  sample_information_region_choose[row.names(test_coordinate)],
  levels = c(phenotype_simu[1], "Background", phenotype_simu[2]))

osmFISH_metadata_cellType <- sample_information_cellType
osmFISH_metadata_region <- sample_information_region
osmFISH_phenotype_simu <- phenotype_simu
usethis::use_data(osmFISH_metadata_cellType, overwrite = TRUE)
usethis::use_data(osmFISH_metadata_region, overwrite = TRUE)
usethis::use_data(osmFISH_phenotype_simu, overwrite = TRUE)

pseudo_bulk_simi <- generate_simulated_bulk_data(
  input_data = sample_information_cellType,
  region_labels = sample_information_region,
  phenotypes = phenotype_simu,
  perturbation_percent = 0.1,
  num_samples = 50,
  mode = "proportion")

pseudo_bulk_df1 <- pseudo_bulk_simi[[1]]
pseudo_bulk_df2 <- pseudo_bulk_simi[[2]]

bulk_decon <- t(as.matrix(cbind(pseudo_bulk_df1, pseudo_bulk_df2)))
bulk_pheno <- rep(c(0, 1), each = 50)
names(bulk_pheno) <- c(colnames(pseudo_bulk_df1), colnames(pseudo_bulk_df2))
family <- "binomial"

osmFISH_bulk_decon <- bulk_decon
osmFISH_bulk_pheno <- bulk_pheno
osmFISH_bulk_coordinate <- test_coordinate
usethis::use_data(osmFISH_bulk_decon, overwrite = TRUE)
usethis::use_data(osmFISH_bulk_pheno, overwrite = TRUE)
usethis::use_data(osmFISH_bulk_coordinate, overwrite = TRUE)
