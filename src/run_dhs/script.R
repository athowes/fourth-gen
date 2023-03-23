#' Uncomment and run the two line below to resume development of this script
# orderly::orderly_develop_start("run_dhs")
# setwd("src/run_dhs")

#' Two-stage cluster design
#' 1. Probability-proportional-to-size sampling of EAs stratified within each district by urban and rural status
#' 2. Systematic sampling of households

data <- readRDS("depends/data.rds")

#' The number of simulated DHS surveys to run
N_sim <- 1

results <- list()
saveRDS(results, "results.rds")
