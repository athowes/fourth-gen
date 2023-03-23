#' Uncomment and run the two line below to resume development of this script
# orderly::orderly_develop_start("run_ea-bsd")
# setwd("src/run_ea-bsd")

#' Two-stage cluster design
#' 1. BSD of EAs with posterior variance acquisition function and distance costs
#' 2. Systematic sampling of households

data <- readRDS("depends/data.rds")

#' The number of simulated DHS surveys to run
N_sim <- 1

results <- list()
saveRDS(results, "results.rds")
