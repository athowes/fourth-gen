#' Uncomment and run the two line below to resume development of this script
# orderly::orderly_develop_start("run_hh-bsd")
# setwd("src/run_hh-bsd")

#' BSD at the household level with posterior variance acquisition function and distance costs

data <- readRDS("depends/data.rds")

#' The number of simulated DHS surveys to run
N_sim <- 1

results <- list()
saveRDS(results, "results.rds")
