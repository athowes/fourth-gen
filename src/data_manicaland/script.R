#' Uncomment and run the two line below to resume development of this script
# orderly::orderly_develop_start("data_manicaland")
# setwd("src/data_manicaland")

#' Load raw data
manicaland <- read_csv("depends/manicaland.csv")

#' Processing steps

#' Save processed data
data <- list()
saveRDS(data, "data.rds")
