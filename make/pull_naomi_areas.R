#' Let's work with Zimbabwe, as it contains the Manicaland province
iso3 <- "ZWE"

#' Name of the report to pull
report <- paste0(tolower(iso3), "_data_areas")
orderly::orderly_pull_archive(report, remote = "naomi2")
