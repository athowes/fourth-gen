#' Set the working directory to the project root
setwd(rprojroot::find_rstudio_root_file())

#' @param i Index of the artefacts to get
archive_to_docs <- function(report, i = 1) {
  #' Artefacts to be moved
  filenames <- yaml::read_yaml(file.path(paste0("src/", report, "/orderly.yml")))$artefacts[[i]]$data$filenames

  #' Latest version in archive
  latest <- orderly::orderly_latest(report)

  #' Copy files over
  files_from <- paste0("archive/", report, "/", latest, "/", filenames)
  files_to <- paste0("docs/", filenames)

  file.copy(from = files_from, to = files_to, overwrite = TRUE)
}

archive_to_docs("docs_paper")
