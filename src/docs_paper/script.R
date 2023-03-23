# orderly::orderly_develop_start("docs_paper")
# setwd("src/docs_paper")

#' Render documents
rmarkdown::render("paper.Rmd")
rmarkdown::render("appendix.Rmd")
