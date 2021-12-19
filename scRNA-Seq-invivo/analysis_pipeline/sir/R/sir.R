#' sir
#'
#' @name sir
#' @docType package
#' @import Seurat
NULL


.onLoad <- function(libname, pkgname) {
  cache_folder <- file.path(dirname(tempdir()), ".rcache")
  dir.create(cache_folder, showWarnings=FALSE)
  cache_file <<- memoise::cache_filesystem(cache_folder)
  useEnsembl <<- memoise::memoise(biomaRt::useEnsembl, cache=cache_file)
  getBM <<- memoise::memoise(biomaRt::getBM, cache=cache_file)
  dir.create("results", showWarnings = FALSE)
  dir.create("processed_data", showWarnings = FALSE)
  invisible(NULL)
}


#' cc.genes_mmusculus
#'
"cc.genes_mmusculus"
