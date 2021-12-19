do_memoize <- function(fun) {
  fun <- substitute(fun)
  if(!eval(substitute(memoise::is.memoized(fun), list(fun=fun))))
    assign(deparse(fun), eval(substitute(memoise::memoize(fun), list(fun=fun))), envir=parent.frame())
}

memoize_if_not <- function(fun) {
  fun <- substitute(fun)
  if(!eval(substitute(memoise::is.memoized(fun), list(fun=fun))))
    eval(substitute(memoise::memoize(fun), list(fun=fun)))
  else
    eval(substitute(fun, list(fun=fun)))
}

#' @export
geom_cell <- function(..., size = 0.2, raster.dpi = 600) ggrastr::geom_point_rast(..., size = size, raster.dpi = raster.dpi)

#' @import magrittr
get_scores <- function(data, features) {
  lapply(features, function(features) {
    colMeans(data[features %>% .[. %in% rownames(data)], ])
  }) %>% as.data.frame()
}

#' @import magrittr
get_scores2 <- function(data, features) {
  lapply(features, function(features) {

    dat <- t(data[features %>% .[. %in% rownames(data)], ])
    res <- tryCatch(irlba::irlba(dat, nv=1),
                    warning = function(w) svd(dat, nu=1, nv=1),
                    error=function(e) svd(dat, nu=1, nv=1))
    with(res, u*sign(mean(v)))
  }) %>% as.data.frame(row.names=colnames(data))
}

ensembl_to_gene_name <- function(ensembl_ids, mart=useEnsembl(biomart="ensembl", dataset="mmusculus_gene_ensembl", GRCh=37)) {
  gene_names <- getBM(filters= "ensembl_gene_id",
                      attributes= c("external_gene_source","external_gene_name","mgi_symbol","ensembl_gene_id"),
                      values=ensembl_ids,
                      mart=mart)
  setDT(gene_names)
  setnames(gene_names, "external_gene_name", "gene_name")
  setkey(gene_names, "ensemble_gene_id")
  gene_names[ensemble_gene_id, external_gene_name]
}

#' convert gene names to ensemble ids using biomarT
#' @param external_gene_names character vector of gene names to convert
#' @import data.table
#' @export
gene_name_to_ensembl <- function(external_gene_names, mart=useEnsembl(biomart="ensembl", dataset="mmusculus_gene_ensembl", GRCh=37)) {
  gene_names <- getBM(filters= "external_gene_name",
                      attributes= c("external_gene_name","ensembl_gene_id"),
                      values=external_gene_names,
                      mart=mart)
  setDT(gene_names)
  gene_names <- gene_names[,.(ensembl_gene_id=paste(unique(ensembl_gene_id), collapse = ", ")), by=external_gene_name]
  setkey(gene_names, "external_gene_name")
  gene_names[external_gene_names, ensembl_gene_id]

}
