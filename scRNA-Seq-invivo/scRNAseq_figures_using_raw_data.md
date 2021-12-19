---
title: scRNAseq analysis plots
author: "Jan Gleixner, Maja Funk"
date: "2021-12-01"
output:
  html_document:
    df_print: paged
    code_folding: hide
    keep_md: true
params:
  aggregated_data_dir: "processed_data"
  processed_data_dir: "processed_data"
  raw_data_dir: "raw_data"
  differential_expression_dir: "processed_data"
  meta_data_dir: "analysis_pipeline"
  results_dir: "results"
---



```
## NULL
```


<details><summary>Show revision, parameter and session details</summary>

<style>
  .main-container {
    max-width: 1658px !important;
  }
</style>
<div data-pagedtable="false">
  <script data-pagedtable-source type="application/json">
{"columns":[{"label":["parameter"],"name":[1],"type":["chr"],"align":["left"]},{"label":["value"],"name":[2],"type":["chr"],"align":["left"]}],"data":[{"1":"aggregated_data_dir","2":"\"processed_data\""},{"1":"processed_data_dir","2":"\"processed_data\""},{"1":"raw_data_dir","2":"\"raw_data\""},{"1":"differential_expression_dir","2":"\"processed_data\""},{"1":"meta_data_dir","2":"\"analysis_pipeline\""},{"1":"results_dir","2":"\"results\""}],"options":{"columns":{"min":{},"max":[10]},"rows":{"min":[10],"max":[10]},"pages":{}}}
  </script>
</div>

```
## Git revision: 4bc7e84
## git diff:
##  M processed_data/DDE_samples_as_replicates.tsv.gz
##  M processed_data/DE_per_cell_type_conditioned_on_bulk.tsv.gz
##  M processed_data/DE_samples_as_replicates.tsv.gz
##  M processed_data/aggregated_expression.tsv.gz
##  M processed_data/dde_one_vs_sum_of_rest_dt.tsv.gz
##  M processed_data/meta_and_reductions.tsv.gz
##  M scRNAseq_figures_using_only_processed_data.Rmd
##  M scRNAseq_figures_using_raw_data.Rmd
## ?? analysis_pipeline/
## ?? results/
## ?? scRNAseq_figures_using_only_processed_data.html
## ?? scRNAseq_figures_using_only_processed_data.md
## ?? scRNAseq_figures_using_only_processed_data_files/
## ?? scRNAseq_figures_using_raw_data.html
## ?? scRNAseq_figures_using_raw_data.md
## ?? scRNAseq_figures_using_raw_data_files/
## Registered S3 method overwritten by 'cli':
##   method     from         
##   print.boxx spatstat.geom
## ─ Session info ───────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
##  setting  value                       
##  version  R version 4.1.0 (2021-05-18)
##  os       CentOS Linux 7 (Core)       
##  system   x86_64, linux-gnu           
##  ui       X11                         
##  language (EN)                        
##  collate  en_US.UTF-8                 
##  ctype    en_US.UTF-8                 
##  tz       Europe/Berlin               
##  date     2021-12-01                  
## 
## ─ Packages ───────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
##  package         * version  date       lib source             
##  abind             1.4-5    2016-07-21 [1] standard (@1.4-5)  
##  assertthat        0.2.1    2019-03-21 [1] CRAN (R 4.1.0)     
##  bslib             0.3.1    2021-10-06 [1] standard (@0.3.1)  
##  cli               3.1.0    2021-10-27 [1] standard (@3.1.0)  
##  cluster           2.1.2    2021-04-17 [2] CRAN (R 4.1.0)     
##  codetools         0.2-18   2020-11-04 [2] CRAN (R 4.1.0)     
##  colorspace        2.0-2    2021-06-24 [1] standard (@2.0-2)  
##  cowplot           1.1.1    2020-12-30 [1] standard (@1.1.1)  
##  crayon            1.4.2    2021-10-29 [1] standard (@1.4.2)  
##  data.table      * 1.14.2   2021-09-27 [1] standard (@1.14.2) 
##  DBI               1.1.1    2021-01-15 [1] standard (@1.1.1)  
##  deldir            1.0-6    2021-10-23 [1] standard (@1.0-6)  
##  digest            0.6.28   2021-09-23 [1] standard (@0.6.28) 
##  dplyr             1.0.7    2021-06-18 [1] standard (@1.0.7)  
##  ellipsis          0.3.2    2021-04-29 [1] standard (@0.3.2)  
##  evaluate          0.14     2019-05-28 [1] standard (@0.14)   
##  fansi             0.5.0    2021-05-25 [1] standard (@0.5.0)  
##  farver            2.1.0    2021-02-28 [1] standard (@2.1.0)  
##  fastmap           1.1.0    2021-01-25 [1] standard (@1.1.0)  
##  fitdistrplus      1.1-6    2021-09-28 [1] standard (@1.1-6)  
##  fs                1.5.0    2020-07-31 [1] standard (@1.5.0)  
##  future            1.23.0   2021-10-31 [1] standard (@1.23.0) 
##  future.apply      1.8.1    2021-08-10 [1] standard (@1.8.1)  
##  gamlss            5.3-4    2021-03-31 [1] standard (@5.3-4)  
##  gamlss.data       6.0-1    2021-03-18 [1] standard (@6.0-1)  
##  gamlss.dist       5.3-2    2021-03-09 [1] standard (@5.3-2)  
##  gamlss.tr         5.1-7    2020-07-13 [1] standard (@5.1-7)  
##  generics          0.1.1    2021-10-25 [1] standard (@0.1.1)  
##  ggforce           0.3.3    2021-03-05 [1] standard (@0.3.3)  
##  ggplot2         * 3.3.5    2021-06-25 [1] standard (@3.3.5)  
##  ggrepel           0.9.1    2021-01-15 [1] standard (@0.9.1)  
##  ggridges          0.5.3    2021-01-08 [1] standard (@0.5.3)  
##  globals           0.14.0   2020-11-22 [1] standard (@0.14.0) 
##  glue              1.5.0    2021-11-07 [1] standard (@1.5.0)  
##  goftest           1.2-3    2021-10-07 [1] standard (@1.2-3)  
##  gridExtra         2.3      2017-09-09 [1] standard (@2.3)    
##  gtable            0.3.0    2019-03-25 [1] standard (@0.3.0)  
##  htmltools         0.5.2    2021-08-25 [1] standard (@0.5.2)  
##  htmlwidgets       1.5.4    2021-09-08 [1] standard (@1.5.4)  
##  httpuv            1.6.3    2021-09-09 [1] standard (@1.6.3)  
##  httr              1.4.2    2020-07-20 [1] standard (@1.4.2)  
##  ica               1.0-2    2018-05-24 [1] standard (@1.0-2)  
##  igraph            1.2.8    2021-11-07 [1] standard (@1.2.8)  
##  irlba             2.3.3    2019-02-05 [1] standard (@2.3.3)  
##  jquerylib         0.1.4    2021-04-26 [1] standard (@0.1.4)  
##  jsonlite          1.7.2    2020-12-09 [1] standard (@1.7.2)  
##  KernSmooth        2.23-20  2021-05-03 [2] CRAN (R 4.1.0)     
##  knitr             1.36     2021-09-29 [1] standard (@1.36)   
##  later             1.3.0    2021-08-18 [1] standard (@1.3.0)  
##  lattice           0.20-45  2021-09-22 [1] standard (@0.20-45)
##  lazyeval          0.2.2    2019-03-15 [1] standard (@0.2.2)  
##  leiden            0.3.9    2021-07-27 [1] CRAN (R 4.1.0)     
##  lifecycle         1.0.1    2021-09-24 [1] standard (@1.0.1)  
##  listenv           0.8.0    2019-12-05 [1] standard (@0.8.0)  
##  lmtest            0.9-39   2021-11-07 [1] standard (@0.9-39) 
##  magrittr          2.0.1    2020-11-17 [1] standard (@2.0.1)  
##  MASS              7.3-54   2021-05-03 [2] CRAN (R 4.1.0)     
##  Matrix            1.3-4    2021-06-01 [2] CRAN (R 4.1.0)     
##  matrixStats       0.61.0   2021-09-17 [1] standard (@0.61.0) 
##  mgcv              1.8-38   2021-10-06 [1] standard (@1.8-38) 
##  mime              0.12     2021-09-28 [1] standard (@0.12)   
##  miniUI            0.1.1.1  2018-05-18 [1] standard (@0.1.1.1)
##  munsell           0.5.0    2018-06-12 [1] standard (@0.5.0)  
##  nlme              3.1-153  2021-09-07 [1] standard (@3.1-153)
##  parallelly        1.28.1   2021-09-09 [1] standard (@1.28.1) 
##  patchwork       * 1.1.1    2020-12-17 [1] standard (@1.1.1)  
##  pbapply           1.5-0    2021-09-16 [1] standard (@1.5-0)  
##  pillar            1.6.4    2021-10-18 [1] standard (@1.6.4)  
##  pkgconfig         2.0.3    2019-09-22 [1] standard (@2.0.3)  
##  plotly            4.10.0   2021-10-09 [1] standard (@4.10.0) 
##  plyr              1.8.6    2020-03-03 [1] standard (@1.8.6)  
##  png               0.1-7    2013-12-03 [1] standard (@0.1-7)  
##  polyclip          1.10-0   2019-03-14 [1] standard (@1.10-0) 
##  promises          1.2.0.1  2021-02-11 [1] standard (@1.2.0.1)
##  purrr             0.3.4    2020-04-17 [1] standard (@0.3.4)  
##  R6                2.5.1    2021-08-19 [1] standard (@2.5.1)  
##  RANN              2.6.1    2019-01-08 [1] standard (@2.6.1)  
##  RColorBrewer      1.1-2    2014-12-07 [1] standard (@1.1-2)  
##  Rcpp              1.0.7    2021-07-07 [1] standard (@1.0.7)  
##  RcppAnnoy         0.0.19   2021-07-30 [1] standard (@0.0.19) 
##  reshape2          1.4.4    2020-04-09 [1] standard (@1.4.4)  
##  reticulate        1.22     2021-09-17 [1] standard (@1.22)   
##  rlang             0.4.12   2021-10-18 [1] standard (@0.4.12) 
##  rmarkdown         2.11     2021-09-14 [1] standard (@2.11)   
##  ROCR              1.0-11   2020-05-02 [1] standard (@1.0-11) 
##  rpart             4.1-15   2019-04-12 [2] CRAN (R 4.1.0)     
##  rprojroot         2.0.2    2020-11-15 [1] standard (@2.0.2)  
##  rstudioapi        0.13     2020-11-12 [1] standard (@0.13)   
##  Rtsne             0.15     2018-11-10 [1] standard (@0.15)   
##  sass              0.4.0    2021-05-12 [1] standard (@0.4.0)  
##  scales            1.1.1    2020-05-11 [1] standard (@1.1.1)  
##  scattermore       0.7      2020-11-24 [1] standard (@0.7)    
##  schelpr         * 0.0.9017 2021-11-30 [1] local              
##  sctransform       0.3.2    2020-12-16 [1] standard (@0.3.2)  
##  sessioninfo       1.1.1    2018-11-05 [2] CRAN (R 4.1.0)     
##  Seurat          * 4.0.5    2021-10-17 [1] standard (@4.0.5)  
##  SeuratObject    * 4.0.3    2021-11-10 [1] standard (@4.0.3)  
##  shiny             1.7.1    2021-10-02 [1] standard (@1.7.1)  
##  spatstat.core     2.3-1    2021-11-02 [1] standard (@2.3-1)  
##  spatstat.data     2.1-0    2021-03-21 [1] standard (@2.1-0)  
##  spatstat.geom     2.3-0    2021-10-09 [1] standard (@2.3-0)  
##  spatstat.sparse   2.0-0    2021-03-16 [1] standard (@2.0-0)  
##  spatstat.utils    2.2-0    2021-06-14 [1] standard (@2.2-0)  
##  stringi           1.7.5    2021-10-04 [1] standard (@1.7.5)  
##  stringr           1.4.0    2019-02-10 [1] standard (@1.4.0)  
##  survival          3.2-13   2021-08-24 [1] standard (@3.2-13) 
##  tensor            1.5      2012-05-05 [1] standard (@1.5)    
##  tibble            3.1.6    2021-11-07 [1] standard (@3.1.6)  
##  tidyr             1.1.4    2021-09-27 [1] standard (@1.1.4)  
##  tidyselect        1.1.1    2021-04-30 [1] standard (@1.1.1)  
##  tweenr            1.0.2    2021-03-23 [1] standard (@1.0.2)  
##  utf8              1.2.2    2021-07-24 [1] standard (@1.2.2)  
##  uwot              0.1.10   2020-12-15 [1] standard (@0.1.10) 
##  vctrs             0.3.8    2021-04-29 [1] standard (@0.3.8)  
##  viridisLite       0.4.0    2021-04-13 [1] standard (@0.4.0)  
##  withr             2.4.2    2021-04-18 [1] standard (@2.4.2)  
##  xfun              0.27     2021-10-18 [1] standard (@0.27)   
##  xtable            1.8-4    2019-04-21 [1] standard (@1.8-4)  
##  yaml              2.2.1    2020-02-01 [1] standard (@2.2.1)  
##  zoo               1.8-9    2021-03-09 [1] standard (@1.8-9)  
## 
## [1] /home/gleixner/R/x86_64-pc-linux-gnu-library/4.1.0b
## [2] /software/r/4.1.0/lib64/R/library
```


</details>

###### load helper functions

```{.r .fold-hide}
sgi_blue    = '#5087C8'
sgi_yellow1 = '#F2EE35'
sgi_yellow2 = '#FED98E'
b110_grey   = '#808080'
b110_grey_light   = '#909090'
b110_transparent_black = ggplot2::alpha('#000000', 0.5)
google_red = '#dd4b39'
google_green = '#0F9D58'
google_yellow = '#F4B400'
google_blue = '#4285F4'
b110_colors <- c(sgi_blue, sgi_yellow1, sgi_yellow2, b110_grey, b110_grey_light, b110_transparent_black, google_red, google_green, google_yellow, google_blue)
legend_dot_size <- 1

maja_red <- alpha('#E64B35CC', 0.8)
maja_blue <- alpha('#4DBBD5CC', 0.8)
maja_gray <- alpha('gray', 0.8)

labeller_set <- function(width = Inf, ...) {
  if (width == Inf) {
    labeller(set = labels_of_set, ...) # returned
  } else {
    wrapped_labels_of_set <- stringr::str_wrap(labels_of_set, width = width)
    names(wrapped_labels_of_set)  <- levels_of_set
    labeller(set = wrapped_labels_of_set, ...)  # returned
  }
}

levels_of_cell_type <- c("pseudo bulk", "dead cells", "Stem",       "TA",       "Enterocyte Progenitor", "Enterocyte",  "EEC",     "Paneth",       "Goblet & Paneth",         "Goblet",       "Tuft"      )
labels_of_cell_type <- c("pseudo bulk", "dead cells", "Stem cells", "TA cells", "EC prog.",              "Enterocytes", "EECs",    "Paneth cells", "Goblet &\nPaneth cells",  "Goblet cells", "Tuft cells")
colors_of_cell_type <- c("#303030",     "#303030",    "#66a61e",    "#1b9e77",  "#418DA6",               "#7570b3",     "#a846a0", "#ce5a02",      "#a6761d",                 "#e6ab02",      "#FE5D26"   )  # old tuft: e7298a #E53D00
names(colors_of_cell_type) <- levels_of_cell_type

scale_x_cell_type <- function(
  name = "cell type",
  breaks = levels_of_cell_type,
  labels = labels_of_cell_type,
  limits = function(x) levels_of_cell_type[levels_of_cell_type %in% x],
  guide = guide_axis(
    #n.dodge=2,
    angle = 45,
    title = element_blank()),
  ...
) scale_x_discrete(name = name, breaks = breaks, labels = labels, guide = guide, limits = limits, ...)

scale_y_cell_type <- function(
  name = "cell type",
  breaks = rev(levels_of_cell_type),
  labels = rev(labels_of_cell_type),
  limits = function(x) rev(levels_of_cell_type[levels_of_cell_type %in% x]),
  guide = guide_axis(
    #n.dodge=2,
    #angle = 45,
    title = element_blank()),
  ...
) scale_y_discrete(name = name, breaks = breaks, labels = labels, guide = guide, limits = limits, ...)

scale_color_cell_type <- function(
  name = "cell type",
  breaks = levels_of_cell_type,
  labels = labels_of_cell_type,
  values = colors_of_cell_type,
  limits = identity,
  na.value = "#808080",
  guide = guide_legend(override.aes = list(size = legend_dot_size)),
  ...
) scale_color_manual(
  name = name,
  breaks = breaks,
  limits = limits,
  labels = labels,
  values = values,
  na.value = na.value,
  guide = guide,
  ...
)

scale_fill_cell_type <- function(
  name = "cell type",
  breaks = levels_of_cell_type,
  labels = labels_of_cell_type,
  values = colors_of_cell_type,
  limits = identity,
  na.value = "#808080",
  guide = guide_legend(override.aes = list(size = legend_dot_size)),
  ...
) scale_fill_manual(
  name = name,
  breaks = breaks,
  limits = limits,
  labels = labels,
  values = values,
  na.value = na.value,
  guide = guide,
  ...
)

scale_x_age <- function(..., limits = c("young", "aged")) scale_x_discrete(..., limits = limits)

age_colors <- c(young = "#91D1C2CC", aged = "#3C5488CC")
age_colors_light <- c(young = "#91D1C2CC", aged = "#8491B4CC")
age_colors_dark <- c(young = "#00A087CC", aged = "#3C5488CC")

scale_color_age <- function(..., limits = names(age_colors), values = age_colors, na.value = "grey", guide = guide_legend(override.aes = list(size = legend_dot_size))) ggplot2::scale_color_manual(..., limits = limits, values = values, na.value = na.value, guide = guide)
#scale_color_age_dark <- function(..., values = age_colors_dark) scale_color_age(..., values = values)

scale_fill_age <- function(..., limits = names(age_colors), values = age_colors, na.value = "grey", guide = guide_legend(override.aes = list(size = legend_dot_size))) ggplot2::scale_fill_manual(..., limits = limits, values = values, na.value = na.value, guide = guide)
scale_fill_age_dark <- function(..., values = age_colors_dark) scale_fill_age(..., values = values)

volcano_plot <- function(dt, alpha = 0.1, highlight = rank(padj)<=20, color = "black", box.padding = 0.25, label_size = 4, ...){
  if(missing(alpha)) message(paste0("FDR <= alpha = ", alpha))
  dt <- copy(dt)
  dt[padj > alpha, sir.class := "not significant"]
  dt[padj <= alpha & -log2FoldChange > 0, sir.class := "up regulated"]
  dt[padj <= alpha & -log2FoldChange < 0, sir.class := "down regulated"]
  labelled_subset_dt <- eval(substitute(copy(dt)[!(highlight), gene := ""]))
  
  dt %>% ggplot(aes(x = -log2FoldChange, y = -log10(pvalue))) +
    #geom_vline(xintercept = 0, linetype = "dotted") +
    geom_vline(xintercept = 0, linetype = "dashed") +
    ggrastr::rasterize(geom_point(data = dt[sir.class == "not significant",],aes(color = sir.class), size=0.5)) +
    geom_point(data = dt[sir.class != "not significant",], aes(color = sir.class), size=0.5) +
    scale_color_manual(breaks = c("down regulated", "up regulated", "not significant"), values = c(maja_blue, maja_red, maja_gray)) +
    guides(color = guide_legend(title = element_blank())) +
    labs(x = "log2fc aged/young") +
    # geom_vline(xintercept = c(0.5, -0.5), linetype = "dotted") + geom_hline(yintercept = 2, linetype = "dotted") +
    if(nrow(labelled_subset_dt)>0) ggrepel::geom_label_repel(aes(label = gene), data = labelled_subset_dt, color = color, min.segment.length = 0, size = label_size, box.padding = box.padding, max.overlaps = Inf, ...) else NULL +
    NULL
}

geom_cell <- function(..., size = 0.2, raster.dpi = 600) ggrastr::geom_point_rast(..., size = size, raster.dpi = raster.dpi)
```

###### load cell meta data

```r
dt <- fread(fs::path(processed_data_dir, "meta_and_reductions.tsv.gz"), sep= "\t")
dt[,cell[1], by = sample_name][, stopifnot(all.equal(order(V1), order(sample_name)))] %>% invisible() # ensure that cell and samples match
```

###### load count data and merge with cell meta data to create seurat object

```r
dirs <- fs::path_dir( fs::dir_ls(raw_data_dir, regexp = "matrix.mtx.gz", recurse = TRUE))
som <- Seurat::Read10X(dirs)
so <- Seurat::CreateSeuratObject(som, meta.data = data.frame(dt, row.names = dt$cell))
so <- so[, !is.na(so$cell)]
```

###### normalize

```r
so %<>% Seurat::NormalizeData()
```

###### marker gene and score tSNE plots

```r
example_marker_genes <- c("Stem cells" = "Olfm4", "Paneth cells" = "Defa24", "Goblet cells" = "Muc2", "Enteroendocrine cells" = "Chgb", "Enterocytes" = "Fabp1", "Tuft cells" = "Dclk1")
cell_type_labels <- c("Stem cells" = "Stem", "Paneth cells" = "Paneth", "Goblet cells" = "Goblet", "Enteroendocrine cells" = "EEC", "Enterocytes" = "Enterocyte", "Tuft cells" = "Tuft")


expression_dt <- as.data.table((exp(Matrix::t(Seurat::GetAssayData(so)[example_marker_genes,]))-1)*1E6/1E4 , keep.rownames = "cell" )
dt <- expression_dt[dt, on = "cell"]

p_gene_s <- lapply(example_marker_genes, function(gene) {
  dt %>% ggplot(aes_string(x = "tSNE_1", y = "tSNE_2", color = gene)) + geom_cell() +  schelpr::scale_color_log1p_scaled(scale=1/100, name = paste0(gene, "\n[CPM]"))
})


p_score_s <- lapply(names(example_marker_genes), function(cell_type) {
  dt %>% ggplot(aes_string(x = "tSNE_1", y = "tSNE_2", color = paste0("`", cell_type, "`"))) + geom_cell() + scale_color_continuous(stringr::str_wrap(paste0(cell_type_labels[cell_type], " score"), 12)) #my_wrap(cell_type))
})

p_name_s <- lapply(names(example_marker_genes), function(cell_type) {
  wrap_elements(plot=grobTree(rectGrob(gp=gpar(col= 'black')), textGrob(cell_type, rot = 90)), clip = TRUE) 
})

wrap_plots(c(p_name_s, p_gene_s, p_score_s), ncol = 3, byrow=FALSE, widths = c(0.13, 1, 1)) + plot_annotation(tag_levels = list(c(letters[seq_along(p_name_s)],character(length(p_name_s)*2))))  
```

<img src="scRNAseq_figures_using_raw_data_files/figure-html/marker gene and score tSNE plots-1.png" width="1600" style="display: block; margin: auto auto auto 0;" />

###### marker gene and score umap plots

```r
example_marker_genes <- c("Stem" = "Olfm4", "Enterocyte.Mature.Proximal" = "Fabp1","Endocrine" = "Chgb", "Paneth" = "Defa24", "Goblet" = "Muc2",  "Tuft" = "Dclk1")
cell_type_labels <- c("Stem" = "Stem", "Paneth" = "Paneth", "Goblet" = "Goblet", "Endocrine" = "EEC", "Enterocyte.Mature.Proximal" = "Enterocyte", "Tuft" = "Tuft")
cell_type_labels2 <- c("Stem" = "Stem cells", "Paneth" = "Paneth cells", "Goblet" = "Goblet cells", "Endocrine" = "EECs", "Enterocyte.Mature.Proximal" = "Enterocytes", "Tuft" = "Tuft cells")

expression_dt <- as.data.table((exp(Matrix::t(Seurat::GetAssayData(so)[example_marker_genes,]))-1)*1E6/1E4 , keep.rownames = "cell" )
dt <- expression_dt[dt, on = "cell"]

p_gene_s <- lapply(example_marker_genes, function(gene) {
  dt %>% ggplot(aes_string(x = "tSNE_1", y = "tSNE_2", color = gene)) + geom_cell() + schelpr::scale_color_log1p_scaled(scale=1/100, name = paste0(gene, "\n[CPM]"))
})
p_score_s <- lapply(names(example_marker_genes), function(cell_type) {
  dt %>% ggplot(aes_string(x = "tSNE_1", y = "tSNE_2", color = paste0("`prediction.score.", cell_type, "`"))) + geom_cell() + scale_color_continuous(stringr::str_wrap(paste0(cell_type_labels[cell_type], " score"), 12))
})

p_name_s <- lapply(names(example_marker_genes), function(cell_type) {
  wrap_elements(plot=grobTree(rectGrob(gp=gpar(col= 'black')), textGrob(cell_type_labels[cell_type], rot = 90)), clip = TRUE) 
})

wrap_plots(c(p_name_s, p_gene_s, p_score_s), ncol = 3, byrow=FALSE, widths = c(0.13, 1, 1)) + plot_annotation(tag_levels = list(c(letters[seq_along(p_name_s)],character(length(p_name_s)*2))))  
```

<img src="scRNAseq_figures_using_raw_data_files/figure-html/marker gene and score umap plots-1.png" width="1600" style="display: block; margin: auto auto auto 0;" />
























