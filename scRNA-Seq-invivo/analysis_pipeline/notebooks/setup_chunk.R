ABS_FIG_WIDTH = 2*2*2*2*2*2*5*5 #2*2*2*5*7*3#2*2*5*5*3*3
PIXEL_RATIO = 1.25
REL_FIG_WIDTH = 2.2
FIG_HEIGHT = 3.5
FIG_WIDTH = FIG_HEIGHT*REL_FIG_WIDTH
DPI = ABS_FIG_WIDTH/FIG_WIDTH*PIXEL_RATIO
INTERACTIVE = is.null(knitr::current_input(dir = TRUE))
ORIGINAL_WD = if(INTERACTIVE) getwd() else fs::path_dir(knitr::current_input(dir = TRUE))
RES_PATH = fs::path_rel(fs::path_abs(fs::path(params$results_dir), start=rprojroot::find_rstudio_root_file()), ORIGINAL_WD)
#FIG_PATH = paste0(fs::path(RES_PATH, 'figures'), "/")

fs::dir_create(RES_PATH)

knitr::opts_knit$set(root.dir = rprojroot::find_rstudio_root_file())

if(INTERACTIVE){

} else {
  # knitr::opts_chunk$set(error = TRUE)
}
knitr::opts_chunk$set(
  #  collapse = TRUE,
  #  comment = "#>",
  fig.align="left",
  out.width = ABS_FIG_WIDTH,
  paged.print = TRUE,
  df_print = "paged",
  fig.width = FIG_WIDTH, # if(is.null(knitr::current_input())) 2*FIG_WIDTH else FIG_WIDTH,
  fig.height = FIG_HEIGHT, #if(is.null(knitr::current_input())) 2*FIG_HEIGHT else FIG_HEIGHT,
  #fig.path = FIG_PATH,
  dev= c("png", "svg", "pdf"),
  dpi=DPI,
  dev.args = list(bg = 'transparent')
)
options(rlang_trace_top_env = rlang::current_env())
options(error = function() {
  print(rlang::trace_back(bottom = sys.frame(-1)), simplify = "none")
  dump_file = fs::path(RES_PATH, "dump.RData")
  save.image(file = dump_file)
  sink()
  cat(paste0("Error. For post-mortem debugging, use: load(\"", fs::path_real(dump_file), '")\n'))
  print(rlang::trace_back(bottom = sys.frame(-1)), simplify = "none")
  q(save="no", status=1, runLast=FALSE)
})

ggplot2::theme_set(ggplot2::theme_classic())
ggplot2::theme_update(aspect.ratio=1)
options(ggplot2.continuous.colour="viridis")
options(ggplot2.continuous.fill="viridis")

donor_colors <- c('#e6194b',  '#3cb44b', '#ffe119', '#4363d8', '#f58231',  '#911eb4', '#46f0f0', '#f032e6', '#bcf60c', '#fabebe', '#008080', '#e6beff', '#9a6324', '#fffac8', '#800000', '#aaffc3', '#808000', '#ffd8b1', '#000075', '#808080', '#ffffff', '#000000')

dicrete_colors <- list(
  c("skyblue", "orange"),
  RColorBrewer::brewer.pal(6, "Dark2"),
  donor_colors
)
options(ggplot2.discrete.colour = dicrete_colors,
        ggplot2.discrete.fill = dicrete_colors)

options("ggrastr.default.dpi" = 600)

#donor  line AUC Gefitinib responder KRAS location UICC_stage
## 1:  D007 D007T          0.86        no  mut   Rectum          3
## 2:  D013 D013T          0.67       yes   wt   Rectum          1
## 3:  D019 D019T          0.95        no  mut  Sigmoid          1
## 4:  D027 D027T          0.46       yes   wt   Rectum          4
## 5:  D046 D046T          0.98        no  mut   Rectum          3
## 6:  D073 D073T          0.35       yes   wt   Rectum          3
## 7:  D080 D080T          0.84        no <NA>   Rectum          3
## 8:  D090 D090T          0.43       yes <NA>  Sigmoid          2
## 9:  D115 D115T          0.45       yes   wt  Sigmoid          4


knitr::knit_hooks$set(fig.path = function(before, options, envir) {
  if(before) paste0("###### ", knitr::opts_current$get()$label, "\n")
})

options(width = 220)
knitr::asis_output(
  "
<style>
  .main-container {
    max-width: 1658px !important;
  }
</style>
")
# print parameters
data.table::data.table(parameter=names(params), value=sapply(params, deparse1))

# and add them to the global env
invisible(list2env(params, globalenv()))


# print git revision and working directory state
cat(paste0("Git revision: ", system2("git", "rev-parse --short HEAD", stdout = TRUE), "\n"))
cat("git diff:\n")
cat(paste0(system2("git", "status -s", stdout = TRUE), collapse = "\n"))

knitr::knit_hooks$set(time_it = local({
  now <- NULL
  function(before, options) {
    if (before) {
      # record the current time before each chunk
      now <<- Sys.time()
    } else {
      # calculate the time difference after a chunk
      res <- difftime(Sys.time(), now)
      # return a character string to show the time
      paste("Time for the chunk", options$label, "to run:", res)
    }
  }})
)

