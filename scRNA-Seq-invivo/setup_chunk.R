ABS_FIG_WIDTH = 2*2*2*2*2*2*5*5 #2*2*2*5*7*3#2*2*5*5*3*3
PIXEL_RATIO = 1.25
REL_FIG_WIDTH = 2.2
FIG_HEIGHT = 3.5
FIG_WIDTH = FIG_HEIGHT*REL_FIG_WIDTH
DPI = ABS_FIG_WIDTH/FIG_WIDTH*PIXEL_RATIO
INTERACTIVE = is.null(knitr::current_input(dir = TRUE))
ORIGINAL_WD = if(INTERACTIVE) getwd() else fs::path_dir(knitr::current_input(dir = TRUE))
RES_PATH = fs::path_rel(fs::path_abs(fs::path(params$results_dir), start=rprojroot::find_rstudio_root_file()), ORIGINAL_WD)

fs::dir_create(RES_PATH)

knitr::opts_knit$set(root.dir = rprojroot::find_rstudio_root_file())

if(INTERACTIVE){
} else {
}
knitr::opts_chunk$set(
  fig.align="left",
  out.width = ABS_FIG_WIDTH,
  paged.print = TRUE,
  df_print = "paged",
  fig.width = FIG_WIDTH, # if(is.null(knitr::current_input())) 2*FIG_WIDTH else FIG_WIDTH,
  fig.height = FIG_HEIGHT, #if(is.null(knitr::current_input())) 2*FIG_HEIGHT else FIG_HEIGHT,
  dev= c("png", "svg", "pdf"),
  dpi=DPI,
  dev.args = list(bg = 'transparent')
)

# change ggplot default theme and colors

ggplot2::theme_set(ggplot2::theme_classic())
ggplot2::theme_update(aspect.ratio=1)
options(ggplot2.continuous.colour="viridis")
options(ggplot2.continuous.fill="viridis")

some_colors <- c('#e6194b',  '#3cb44b', '#ffe119', '#4363d8', '#f58231',  '#911eb4', '#46f0f0', '#f032e6', '#bcf60c', '#fabebe', '#008080', '#e6beff', '#9a6324', '#fffac8', '#800000', '#aaffc3', '#808000', '#ffd8b1', '#000075', '#808080', '#ffffff', '#000000')

dicrete_colors <- list(
  c("skyblue", "orange"),
  RColorBrewer::brewer.pal(6, "Dark2"),
  some_colors
)
options(ggplot2.discrete.colour = dicrete_colors,
        ggplot2.discrete.fill = dicrete_colors)


options("ggrastr.default.dpi" = 600)


# to print the chucnk label in front of each chunk
knitr::knit_hooks$set(fig.path = function(before, options, envir) {
  if(before) paste0("###### ", knitr::opts_current$get()$label, "\n")
})

knitr::asis_output(
"
<details><summary>Show revision, parameter and session details</summary>
"
)

# use more screen area in output
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

# enable timing of chunks
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

sessioninfo::session_info()

knitr::asis_output(
  "
</details>
"
)
