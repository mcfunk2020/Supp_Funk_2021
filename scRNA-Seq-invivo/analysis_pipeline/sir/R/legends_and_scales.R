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
# scales::show_col(b110_colors)

#' @export
levels_of_set <- c(NA, "t2d_candidates", "crc_candidates", "cancer_syndromes", "crohns_candidates", "targeting_negative_controls", "essential", "selected_strong_RNA_phenotypes")

#' @export
labels_of_set <- c(NA, "type II diabetes", "colorectal cancer", "cancer syndromes", "Crohn's disease", "targeting negative controls", "essential genes", "expression effect controls")
names(labels_of_set) <- levels_of_set

extended_dark <- c("#303030", "#1b9e77", "#ce5a02", "#7570b3", "#66a61e", "#e7298a", "#a6761d", "#e6ab02", "#72313e", "#666666")

#' @export
#' @import ggplot2
scale_set <- function(...) scale_color_manual(
  breaks = levels_of_set,
  labels = labels_of_set,
  values = c("#808080", '#ffe119', '#4363d8', '#f58231',  '#3cb44b',  '#e6194b', '#911eb4', '#46f0f0', '#f032e6'),
  na.value = "#808080",
  ...
)

#' @export
#' @import ggplot2
labeller_set <- function(width = Inf, ...) {
  if (width == Inf) {
    labeller(set = labels_of_set, ...) # returned
  } else {
    wrapped_labels_of_set <- stringr::str_wrap(labels_of_set, width = width)
    names(wrapped_labels_of_set)  <- levels_of_set
    labeller(set = wrapped_labels_of_set, ...)  # returned
  }
}

#' @export
levels_of_cell_type <- c("pseudo bulk", "dead cells", "Stem",       "TA",       "Enterocyte Progenitor", "Enterocyte",  "EEC",     "Paneth",       "Goblet & Paneth",         "Goblet",       "Tuft"      )

#' @export
labels_of_cell_type <- c("pseudo bulk", "dead cells", "Stem cells", "TA cells", "EC prog.",              "Enterocytes", "EECs",    "Paneth cells", "Goblet &\nPaneth cells",  "Goblet cells", "Tuft cells")

#' @export
colors_of_cell_type <- c("#303030",     "#303030",    "#66a61e",    "#1b9e77",  "#418DA6",               "#7570b3",     "#a846a0", "#ce5a02",      "#a6761d",                 "#e6ab02",      "#FE5D26"   )  # old tuft: e7298a #E53D00

names(colors_of_cell_type) <- levels_of_cell_type

#' @export
#' @import ggplot2
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

#' @export
#' @import ggplot2
scale_y_cell_type <- function(
  name = "cell type",
  breaks = levels_of_cell_type,
  labels = labels_of_cell_type,
  limits = function(x) levels_of_cell_type[levels_of_cell_type %in% x],
  guide = guide_axis(
    #n.dodge=2,
    #angle = 45,
    title = element_blank()),
  ...
) scale_y_discrete(name = name, breaks = breaks, labels = labels, guide = guide, limits = limits, ...)

#' @export
#' @import ggplot2
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

#' @export
#' @import ggplot2
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

#' @export
#' @import ggplot2
scale_x_age <- function(..., limits = c("young", "aged")) scale_x_discrete(..., limits = limits)


age_colors <- c(young = "#91D1C2CC", aged = "#3C5488CC")
#' @export
#' @import ggplot2
scale_color_age <- function(..., limits = names(age_colors), values = age_colors, na.value = "grey", guide = guide_legend(override.aes = list(size = legend_dot_size))) ggplot2::scale_color_manual(..., limits = limits, values = values, na.value = na.value, guide = guide)
#' @export
#' @import ggplot2
scale_fill_age <- function(..., limits = names(age_colors), values = age_colors, na.value = "grey", guide = guide_legend(override.aes = list(size = legend_dot_size))) ggplot2::scale_fill_manual(..., limits = limits, values = values, na.value = na.value, guide = guide)


#' @export
#' @import ggplot2
scale_set <- function(...) scale_color_manual(
  breaks = levels_of_set,
  labels = labels_of_set,
  values = c("#808080", '#ffe119', '#4363d8', '#f58231',  '#3cb44b',  '#e6194b', '#911eb4', '#46f0f0', '#f032e6'),
  ...
)
