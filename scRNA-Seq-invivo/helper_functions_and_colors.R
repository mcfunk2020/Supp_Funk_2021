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
  breaks = levels_of_cell_type,
  labels = labels_of_cell_type,
  limits = function(x) levels_of_cell_type[levels_of_cell_type %in% x],
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

scale_x_age <- function(..., limits = c("young", "aged")) scale_x_discrete(..., limits = limits)

age_colors <- c(young = "#91D1C2CC", aged = "#3C5488CC")

scale_color_age <- function(..., limits = names(age_colors), values = age_colors, na.value = "grey", guide = guide_legend(override.aes = list(size = legend_dot_size))) ggplot2::scale_color_manual(..., limits = limits, values = values, na.value = na.value, guide = guide)

scale_fill_age <- function(..., limits = names(age_colors), values = age_colors, na.value = "grey", guide = guide_legend(override.aes = list(size = legend_dot_size))) ggplot2::scale_fill_manual(..., limits = limits, values = values, na.value = na.value, guide = guide)

scale_set <- function(...) scale_color_manual(
  breaks = levels_of_set,
  labels = labels_of_set,
  values = c("#808080", '#ffe119', '#4363d8', '#f58231',  '#3cb44b',  '#e6194b', '#911eb4', '#46f0f0', '#f032e6'),
  ...
)

volcano_plot <- function(dt, alpha = 0.1, highlight = rank(padj)<=20, color = "black", box.padding = 0.25, label_size = 4, ...){
  if(missing(alpha)) message(paste0("FDR <= alpha = ", alpha))
  labelled_subset_dt <- eval(substitute(copy(dt)[!(highlight), gene := ""]))
  dt <- copy(dt)
  dt[padj > alpha, sir.class := "not significant"]
  dt[padj <= alpha & -log2FoldChange > 0, sir.class := "up regulated"]
  dt[padj <= alpha & -log2FoldChange < 0, sir.class := "down regulated"]
  
  dt %>% ggplot(aes(x = -log2FoldChange, y = -log10(pvalue))) +
    #geom_vline(xintercept = 0, linetype = "dotted") +
    geom_point(aes(color = sir.class)) +
    scale_color_manual(breaks = c("down regulated", "up regulated", "not significant"), values = c(maja_blue, maja_red, maja_gray)) +
    guides(color = guide_legend(title = element_blank())) +
    labs(x = "log2fc aged/young") +
    if(nrow(labelled_subset_dt)>0) ggrepel::geom_text_repel(aes(label = gene), data = labelled_subset_dt, color = color, min.segment.length = 0, size = label_size, box.padding = box.padding, max.overlaps = Inf, ...) else NULL +
    NULL
}

geom_cell <- function(..., size = 0.2, raster.dpi = 600) ggrastr::geom_point_rast(..., size = size, raster.dpi = raster.dpi)
