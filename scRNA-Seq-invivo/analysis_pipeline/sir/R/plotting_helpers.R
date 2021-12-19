maja_red <- alpha('#E64B35CC', 0.8)
maja_blue <- alpha('#4DBBD5CC', 0.8)
maja_gray <- alpha('gray', 0.8)

#' @export
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
