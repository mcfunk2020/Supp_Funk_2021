library(analysismaker)

analysis <- new_analysis() %>%
  add_notebook(haber_so_dir = "0_prepare_Haber_data.Rmd") %>%
  add_notebook(qc_results_dir = "0_qc_filtering.Rmd") %>%
  add_notebook("0x1_unbiased_dimensionality_reduction.Rmd", seurat_object_dir = qc_results_dir) %>%
  add_notebook("0x1_integration_with_Haber.Rmd", processed_data_dir=qc_results_dir) %>%
  add_notebook("0x1_qc_stats.Rmd") %>%
  add_notebook(scaled_data_dir = "1_dimensionality_reduction.Rmd", processed_data_dir=qc_results_dir) %>%
  add_notebook(seurat_object_dir = "2_cell_type_assignment.Rmd", so_dir = scaled_data_dir, haber_so_dir) %>%
  add_notebook("2x1_stemness.Rmd", processed_data_dir = seurat_object_dir) %>%
  add_notebook("2x1_subcluster_analysis.Rmd", processed_data_dir = seurat_object_dir) %>%
  add_notebook("2x1_plot_features.Rmd", processed_data_dir = seurat_object_dir) %>%
  add_notebook("2x1_plot_marker_genes.Rmd", processed_data_dir = seurat_object_dir) %>%
  add_notebook("2x1_plot_meta.Rmd", processed_data_dir = seurat_object_dir) %>%
  add_notebook("2x1_differenatial_composition.Rmd", processed_data_dir = seurat_object_dir) %>%
  add_notebook(aggregated_data_dir = "3_aggregation.Rmd", processed_data_dir = seurat_object_dir) %>%
  add_notebook("3x1_agg_plots.Rmd", aggregated_data_dir) %>%
  add_notebook("3x1_aggregation_viz.Rmd", aggregated_data_dir, processed_data_dir = seurat_object_dir) %>%
  add_notebook("3x1_sample_pca_embedding.Rmd", aggregated_data_dir) %>%
  add_notebook(differential_expression_dir = "4_differential_expression.Rmd", aggregated_data_dir) %>%
  add_notebook(dde_dir = "4_differential_differential_expression.Rmd", aggregated_data_dir) %>%
  add_notebook(ddeII_dir = "4_differential_differential_expression_II.Rmd", aggregated_data_dir) %>%
  add_notebook("4x1_differential_differential_expression_II_plots.Rmd", aggregated_data_dir, de_results_dir = differential_expression_dir, dde_results_dir=ddeII_dir) %>%
  add_notebook("4x1_cell_type_specific_de_plots.Rmd", aggregated_data_dir, differential_expression_dir) %>%
  add_notebook("4x1_plot_differential_differential_expression.Rmd", aggregated_data_dir, dde_dir) %>%
  add_notebook("9_paper_figure.Rmd", processed_data_dir = seurat_object_dir, aggregated_data_dir, differential_expression_dir) %>%
  add_notebook("4x1_differential_expression_comparison_with_bulk.Rmd", single_cell_DE_data_dir = differential_expression_dir) %>%
  identity()

analysis_in_vivo <-  analysis %>% bind_parameters(
  data_dir = "raw_data",
  subset_10x_samples =  c("18.1", "18.2", "18.3"),
  model_formula = ~ condition + Tx_run_number,
  only_untreated = TRUE,
  merge_goblet_and_paneth = FALSE,
  bulk_data_file = "raw_data/DE_young_aged_bulk_in_vivo.tsv",
  selected_cell_types = c("Paneth", "Goblet")
)
make_makefile(analysis_in_vivo)
