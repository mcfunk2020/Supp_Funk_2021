analysis_joint17: results_human/analysis_joint17/4x1_differential_expression_comparison_organoid_in_vivo/4x1_differential_expression_comparison_organoid_in_vivo.html results/4x1_differential_expression_comparison_organoid_in_vivo/36aa3e20/4x1_differential_expression_comparison_organoid_in_vivo.html

clean_analysis_joint17: clean_results/4x1_differential_expression_comparison_organoid_in_vivo/36aa3e20

clean_results/4x1_differential_expression_comparison_organoid_in_vivo/36aa3e20:
	-rm -rf "results/4x1_differential_expression_comparison_organoid_in_vivo/36aa3e20"
	-rm -rf "results_human/analysis_joint17/4x1_differential_expression_comparison_organoid_in_vivo"

results/4x1_differential_expression_comparison_organoid_in_vivo/36aa3e20/4x1_differential_expression_comparison_organoid_in_vivo.html: notebooks/4x1_differential_expression_comparison_organoid_in_vivo.Rmd
	-rm -rf "results/4x1_differential_expression_comparison_organoid_in_vivo/36aa3e20"
	-mkdir -p "results/4x1_differential_expression_comparison_organoid_in_vivo/36aa3e20"
	Rscript -e 'rmarkdown::render(input = structure("notebooks/4x1_differential_expression_comparison_organoid_in_vivo.Rmd", class = c("fs_path", "character")), output_file = structure("results/4x1_differential_expression_comparison_organoid_in_vivo/36aa3e20/4x1_differential_expression_comparison_organoid_in_vivo.html", class = c("fs_path", "character")), output_dir = structure("results/4x1_differential_expression_comparison_organoid_in_vivo/36aa3e20", class = c("fs_path", "character")), params = list(bulk_data_file_organoid = "raw_data/DE_young_aged_bulk_organoid.tsv",     bulk_data_file_in_vivo = "raw_data/DE_young_aged_bulk_in_vivo.tsv",     single_cell_data_file_organoid = "results_human/analysis_17only/4_differential_expression/DE_samples_as_replicates.tsv",     single_cell_data_file_in_vivo = "results_human/analysis_in_vivo/4_differential_expression/DE_samples_as_replicates.tsv",     results_dir = structure("results/4x1_differential_expression_comparison_organoid_in_vivo/36aa3e20", class = c("fs_path",     "character"))), intermediates_dir = tempdir(), clean = FALSE)'

results_human/analysis_joint17/4x1_differential_expression_comparison_organoid_in_vivo/4x1_differential_expression_comparison_organoid_in_vivo.html: results/4x1_differential_expression_comparison_organoid_in_vivo/36aa3e20/4x1_differential_expression_comparison_organoid_in_vivo.html
	-rm -rf "results_human/analysis_joint17/4x1_differential_expression_comparison_organoid_in_vivo"
	mkdir -p "results_human/analysis_joint17"
	-ln -s "$$(realpath --relative-to="results_human/analysis_joint17" "results/4x1_differential_expression_comparison_organoid_in_vivo/36aa3e20")" "results_human/analysis_joint17/4x1_differential_expression_comparison_organoid_in_vivo"
