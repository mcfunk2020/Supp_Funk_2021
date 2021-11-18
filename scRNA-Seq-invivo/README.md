# Analysis of single cell RNA sequencing data of Funk et al. 2021
This analysis is separated into three independent parts:
1. (analysis_pipeline) - An R project (including a makefile and a helper package) that allows reproduction of the generation of processed data (like `cell_meta.tsv.gz` available on GEO (####)) from the cell ranger outputs available on GEO (####) [requires cell ranger output in `raw_data`]
2.  (scRNAseq_figures_using_raw_data.Rmd) - An R markdown notebook to generate figures ###-### [requires cell ranger output in `raw data` and `cells_meta.tsv.gz` etc. in `processed data`]
3.  [scRNAseq_figures_using_only_processed_data.Rmd](scRNAseq_figures_using_only_processed_data.Rmd) - An R markdown notebook to generate figures ###-### [requires only `processed data` (either downloaded from GEO or generated with 1)]

