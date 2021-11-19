# Analysis of single cell RNA sequencing data of Funk et al. 2021
This analysis is separated into three independent parts:
1. [analysis_pipeline](analysis_pipeline) - An R project (including a makefile and a helper package) that allows reproduction of [`processed_data`](processed_data) from the cell ranger outputs available on GEO (####) 
  <br>[requires cell ranger output in `raw_data`, reproduces [`processed data`](processed_data)]
3.  [scRNAseq_figures_using_raw_data.Rmd](scRNAseq_figures_using_raw_data.Rmd) - An R markdown notebook to generate figures ###-### <br>
  [requires cell ranger output in `raw data` and [`processed data`](processed_data)]
5.  [scRNAseq_figures_using_only_processed_data.Rmd](scRNAseq_figures_using_only_processed_data.Rmd) - An R markdown notebook to generate figures ###-### <br> 
    [requires only [`processed data`](processed_data)]

