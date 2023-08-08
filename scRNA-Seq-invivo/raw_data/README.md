This directory is meant to hold cell ranger output.
```
$ tree
.
├── 10X_samples_meta_data.tsv
├── Tx18_1_1_Y
│   └── outs
│       └── filtered_feature_bc_matrix
│           ├── barcodes.tsv.gz
│           ├── features.tsv.gz
│           └── matrix.mtx.gz
├── Tx18_1_2_A
│   └── outs
│       └── filtered_feature_bc_matrix
│           ├── barcodes.tsv.gz
│           ├── features.tsv.gz
│           └── matrix.mtx.gz
├── Tx18_2_1_Y
│   └── outs
│       └── filtered_feature_bc_matrix
│           ├── barcodes.tsv.gz
│           ├── features.tsv.gz
│           └── matrix.mtx.gz
├── Tx18_2_2_A
│   └── outs
│       └── filtered_feature_bc_matrix
│           ├── barcodes.tsv.gz
│           ├── features.tsv.gz
│           └── matrix.mtx.gz
├── Tx18_3_1_Y
│   └── outs
│       └── filtered_feature_bc_matrix
│           ├── barcodes.tsv.gz
│           ├── features.tsv.gz
│           └── matrix.mtx.gz
└── Tx18_3_2_A
└── outs
└── filtered_feature_bc_matrix
├── barcodes.tsv.gz
├── features.tsv.gz
└── matrix.mtx.gz
18 directories, 19 files 

```
To populate it, you might either download it from GEO ([GSE190286](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?&acc=GSE190286)) or download the raw resquencing data  ([GSE190286](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?&acc=GSE190286)) and re-run cellranger.
We ran cellranger version 3.0.1 with parameters `--transcriptome=mm10-1.2.0` and `--expect-cells=5000`.

