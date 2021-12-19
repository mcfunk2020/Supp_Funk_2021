library(data.table)
library(schelpr)

fread("https://docs.google.com/spreadsheets/u/1/d/1DTv_x4kSwWJtnq_0DtkWxTnwSYA-s8Vg8chjC0hZgv4/export?format=csv&id=1DTv_x4kSwWJtnq_0DtkWxTnwSYA-s8Vg8chjC0hZgv4&gid=1799628969") %>%
  fwrite("SI_cell_type_markers.tsv", sep="\t")
fread("https://docs.google.com/spreadsheets/u/1/d/1DTv_x4kSwWJtnq_0DtkWxTnwSYA-s8Vg8chjC0hZgv4/export?format=tsv&id=1DTv_x4kSwWJtnq_0DtkWxTnwSYA-s8Vg8chjC0hZgv4&gid=1055237546") %>%
  fwrite("SI_feature_gene_lists.tsv", sep="\t")
fread("https://docs.google.com/spreadsheets/d/1AQvsALrlfY8aEZS-446gLSIpNFITZ-sVHcacyTexGa4/export?format=tsv&id=1AQvsALrlfY8aEZS-446gLSIpNFITZ-sVHcacyTexGa4&gid=0") %>%
  fwrite("Sample meta - Sheet1.tsv", sep="\t")
