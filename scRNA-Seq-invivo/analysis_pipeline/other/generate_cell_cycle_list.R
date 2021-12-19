library(tidyverse)
library(biomaRt)
library(Seurat)

mart1 = useMart("ensembl", dataset="hsapiens_gene_ensembl")
mart2 = useMart("ensembl", dataset="mmusculus_gene_ensembl")


# human / mouse
orthologues <-
  getLDS(attributes=c("external_gene_name"),
         filters="external_gene_name", values=c(Seurat::cc.genes$s.genes, Seurat::cc.genes$g2m.genes), mart=mart1,
         attributesL=c("external_gene_name"), martL=mart2) %>% as_tibble
s.genes <- orthologues %>% filter(Gene.name %in% cc.genes$s.genes) %>% with(Gene.name.1)
g2m.genes <- orthologues %>% filter(Gene.name %in% cc.genes$g2m.genes) %>% with(Gene.name.1)
cell_cycle_genes <- c(s.genes, g2m.genes)
cc.genes_mmusculus <- list(s.genes=s.genes, g2m.genes=g2m.genes)

usethis::use_data(cc.genes_mmusculus, cc.genes_mmusculus, pkg="sir")
