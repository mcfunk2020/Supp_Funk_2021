library(data.table)

m_dt <- as.data.table(msigdbr::msigdbr(species = "Mus musculus", category = "H") )
inflammaging_gene_sets <- c("HALLMARK_INTERFERON_ALPHA_RESPONSE", "HALLMARK_INTERFERON_GAMMA_RESPONSE",
                            "HALLMARK_ALLOGRAFT_REJECTION", "HALLMARK_TNFA_SIGNALING_VIA_NFKB",
                            "HALLMARK_INFLAMMATORY_RESPONSE", "HALLMARK_IL6_JAK_STAT3_SIGNALING")
m_dt <- m_dt[gs_name%in% inflammaging_gene_sets]

fwrite(m_dt, "inflammaging_signature.tsv")
