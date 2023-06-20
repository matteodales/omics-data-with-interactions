library('clusterProfiler')
library('org.Hs.eg.db')
library('ReactomePA')
library('gson')



## This script generates enrichment analysis results (Figures 13 and 14)



library("msigdbr")
msig_h <- msigdbr(species = "Homo sapiens", category = "H") %>%
  dplyr::select(gs_name, entrez_gene) %>%
  dplyr::rename(ont = gs_name, gene = entrez_gene)
msig_h


# selected_betas_coop contains the number of times each feature
# was included in the model for the 10 iterations


counts_beta = selected_betas_coop_axi
#counts_beta = selected_betas_coop_nilo
gene_list_tot = names(counts_beta)
id_list_tot = bitr(gene_list_tot, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")

min_counts = 2
gene_list = names(counts_beta[counts_beta>=min_counts])
id_list = bitr(gene_list, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")

enrich_hallmark_results = enricher(gene = id_list$ENTREZID, universe = id_list_tot$ENTREZID, TERM2GENE = msig_h, pvalueCutoff = 1000)
enrichplot::dotplot(enrich_hallmark_results)



