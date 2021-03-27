colnames(gl_scor_m2) <-  colnames(lipids_for_simplecor)


library(cluster) #agnes
library(ggplot2)
library(ggpubr)
library(dendextend) #entanglement, tanglegram
library(clusterProfiler)
library(ComplexHeatmap)
library(org.Hs.eg.db)
library(enrichplot)


gl_scor_m2 <- as.data.frame(t(gl_scor_m2[-c(1680)]))


hc <- agnes(gl_scor_m2, method = "ward")
hc$ac

pltree(hc)

hc_comp <- agnes(gl_scor_m2, method = "complete")
hc_comp$ac
pltree(hc_comp)

ggarrange(pltree(hc_comp), pltree(hc))

dend1 <- as.dendrogram (hc)
dend2 <- as.dendrogram (hc_comp)
tanglegram(dend1, dend2) #супер долго

entanglement(dend1, dend2)

d <- dist(gl_scor_m2, method = "euclidean")

dend3 <- as.dendrogram (hclust(d))
dend4 <- as.dendrogram (hclust(d,method = "ward.D"))

entanglement(dend3, dend4)
entanglement(dend1, dend3)


plot(hclust(d))
plot(hclust(d, method = "ward.D"))
?hclust

my_clust <- hclust(d, method = "ward.D2")
genes_clusters <- cutree(my_clust, k = 2)

Heatmap(t(gl_scor_m2), column_split = genes_clusters, show_column_names = F, row_title = "HS lipids", 
        column_title = "HS genes", heatmap_legend_param = list(title = "correlation coefficient"), 
        column_gap = unit(3, "mm"))

length(genes_clusters[genes_clusters == 1])
length(genes_clusters[genes_clusters == 2])

colnames(gl_scor_m2)
ego_CC <- enrichGO(gene         = names(genes_clusters[genes_clusters == 1]),
                 OrgDb         = org.Hs.eg.db,
                 keyType       = 'ENSEMBL',
                 ont           = "CC",
                 pAdjustMethod = "BH",
                 pvalueCutoff  = 0.05,
                 qvalueCutoff  = 0.05)

barplot(ego_CC) + ggtitle("Gene cluster 1 CC ontology")
dotplot(ego_CC, x="count") + ggtitle("Gene cluster 1 CC ontology")

#для набора 1 нет обогащенных BP на катоф 0.1
#для MF получается только cadherin binging

ego_2 <- enrichGO(gene         = names(genes_clusters[genes_clusters == 2]),
                OrgDb         = org.Hs.eg.db,
                keyType       = 'ENSEMBL',
                ont           = "BP",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.05,
                qvalueCutoff  = 0.05)

ego_2_СС <- enrichGO(gene         = names(genes_clusters[genes_clusters == 2]),
                  OrgDb         = org.Hs.eg.db,
                  keyType       = 'ENSEMBL',
                  ont           = "CC",
                  pAdjustMethod = "BH",
                  pvalueCutoff  = 0.05,
                  qvalueCutoff  = 0.05)

ego_2_MF <- enrichGO(gene         = names(genes_clusters[genes_clusters == 2]),
                     OrgDb         = org.Hs.eg.db,
                     keyType       = 'ENSEMBL',
                     ont           = "MF",
                     pAdjustMethod = "BH",
                     pvalueCutoff  = 0.05,
                     qvalueCutoff  = 0.05)
barplot(ego_2)
dotplot(ego_2_MF) + ggtitle("Gene cluster 2 MF ontology")
ego

rm(enriched, figure, g_0.75, bp, c, genes_for_six_func)

kk_2 <- enrichKEGG(gene         = ENTREZID_2$ENTREZID,
                 organism     = 'hsa',
                 pvalueCutoff = 0.05)

ENTREZID_2 <- bitr(geneID = names(genes_clusters[genes_clusters == 2]), fromType = "ENSEMBL", toType = "ENTREZID",
                   OrgDb="org.Hs.eg.db")

barplot(kk_2)
dotplot(kk_2) + ggtitle("Gene cluster 2 enriched KEGG pathways")

ENTREZID_1 <- bitr(geneID = names(genes_clusters[genes_clusters == 1]), fromType = "ENSEMBL", toType = "ENTREZID",
                   OrgDb="org.Hs.eg.db")


kk_1 <- enrichKEGG(gene         = ENTREZID_1$ENTREZID,
                   organism     = 'hsa',
                   pvalueCutoff = 0.05)
barplot(kk_1)

cnetplot(kk_2, circular = TRUE, colorEdge = TRUE)
cnetplot(kk_2, node_label="category")
emapplot(kk_2)
emapplot(ego_2)

allgenes <- colnames(h_genes)

ego_2_bg <- enrichGO(gene         = names(genes_clusters[genes_clusters == 2]),
                  OrgDb         = org.Hs.eg.db,
                  keyType       = 'ENSEMBL',
                  ont           = "ALL",
                  pAdjustMethod = "BH",
                  pvalueCutoff  = 0.05,
                  qvalueCutoff  = 0.05, 
                  universe = allgenes)

ego_2_BP_bg <- enrichGO(gene         = names(genes_clusters[genes_clusters == 2]),
                     OrgDb         = org.Hs.eg.db,
                     keyType       = 'ENSEMBL',
                     ont           = "BP",
                     pAdjustMethod = "BH",
                     pvalueCutoff  = 0.05,
                     qvalueCutoff  = 0.05,
                     universe = allgenes)

ego_2_СС_bg <- enrichGO(gene         = names(genes_clusters[genes_clusters == 2]),
                        OrgDb         = org.Hs.eg.db,
                        keyType       = 'ENSEMBL',
                        ont           = "CC",
                        pAdjustMethod = "BH",
                        pvalueCutoff  = 0.05,
                        qvalueCutoff  = 0.05,
                        universe = allgenes)

ego_2_MF_bg <- enrichGO(gene         = names(genes_clusters[genes_clusters == 2]),
                     OrgDb         = org.Hs.eg.db,
                     keyType       = 'ENSEMBL',
                     ont           = "MF",
                     pAdjustMethod = "BH",
                     pvalueCutoff  = 0.05,
                     qvalueCutoff  = 0.05,
                     universe = allgenes)

dotplot(ego_2_MF_bg, showCategory=20)


ego_CC_bg <- enrichGO(gene         = names(genes_clusters[genes_clusters == 1]),
                   OrgDb         = org.Hs.eg.db,
                   keyType       = 'ENSEMBL',
                   ont           = "CC",
                   pAdjustMethod = "BH",
                   pvalueCutoff  = 0.05,
                   qvalueCutoff  = 0.05,
                   universe = allgenes)


ego_MF_bg <- enrichGO(gene         = names(genes_clusters[genes_clusters == 1]),
                      OrgDb         = org.Hs.eg.db,
                      keyType       = 'ENSEMBL',
                      ont           = "MF",
                      pAdjustMethod = "BH",
                      pvalueCutoff  = 0.05,
                      qvalueCutoff  = 0.05,
                      universe = allgenes)

ego_BP_bg <- enrichGO(gene         = names(genes_clusters[genes_clusters == 1]),
                      OrgDb         = org.Hs.eg.db,
                      keyType       = 'ENSEMBL',
                      ont           = "BP",
                      pAdjustMethod = "BH",
                      pvalueCutoff  = 0.05,
                      qvalueCutoff  = 0.05,
                      universe = allgenes)


dotplot(ego_CC_bg, showCategory=20) +  ggtitle("Gene cluster 1 CC ontology")


ENTREZID_all <- bitr(geneID = allgenes, fromType = "ENSEMBL", toType = "ENTREZID",
                   OrgDb="org.Hs.eg.db")


kk_1_bg <- enrichKEGG(gene         = ENTREZID_1$ENTREZID,
                   organism     = 'hsa',
                   pvalueCutoff = 0.05,
                   universe = ENTREZID_all$ENTREZID)

kk_2_bg <- enrichKEGG(gene         = ENTREZID_2$ENTREZID,
                   organism     = 'hsa',
                   pvalueCutoff = 0.05,
                   universe = ENTREZID_all$ENTREZID)

dotplot(ego_2_bg, showCategory=20) +  ggtitle("Gene cluster 2 enriched KEGG pathways")


dotplot(kk_2_bg, showCategory=20)
cnetplot(kk_2, circular = TRUE, colorEdge = TRUE,  showCategory=10)
cnetplot(kk_2,  showCategory=5)

emapplot(kk_2_bg)
emapplot(ego_2_СС_bg)







rm(g_0.5, g_0.55, g_0.6, g_0.65, g_0.7, g_0.75, g_0.8, g_ENSG00000170271, ggene_1, ggene_10, 
   ggene_11, ggene_12, ggene_13, ggene_14, ggene_15, ggene_2, ggene_3, ggene_4, ggene_5, ggene_6, ggene_7,
   ggene_8, ggene_9)
rm(bp, c,enriched, figure, genes_for_six_func, l, lipids_for_six_func, m, names_table, normalazid_freq_0.65,
   psize_table, size_table, thresh_65, thresh_65_vs_all, winners_df, winners_list, x)
rm(Column_renamer, Greper_b, Greper_c, Greper_h, Greper_m, Human_specific_features, Multiple_Comparisons, 
   Substringer, try_func)
rm(winners, z, uniqe_genes, unique_lipids,threshold, sizes, size, p, i, lipids_genes)
rm(winers,lipid_genes_ENSEMBLE,goal, dbs,counter, correlation_num, corr_genes)
rm(lipid_genes_th)
rm(Genes_lipids_correlation_simple)
rm(p_gl_scor_m1, p_gl_scor_m1_data, p_gl_scor_m1_vs_all, p_gl_scor_m1_vs_all_data)
rm(lipids, lipids_for_simplecor, heatmap_matrix, gl_scor_m1, gl_scor_m1_data, gl_scor_m1_vs_all, 
   gl_scor_m1_vs_all_data)
rm(genes, genes_for_simplecor, genes_lipids_simplecor_m1, genes_lipids_simplecor_m1_vs_all, genes_lipids_simplecor_m2)


library(psych)
pairs.panels(human_lipids[,c(39,332,62,225)],
             hist.col = "lightblue")

