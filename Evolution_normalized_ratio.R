####################################################################
#Function for correlation (and p-values) counting between every pair 
#of gene and lipid in two dfs

#INPUT:
#df_lipids = n lipids columns, m regions  rows
#df_genes = l genes columns, m regions rows

#OUTPUT:
#df with n*2 rows (first n - correlations, second - p-values)
#rows - lipids, columns - genes => (n*2,l)
#rows aren't named 
####################################################################
Genes_lipids_correlation_simple <- function(df_lipids, df_genes){
  
  pvalue_df <- data.frame()
  corr_df <- data.frame()
  
  for (i in 1:dim(df_lipids)[2]){
    cur_lipid <- c()
    pcur_lipid <- c()
    for (j in 1:dim(df_genes)[2]){
      cor_data <- cor.test(df_lipids[,i], df_genes[,j])
      cur_lipid <- c(cur_lipid, cor_data$estimate)
      pcur_lipid <- c(pcur_lipid, cor_data$p.value)
    }
    pvalue_df <- rbind(pvalue_df, pcur_lipid)
    corr_df <- rbind(corr_df, cur_lipid)
  }
  
  colnames(corr_df) <- colnames(df_genes)
  colnames(pvalue_df) <- colnames(df_genes)
  result_table <- rbind(corr_df, pvalue_df)
  return(result_table)
}  

####################################################################
#Function for counting HS ratio (genes or lipids)
#(|h-b| + |h-c| + |h-m|)/(|b-c| + |b-m| + |c-m|) with evolutionary distance between species

#INPUT: 4 df in the following order
#Human_data, Bonobo_data, Chimpanzee_data, Macaque_data

#OUTPOT: ratio df
####################################################################

Ratio_counter <- function(h_df, b_df, c_df, m_df){
  diff1 <- abs(h_df - b_df) + abs(h_df - c_df) + (6/25)*abs(h_df - m_df)  #не делю на три потому что потом все равно беру отношение
  diff2 <- (2/25)*abs(b_df - m_df) + abs(b_df - c_df) + (2/25)*abs(c_df - m_df) 
  return(diff1/diff2)
}


load("Desktop/Khrameeva_2019/Coding/DATA/Basic.RData")

library(ggplot2)
library(reshape2)
library(clusterProfiler)

lipids_ratio <- Ratio_counter(human_lipids, bonobo_lipids, chimp_lipids, macaque_lipids)
genes_ratio <- Ratio_counter(h_genes, b_genes, c_genes, m_genes)

row_names <- paste(as.character(seq(1,33)),  substr(rownames(b_genes),4, nchar(rownames(b_genes))))
rownames(lipids_ratio) <- row_names
rownames(genes_ratio) <- row_names

apply(lipids_ratio, MARGIN = 1, FUN = mean)
apply(genes_ratio, MARGIN = 1, FUN = mean)

#82 значения бесконечности, нужно убрать эти гены 
dim(genes_ratio)[1] * dim(genes_ratio)[2] - sum(apply(genes_ratio, MARGIN = 2, FUN = is.finite))
# а что это за гены и есть ли они среди HS
. <- apply(genes_ratio,  MARGIN = 2, function(x) all(is.finite(x)) == T) 
. <- colnames(genes_ratio[. == F])
intersect(., HS_genes_relaxed_3) # пересечений нет, точно убираем 

genes_ratio_checked <- genes_ratio[, -which(colnames(genes_ratio) %in% .)]

#все починилось!
apply(genes_ratio_checked, MARGIN = 1, FUN = mean)

genes_ratio_checked_log = log(genes_ratio_checked)
lipids_ratio_log = log(lipids_ratio)
ComplexHeatmap::Heatmap(genes_ratio_checked_log, use_raster =F, show_column_names = F,row_split = 3,
                        column_split = 3,  clustering_method_rows = "ward.D2", clustering_method_columns = "ward.D2",
                        row_title = "Brain region", column_title = "Genes", name = "log ratio")

ComplexHeatmap::Heatmap(lipids_ratio_log, use_raster =F, show_column_names = F,row_split = 3,
                        column_split = 3,  clustering_method_rows = "ward.D2", clustering_method_columns = "ward.D2",
                        row_title = "Brain region", column_title = "Lipids", name = "log ratio")

#то же на человеко-специфичных товарищах
ComplexHeatmap::Heatmap(lipids_ratio_log[, which(names(lipids_ratio_log) %in% HS_lipids_medium_2)], use_raster =F, show_column_names = F,row_split = 3,
                        column_split = 3,  clustering_method_rows = "ward.D2", clustering_method_columns = "ward.D2",
                        row_title = "Brain region", column_title = "Lipids", name = "log ratio")

ComplexHeatmap::Heatmap(genes_ratio_checked_log[, which(names(genes_ratio_checked_log) %in% HS_genes_medium_2)], use_raster =F, show_column_names = F,row_split = 3,
                        column_split = 3,  clustering_method_rows = "ward.D2", clustering_method_columns = "ward.D2",
                        row_title = "Brain region", column_title = "Genes", name = "log ratio")

#рейтинг регионов 
region_rating <- data.frame("By_lipids_ratio_mean" = apply(lipids_ratio_log, MARGIN = 1, FUN = mean),
                            "By_genes_ratio_mean" = apply(genes_ratio_checked_log, MARGIN = 1, FUN = mean))

region_rating <- region_rating[order(region_rating$By_lipids_ratio_mean, decreasing = T),]

#рейтинг липидов против рейтинга генов 
ggplot(region_rating, aes(By_lipids_ratio_mean, By_genes_ratio_mean)) + geom_point() + geom_smooth(method='lm', formula= y~x) +
  labs(x = "Mean region lipids log ratio", 
       y = "Mean region genes log ratio", 
       title = paste("Pearson correlation coefficient",
                     round(cor(region_rating$By_lipids_ratio_mean, region_rating$By_genes_ratio_mean),2),
                     sep = ":\n"))

. <- lm( region_rating$By_genes_ratio_mean ~ region_rating$By_lipids_ratio_mean)
summary(.)$r.squared

region_rating <- region_rating[order(region_rating$By_lipids_ratio_mean, decreasing = T),]
region_rating$lipid_num = seq(1:33)
region_rating <- region_rating[order(region_rating$By_genes_ratio_mean, decreasing = T),]
region_rating$genes_num = seq(1:33)

ggplot(region_rating, aes(lipid_num, genes_num)) + geom_point() + geom_smooth(method='lm', formula= y~x) +
  labs(x = "Region Rank by mean lipid ratio", 
       y = "Region Rank by mean gene ratio", 
       title = paste("Spearman correlation coefficient",
                     round(cor(region_rating$lipid_num, region_rating$genes_num),2),
                     sep = ":\n"))


region_rating$By_HSlipids_ratio_mean = 
  apply(lipids_ratio_log[, which(names(lipids_ratio_log) %in% HS_lipids_medium_2)], MARGIN = 1, FUN = mean)

region_rating$By_HSgenes_ratio_mean = 
  apply(genes_ratio_checked_log[, which(names(genes_ratio_checked_log) %in% HS_genes_medium_2)], MARGIN = 1, FUN = mean)


corr.test(region_rating$By_HSlipids_ratio_mean, region_rating$By_HSgenes_ratio_mean)

ggplot(region_rating, aes(By_HSlipids_ratio_mean, By_HSgenes_ratio_mean)) + geom_point() + geom_smooth(method='lm', formula= y~x) +
  labs(x = "Mean region HS lipids log ratio", 
       y = "Mean region HS genes log ratio", 
       title = paste(paste("Pearson correlation coefficient",
                     round(cor(region_rating$By_HSlipids_ratio_mean, region_rating$By_HSgenes_ratio_mean),2),
                     sep = ": "), "P-value: 0.04", sep = "\n")) + geom_abline(slope = 1, intercept = 0, color = "red")


region_rating <- region_rating[order(region_rating$By_HSlipids_ratio_mean, decreasing = T),]
region_rating$HSlipid_num = seq(1:33)
region_rating <- region_rating[order(region_rating$By_HSgenes_ratio_mean, decreasing = T),]
region_rating$HSgenes_num = seq(1:33)

ggplot(region_rating, aes(HSlipid_num, HSgenes_num)) + geom_point() + geom_smooth(method='lm', formula= y~x) +
  labs(x = "Region Rank by mean HS lipid ratio", 
       y = "Region Rank by mean HS gene ratio", 
       title = paste("Spearman correlation coefficient",
                     round(cor(region_rating$HSlipid_num, region_rating$HSgenes_num),2),
                     sep = ":\n"))

corr.test(region_rating$By_HSlipids_ratio_mean, region_rating$By_HSgenes_ratio_mean, method = "spearman")

#финальная таблица с рейтингами 
region_rating_ready_to_use <- data.frame(
  "lipids_ratio_mean" =  rownames(region_rating[order(region_rating$By_lipids_ratio_mean, decreasing = T),]), 
  "genes_ratio_mean" =  rownames(region_rating[order(region_rating$By_genes_ratio_mean, decreasing = T),]),
  "HSlipids_ratio_mean" =  rownames(region_rating[order(region_rating$By_HSlipids_ratio_mean, decreasing = T),]),
  "HSgenes_ratio_mean" =  rownames(region_rating[order(region_rating$By_HSgenes_ratio_mean, decreasing = T),]))


#графики с визуализацией рейтингов 
ggplot(region_rating[order(region_rating$By_lipids_ratio_mean, decreasing = T),],
       aes(x = reorder(rownames(region_rating[order(region_rating$By_lipids_ratio_mean, decreasing = T),]), -By_lipids_ratio_mean),
           y = By_lipids_ratio_mean, fill = By_lipids_ratio_mean)) +
  geom_col() + theme(axis.text.x=element_text(angle=90, hjust=1)) +theme(legend.position = "none") + 
  labs(x = "Brain regions", y = "Mean lipids log(ratio)", 
       title = "Region rating by lipid ratio")

ggplot(region_rating[order(region_rating$By_genes_ratio_mean, decreasing = T),],
       aes(x = reorder(rownames(region_rating[order(region_rating$By_genes_ratio_mean, decreasing = T),]), -By_genes_ratio_mean),
           y = By_genes_ratio_mean, fill = By_genes_ratio_mean)) +
  geom_col() + theme(axis.text.x=element_text(angle=90, hjust=1)) +theme(legend.position = "none") + 
  labs(x = "Brain regions", y = "Mean genes log(ratio)", 
       title = "Region rating by genes ratio")


lipids_ratio_log <- as.data.frame(t(lipids_ratio_log))
genes_ratio_checked_log <- as.data.frame(t(genes_ratio_checked_log))

reshape2::melt(lipids_ratio_log) %>% 
  ggplot(aes(variable, value, fill = variable)) + geom_boxplot(notch = T) + 
  theme(axis.text.x=element_text(angle=90, hjust=1), legend.position = "none") + 
  labs(x = "", y = "Ratio", title = "Distribution of lipid log(ratio) in brain regions")

reshape2::melt(genes_ratio_checked_log) %>% 
  ggplot(aes(variable, value, fill = variable)) + geom_boxplot(notch = T) + 
  theme(axis.text.x=element_text(angle=90, hjust=1), legend.position = "none") + 
  labs(x = "", y = "Ratio", title = "Distribution of genes log(ratio) in brain regions")

#по чловеко-специфичным
genes_ratio_checked_log[, which(names(genes_ratio_checked_log) %in% HS_genes_medium_2)] %>% t %>% 
  melt %>% 
  ggplot(aes(Var2, value, fill = Var2)) + geom_boxplot() + 
  theme(axis.text.x=element_text(angle=90, hjust=1), legend.position = "none") + 
  labs(x = "", y = "Ratio", title = "Distribution of genes log(ratio) in brain regions")


##########################################################################################
#корреляции на основе ratio


lipids_for_simplecor <- dplyr::select(lipids_ratio, HS_lipids_medium_2)
genes_for_simplecor <- dplyr::select(genes_ratio, HS_genes_medium_2)


ratio_corr_m2 <- Genes_lipids_correlation_simple(lipids_for_simplecor,genes_for_simplecor)

p_gl_scor_m2 <- ratio_corr_m2[33:64,]
gl_scor_m2 <- ratio_corr_m2[1:32,]
gl_scor_m2$lipid <- colnames(lipids_for_simplecor)
p_gl_scor_m2$lipid <- colnames(lipids_for_simplecor)

gl_scor_m2_data <- melt(gl_scor_m2, id = "lipid")
p_gl_scor_m2_data <- melt(p_gl_scor_m2, id = "lipid")
colnames(gl_scor_m2_data) <- c("lipid", "gene", "correlation")
gl_scor_m2_data$p.value <- p_gl_scor_m2_data$value
gl_scor_m2_data$p.value <-p.adjust(gl_scor_m2_data$p.value, "BH")

#320 нулей #160 колонок
sum(apply(gl_scor_m2_data, 1, is.na))

gl_scor_m2_data <- na.omit(gl_scor_m2_data)

length(unique(gl_scor_m2_data$gene[gl_scor_m2_data$p.value < 0.05]))
length(gl_scor_m2_data$gene[gl_scor_m2_data$p.value < 0.05])
length(unique(gl_scor_m2_data$lipid[gl_scor_m2_data$p.value < 0.05]))
length(gl_scor_m2_data$gene)

mean(abs(gl_scor_m2_data$correlation[gl_scor_m2_data$p.value < 0.05]))
max(abs(gl_scor_m2_data$correlation[gl_scor_m2_data$p.value < 0.05]))
min(abs(gl_scor_m2_data$correlation[gl_scor_m2_data$p.value < 0.05]))
median(abs(gl_scor_m2_data$correlation[gl_scor_m2_data$p.value < 0.05]))



ggarrange(ggplot(gl_scor_m2_data, aes(correlation)) +
            geom_histogram(fill = "lightblue", color = "black") +
            labs(x = "Pearson correlation coefficient",y = "Count", 
                 title = "Correlation coefficient between HS-lipids and HS-genes (Medium 2)"),
          ggplot(gl_scor_m2_data[gl_scor_m2_data$p.value < 0.05,], aes(correlation)) +
            geom_histogram(fill = "pink", color = "black")+ 
            labs(x = "Pearson correlation coefficient",y = "Count", 
                 title = "Correlation coefficient between HS-lipids and HS-genes (Medium 2)
                 with p-values < 0.05"), nrow = 2)

gl_scor_m2_data$lipid <-  factor(gl_scor_m2_data$lipid,
                                 levels = c("HexCer(d41:0)", "HexCer(d42:0)", "HexCer(d43:2)","HexCer(t40:2)", "HexCer(t41:2)", "HexCer(t42:2)",
                                            "TG(58:8)", "PG(38:1)", "Cer(m40:1)", "PE(36:2)", "PE(P-36:5)", "PE(P-36:6)",
                                            "MG(20:3)", "MG(18:2)", "MG(16:2)", "MG(22:3)", "PC(34:3)", "PC(36:4)", 
                                            "PC(40:8)", "PC(42:10)", "SulfoHexCer(d41:1)","SulfoHexCer(d42:1)","SulfoHexCer(d43:1)",
                                            "SulfoHexCer(d43:2)","SulfoHexCer(t40:1)","SulfoHexCer(t41:1)","SulfoHexCer(t41:2)",
                                            "SulfoHexCer(t42:1)","SulfoHexCer(t43:1)","SulfoHexCer(t43:2)"))

gl_scor_m2_data$lipid_class <- unlist(strsplit(gl_scor_m2_data$lipid, split = "(",fixed = T))[seq(1, 107136,2)]

#гистограмма для каждого липида
ggplot(gl_scor_m2_data[gl_scor_m2_data$p.value < 0.05,], aes(correlation, color = lipid_class)) + 
  geom_freqpoly() + facet_wrap(~lipid) +  theme(legend.position = "none") + 
  labs(x = "Correlation coefficient", y = "Number of correlation coefficient", 
       title = "Distribution of significant correlation coefficients between  genes and  lipids ratio by lipid")


ggplot(gl_scor_m2_data[gl_scor_m2_data$p.value < 0.05,], aes(correlation, fill = lipid_class)) + 
  geom_histogram(binwidth = 0.02) + facet_wrap(~lipid_class) + theme(legend.position = "none") +
  labs(x = "", y = "Correlation coefficient", title = "Distribution of ratio correlation coefficient by lipid class")

#вайолин плоты с боксплотами внутри по классам
ggplot(gl_scor_m2_data, aes(x = lipid_class, y = correlation, fill = lipid_class)) + 
  geom_violin(width = 1.4)+ theme(legend.position = "none") + geom_boxplot(alpha=0.2, width = 0.1, fill = "gray",
                                                                           outlier.size = 3, color = "blue") +
  labs(x = "", y = "Correlation coefficient", title = "Distribution of ratio correlation coefficient by lipid class")

ggplot(gl_scor_m2_data[gl_scor_m2_data$p.value < 0.05,], aes(x = lipid_class, y = correlation, fill = lipid_class)) + 
  geom_violin(width = 1.4)+ theme(legend.position = "none") + geom_boxplot(alpha=0.2, width = 0.1, fill = "gray",
                                                                           outlier.size = 3, color = "blue") +
  labs(x = "", y = "Correlation coefficient", title = "Distribution of ratio significant correlation coefficient by lipid class")


#количество коэффициентов корреляций по классам
ggplot(gl_scor_m2_data[gl_scor_m2_data$p.value < 0.05,], aes(lipid_class, fill = lipid_class)) + geom_bar() +
  labs(title = "Number of ratio correlation coefficient with p-value < 0.05", 
       x = "", y = "Number of correlation coefficient") + theme(axis.text.x=element_text(angle=45, hjust=1),
                                                                legend.position = "none") 


#отдельно для каждого липида 
ggplot(gl_scor_m2_data[gl_scor_m2_data$p.value < 0.05,], aes(lipid, fill = lipid_class)) + geom_bar() +
  labs(title = "Number of correlation coefficient with p-value < 0.05 for each lipid", x = "", 
       y = "Number of correlation coefficient") +
  theme(axis.text.x=element_text(angle=60, hjust=1), legend.position = "none") 

ggplot(gl_scor_m2_data[gl_scor_m2_data$p.value < 0.05,], aes(x = lipid, y = abs(correlation), fill = lipid_class)) +
  geom_boxplot() +
  labs(title = "Correlation coefficients with p-value < 0.05 for each lipid", x = "", y = "") +
  theme(axis.text.x=element_text(angle=60, hjust=1), legend.position = "none") 

rownames(gl_scor_m2) = gl_scor_m2$lipid
ComplexHeatmap::Heatmap(t(na.omit(t(gl_scor_m2[,1:1679]))), show_column_names = F, row_split = 3,
                        column_split = 3,  clustering_method_rows = "ward.D2", clustering_method_columns = "ward.D2", 
                        name = "Pearson\ncorrelation\ncoefficient", column_title = "HS genes", 
                        row_title = "HS lipids")


lipids_genes <- read.csv("/Users/ann/Desktop/Khrameeva_2019/Coding/DATA/genes_associated_with_lipids/lipid_enzymes.txt", sep = "\t", header = F)
lipids_genes <- bitr(lipids_genes$V1, fromType = "ENTREZID", toType = "ENSEMBL", OrgDb="org.Hs.eg.db")

length(intersect(lipids_genes$ENSEMBL,gl_scor_m2_data$gene[gl_scor_m2_data$p.value < 0.05]))


winners <- intersect(lipids_genes$ENSEMBL,unique(gl_scor_m2_data$gene[gl_scor_m2_data$p.value < 0.05]))
lipids_genes_data <- gl_scor_m2_data[gl_scor_m2_data$gene %in% winners & gl_scor_m2_data$p.value < 0.05,]
lipids_genes_data <- lipids_genes_data[order(lipids_genes_data$correlation, decreasing = T),]

dim(lipids_genes_data)
unique(lipids_genes_data$lipid)
range(lipids_genes_data$correlation)


ggplot(lipids_genes_data, aes(lipid_class, fill = lipid_class)) + geom_bar() + 
  labs(x = "", y = "Number of sufficient correlation", title = "Number of lipid ratio correlation coefficient with lipid genes in lipid classes") +
  theme(legend.position = "none")

ggplot(lipids_genes_data, aes(lipid_class, correlation, fill = lipid_class)) + geom_violin(width = 1) + 
  geom_boxplot(varwidth = T, alpha = 0.3) + theme(legend.position = "none") + 
  labs(x = "", title = "Correlation coefficients of ratio of lipids correlating with ratio of lipid genes")

lipids_genes_data$lipid <- factor(lipids_genes_data$lipid , levels = unique(lipids_genes_data$lipid))

ggplot(lipids_genes_data, aes(lipid, correlation, color = lipid_class)) + geom_count(size = 6) +
  theme(axis.text.x=element_text(angle=60, hjust=1)) +
  labs(x = "", y = "Correlation coefficient",
       title = "Correlation coefficients of ratio of lipids correlating with ratio of lipid genes")


#посмотрим на GO-terms генов, связанных с липидами при помощи биомарта 
library(biomaRt)
mart <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
results <- getBM(attributes = c("ensembl_gene_id", "go_id","name_1006"), filters = c("ensembl_gene_id"),
                 values = winners, mart = mart)

listAttributes(mart)
listFilters(mart)

for (i in lipids_genes_data$gene[10 : 23]){
  print("-------------------------------------------------------------------------------------------")
  print(paste("Gene", i))
  print("-------------------------------------------------------------------------------------------")
  print(results[results$ensembl_gene_id == i,]
)
}

#а теперь на гены с супер-высокими корреляциями (можно сделать энричмент анализ!)

#набор генов с высокими значимыми корреляциями 
ego <- enrichGO(gene         = unique(gl_scor_m2_data$gene[gl_scor_m2_data$p.value < 0.05 & gl_scor_m2_data$correlation > 0.7]),
                     OrgDb         = org.Hs.eg.db,
                     keyType       = 'ENSEMBL',
                     universe      = colnames(h_genes),
                     ont           = "BP",
                     pAdjustMethod = "BH",
                     pvalueCutoff  = 0.05,
                     qvalueCutoff  = 0.05)
dotplot(ego, showCategory=30) + ggtitle( "Enriched GO BP terms in set of genes with ratio\ncorrelation coefficients with lipids > 0.7")
emapplot(ego)
length(unique(gl_scor_m2_data$gene[gl_scor_m2_data$p.value < 0.05 & gl_scor_m2_data$correlation > 0.7]))

#набор генов со значимыми корреляциями
ego_2 <- enrichGO(gene         = unique(gl_scor_m2_data$gene[gl_scor_m2_data$p.value < 0.05]),
                OrgDb         = org.Hs.eg.db,
                keyType       = 'ENSEMBL',
                universe      = colnames(h_genes),
                ont           = "BP",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.05,
                qvalueCutoff  = 0.05)

dotplot(ego_2, showCategory=30) + ggtitle("Enriched GO BP terms in set of genes with sugnificant ratio\ncorrelation coefficients with lipids")
length(unique(gl_scor_m2_data$gene[gl_scor_m2_data$p.value < 0.05]))

#пусто...
ego_3 <- enrichGO(gene         = unique(gl_scor_m2_data$gene[gl_scor_m2_data$p.value < 0.05 & gl_scor_m2_data$correlation > 0.85]),
                  OrgDb         = org.Hs.eg.db,
                  keyType       = 'ENSEMBL',
                  universe      = colnames(h_genes),
                  ont           = "BP",
                  pAdjustMethod = "BH",
                  pvalueCutoff  = 0.05,
                  qvalueCutoff  = 0.05)

dotplot(ego_3, showCategory=30) 

for_kegg <- bitr(unique(gl_scor_m2_data$gene[gl_scor_m2_data$p.value < 0.05 & gl_scor_m2_data$correlation > 0.7]), fromType = "ENSEMBL", toType = "ENTREZID",
                 OrgDb="org.Hs.eg.db")

#кег по нулям для обоих случаев....
kk <- enrichKEGG(gene         = for_kegg$ENTREZID,
                 organism     = 'hsa',
                 pvalueCutoff = 0.05,
                 universe      = colnames(h_genes))



gl_scor_m2_data <- gl_scor_m2_data[order(gl_scor_m2_data$correlation, decreasing = T),]
View(gl_scor_m2_data[gl_scor_m2_data$p.value < 0.05,])

. <- gl_scor_m2_data$gene[gl_scor_m2_data$p.value < 0.05]
. <- .[1:23]

results <- getBM(attributes = c("ensembl_gene_id", "go_id","name_1006"), filters = c("ensembl_gene_id"),
                 values = ., mart = mart)
for (i in .){
  print("-------------------------------------------------------------------------------------------")
  print(paste("Gene", i))
  print("-------------------------------------------------------------------------------------------")
  print(results[results$ensembl_gene_id == i,]
  )
}

