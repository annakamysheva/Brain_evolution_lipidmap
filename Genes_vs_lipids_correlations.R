####################################################################
#Function for correlation (and p-values) counting between every pair 
#of gene and lipid in two dfs

#INPUT:
#df_lipids = n lipids columns, m regions (originally 132 h + b + c + m) rows
#df_genes = l genes columns, m regions (originally 132 h + b + c + m) rows

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
  


###################################################################################################################
#Анализ содержимого таблицы с генами, функционально связанными с липидами (lipid_enzymes.txt)
####################################################################################################################
library(clusterProfiler)

#наши данные в ENSEMBLE, таблица в ENTREZID - переведем и посмотрим есть ли пересечения
x <- bitr(HS_genes_medium_2, fromType = "ENSEMBL", toType = "ENTREZID", OrgDb="org.Hs.eg.db")

#нашлись все гены, но почему-то есть дупликаты????
setdiff(x$ENSEMBL, HS_genes_medium_2)
x$ENSEMBL[duplicated(x$ENSEMBL)]

#посмотрим, есть ли в моем списке вообще гены, связанные с липидами (табличка от Екатерины)
lipids_genes <- read.csv("Desktop/Khrameeva_2019/Coding/DATA/genes_associated_with_lipids/lipid_enzymes.txt", sep = "\t", header = F)


#всего 0.1% генов по синтезу липидов оказались HS M2... не густо 
length(intersect(x$ENTREZID, lipids_genes$V1))/length(lipids_genes$V1)

#а если не M2? - сделаем то же самое для всех
# до 0.2% ID не перевелось - не очень много, 2 гена из каждого Relaxed (вектор loss)

HS_ENTREZID <- list()
lipid_percent <- c()
loss <- c()
num_HS <- c()
num_HS_and_lipids <- c()
for (i in colnames(HS_genes_data)){
  x <- na.omit(HS_genes_data[i])
  num_HS[i] <- dim(x)[1]
  y <- bitr(unlist(x), fromType = "ENSEMBL", toType = "ENTREZID", OrgDb="org.Hs.eg.db")
  HS_ENTREZID[[i]]  <- y$ENTREZID
  loss[i] = dim(x)[1] - length(unique(y$ENSEMBL))
  lipid_percent[i] = length(intersect(y$ENTREZID, lipids_genes$V1))/length(lipids_genes$V1) 
  num_HS_and_lipids[i] = length(intersect(y$ENTREZID, lipids_genes$V1))
}

# 0.12*0.01 * length(rownames(genes))  =  13 генов не перевелись 
y <- bitr(rownames(genes), fromType = "ENSEMBL", toType = "ENTREZID", OrgDb="org.Hs.eg.db")
HS_ENTREZID[["All"]]  <- y$ENTREZID
loss["All"] = length(rownames(genes)) - length(unique(y$ENSEMBL))
lipid_percent["All"] = length(intersect(y$ENTREZID, lipids_genes$V1))/length(lipids_genes$V1)
#в нашем списке всего 63% липидных генов

#75 дупликатов??? (0.007)
length(y$ENSEMBL[duplicated(y$ENSEMBL)])

lipid_percent_df <- as.data.frame(lipid_percent)
lipid_percent_df$Metric_type <- rownames(lipid_percent_df)

lipid_percent_df$Metric_type <- factor(lipid_percent_df$Metric_type, 
                                       levels = c("Strict_1" ,"Strict_2",  "Strict_3", 
                                                  "Medium_1", "Medium_2",  "Medium_3",
                                                  "Relaxed_1", "Relaxed_2", "Relaxed_3","All"))

lipid_percent_df$lipid_percent <-  round(lipid_percent_df$lipid_percent,2)

library(ggplot2)
# график доли HS среди липидных генов 
g1 <- ggplot(lipid_percent_df, aes(x = Metric_type, y = lipid_percent, fill = Metric_type)) +
  geom_col() + theme(legend.position = "none") + 
  geom_text(aes(label=round(lipid_percent*909)), position=position_dodge(width=0.9), vjust=-0.25) + 
  labs(x = "", y = "Share", title = "Share of HS genes in genes assosiated with lipid metabolism")


df2 <- lipid_percent_df
colnames(df2) <- c("Share", colnames(df2)[2])
num_HS["All"] = 11176
num_HS_and_lipids["All"] <-  length(intersect(y$ENTREZID, lipids_genes$V1))
df2$Number_of_HS <- num_HS
df2$Number_of_HS_and_lipids <- num_HS_and_lipids
#df2$Number_of_HS_and_lipids <- df2$Number_of_HS - num_HS_and_lipids
df2$Number_of_HS <- df2$Number_of_HS - df2$Number_of_HS_and_lipids 

df2.0 <- reshape2::melt(df2[1:nrow(df2)-1,],id.vars = colnames(df2)[c(1,2)])
df2.0$variable <- factor(df2.0$variable, levels = levels(df2.0$variable)[seq(2,1,-1)])

# график доли липидных генов среди HS  
ggplot(df2.0, aes(x = Metric_type, y = value, fill = variable)) +
  geom_col(position = "stack") + theme(legend.position = "bottom") +
  labs(x = "", y = "Number of Genes", title = "Share of genes assosiated with lipid metabolism in HS genes") +
  scale_fill_discrete(name = "Involvment in lipid metabolism", labels = c("Gene Involved", "Gene not involved"))

#перейдем к корреляциям 
####################################################################
#Data preparation
####################################################################
lipids_for_simplecor <- rbind(human_lipids, bonobo_lipids, chimp_lipids, macaque_lipids)
genes_for_simplecor <- rbind(h_genes, b_genes, c_genes, m_genes)

library(dplyr)
lipids_for_simplecor <- select(lipids_for_simplecor, HS_lipids_medium_2) 
genes_for_simplecor <- select(genes_for_simplecor, HS_genes_medium_2)

genes_lipids_simplecor_m2 <- Genes_lipids_correlation_simple(lipids_for_simplecor,genes_for_simplecor)
genes_lipids_simplecor_m1 <- Genes_lipids_correlation_simple(lipids_for_simplecor,genes_for_simplecor)

p_gl_scor_m1 <- genes_lipids_simplecor_m1[27:52,]
gl_scor_m1 <- genes_lipids_simplecor_m1[1:26,]
gl_scor_m1$lipid <- colnames(lipids_for_simplecor)
p_gl_scor_m1$lipid <- colnames(lipids_for_simplecor)

#для m2
p_gl_scor_m2 <- genes_lipids_simplecor_m2[33:64,]
gl_scor_m2 <- genes_lipids_simplecor_m2[1:32,]
gl_scor_m2$lipid <- colnames(lipids_for_simplecor)
p_gl_scor_m2$lipid <- colnames(lipids_for_simplecor)



library(reshape2)
gl_scor_m1_data <- melt(gl_scor_m1, id = "lipid")
p_gl_scor_m1_data <- melt(p_gl_scor_m1, id = "lipid")
colnames(gl_scor_m1_data) <- c("lipid", "gene", "correlation")
gl_scor_m1_data$p.value <- p_gl_scor_m1_data$value
gl_scor_m1_data$p.value <-p.adjust(gl_scor_m1_data$p.value, "BH")

#для m2
gl_scor_m2_data <- melt(gl_scor_m2, id = "lipid")
p_gl_scor_m2_data <- melt(p_gl_scor_m2, id = "lipid")
colnames(gl_scor_m2_data) <- c("lipid", "gene", "correlation")
gl_scor_m2_data$p.value <- p_gl_scor_m2_data$value
gl_scor_m2_data$p.value <-p.adjust(gl_scor_m2_data$p.value, "BH")


####################################################################
#Визуализация
####################################################################
library(ggplot2)
library(ggpubr)
#гистограммы распрделения коэффициентов корреляций (p-value < 0.05 и все)
ggarrange(ggplot(gl_scor_m2_data, aes(correlation)) +
  geom_histogram(fill = "lightblue", color = "black") +
  labs(x = "Pearson correlation coefficient",y = "Count", 
       title = "Correlation coefficient between HS-lipids and HS-genes (Medium 2)"),
  ggplot(gl_scor_m2_data[gl_scor_m2_data$p.value < 0.05,], aes(correlation)) +
  geom_histogram(fill = "pink", color = "black")+ 
  labs(x = "Pearson correlation coefficient",y = "Count", 
       title = "Correlation coefficient between HS-lipids and HS-genes (Medium 2)
with p-values < 0.05"), nrow = 2)

#денсити плот с отдельными липидами
z <-levels(gl_scor_m2_data$lipid)
gl_scor_m2_data$lipid <-  factor(gl_scor_m2_data$lipid, levels = c(z[3:8], z[1:2], z[32], z[9:31] ))
z <-levels(gl_scor_m2_data$lipid)
gl_scor_m2_data$lipid <-  factor(gl_scor_m2_data$lipid, levels = c(z[1:6], z[9], z[7:8], z[21], z[10:20], z[22:32]))

ggplot(gl_scor_m2_data[gl_scor_m2_data$p.value < 0.05,], aes(correlation, fill = lipid_class)) + 
  geom_density() + facet_wrap(~lipid, ncol = 7) +  theme(legend.position = "none") + 
  labs(x = "Correlation coefficient", y = "Number of correlation coefficient", 
       title = "Distribution of correlation coefficient between HS genes and particular HS lipid (p-value < 0.05)")

#добавим инфо про класс липида 
gl_scor_m2_data$lipid_class <- unlist(strsplit(gl_scor_m2_data$lipid, split = "(",fixed = T))[seq(1,107456,2)]

ggplot(gl_scor_m2_data[gl_scor_m2_data$p.value < 0.05,], aes(correlation, fill = lipid_class)) + 
  geom_histogram(binwidth = 0.05) + facet_wrap(~lipid_class) + theme(legend.position = "none") +
  labs(x = "", y = "Correlation coefficient", title = "Distribution of correlation coefficient by lipid class")

#вайолин плоты с боксплотами внутри по классам
ggplot(gl_scor_m2_data, aes(x = lipid_class, y = correlation, fill = lipid_class)) + 
  geom_violin(width = 1.4)+ theme(legend.position = "none") + geom_boxplot(alpha=0.2, width = 0.1, fill = "gray",
                                                                           outlier.size = 3, color = "blue") +
  labs(x = "", y = "Correlation coefficient", title = "Distribution of correlation coefficient by lipid class")

#количество коэффициентов корреляций по классам
ggplot(gl_scor_m2_data[gl_scor_m2_data$p.value < 0.05,], aes(lipid_class, fill = lipid_class)) + geom_bar() +
  labs(title = "Number of correlation coefficient with p-value < 0.05", 
       x = "", y = "Number of correlation coefficient") + theme(axis.text.x=element_text(angle=45, hjust=1),
                                                                legend.position = "none") 
  

#отдельно для каждого липида 
ggplot(gl_scor_m2_data[gl_scor_m2_data$p.value < 0.05,], aes(lipid, fill = lipid_class)) + geom_bar() +
  labs(title = "Number of correlation coefficient with p-value < 0.05 for each lipid", x = "", 
       y = "Number of correlation coefficient") +
  theme(axis.text.x=element_text(angle=60, hjust=1), legend.position = "none") 

ggplot(gl_scor_m2_data[gl_scor_m2_data$p.value < 0.05,], aes(x = lipid, y = abs(correlation), fill = lipid_class)) +
  geom_violin() +
  labs(title = "Correlation coefficients with p-value < 0.05 for each lipid", x = "", y = "") +
  theme(axis.text.x=element_text(angle=60, hjust=1), legend.position = "none") 

length(unique(gl_scor_m2_data$gene[gl_scor_m2_data$p.value < 0.05]))
length(gl_scor_m2_data$gene[gl_scor_m2_data$p.value < 0.05])
length(unique(gl_scor_m2_data$lipid[gl_scor_m2_data$p.value < 0.05]))
length(gl_scor_m2_data$gene)

mean(abs(gl_scor_m2_data$correlation[gl_scor_m2_data$p.value < 0.05]))
max(abs(gl_scor_m2_data$correlation[gl_scor_m2_data$p.value < 0.05]))
min(abs(gl_scor_m2_data$correlation[gl_scor_m2_data$p.value < 0.05]))

threshold = 0.65

thresh_65 <- gl_scor_m2_data[abs(gl_scor_m2_data$correlation) > threshold,]
dim(gl_scor_m2_data[abs(gl_scor_m1_data$correlation) > threshold,])
length(unique(gl_scor_m2_data$gene[abs(gl_scor_m2_data$correlation) > threshold]))
length(unique(gl_scor_m2_data$lipid[abs(gl_scor_m2_data$correlation) > threshold]))

#гены, перешагнувшие через порог
winers <- unique(gl_scor_m2_data$gene[abs(gl_scor_m2_data$correlation) > threshold])
#переведем табличку в ENSEMBL
lipid_genes_ENSEMBLE <- bitr(lipids_genes$V1, fromType = "ENTREZID", toType = "ENSEMBL", OrgDb="org.Hs.eg.db")

#2 гена не перевелись, 57 дупликатов
setdiff(lipids_genes$V1,lipid_genes_ENSEMBLE$ENTREZID)
unique(lipid_genes_ENSEMBLE$ENTREZID[duplicated(lipid_genes_ENSEMBLE$ENTREZID)])
 

goal <- intersect(lipid_genes_ENSEMBLE$ENSEMBL,winers)
intersect(x$ENTREZID, lipids_genes$V1) #- 15 генов для трешхолда 0.65


#сделаем табличку для разных трешхолдов
uniqe_genes <- c()
unique_lipids <- c()
correlation_num <- c()
lipid_genes_th <- c()
winners_list <- list()
for (threshold in seq(0.5, 0.8, 0.05)){
  x <- gl_scor_m2_data[abs(gl_scor_m2_data$correlation) > threshold,]
  correlation_num <-append(correlation_num, dim(gl_scor_m2_data[abs(gl_scor_m1_data$correlation) > threshold,])[1])
  uniqe_genes <-append(uniqe_genes,length(unique(gl_scor_m2_data$gene[abs(gl_scor_m2_data$correlation) > threshold])))
  unique_lipids <-c(unique_lipids,length(unique(gl_scor_m2_data$lipid[abs(gl_scor_m2_data$correlation) > threshold])))
  winners_list[[as.character(threshold)]] <- x
  lipid_genes_th <- c(lipid_genes_th,length(intersect(lipid_genes_ENSEMBLE$ENSEMBL, 
                                               unique(gl_scor_m2_data$gene[abs(gl_scor_m2_data$correlation) > threshold]))))
}

winners_df <- data.frame("Threshold" = seq(0.5, 0.8, 0.05),"Number_of_correlation" = correlation_num, 
                        "Number_of_genes_assosiated_with_lipids" = lipid_genes_th,
                        "Unique_genes" = uniqe_genes, "Unique_lipids" = unique_lipids)

#как меняется распределение по классам с устрожением критерия
colors_for_lipid_class <- c("SulfoHexCer"= "#00BFC4", "HexCer"= "#F8766D", "Cer"= "#00A9FF" , "PE" = "#7CAE00", "PC" = "#C77CFF",
                            "MG"= "#FF61CC", "TG" = "#00BE67", "PG"  ="#00C198", "DG"= "#F564E3", "SM" = "#619CFF", "LPE" = "#DE8C00")

for (i in seq(0.5, 0.8, 0.05)){
  assign (paste("g", as.character(i), sep = "_"),
          ggplot(winners_list[[as.character(i)]], aes(lipid_class, fill = lipid_class)) + geom_bar() +
            theme(legend.position = "none") +  scale_fill_manual(values=colors_for_lipid_class) + 
            labs(x= "",title = paste("Number of correlation coefficient >", as.character(i)),
                 y = "Count"))
}

ggarrange(g_0.5, g_0.55, g_0.6, g_0.65, g_0.7, g_0.75, g_0.8)

#Сколько из победителей пересекаются с более строгими наборами HS - 
#насколько устойчивым будет результат если взять другую метрику
intersect(unique(winners_list["0.65"]$`0.65`$lipid), HS_lipids_strict_1) #15  из 15 (всего 26)
intersect(unique(winners_list["0.65"]$`0.65`$lipid), HS_lipids_strict_2) #17 из 19
intersect(unique(winners_list["0.65"]$`0.65`$lipid), HS_lipids_strict_3) #17 из 23
intersect(unique(winners_list["0.65"]$`0.65`$lipid), HS_lipids_medium_1) #24 из 26

intersect(unique(winners_list["0.75"]$`0.75`$lipid), HS_lipids_strict_1) #7  из 15  (всего 13)
intersect(unique(winners_list["0.75"]$`0.75`$lipid), HS_lipids_strict_2) #7 из 19
intersect(unique(winners_list["0.75"]$`0.75`$lipid), HS_lipids_strict_3) #7 из 23
intersect(unique(winners_list["0.75"]$`0.75`$lipid), HS_lipids_medium_1) #13 из 26

intersect(unique(winners_list["0.8"]$`0.8`$lipid), HS_lipids_strict_1) #3  из 15  (всего 6)
intersect(unique(winners_list["0.8"]$`0.8`$lipid), HS_lipids_strict_2) #3 из 19
intersect(unique(winners_list["0.8"]$`0.8`$lipid), HS_lipids_strict_3) #3 из 23
intersect(unique(winners_list["0.8"]$`0.8`$lipid), HS_lipids_medium_1) #6 из 26

#нинада
#ggplot(thresh_65[thresh_65$p.value< 0.05,], aes(x = lipid_class, y = correlation, fill = lipid_class)) + geom_violin() +
#  labs(title = "Correlation coefficients > 0.65")


ggplot(thresh_65[thresh_65$p.value < 0.05,], aes(lipid_class, fill = lipid_class)) + geom_bar() +
  labs(title = "Number of correlation coefficient  > 0.65
with p-value < 0.05 colored by lipid class",  x = "", y = "Number of correlation coefficient") + 
  theme(axis.text.x=element_text(angle=30, hjust=1), legend.position = "none") 

ggplot(thresh_65[thresh_65$p.value < 0.05,], aes(lipid, fill = lipid_class)) + geom_bar() +
  labs(title = "Number of correlation coefficient > 0.65 with p-value < 0.05 
colored by lipid class", x = "", y = "Number of correlation coefficient") + 
  theme(axis.text.x=element_text(angle=60, hjust=1), legend.position = "none") 

#посчитаем среднюю корреляцию для липида в классе (то же, но с нормализацией, то бишь среднее)
normalazid_freq_0.65 <- as.data.frame(table(thresh_65$lipid_class)/table(unlist(strsplit(HS_lipids_medium_2,
                                                                                         "(", fixed = T))[seq(1,64,2)])[1:7])
ggplot(normalazid_freq_0.65, aes(Var1, Freq, fill = Var1)) + geom_col() +
  labs(title = "Mean correlation coefficient > 0.65 with p-value < 0.05 
in lipid class", x = "", y = "Number of correlation coefficient") + 
  theme(axis.text.x=element_text(angle=30, hjust=1), legend.position = "none") 

#помотрим на самые частовстречающиеся гены 
l <- as.data.frame(table(thresh_65$gene))
l <- l[order(l$Freq, decreasing = T),]
dim(l[l$Var1 %in% goal,])
m <-  l[l$Var1 %in% goal,]

l[7,]
ggplot(l[l$Freq > 0,], aes(Freq)) + geom_bar(fill = "skyblue") +
  labs(x = "Gene frequency", y = "Number of genes", 
       title = "Distribution of gene frequency among genes\nwith correlation with HS lipid > 0.65")

#посмотрим на гены с самыми высокими корреляциями 
thresh_65 <- thresh_65[order(abs(thresh_65$correlation), decreasing = T),]

#посмотрим на то, с какими липидами коррелируют липидные гены у которых корреляция была высокой (для 15) 
counter = 1
for (i in goal){
  assign(paste("ggene", counter, sep ="_"),
         ggplot(thresh_65[thresh_65$gene == i,], aes(lipid_class, fill = lipid_class)) + geom_bar() + 
           scale_fill_manual(values=colors_for_lipid_class) + ylim(0,11) +
           labs(x = '', y = "Count", title =  i) + theme(legend.position = "none", 
                                                         axis.text.x=element_text(angle=30, hjust=1)))
  counter = counter + 1
}


figure <- ggarrange(ggene_3, ggene_4, ggene_5, ggene_6, ggene_9, ggene_10, 
                    ggene_12, ggene_13, ggene_15, ggene_14, ggene_1, ggene_7, ggene_8, ggene_2, ggene_11)

annotate_figure(figure, top = text_grob("Number of correlation coefficint > 0.65 for lipid genes", size = 18))

#изготовление радужного хитмэпа корреляций 
heatmap_matrix <- as.matrix(gl_scor_m2[,-1680])
rownames(heatmap_matrix) <-gl_scor_m2$lipid

heatmap_matrix <- abs(heatmap_matrix)

library(ComplexHeatmap)

Heatmap(heatmap_matrix, name = "Correlation\ncoefficient", col = rev(rainbow(10)),
        show_column_names = FALSE, show_row_names = T, clustering_method_rows = "single", 
        clustering_method_columns = "ward.D2", 
        column_title = "Correlation between HS genes and lipids (Medium 2)", 
        row_title =  "HS lipids")



BiocManager::install("clusterProfiler")

corr_genes <- as.character(unique(gl_scor_m1_data$gene[abs(gl_scor_m1_data$correlation) > threshold]))
#library(enrichR)
#dbs <- listEnrichrDbs()
#dbs <- c("GO_Molecular_Function_2017", "GO_Cellular_Component_2017", "GO_Biological_Process_2017","KEGG_2019_Human")
#enriched <- enrichr(corr_genes, dbs)
#printEnrich(enriched, "output.txt" , sep = "\t", columns = c(1:9))
#bp <- enriched[["GO_Biological_Process_2017"]]

genes_lipids_simplecor_m1_vs_all <- genes_lipids_simplecor_m1

p_gl_scor_m1_vs_all <- genes_lipids_simplecor_m1_vs_all[95:188,]
gl_scor_m1_vs_all <- genes_lipids_simplecor_m1_vs_all[1:94,]
gl_scor_m1_vs_all$lipid <- colnames(lipids_for_simplecor)
p_gl_scor_m1_vs_all$lipid <- colnames(lipids_for_simplecor)

gl_scor_m1_vs_all_data <- melt(gl_scor_m1_vs_all, id = "lipid")
p_gl_scor_m1_vs_all_data <- melt(p_gl_scor_m1_vs_all, id = "lipid")
colnames(gl_scor_m1_vs_all_data) <- c("lipid", "gene", "correlation")
gl_scor_m1_vs_all_data$p.value <- p_gl_scor_m1_vs_all_data$value
gl_scor_m1_vs_all_data$p.value <-p.adjust(gl_scor_m1_vs_all_data$p.value, "BH")

length(unique(gl_scor_m1_vs_all_data$gene[gl_scor_m1_vs_all_data$p.value < 0.05]))
length(gl_scor_m1_vs_all_data$gene[gl_scor_m1_vs_all_data$p.value < 0.05])
length(unique(gl_scor_m1_vs_all_data$lipid[gl_scor_m1_vs_all_data$p.value < 0.05]))
length(gl_scor_m1_vs_all_data$gene)

mean(abs(gl_scor_m1_vs_all_data$correlation))
max(abs(gl_scor_m1_vs_all_data$correlation))
min(abs(gl_scor_m1_vs_all_data$correlation))

ggplot(gl_scor_m1_vs_all_data, aes(correlation)) + geom_histogram(fill = "lightblue", color = "black", binwidth = 0.05) + 
  labs(x = "Pearson correlation coefficient",y = "Count", 
       title = "Correlation coefficient between HS-lipids (Relaxed 3) and all genes")

ggplot(gl_scor_m1_vs_all_data[gl_scor_m1_vs_all_data$p.value < 0.05,], aes(correlation)) + geom_histogram(fill = "pink", color = "black")+
  labs(x = "Pearson correlation coefficient", y = "Count",
       title = "Correlation coefficient between HS-lipids (Relaxed 3) and  all genes \nwith p-values < 0.05")

ggplot(gl_scor_m1_vs_all_data[gl_scor_m1_vs_all_data$p.value < 0.05,], aes(correlation, fill = lipid)) + geom_histogram() + 
  facet_wrap(~lipid)

gl_scor_m1_vs_all_data$lipid_class <- unlist(strsplit(gl_scor_m1_vs_all_data$lipid, split = "(",fixed = T))[seq(1,2101088,2)]
ggplot(gl_scor_m1_vs_all_data[gl_scor_m1_vs_all_data$p.value < 0.05,], aes(correlation, fill = lipid_class)) + geom_histogram(binwidth = 0.05) + facet_wrap(~lipid_class) 

ggplot(gl_scor_m1_vs_all_data[gl_scor_m1_vs_all_data$p.value < 0.05,], aes(lipid_class, fill = lipid_class)) + geom_bar() +
  labs(title = "Number of correlation coefficient with p-value < 0.05")

ggplot(gl_scor_m1_vs_all_data[gl_scor_m1_vs_all_data$p.value < 0.05,], aes(lipid, fill = lipid_class)) + geom_bar() +
  labs(title = "Number of correlation coefficient with p-value < 0.05")

thresh_65_vs_all <- gl_scor_m1_vs_all_data[abs(gl_scor_m1_vs_all_data$correlation) > threshold,]
dim(gl_scor_m1_vs_all_data[abs(gl_scor_m1_vs_all_data$correlation) > threshold,])
length(unique(gl_scor_m1_vs_all_data$gene[abs(gl_scor_m1_vs_all_data$correlation) > threshold]))
length(unique(gl_scor_m1_vs_all_data$lipid[abs(gl_scor_m1_vs_all_data$correlation) > threshold]))

ggplot(thresh_65_vs_all[thresh_65_vs_all$p.value< 0.05,], aes(x = lipid_class, y = correlation, fill = lipid_class)) + geom_violin() +
  labs(title = "Correlation coefficients > 0.65")

ggplot(gl_scor_m1_vs_all_data, aes(x = lipid_class, y = correlation, fill = lipid_class)) + geom_violin() 

ggplot(thresh_65_vs_all[thresh_65_vs_all$p.value< 0.05,], aes(lipid_class, fill = lipid_class)) + geom_bar() +
  labs(title = "Number of correlation coefficient >0.65 (<-0.65) with p-value < 0.05")

ggplot(thresh_65_vs_all[thresh_65_vs_all$p.value< 0.05,], aes(lipid, fill = lipid_class)) + geom_bar() +
  labs(title = "Number of correlation coefficient  >0.65 (<-0.65) with p-value < 0.05")

ggplot(thresh_65_vs_all[thresh_65_vs_all$p.value< 0.05,], aes(lipid, fill = lipid_class)) + geom_bar() +
  labs(title = "Number of correlation coefficient with p-value < 0.05")








