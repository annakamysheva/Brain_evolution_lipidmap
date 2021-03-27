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
#(|h-b| + |h-c| + |h-m|)/(|b-c| + |b-m| + |c-m|)

#INPUT: 4 df in the following order
#Human_data, Bonobo_data, Chimpanzee_data, Macaque_data

#OUTPOT: ratio df
####################################################################

Ratio_counter <- function(h_df, b_df, c_df, m_df){
  diff1 <- abs(h_df - b_df) + abs(h_df - c_df) + abs(h_df - m_df)  #не делю на три потому что потом все равно беру отношение
  diff2 <- abs(b_df - m_df) + abs(b_df - c_df) + abs(c_df - m_df) 
  return(diff1/diff2)
}

#Helping function for vector elements composition
Composition <- function(vect){
  res = 1
  for (i in vect){
    res = res*i
  }
  return(res)
}

#Function for counting distribution mode
Getmode <- function(v) {
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}


library(ggplot2)
library(reshape2)
library(ComplexHeatmap)
library(dplyr)
library(tidyr)
library(clusterProfiler)
library(ggpubr)

lipids_ratio <- Ratio_counter(human_lipids, bonobo_lipids, chimp_lipids, macaque_lipids)
genes_ratio <- Ratio_counter(h_genes, b_genes, c_genes, m_genes)

row_names <- paste(as.character(seq(1,33)),  substr(rownames(b_genes),4, nchar(rownames(b_genes))))
rownames(lipids_ratio) <- row_names
rownames(genes_ratio) <- row_names


############################################################################################
#HS regions based on ratio
############################################################################################

region_ratio_sum <- as.data.frame(sort(apply(X = lipids_ratio,FUN =  sum,MARGIN =1)))
region_ratio_comp <- as.data.frame(sort(apply(X = lipids_ratio,FUN =  Composition,MARGIN =1)))
colnames(region_ratio_comp) <- c("Composition")
colnames(region_ratio_sum) <- c("Sum")

region_ratio_sum$region = rownames(region_ratio_sum)
region_ratio_sum$region <- factor(region_ratio_sum$region, levels = rownames(region_ratio_sum))
levels(region_ratio_sum$region)

region_ratio_comp$region = rownames(region_ratio_comp)
region_ratio_comp$region <- factor(region_ratio_sum$region, levels = rownames(region_ratio_sum))


ggplot(region_ratio_sum, aes(x = region, y = Sum)) + geom_col(stat = "identity", fill = "skyblue") +
  theme(axis.text.x=element_text(angle=90, hjust=1)) + labs(x = "", y = "Ratio Sum", 
                                                            title = "Total sum of lipid ratio in brain regions")


HS_regions_rating <- data.frame("Sum rating" = region_ratio_sum$region[seq(33,1,-1)],
                                "Composition rating" = region_ratio_comp$region[seq(33,1,-1)],
                                "HS position" = seq(1,33), check.names = F)


melt(lipids_ratio) %>% 
  ggplot(aes(variable, value, fill = variable)) + geom_boxplot() + 
  theme(axis.text.x=element_text(angle=90, hjust=1), legend.position = "none") + 
  labs(x = "", y = "Ratio", title = "Distribution of ratio in brain regions")

melt(as.data.frame(t(1/lipids_ratio))) %>% 
  ggplot(aes(variable, value, fill = variable)) + geom_boxplot(notch = T) + 
  theme(axis.text.x=element_text(angle=90, hjust=1), legend.position = "none") + 
  labs(x = "", y = "Ratio", title = "Distribution of 1/ratio in brain regions")

min(melt(as.data.frame(t(lipids_ratio)))[2])
max(melt(as.data.frame(t(lipids_ratio)))[2])
median(unlist((melt(as.data.frame(t(lipids_ratio)))[2])))
mean(unlist((melt(as.data.frame(t(lipids_ratio)))[2])))


HS_regions_rating$mean = rownames(as.data.frame(sort(apply(as.data.frame(t(lipids_ratio)), MARGIN = 2, FUN = mean))))[seq(33,1,-1)]
HS_regions_rating$median = rownames(as.data.frame(sort(apply(as.data.frame(t(lipids_ratio)), MARGIN = 2, FUN = median))))[seq(33,1,-1)]
HS_regions_rating$median = rownames(as.data.frame(sort(apply(as.data.frame(t(lipids_ratio)), MARGIN = 2, FUN = Getmode))))[seq(33,1,-1)]

lipids_ratio <- as.data.frame(t(lipids_ratio))
genes_ratio <- as.data.frame(t(genes_ratio))

outlier_list <- list()
for (i in colnames(lipids_ratio)){
  x = unlist(lipids_ratio[i])
  qnt <- quantile(x, probs=c(.25, .75))
  H <- 1.5 * IQR(x)
  nam <-rownames(lipids_ratio)
  outlier_list[[i]] <- c(nam[which(x < (qnt[1] - H))],   nam[which(x > (qnt[2] + H))])
}

#для обратной ситуации
outlier_list_new <- list()
for (i in colnames(lipids_ratio)){
  x = unlist(1/lipids_ratio[i])
  qnt <- quantile(x, probs=c(.25, .75))
  H <- 1.5 * IQR(x)
  nam <-rownames(lipids_ratio)
  outlier_list_new[[i]] <- c(nam[which(x < (qnt[1] - H))],   nam[which(x > (qnt[2] + H))])
}

HS_intersect <- c()
for (i in colnames(lipids_ratio)){
  HS_intersect <- c(HS_intersect, round(100*length(intersect(outlier_list[[i]], HS_lipids_medium_2))/length(outlier_list[[i]]),2))
                    
}

data.frame("HS" = HS_intersect, "regions" = colnames(lipids_ratio)) %>%
  ggplot(aes(x = regions, y = HS, fill = regions)) + geom_col() +
  theme(axis.text.x=element_text(angle=90, hjust=1), legend.position = "none") + 
  labs(x = "", y = "Percent", title = "Percent of outliers which turned out to be HS (M2)")


############################################################################################
#lipids-genes correlations ratio-based.Visualisation.
############################################################################################
ratio_corr <- Genes_lipids_correlation_simple(lipids_ratio, genes_ratio)

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


#почему-то 5 генов не посчитались ??? 
weird_genes <- unique(gl_scor_m2_data$gene[is.na(gl_scor_m2_data$correlation)])
gl_scor_m2[weird_genes]

genes_ratio[weird_genes] #содержат бесконечности, океюшки
#рип 5 генов
gl_scor_m2_data  <- drop_na(gl_scor_m2_data )
gl_scor_m2_data$lipid_class <- unlist(strsplit(gl_scor_m2_data$lipid, split = "(",fixed = T))[seq(1, 107136,2)]
gl_scor_m2_data  <- drop_na(gl_scor_m2_data, lipid) #доудаляем оставшееся в колонке с липидами

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

gl_scor_m2_data$lipid <-  factor(gl_scor_m2_data$lipid,
                                 levels = c("HexCer(d41:0)", "HexCer(d42:0)", "HexCer(d43:2)","HexCer(t40:2)", "HexCer(t41:2)", "HexCer(t42:2)",
                                            "TG(58:8)", "PG(38:1)", "Cer(m40:1)", "PE(36:2)", "PE(P-36:5)", "PE(P-36:6)",
                                            "MG(20:3)", "MG(18:2)", "MG(16:2)", "MG(22:3)", "PC(34:3)", "PC(36:4)", 
                                            "PC(40:8)", "PC(42:10)", "SulfoHexCer(d41:1)","SulfoHexCer(d42:1)","SulfoHexCer(d43:1)",
                                            "SulfoHexCer(d43:2)","SulfoHexCer(t40:1)","SulfoHexCer(t41:1)","SulfoHexCer(t41:2)",
                                            "SulfoHexCer(t42:1)","SulfoHexCer(t43:1)","SulfoHexCer(t43:2)"))


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

  
length(unique(gl_scor_m2_data$gene[gl_scor_m2_data$p.value < 0.05]))
length(gl_scor_m2_data$gene[gl_scor_m2_data$p.value < 0.05])
length(unique(gl_scor_m2_data$lipid[gl_scor_m2_data$p.value < 0.05]))
length(gl_scor_m2_data$gene)

mean(abs(gl_scor_m2_data$correlation[gl_scor_m2_data$p.value < 0.05]))
max(abs(gl_scor_m2_data$correlation[gl_scor_m2_data$p.value < 0.05]))
min(abs(gl_scor_m2_data$correlation[gl_scor_m2_data$p.value < 0.05]))
median(abs(gl_scor_m2_data$correlation[gl_scor_m2_data$p.value < 0.05]))


############################################################################################
#lipids-genes correlations ratio-based. Gene finding.
############################################################################################
lipids_genes <- read.csv("/Users/ann/Desktop/Khrameeva_2019/Coding/DATA/genes_associated_with_lipids/lipid_enzymes.txt", sep = "\t", header = F)
lipids_genes <- bitr(lipids_genes$V1, fromType = "ENTREZID", toType = "ENSEMBL", OrgDb="org.Hs.eg.db")

length(intersect(lipids_genes$ENSEMBL,unique(gl_scor_m2_data$gene[gl_scor_m2_data$p.value < 0.05])))


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

ggplot(lipids_genes_data, aes(lipid, correlation, color = lipid_class)) + geom_count(size = 7) +
  theme(axis.text.x=element_text(angle=60, hjust=1)) +
  labs(x = "", y = "Correlation coefficient",
       title = "Correlation coefficients of ratio of lipids correlating with ratio of lipid genes")

gl_scor_m2_data <- gl_scor_m2_data[order(gl_scor_m2_data$correlation, decreasing = T),]
View(gl_scor_m2_data[gl_scor_m2_data$p.value < 0.05,])
