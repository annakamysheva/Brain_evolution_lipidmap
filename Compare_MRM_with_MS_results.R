tr <- function(df){
  #' Transponating function
  #' 
  #' Returns transponated dataframe
  return(as.data.frame(t(df)))
}

rename_MRM <- function(MRM_df){
  #' Function for renaming MRM dataframe
  #' @description "SHexCer_40:2;2_d18:1_result_area" -> "SHexCer d40:2"
  #' @details based on MRM_targetlist table from "Desktop/Khrameeva_2019/Coding/DATA_actual/MRM/MRM_raw_tables/MRM_human_TARGETLIST.xlsx"
  #' can be rid with read_excel from library(readxl)
  return(sapply(colnames(MRM_df), function(x) MRM_targetlist$`Lipid species ms`[MRM_targetlist$peak == x]))
}

rename_MS <- function(MS_df){
  #' Function for renaming MS dataframe
  #' @description "negFT21240" -> "SHexCer d40:2"
  #' @details based on lipids_anno table from "Desktop/Khrameeva_2019/Coding/DATA_actual/MS_anno_new/MS_new_anno_raw/MS_new_anno_lipids_annotation.csv"
  return(sapply(colnames(MS_df), function(x) lipids_anno$Bulk.structure[lipids_anno$index_newtable_4sp == x]))
}

get_lipid_class <- function(name_vec){
  #' Function for getting lipid class
  #' @description "SHexCer d40:2" -> "SHexCer"
  #' @details splitting by space
  return(unlist(strsplit(name_vec, " ", fixed = TRUE))[seq(1, length(name_vec)*2, by = 2)])
}

hm_func <- function(matrix, col = "RdBu", title = ""){
  #' Lipid datasets correlation heatmap
  #'@description big heatmap magic exectly for correlations between MS and MRM
  hp <- ComplexHeatmap::Heatmap(matrix, 
                row_order = order(as.numeric(gsub("row", "", rownames(matrix)))),
                column_order = order(as.numeric(gsub("column", "", colnames(matrix)))), 
                row_split = get_lipid_class(rownames(matrix)), 
                column_split = get_lipid_class(colnames(matrix)), 
                border = TRUE, 
                col = rev(
                  RColorBrewer::brewer.pal(8, col)),  
                column_title_gp = gpar(fill = MRM_colors, font = 1:3),  
                row_title_gp = gpar(fill = MRM_colors, font = 1:3), 
                row_gap = unit(0, "mm"), column_gap= unit(0, "mm"), #setting size of gaps between groups
                show_row_names = FALSE, show_column_names = F, 
                cell_fun = function(j, i, x, y, width, height, fill) {
                  if(i == j)
                    grid.points(x, y, pch = 1, size = unit(1, "mm"))
                }
                
  )
  return(hp)
}

box_plot <- function(mtx, thresh1 = 0.7, thresh2 = 0.5, title = ""){
  #' Plot Boxplot of correlation coefficient by lipid
  #' @description Colored by lipid class 
  #' Red dots display correlation of lipid with itself 
  #' Line display threshhold
  df <-reshape2::melt(mtx)
  colnames(df) <- c("Lipid1", "Lipid2", "Correlation_coefficient")
  df$Lipid1 <- as.character(df$Lipid1)
  df$Lipid_class <- get_lipid_class(df$Lipid1)
  cors <- diag(mtx)
  cors_NA <- c()
  for (i in cors){
    cors_NA <- c(cors_NA, i, rep(NA, dim(mtx)[1] - 1)) 
  }
  df$cors_NA <- cors_NA
  ggplot(df, aes(x = Lipid1, y = Correlation_coefficient, fill = Lipid_class, color = Lipid_class)) + geom_boxplot() +
      geom_point(aes(x = Lipid2, y = cors_NA), color = "red", shape = 19) + 
      geom_hline(yintercept= thresh1, linetype="dashed")+ 
      geom_hline(yintercept= thresh2, linetype="dashed", color = "red")+
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
      labs(title = paste("Correlation coefficient distribution between MS new anno and MRM in", title), 
          x = "Lipid species", y = "Pearson Correlation Coefficient")
}

box_plot_by_class <- function(mtx, title =""){
  df <- data.frame("Lipid_class" = get_lipid_class(colnames(mtx)), 
                   "Correlation_coefficient" = diag(mtx))
  print(df)
  ggplot(df, aes(x = Lipid_class, y = Correlation_coefficient, fill = Lipid_class)) +
    geom_violin() + 
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
    labs(title = title, 
         x = "Lipid class", y = "Pearson Correlation Coefficients") 
}

#################################################
#Data loading
#################################################
#Full enviromant 
load("/Users/ann/Desktop/Khrameeva_2019/Coding/MS_new_anno/VS_MRM/MS_new_anno_VS_MRM.RData")

#by step
load("/Users/ann/Desktop/Khrameeva_2019/Coding/DATA_actual/MRM/MRM_processed/MRM_preprocessed.RData")
load("/Users/ann/Desktop/Khrameeva_2019/Coding/DATA_actual/MS_anno_new/MS_new_anno_preprocessed/MS_new_anno_preprocessed.RData")
c_MRM_33 <- c_MRM_33_no_CHA
rm(b_MRM_33_no_BB, c_MRM_33_no_CHA)
lipids_anno <-read.csv("Desktop/Khrameeva_2019/Coding/DATA_actual/MS_anno_new/MS_new_anno_raw/MS_new_anno_lipids_annotation.csv")
MRM_targetlist <- read_excel("Desktop/Khrameeva_2019/Coding/DATA_actual/MRM/MRM_raw_tables/MRM_human_TARGETLIST.xlsx")

library(docstring)
library(readxl)
library(dplyr) #select
library(ggplot2)
library(ggpubr)
library(ComplexHeatmap)
library(RColorBrewer)
library(reshape2)

setwd("Desktop/Khrameeva_2019/Coding/MS_new_anno/VS_MRM/")

MRM_colors <- c("Cer" = "#F8766D", "Cholesterol" = "#E18A00", "DG" = "#BE9C00", "HexCer" = "#8CAB00", 
                "LPC" = "#24B700", "LPE" = "#00BE70", "MG" = "#00C1AB", "PC" = "#00BBDA", "PE" = "#00ACFC", "PG" = "#8B93FF", 
                "SM" = "#616bff", "SHexCer" = "#F962DD","SulfoHexCer" = "#F962DD", "TG" = "#FF65AC", 
                "PS" = "#D575FE", "CE" = "#7a83ff", "PI" = "#f64d41")
#################################################
#preprocessing
#################################################
all_datasets <-  list(h_MS, b_MS, c_MS, m_MS, h_MRM_33, b_MRM_33, c_MRM_33, m_MRM_33)
all_datasets_names <- c("h_MS", "b_MS", "c_MS", "m_MS", "h_MRM_33", "b_MRM_33", "c_MRM_33", "m_MRM_33")
names(all_datasets) <-  all_datasets_names

#sorting columsn in MS (they are in different order in comparison with with MRM)
#lets make an alphabetic order 
for (i in all_datasets_names[1:4]){
  . <- all_datasets[[i]]
  . <- .[,order(colnames(.))]
  assign(i, .)
}

#transponating 
all_datasets <-  list(h_MS, b_MS, c_MS, m_MS, h_MRM_33, b_MRM_33, c_MRM_33, m_MRM_33)
names(all_datasets) <-  all_datasets_names
for (i in all_datasets_names){
  assign(i, tr(all_datasets[[i]]))
}

#renaming 
#MRM 
all_datasets <-  list(h_MS, b_MS, c_MS, m_MS, h_MRM_33, b_MRM_33, c_MRM_33, m_MRM_33)
names(all_datasets) <-  all_datasets_names
for (i in all_datasets_names[5:8]){
  colnames(all_datasets[[i]]) <- rename_MRM(all_datasets[[i]])
  assign(i, all_datasets[[i]])
}

#MS 
for (i in all_datasets_names[1:4]){
  colnames(all_datasets[[i]]) <- rename_MS(all_datasets[[i]])
  assign(i, all_datasets[[i]])
}

MRM_4sp = rbind(h_MRM_33, b_MRM_33, c_MRM_33, m_MRM_33)
MS_4sp = rbind(h_MS, b_MS, c_MS, m_MS)

##############################################################
#Correlation between MS and MRM 
###############################################################

#1) Common lipids findig
MS_MRM_common <- intersect(colnames(h_MS), colnames(h_MRM_33)) #151 common all unique

#2) Make dfs with only common lipids 
#add MRM_4sp, MS_4sp
all_datasets <-  list(h_MS, b_MS, c_MS, m_MS, h_MRM_33, b_MRM_33, c_MRM_33, m_MRM_33, MS_4sp, MRM_4sp)
all_datasets_names <- c("h_MS", "b_MS", "c_MS", "m_MS", "h_MRM_33", "b_MRM_33", "c_MRM_33", "m_MRM_33", "MS_4sp", "MRM_4sp")
names(all_datasets) <-  all_datasets_names

for (name in all_datasets_names){
  assign(name, select(all_datasets[[name]], MS_MRM_common))
}

#3)Correlation counting
cor_list <- list()
all_datasets <-  list(h_MS, b_MS, c_MS, m_MS,MS_4sp, h_MRM_33, b_MRM_33, c_MRM_33, m_MRM_33,MRM_4sp)

for (i in 1:5){
  cor_list[[i]] <- cor(all_datasets[[i]], all_datasets[[i+5]])
}

#diagonal_elements 
diag_elems <- list()
for (i in 1:5) diag_elems[[i]] <- diag(cor_list[[i]])

#4)Легкая статистика 
for (i in 1:5){
  print(paste("****************", graph_names[i], "*****************", sep = ""))
  print(paste("Range:", range(diag_elems[[i]])[1], range(diag_elems[[i]])[2]))
  print(paste("Mean:",mean(diag_elems[[i]])))
  print(paste("Median:",median(diag_elems[[i]])))
}
##################################################################
#Correlation plotting 
##################################################################
#1) Correlation Hists 
pdf(file = "MS_new_anno_VS_MRM_corr_hists.pdf",   # The directory you want to save the file in
    width = 10.5, # The width of the plot in inches
    height = 6)

graph_names <- c("HUMAN", 
                 "BONOBO", 
                 "CHIMPANZEE",
                 "MACAQUE", 
                 "All species together")

names(graph_names) <-  1:5

colors <- c("#24B700","#F8766D", "#FF65AC", "#F962DD", "#00BBDA")
names(colors) <-  1:5

plotlist <-  list()
for(i in 1:5){ 
  . <- as.data.frame(diag_elems[[i]])
  colnames(.) <- "x"
  plotlist[[i]] <- ggplot(., aes(x=x)) + 
    geom_histogram(bins = 30, fill = colors[i]) +
    labs(title = graph_names[i], x = "Pearson correlation coefficient", y = "Frequency") +
    theme(legend.position = "none")
}

cor_hists <- ggarrange(plotlist[[5]], ggarrange(plotlist = plotlist[c(1,2,3,4)]))

annotate_figure(cor_hists,
                top = text_grob("Correlations between same lipids in MS and MRM datasets", size = 16, face = "bold"))

dev.off()

#2) Heatmaps

for (i in 1:5){
  pdf(file = paste("MS_new_anno_VS_MRM_corr_heatmap_", graph_names[i] , ".pdf", sep = ""),   
      width = 15, # The width of the plot in inches
      height = 15)
  
  . <- hm_func(cor_list[[i]])
  draw(.,column_title = paste("Correlations between MS new anno and MRM common lipids in",  graph_names[i]), column_title_gp = gpar(fontsize = 16))
  
  dev.off()
}

#3) Boxplots 
box_plot(cor_list[[i]])

for (i in 1:5){
  pdf(file = paste("MS_new_anno_VS_MRM_corr_boxplot_", graph_names[i] , ".pdf", sep = ""),   
      width = 15, # The width of the plot in inches
      height = 8)
  
  plot(box_plot(cor_list[[i]], title = graph_names[i]))
  dev.off()
}


#by lipid class 
for (i in 1:5){
  pdf(file = paste("MS_new_anno_VS_MRM_corr_boxplot_byclass", graph_names[i], ".pdf", sep = ""),   
      width = 8, # The width of the plot in inches
      height = 6)
  
  plot(box_plot_by_class(cor_list[[i]], title = graph_names[i]))
  dev.off()
}

#all in one picture
pdf(file = "MS_new_anno_VS_MRM_corr_boxplot_byclass_common_pic2.pdf",   
    width = 20,
    height = 10)

plotlist <- list()
for (i in 1:5){
  plotlist[[i]] <- box_plot_by_class(cor_list[[i]], title = graph_names[i])
}

cor_hists <- ggarrange(plotlist[[5]], ggarrange(plotlist = plotlist[c(1,2,3,4)], legend = "none"), common.legend = T, legend = "right")

annotate_figure(cor_hists,
                top = text_grob("Correlation coefficients distribution between same lipids\nin MS new anno and MRM datasets by lipid class", size = 16, face = "bold"))

dev.off()



#######################################################################
#ANOVA Human-specificity
#######################################################################



