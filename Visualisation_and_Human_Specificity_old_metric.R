############################
#new data analysis and full pipeline 


#######################################################################
#Functions
#######################################################################

name_convertator <- function(MS_name_vec, brackets = T){
  #converts names like "posFT25701" to "DG 36:3" if brackets = False, or "DG(36:3)"
  names_vec <- unlist(lapply(MS_name_vec, function(x) lipids_anno$Bulk.structure[lipids_anno$index_newtable_4sp == x]))
  if (brackets){
    names_vec <- sub(" ", "\\(", names_vec)
    names_vec <- unlist(lapply(names_vec, function(x) paste(x, ")", sep = "")))
  }
  return(names_vec)
}

lipid_class <- function(names_vec, class_names = T, brackets = F){
  #class_names = T if name like "DG 36:3" or "DG(36:3)" and = F if name like "posFT25701"
  #brackets = T if names contains brackets like "DG(36:3)"
  if (class_names & !brackets){
    return(unlist(sapply(names_vec, function(x) lipids_anno$Lipid.class[lipids_anno$Bulk.structure ==  x])))
  } else if  (class_names & brackets){
    return(sapply(strsplit(names_vec, split = "(", fixed = T), "[[", 1))
  } else{
    return(sapply(names_vec, function(x) lipids_anno$Lipid.class[lipids_anno$index_newtable_4sp ==  x]))
  }
}

species_getter <- function(weird_names){
  #work for things like "MS146_MA" - gives sp letter
  result <- weird_names %>% strsplit("_") %>% sapply("[[", 2) %>% substr(1,1)
  return(result)
}

brain_getter <- function(weird_names){
  #work for things like "MS146_MA" - gives sample "HA" and others
  result <- weird_names %>% strsplit("_") %>% sapply("[[", 2)
  return(result)
}

region_getter <- function(weird_names){
  #work for things like "MS146_MA" - gives region with 1 to 75 number
  weird_names <- weird_names %>% strsplit("_") %>% sapply("[[", 1)
  return(sapply(weird_names, function(x) sample_anno$Region[sample_anno$Weird_flex2 == x]))
}

dif_perplex_count <- function(df, 
                              perplex_vec = c(2, 5, 10, 20, 30, 40, 50, 70, 100), 
                              perplexity_col = T){
  
  #count Tsne coordinate df with  different perplexity 
  #perpexity col - add column with vector indicating perplexity meanings
  
  dif_perplex = data.frame()
  
  for (i in perplex_vec){
    tsne <- Rtsne(df, dims = 2, perplexity=i, verbose=TRUE, max_iter = 500)
    tsne <- as.data.frame(tsne$Y)
    dif_perplex <- rbind(dif_perplex, tsne)
  }
  
  colnames(dif_perplex) <-  c("X", "Y")
  
  if (perplexity_col == T){
    dif_perplex$perplexity <- rep(paste("perplexity: ", as.character(perplex_vec)), 
                                  each = dim(df)[1])
    dif_perplex$perplexity <- factor(dif_perplex$perplexity, levels = paste("perplexity: ", as.character(perplex_vec)))
  }
  return(dif_perplex)
}


tSNE_plot <- function(tSNE_df, colorby, title = "", legend_title = "", facet_w = T, point_size = 1.5){
  #plot tSNE with different perpexities (facet_w = T)
  #or just single tSNE 
  
  if (facet_w == T){
    ggplot(tSNE_df, aes(x = X, y = Y, color = colorby)) + 
      geom_point(size = 1.5) + facet_wrap(~perplexity, scales = "free") +
      labs(x = "Coordinate 1", y = "Coordinat 2", title = title, colour = legend_title)
  } else{
    ggplot(tSNE_df, aes(x = X, y = Y, color = colorby)) + 
      geom_point(size = point_size) +
      labs(x = "Coordinate 1", y = "Coordinat 2", title = title, colour = legend_title)
  }
}

one_tSNE_plot <- function(df, colorby, title = "", legend_title = "", point_size = 1.5){
  #for single tSNE 
  tsne <- Rtsne(df, dims = 2, perplexity=30, verbose=TRUE, max_iter = 500)
  tsne <- as.data.frame(tsne$Y)
  colnames(tsne) <-  c("X", "Y")
  ggplot(tsne, aes(x = X, y = Y, color = colorby)) + 
    geom_point(size = point_size) +
    labs(x = "Coordinate 1", y = "Coordinat 2", title = title, colour = legend_title)
}


PCA_plot <- function(df, colorby, title = "", legend_title = "", point_size = 1.5){
  pca <- prcomp(df, center = TRUE,scale. = TRUE)
  pca <- pca$x 
  pca <- as.data.frame(pca[1:nrow(pca),1:2]) 
  ggplot(pca, aes(x = PC1, y = PC2, color = colorby)) + geom_point(size = point_size) +
    labs(title = title, colour = legend_title)
}

mean_region_counter <-  function(sp_letter, 
                                 all_reg_names = unique(region_getter(colnames(MS_a2_norm))), 
                                 lipids_df = MS_a2_norm, 
                                 anno_df = sample_anno,
                                 anno_name_col = "actual_name", 
                                 anno_sp_col = "species", 
                                 anno_reg_col = "Region"){
  
  #sp_letter - which species is need to be average "C" or "H" or "B" or "M"
  #all_reg_names - 33 reg names 
  #lipids_df - df with lipids contentrations
  #anno_df - df with sample annotations 
  #anno_name_col - name of column from anno_df whicth contains same as colnames(lipids_df)
  #anno_sp_col - name of column from anno_df whicth contains species in form of one letter - "C" or "H" or "B" or "M"
  #anno_reg_col - name of column from anno_df whicth contains region names 
  
  df <- data.frame(matrix(NA, nrow = dim(lipids_df)[1], ncol =0))
  
  for (i in all_reg_names){
    samples_names <- anno_df[anno_name_col][anno_df[anno_sp_col] == sp_letter & anno_df[anno_reg_col] == i]
    samples_names <- intersect(samples_names, colnames(lipids_df))
    df <- cbind(df, apply(select(lipids_df, samples_names), MARGIN = 1, mean))
  }
  colnames(df) <- all_reg_names
  return(df)
}

Multiple_Comparisons <- function(x, method = "BH"){
  return(p.adjust(x, method, length(x)))
}


Human_specific_features <- function(df){
  #df = H(regions_num, lipids_num) + B(regions_num, lipids_num) + C(regions_num, lipids_num) + M(regions_num, lipids_num)
  out_table <- data.frame("Затравка" = c(1:9))
  len <-  length(df)
  name <- colnames(df)[1:(len/4)]
  for (i in seq(1,(len/4))){
    for_species <- df[,i]
    for (j in seq(i+(len/4), len, (len/4))){
      for_species <- cbind(for_species, df[,j])
    }
    colnames(for_species) <- seq(1:4)
    rownames(for_species) <- seq(1:33)
    corr <- corr.test(for_species)
    corr_dist <- 1/abs(c(corr$r[1,2], corr$r[1,3], corr$r[1,4], corr$r[2,3], corr$r[2,4], corr$r[3,4]))
    ans1 <- ifelse(min(corr_dist[1],corr_dist[2], corr_dist[3]) > 
                     max(corr_dist[4], corr_dist[5], corr_dist[6]), 1, 0)
    ans2 <- ifelse(sort(c(corr_dist[1],corr_dist[2], corr_dist[3]))[2] > 
                     max(corr_dist[4], corr_dist[5], corr_dist[6]), 1, 0)
    ans3 <- ifelse(max(corr_dist[1],corr_dist[2], corr_dist[3]) > 
                     max(corr_dist[4], corr_dist[5], corr_dist[6]), 1, 0)
    #without multiple comparison correction
    pvalue <- c(corr$p[2,1], corr$p[3,1], corr$p[4,1], corr$p[3,2], corr$p[4,2], corr$p[4,3])
    #with multiple comparison correction
    #pvalue <- c(corr$p[1,2], corr$p[1,3], corr$p[1,4], corr$p[2,3], corr$p[2,4], corr$p[3,4])
    res <- c(ans1, ans2, ans3, pvalue)
    out_table[name[i]] <- res
  }
  rownames(out_table) <- c("Strict_HS", "Medium_HS", "Relaxed_HS",
                           "HB_p.value", "HC_p.value", "HM_p.value", 
                           "BC_p.value", "BM_p.value", "CM_p.value")
  return(out_table[,2:ncol(out_table)])
}


HS_wich_pass_criteria <- function(lipids_names, HS_list){
  #func to help HS_lipids_lists with df creation 
  #lipids_names - names of all lipids which was tested to HS 
  #HS_list - list of lipids passes criteria 
  return(ifelse(lipids_names %in% HS_list, 1, 0))
}

HS_lipids_lists <- function(HS_df = HS_lipids, p = 0.05, need_df_type1 = T, need_df_type2 = F){
  
  #func count HS based on p-value 
  #HS_df - df from Human_specific_features function result 
  #p - p_vlaue 
  #need_df_type1 - will create df with lipid names as rows and 9 criteria columns with 0 - lipid doesn't pass criteria,
  #1 - lipid passes criteria
  #need_df_type1 - will create df with nine columns of criteria, containing list of lipids passed criteria and NA 
  
  HS_lipids_s1 <<- rownames(subset(as.data.frame(t(HS_df)), Strict_HS == 1 & HB_p.value < p & HC_p.value < p & HM_p.value < p & BC_p.value < p & BM_p.value < p & CM_p.value < p))
  HS_lipids_r1 <<- rownames(subset(as.data.frame(t(HS_df)), Relaxed_HS == 1 & HB_p.value < p & HC_p.value < p & HM_p.value < p & BC_p.value < p & BM_p.value < p & CM_p.value < p))
  HS_lipids_m1 <<- rownames(subset(as.data.frame(t(HS_df)), Medium_HS == 1 & HB_p.value < p & HC_p.value < p & HM_p.value < p & BC_p.value < p & BM_p.value < p & CM_p.value < p))
  
  HS_lipids_s2 <<- rownames(subset(as.data.frame(t(HS_df)), Strict_HS == 1 &  (HB_p.value + HC_p.value + HM_p.value + BC_p.value + BM_p.value + CM_p.value)/6 < p))
  HS_lipids_r2 <<- rownames(subset(as.data.frame(t(HS_df)), Relaxed_HS == 1 &  (HB_p.value + HC_p.value + HM_p.value + BC_p.value + BM_p.value + CM_p.value)/6 < p))
  HS_lipids_m2 <<- rownames(subset(as.data.frame(t(HS_df)), Medium_HS == 1 &  (HB_p.value + HC_p.value + HM_p.value + BC_p.value + BM_p.value + CM_p.value)/6 < p))
  
  HS_lipids_s3 <<- rownames(subset(as.data.frame(t(HS_df)), Strict_HS == 1 &  (HB_p.value * HC_p.value * HM_p.value * BC_p.value * BM_p.value * CM_p.value) < p**6))
  HS_lipids_r3 <<- rownames(subset(as.data.frame(t(HS_df)), Relaxed_HS == 1 &  (HB_p.value * HC_p.value * HM_p.value * BC_p.value * BM_p.value * CM_p.value) < p**6))
  HS_lipids_m3 <<- rownames(subset(as.data.frame(t(HS_df)), Medium_HS == 1 &  (HB_p.value * HC_p.value * HM_p.value * BC_p.value * BM_p.value * CM_p.value) < p**6))
  
  if (need_df_type1){
    HS_lipids1 <<- data.frame(HS_wich_pass_criteria(colnames(HS_lipids), HS_lipids_s1), 
                                HS_wich_pass_criteria(colnames(HS_lipids), HS_lipids_s2), 
                                HS_wich_pass_criteria(colnames(HS_lipids), HS_lipids_s3), 
                                HS_wich_pass_criteria(colnames(HS_lipids), HS_lipids_m1), 
                                HS_wich_pass_criteria(colnames(HS_lipids), HS_lipids_m2), 
                                HS_wich_pass_criteria(colnames(HS_lipids), HS_lipids_m3), 
                                HS_wich_pass_criteria(colnames(HS_lipids), HS_lipids_r1), 
                                HS_wich_pass_criteria(colnames(HS_lipids), HS_lipids_r2), 
                                HS_wich_pass_criteria(colnames(HS_lipids), HS_lipids_r3))
    
    rownames(HS_lipids1) <<- colnames(HS_lipids)
    colnames(HS_lipids1) <<- c("Strict_1", "Strict_2", "Strict_3", "Medium_1", "Medium_2", "Medium_3",
                                 "Relaxed_1", "Relaxed_2", "Relaxed_3")
    
  }
  
  if (need_df_type2){
    size = length(HS_lipids_r3)
    HS_lipids2 <<- data.frame(c(HS_lipids_s1, rep(NA, size - length(HS_lipids_s1))),
                                c(HS_lipids_s2, rep(NA, size - length(HS_lipids_s2))),
                                c(HS_lipids_s3, rep(NA, size - length(HS_lipids_s3))),
                                c(HS_lipids_m1, rep(NA, size - length(HS_lipids_m1))),
                                c(HS_lipids_m2, rep(NA, size - length(HS_lipids_m2))),
                                c(HS_lipids_m3, rep(NA, size - length(HS_lipids_m3))),
                                c(HS_lipids_r1, rep(NA, size - length(HS_lipids_r1))),
                                c(HS_lipids_r2, rep(NA, size - length(HS_lipids_r2))),
                                c(HS_lipids_r3, rep(NA, size - length(HS_lipids_r3))))
    
    colnames(HS_lipids2) <<- c("Strict_1", "Strict_2", "Strict_3", "Medium_1", "Medium_2", "Medium_3",
                                 "Relaxed_1", "Relaxed_2", "Relaxed_3")
    
  }
}

Stacked_bar_plot_HS <- function(HS_lipids2, all_lipids_names, title = ""){
  
  melted <- HS_lipids2 %>% melt( measure.vars  = colnames(HS_lipids2)) %>% na.omit()
  colnames(melted) <- c("Criteria", "Lipid_name")
  
  
  melted <- rbind(melted, data.frame("Criteria" = rep("MS_data", length(all_lipids_names)), 
                                     "Lipid_name" = all_lipids_names))
  melted$Class <- lipid_class(melted$Lipid_name, class_names = F)
  melted$Criteria <- relevel(melted$Criteria, "MS_data")
  
  ggplot(melted, aes(x = Criteria, fill = Class)) + geom_bar(position = "fill") +
    scale_fill_manual(values = Colors) + 
    labs(x = "Criterion", 
         y = "Share of lipid class", 
         title = title) 
}

#######################################################################
#packages
#######################################################################

library(magrittr)
library(Rtsne)
library(ggplot2)
library(dplyr)
library (psych)
library(reshape2)

#######################################################################
#Design
#######################################################################

Colors_old <- c("Cer" = "#F8766D", "Cholesterol" = "#E18A00", "DG" = "#BE9C00", "HexCer" = "#8CAB00", 
                "LPC" = "#24B700", "LPE" = "#00BE70", "MG" = "#00C1AB", "PC" = "#00BBDA", "PE" = "#00ACFC", "PG" = "#8B93FF", 
                "SM" = "#D575FE", "SulfoHexCer" = "#F962DD", "TG" = "#FF65AC", 
                "PS" = "#616bff", "CE" = "#f64d41", "PI" = "#7d83ff", "CAR" = "#FF6C90", "FA" = "#00C1AB", "PC_O" = "#00BB6A",
                "PC_P" = "#00BB2A", "PE_P" = "#00AC2C", "PE_O" = "#00AC6C")

class_order <- c("CAR",  "CE", "Cer","Cholesterol",  "DG", "FA", "HexCer" , "LPC" , "LPE",
                 "PC","PC_O", "PC_P",  "PE","PE_O", "PE_P", "PG" , "PI" , "PS" , "SM",  "SulfoHexCer", "TG")

Colors <- c("CAR" = "#FF61C9",  "CE" = "#FF689E", "Cer" = "#F8766D","Cholesterol" = "#F37B59",
                       "DG" = "#ED8141", "FA" = "#E09000", "HexCer" = "#BB9D00", "LPC" = "#A3A500" , "LPE" = "#85AD00",
                       "PC" = "#5BB300","PC_O" = "#00B81F", "PC_P" = "#00BE6C",  
                       "PE" = "#00C1AA","PE_O" = "#00BBDB", "PE_P" = "#00B0F6",
                       "PG" = "#529EFF" , "PI" = "#9590FF", "PS" = "#BF80FF", 
                       "SM" = "#DC71FA",  "SulfoHexCer" = "#F763E0", "TG" = "#FF65AE")

#######################################################################
#Data reading
#######################################################################
setwd("Desktop/Khrameeva_2019/Coding/MS_new_anno/")

#new MS data with normalization by sample (bybrain) 
MS_a2_norm <-  read.csv("http://arcuda.skoltech.ru/~khrameeva/brainmap/newTL/tables/newTL.33regions_norm.txt.bybrain", sep = "\t", check.names = F)
#lipids annotation
lipids_anno <- read.csv("http://arcuda.skoltech.ru/~khrameeva/brainmap/newTL/tables/newTL_4sp.txt", sep = "\t")
#names of lipids have duplicates, therefore will use original rownames for better statistics

#brain sample annotation 
sample_anno <-  read.csv("/Users/ann/Desktop/Khrameeva_2019/Coding/DATA/region_names_4species.txt", sep = "\t", header = F, col.names = c("Num", "Full_name","Exp_type", "Weird_flex1", "Weird_flex2","Sample", "Sammple_ID",  "Region", "Region_type"))

anova_res <- read.csv("http://arcuda.skoltech.ru/~khrameeva/brainmap/newTL/tables/ANOVA_regions_and_species.txt", sep = "\t", 
                      header = F, col.names = c("Lipid_name", "Anova_p_value"))

#number of samples
colnames(MS_a2_norm) %>% strsplit("_") %>% sapply("[[", 1) %>%
  sapply(function(x) sample_anno$Sample[sample_anno$Weird_flex2 == x]) %>%table(.)
#       BA   BB   BC  CHA  CHB  CHC  CHD   HA   HB   HC   HD   MA   MB   MC   QC wash 
#   0   34   29   33    0   33   32   33   33   33   35   33   31   31   32    0    0 

#######################################################################
#Visualizing (tSNE + PCA)
#######################################################################

#MS_anno2_tSNE_sp
tSNE_plot(dif_perplex_count(t(MS_a2_norm)), colorby = rep(species_getter(colnames(MS_a2_norm)), times = 9), point_size = 1, 
          legend_title = "Species",title = "MS new annotation tSNE by species")
#MS_anno2_tSNE_reg
tSNE_plot(dif_perplex_count(t(MS_a2_norm)), colorby = rep(region_getter(colnames(MS_a2_norm)), times = 9),  point_size = 1,
          legend_title = "Region",title = "MS new annotation tSNE by region")
#MS_anno2_tSNE_brain
tSNE_plot(dif_perplex_count(t(MS_a2_norm)), colorby = rep(brain_getter(colnames(MS_a2_norm)), times = 9),  point_size = 1,
          legend_title = "Brain",title = "MS new annotation tSNE by brain sample")

#MS_anno2_PCA_sp
PCA_plot(t(MS_a2_norm), colorby = species_getter(colnames(MS_a2_norm)), 
         legend_title = "Species",title = "MS new annotation PCA by species")
#MS_anno2_PCA_reg
PCA_plot(t(MS_a2_norm), colorby = region_getter(colnames(MS_a2_norm)), 
         legend_title = "Region",title = "MS new annotation PCA by region")
#MS_anno2_PCA_brain
PCA_plot(t(MS_a2_norm), colorby = brain_getter(colnames(MS_a2_norm)), 
         legend_title = "Brain",title = "MS new annotation PCA by brain sample")


#######################################################################
#Human-specificity
#######################################################################
#first need to add column with column_names to sample_annotation to count mean 
sample_anno$actual_name <- paste(sample_anno$Weird_flex2, sample_anno$Sample, sep = "_")
sample_anno$species <- substr(sample_anno$Sample, 1,1)

#counting mean (all params set as default in mean_region_counter)
b_MS <- mean_region_counter("B")
h_MS <- mean_region_counter("H")
c_MS <- mean_region_counter("C")
m_MS <- mean_region_counter("M")
#cheking for success
sum(is.na.data.frame(c_MS))

#HS_df creating 
HS_lipids <- Human_specific_features(as.data.frame(t(rbind(h_MS, b_MS, c_MS, m_MS))))

#Multiple comparison
for (i in 4:9){
  HS_lipids[i,] <- Multiple_Comparisons(HS_lipids[i,])
}

#readable HS_df creating (HS_lipids1, HS_lipids2 and HS vectors)
HS_lipids_lists(HS_lipids, need_df_type1 = T, need_df_type2 = T)
write.csv(HS_lipids1, "MS_anno2_HS_lipids.csv")
write.csv(HS_lipids, "MS_anno2_HS_lipids_raw.csv")


#######################################################################
#Human-specificity visualisation 
#######################################################################
#unique classes in this df
sort(unique(lipid_class(colnames(HS_lipids), class_names = F)))

#MS_HS_barplot
Stacked_bar_plot_HS(HS_lipids2, rownames(MS_a2_norm),
                    title = "Classes representation among  human-specific lipids in MS anno2 dataset")


#######################################################################
#Human-specificity within ANOVA results
#######################################################################
anova_specific <- as.character(anova_res$Lipid_name[anova_res$Anova_p_value < 0.05])

HS_lipids <- as.data.frame(t(HS_lipids))
anova_HS <- data.frame("Strict" = c(intersect(rownames(HS_lipids[HS_lipids$Strict_HS == 1,]),anova_specific) , rep(NA, 102 - 23)), 
                       "Medium" = c(intersect(rownames(HS_lipids[HS_lipids$Medium_HS == 1,]),anova_specific), rep(NA, 102 - 46)), 
                       "Relaxed" = intersect(rownames(HS_lipids[HS_lipids$Relaxed_HS == 1,]),anova_specific))
#MS_HS_barplot_ANOVA
Stacked_bar_plot_HS(anova_HS, anova_specific, "Classes representation among  human-specific lipids\nin species-specific lipids (ANOVA)")

write.csv(anova_HS, "MS_anno2_HS_lipids_with_ANOVA_filter.csv")






