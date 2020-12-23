OB heatmaps
================
kanazian
2020-12-23

# Goal:

Heatmaps of OR gene expression across single dimension positions
annotated with features, examination of gene covariance and calculation
of glomerulus symmetry

# Setup

Load packages and functions

``` r
#knitr::opts_chunk$set(warning=F)
#packages
library(heatmaply)
library(reshape2)
library(gtools)
library(patchwork)
library(plotly)
library(tidyverse)
library(cowplot)
library(uwot)

#functions
#x is a vector of tpm values for an OR across samples
NormalizeMaxMin <- function(x) {
  #x needs to be pseudocounted
  (x-min(x))/(max(x)-min(x))
} #end NormalizeMaxMin

GeoMean <- function(x, na.rm=TRUE) {
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
} #end GeoMean

#given raw TPM data, normalize values to between 0 and 1 (inclusive)
NormTPM <- function(df_tpm) {
  x_vals <- df_tpm %>% select(-slice) %>% t()
  x_vals <- x_vals + 0.25
  
  x_norm <- matrix(nrow = nrow(x_vals), ncol = ncol(x_vals))
  for (i in 1:nrow(x_vals)) {
    if (min(x_vals[i,]) - max(x_vals[i,]) == 0) {
      x_norm[i,] <- rep(0, length(x_vals[i,]))
    } else {
    x_norm[i,] <- NormalizeMaxMin(x_vals[i,])
    } #endif
  } #endfor
  
  colnames(x_norm) <- df_tpm$slice
  rownames(x_norm) <- rownames(x_vals)
  return(x_norm)
} # end NormTPM

#given normalized data and a list with "gene" and "sortrank" column, sort the normalized data
#list must have ORs in "gene" column and rank in "sortrank" column
SortByList <- function(df_norm, list) {
  df_tib <- as_tibble(df_norm)
  df_tib$gene <- rownames(df_norm)
  df_info <- left_join(df_tib, list, by = "gene") %>% arrange(sortrank)
  df_matrix <- t(t(df_info %>% select(-gene, -sortrank, -wavgs)))
  rownames(df_matrix) <- df_info$gene
  return(df_matrix)
} #end SortByList

#given tpm data, apply factors based on kallisto quantified OMP values
Ompnormify <- function(df_tpm) {
  soi <- df_tpm %>% select(-slice) %>% t()
  soi <- soi + 0.25
  for (i in 1:ncol(soi)) {
    soi[,i] = soi[,i] / ompNorms[i]
  } #endfor
  soi_norm <- matrix(nrow = nrow(soi), ncol = ncol(soi))
  for (i in 1:nrow(soi)) {
    soi_norm[i,] <- NormalizeMaxMin(soi[i,])
  } #endfor
  colnames(soi_norm) <- df_tpm$slice
  rownames(soi_norm) <- rownames(soi)
  return(soi_norm)
} #end Ompnormify

#given tpm data, apply factors based on OB shape from MRI model
Voxnormify <- function(df_tpm, voxweights) {
  soi <- df_tpm %>% select(-slice) %>% t()
  soi <- soi + 0.25
  for (i in 1:ncol(soi)) {
    soi[,i] <- soi[,i] * voxweights$weight[i]
  } #endfor
  soi_norm <- matrix(nrow = nrow(soi), ncol = ncol(soi))
  for (i in 1:nrow(soi)) {
    soi_norm[i,] <- NormalizeMaxMin(soi[i,])
  } #endfor
  colnames(soi_norm) <- df_tpm$slice
  ornames <- rownames(soi)
  rownames(soi_norm) <- ornames
  return(soi_norm)
} #end Voxnormify

normNormDF4 <- function(df1, df2, df3, df4) {
  df_merge <- matrix(nrow = nrow(df1), ncol = ncol(df1))
  for (i in 1:ncol(df1)) {
    df_merge[,i] <- df1[,i] + df2[,i] + df3[,i] + df4[,i]
  } #endfor
  
  #normalize merged normalized values
  merge_vals <- df_merge
  merge_norm <- matrix(nrow = nrow(merge_vals), ncol = ncol(merge_vals))
  for (i in 1:nrow(merge_vals)) {
    merge_norm[i,] <- NormalizeMaxMin(merge_vals[i,])
  } #endfor
  
  #name rows and columns
  colnames(merge_norm) <- colnames(df1)
  rownames(merge_norm) <- rownames(df1)
  return(merge_norm)
} #end normNormDF4

#need to learn how to do ... args
normNormDF3 <- function(df1, df2, df3) {
  df_merge <- matrix(nrow = nrow(df1), ncol = ncol(df1))
  for (i in 1:ncol(df1)) {
    df_merge[,i] <- df1[,i] + df2[,i] + df3[,i]
  } #endfor
  
  #normalize merged normalized values
  merge_vals <- df_merge
  merge_norm <- matrix(nrow = nrow(merge_vals), ncol = ncol(merge_vals))
  for (i in 1:nrow(merge_vals)) {
    merge_norm[i,] <- NormalizeMaxMin(merge_vals[i,])
  } #endfor
  
  #name rows and columns
  colnames(merge_norm) <- colnames(df1)
  ornames <- rownames(df1)
  rownames(merge_norm) <- ornames
  return(merge_norm)
} #end normNormDF3

#given a set of normalized values, merge them into 1 df, renormalize, sort, return sortlist
MakeSort4 <- function(df1, df2, df3, df4, method = "posmean") {
  merge_norm <- normNormDF4(df1, df2, df3, df4)
  merge_bin <- merge_norm
  ornames <- rownames(df1)
  
  if (method == "kzsort") {
    #sort based on position of maximum expression section and 
    #distance to the mean of positions for the three highest expression sections
    #higher distances being placed further away, no management for ties.
    #create a maxsec variable
    #anything lower than 1 becomes 0
    merge_bin[merge_bin < 1] <- 0
    maxsec <- vector(mode = "numeric", length = nrow(merge_bin))
    for (i in 1:nrow(merge_bin)) {
      for (s in 1:ncol(merge_bin)) {
        if (merge_bin[i,s] == 1) {
          maxsec[i] <- s
        } #endif
      } #endfor
    } #endfor
    
    #find avg position of top 3 sections
    avgpos3 <- vector(mode = "numeric", length = nrow(merge_bin))
    for (i in 1:nrow(merge_bin)) {
      avgpos3[i] <- mean(which(min_rank(desc(merge_norm[i,])) <= 3))
    } #endfor
    
    max2avg <- vector(mode = "numeric", length = nrow(merge_bin))
    for (i in 1:length(maxsec)) {
      max2avg[i] <- abs(maxsec[i] - avgpos3[i])
    } #endfor
    
    sorttable <- tibble(gene = ornames, maxsec, avgpos3, max2avg) %>% 
      arrange(desc(maxsec, max2avg)) %>% 
      mutate(sortrank = 1:length(ornames))
    
  } else if (method == "posmean") {
    #weighted avg based on python code for "find_1_glom"
    #weighted_avg = sum(x*y for x,y in zip(indexes,vals))/sum(vals) 
    #zip(c(1,2,3), c("a","b","c")) = (1,"a"), (2,"b"), (3,"c")
    secs <- c(1:23)
    wavgs <- vector(mode = "numeric", length = nrow(merge_bin))
    for (i in 1:nrow(merge_norm)) {
      vals <- merge_norm[i,]
      sumvals <- sum(vals)
      prods <- vector(mode = "numeric", length = length(vals))
      for (s in 1:length(vals)) {
        prods[s] <- vals[s] * secs[s]
      }
      wavgs[i] <- sum(prods)/sumvals
    }
    
    sorttable <- tibble(gene = ornames, wavgs) %>%
      arrange(desc(wavgs)) %>%
      mutate(sortrank = 1:length(ornames))
    
  } else if (method == "123") {
    #which secs are ranked 1,2,3
    onesec <- vector(mode = "numeric", length = nrow(merge_norm))
    for (i in 1:nrow(merge_norm)) {
      onesec[i] <- which(min_rank(desc(merge_norm[i,])) == 2)
    } #endfor
    twosec <- vector(mode = "numeric", length = nrow(merge_norm))
    for (i in 1:nrow(merge_norm)) {
      twosec[i] <- which(min_rank(desc(merge_norm[i,])) == 2)
    } #endfor
    threesec <- vector(mode = "numeric", length = nrow(merge_norm))
    for (i in 1:nrow(merge_norm)) {
      threesec[i] <- which(min_rank(desc(merge_norm[i,])) == 3)
    } #endfor
    
    sorttable <- tibble(gene = ornames, onesec, twosec, threesec) %>%
      arrange(desc(onesec, twosec, threesec)) %>%
      mutate(sortrank = 1:length(ornames))
  } #endif
  return(sorttable)
} #end MakeSort4

#given a set of normalized values, merge them into 1 df, renormalize, sort, return sortlist
MakeSort3 <- function(df1, df2, df3, method = "posmean") {
  merge_norm <- normNormDF3(df1, df2, df3)
  merge_bin <- merge_norm
  ornames <- rownames(df1)
  
  if (method == "kzsort") {
    #sort based on position of maximum expression section and 
    #distance to the mean of positions for the three highest expression sections
    #higher distances being placed further away, no management for ties.
    #create a maxsec variable
    #anything lower than 1 becomes 0
    merge_bin[merge_bin < 1] <- 0
    maxsec <- vector(mode = "numeric", length = nrow(merge_bin))
    for (i in 1:nrow(merge_bin)) {
      for (s in 1:ncol(merge_bin)) {
        if (merge_bin[i,s] == 1) {
          maxsec[i] <- s
        } #endif
      } #endfor
    } #endfor
    
    #find avg position of top 3 sections
    avgpos3 <- vector(mode = "numeric", length = nrow(merge_bin))
    for (i in 1:nrow(merge_bin)) {
      avgpos3[i] <- mean(which(min_rank(desc(merge_norm[i,])) <= 3))
    } #endfor
    
    max2avg <- vector(mode = "numeric", length = nrow(merge_bin))
    for (i in 1:length(maxsec)) {
      max2avg[i] <- abs(maxsec[i] - avgpos3[i])
    } #endfor
    
    sorttable <- tibble(gene = ornames, maxsec, avgpos3, max2avg) %>% 
      arrange(desc(maxsec, max2avg)) %>% 
      mutate(sortrank = 1:length(ornames))
    
  } else if (method == "posmean") {
    #weighted avg based on python code for "find_1_glom"
    #weighted_avg = sum(x*y for x,y in zip(indexes,vals))/sum(vals) 
    #zip(c(1,2,3), c("a","b","c")) = (1,"a"), (2,"b"), (3,"c")
    secs <- c(1:23)
    wavgs <- vector(mode = "numeric", length = nrow(merge_bin))
    for (i in 1:nrow(merge_norm)) {
      vals <- merge_norm[i,]
      sumvals <- sum(vals)
      prods <- vector(mode = "numeric", length = length(vals))
      for (s in 1:length(vals)) {
        prods[s] <- vals[s] * secs[s]
      } #endfor
      wavgs[i] <- sum(prods)/sumvals
    } #endfor
    
    sorttable <- tibble(gene = ornames, wavgs) %>%
      arrange(desc(wavgs)) %>%
      mutate(sortrank = 1:length(ornames))
    
    
  } else if (method == "123") {
    #which secs are ranked 1,2,3
    onesec <- vector(mode = "numeric", length = nrow(merge_norm))
    for (i in 1:nrow(merge_norm)) {
      onesec[i] <- which(min_rank(desc(merge_norm[i,])) == 2)
    } #endfor
    twosec <- vector(mode = "numeric", length = nrow(merge_norm))
    for (i in 1:nrow(merge_norm)) {
      twosec[i] <- which(min_rank(desc(merge_norm[i,])) == 2)
    } #endfor
    threesec <- vector(mode = "numeric", length = nrow(merge_norm))
    for (i in 1:nrow(merge_norm)) {
      threesec[i] <- which(min_rank(desc(merge_norm[i,])) == 3)
    } #endfor
    
    sorttable <- tibble(gene = ornames, onesec, twosec, threesec) %>%
      arrange(desc(onesec, twosec, threesec)) %>%
      mutate(sortrank = 1:length(ornames))
  } #endif
  return(sorttable)
} #end MakeSort3

#given a normalized matrix, output old style Black(0) to Red(1) heatmaps
MakeBlackRedHeatmap <- function(df_norm, title = NA) {
  df_df <- as.data.frame(df_norm)
  df_df$target = rownames(df_df)
  df_melt <- melt(df_df, id.vars ="target")
  names(df_melt)[2:3] <- c("section","normTPM")
  #make genelist into character vector
  df_melt$target <- as.character(df_melt$target)
  df_melt$target <- factor(df_melt$target, levels=unique(df_melt$target))
  
  if (is.na(title)) {
    #HEAT DA MAP generic
    plot <- ggplot(df_melt, aes(section, target, height = 3)) + 
      geom_tile(aes(fill = normTPM), color = "grey") +
      scale_fill_gradient2(low = "black", mid = "red", high = "red3", midpoint = 0.65) +
      ylab("Olfr Genes") +
      theme(legend.title = element_text(size = 11, vjust = 1),
            legend.text = element_text(size = 7),
            legend.position = "none", #none if you dont want expression level legend
            plot.title = element_text(size = 24),
            axis.title.x = element_text(size = 0, vjust = 50, face="bold"),
            axis.title.y = element_text(size = 0, vjust = 0, face="bold"),
            axis.text.x = element_text(angle = 330, hjust = 0, vjust = 2, size = 0),
            axis.text.y = element_text(size = 0),
            axis.ticks = element_blank()) +
      labs(fill = "Expression Level")
    return(plot)
  } else {
    #HEAT DA MAP generic
    plot <- ggplot(df_melt, aes(section, target, height = 3)) + 
      geom_tile(aes(fill = normTPM), color = "grey") +
      scale_fill_gradient2(low = "black", mid = "red", high = "red3", midpoint = 0.65) +
      ylab("Olfr Genes") +
      theme(legend.title = element_text(size = 11, vjust = 1),
            legend.text = element_text(size = 7),
            legend.position = "none", #none if you dont want expression level legend
            plot.title = element_text(size = 24),
            axis.title.x = element_text(size = 0, vjust = 50, face="bold"),
            axis.title.y = element_text(size = 0, vjust = 0, face="bold"),
            axis.text.x = element_text(angle = 330, hjust = 0, vjust = 2, size = 0),
            axis.text.y = element_text(size = 0),
            axis.ticks = element_blank()) +
      labs(fill = "Expression Level") + 
      ggtitle(title)
    return(plot)
  } #endif
} #end MakeBlackRedHeatmap

#Make Black Red Heatmaps based on branch subsets
MakeBLRHtree3 <- function(df1, df2, df3, tree, title = NA) {
  treeset1 <- df1[which(rownames(df1) %in% tree),]
  treeset2 <- df2[which(rownames(df2) %in% tree),]
  treeset3 <- df3[which(rownames(df3) %in% tree),]
  
  treesort <- MakeSort3(treeset1, treeset2, treeset3, method = "posmean")
  
  merge <- normNormDF3(treeset1, treeset2, treeset3)
  mergesort <- SortByList(merge, treesort)
  
  p <- MakeBlackRedHeatmap(mergesort, title)
  return(p)
} #end MakeBLRHtree3

MakeBLRHtree4 <- function(df1, df2, df3, df4, tree, title = NA, out = "heat") {
  treeset1 <- df1[which(rownames(df1) %in% tree),]
  treeset2 <- df2[which(rownames(df2) %in% tree),]
  treeset3 <- df3[which(rownames(df3) %in% tree),]
  treeset4 <- df4[which(rownames(df4) %in% tree),]
  
  treesort <- MakeSort4(treeset1, treeset2, treeset3, treeset4, method = "posmean")
  
  merge <- normNormDF4(treeset1, treeset2, treeset3, treeset4)
  mergesort <- SortByList(merge, treesort)
  
  p <- MakeBlackRedHeatmap(mergesort, title)
  if (out == "heat") {
    return(p)
  } else if (out == "names") {
  return(treesort)
  } #endif
}#end MakeBLRHtree4

#given 2 matrixes, compute spearman correlation of rows (OR vs OR) 
SpearMatrix <- function(df1, df2) {
  spear_vec <- vector(length = nrow(df1), mode = "numeric")
  for (i in 1:nrow(df1)) {
    spear_vec[i] <- cor(df1[i,], df2[i,], method = "spearman")
  } #endfor
  return(spear_vec)
} #end SpearMatrix

#apply SpearMatrix to multiple matrices(just 3 or 4)
#probably need to rewrite using permutation/combination
MultiSpearMatrix <- function(...){
  arguments <- list(...)
  if (length(arguments) == 3) {
    mat3 <- matrix(ncol = 3, nrow = nrow(arguments[[1]]))
    colnames(mat3) <- c("1v2", "1v3", "2v3")
    ornames <- rownames(arguments[[1]])
    rownames(mat3) <- ornames
    mat3[,1] <- SpearMatrix(arguments[[1]], arguments[[2]])
    mat3[,2] <- SpearMatrix(arguments[[1]], arguments[[3]])
    mat3[,3] <- SpearMatrix(arguments[[2]], arguments[[3]])
    mat3
  } else if(length(arguments) == 4) {
    mat4 <- matrix(ncol = 6, nrow = nrow(arguments[[1]]))
    colnames(mat4) <- c("1v2", "1v3", "1v4", "2v3", "2v4", "3v4")
    ornames <- rownames(arguments[[1]])
    rownames(mat4) <- ornames
    mat4[,1] <- SpearMatrix(arguments[[1]], arguments[[2]])
    mat4[,2] <- SpearMatrix(arguments[[1]], arguments[[3]])
    mat4[,3] <- SpearMatrix(arguments[[1]], arguments[[4]])
    mat4[,4] <- SpearMatrix(arguments[[2]], arguments[[3]])
    mat4[,5] <- SpearMatrix(arguments[[2]], arguments[[4]])
    mat4[,6] <- SpearMatrix(arguments[[3]], arguments[[4]])
    mat4
  } #endif
} #end MultiSpearMatrix

#Find an Anterior Peak and a Posterior Peak
#note that the rank df solves ties using "random" in order to prevent passing a vector when searching for rank X.
AntPeakPostPeak <- function(df_norm) {
  #rank df
  rank_df <- matrix(nrow = nrow(df_norm), ncol = ncol(df_norm))
  for (i in 1:nrow(df_norm)) {
    rank_df[i,] <- rank(desc(df_norm[i,]), ties.method = "random")
  } #endfor
  
  #find max sec
  maxsec <- vector(mode = "numeric", length = nrow(rank_df))
  maxval <- vector(mode = "numeric", length = nrow(rank_df))
  for (i in 1:nrow(rank_df)) {
    for (s in 1:ncol(rank_df)) {
      if (rank_df[i,s] == min(rank_df[i,])) {
        maxsec[i] <- s
        maxval[i] <- df_norm[i,s]
      } #endif
    } #endfor
  } #endfor
  
  #find next
  nextsec <- rep(NA, length(maxsec))
  nextval <- vector(mode = "numeric", length = nrow(rank_df))
  nextmany <- vector(mode = "numeric", length = nrow(rank_df))
  maxsecin <- ncol(rank_df)
  for (i in 1:length(maxsec)) {
    #set findrank
    findrank <- 2
    #set neighbors
    if (maxsec[i] == 1) {
      min_neigh <- 1
      max_neigh <- 2
    } else if (maxsec[i] == maxsecin) {
      min_neigh <- maxsecin - 1
      max_neigh <- maxsecin
    } else {
      min_neigh <- maxsec[i] - 1
      max_neigh <- maxsec[i] + 1
    } #endif
    
    #while we have not found nextsec
    while (is.na(nextsec[i])) {
      potentialsec <- which(rank_df[i,] == findrank)
      #if potentialsec is a neighbor, bump findrank and neighbor
      if (between(potentialsec, min_neigh, max_neigh)) {
        findrank <- findrank + 1
        #bump neighbor
        if (potentialsec == min_neigh) {
          min_neigh <- min_neigh - 1
        } else if (potentialsec == max_neigh) {
          max_neigh <- max_neigh + 1
        } #endif
      } else if (between(potentialsec, min_neigh, max_neigh) == FALSE) {
        nextsec[i] <- potentialsec
        nextval[i] <- df_norm[i, potentialsec]
      } else {
        nextsec[i] <- -1
        nextval[i] <- -1
      } #endif
    } #endwhile
  } #endfor
  
  #notate Ant and Post Peaks
  antpeak <- vector(mode = "numeric", length = nrow(rank_df))
  postpeak <- vector(mode = "numeric", length = nrow(rank_df))
  antval <- vector(mode = "numeric", length = nrow(rank_df))
  postval <- vector(mode = "numeric", length = nrow(rank_df))
  for (i in 1:nrow(rank_df)) {
    if (maxsec[i] < nextsec[i]) {
      antpeak[i] <- maxsec[i]
      antval[i] <- maxval[i]
      postpeak[i] <- nextsec[i]
      postval[i] <- nextval[i]
    } else {
      antpeak[i] <- nextsec[i]
      antval[i] <- nextval[i]
      postpeak[i] <- maxsec[i]
      postval[i] <- maxval[i]
    } #endif
  } #endfor
  
  #output
  gene <- rownames(df_norm)
  record <- tibble(gene, antpeak, postpeak, antval, postval)
  return(record)
} #end AntPeakPostPeak

#join APPP output, names should be like "_mX"
MultiAPPP <- function(df1, df2, df3, df4, name1, name2, name3, name4) {
  rec1 <- AntPeakPostPeak(df1)
  rec2 <- AntPeakPostPeak(df2)
  rec3 <- AntPeakPostPeak(df3)
  rec4 <- AntPeakPostPeak(df4)
  join1 <- left_join(rec1, rec2, by = "gene", suffix = c(name1, name2))
  join2 <- left_join(rec3, rec4, by = "gene", suffix = c(name3, name4))
  endjoin <- left_join(join1, join2, by = "gene")
  return(endjoin)
} #end MultiAPPP

#using APPP, take a mean position of peaks
SymLiner <- function(ap, vd, ml, out = "plot") {
  ap_appp <- AntPeakPostPeak(ap)
  vd_appp <- AntPeakPostPeak(vd) %>% 
    rename(venpeak = antpeak,
           dorpeak = postpeak,
           venval = antval,
           dorval = postval)
  ml_appp <- AntPeakPostPeak(ml) %>% 
    rename(medpeak = antpeak,
           latpeak = postpeak,
           medval = antval,
           latval = postval)
  allpeaks <- left_join(ap_appp, vd_appp, by = "gene") %>% 
    left_join(ml_appp, by = "gene")
  
  peakdifs <- allpeaks %>%
    rowwise() %>%
    mutate(apdif = mean(c(antpeak, postpeak)),
           vddif = mean(c(venpeak, dorpeak)),
           mldif = mean(c(medpeak, latpeak))) %>%
    ungroup()
  
  if (out == "plot") {
    plot <- ggplot(peakdifs) +
      geom_point(aes(apdif, mldif, alpha = 0.1)) +
      geom_smooth(aes(apdif, mldif), 
                  method = "lm", formula = y~x, color = "blue") +
      geom_smooth(aes(apdif, mldif), 
                  method = "loess", formula = y~x, 
                  color = "red", se = F) +
      theme(legend.position = "none")
    return(plot)
  } else if (out == "fit") {
      fit <- coef(lm(peakdifs$vddif ~ peakdifs$apdif))
      fitout <- summary(fit)
      return(fit)
  }  else {
    return(peakdifs)
  } #endif
} #end SymLiner

#Make triple heatmap of class, tanzone, and matsunami diff OE
PlotFeatures <- function(sortlistin, featureOut = "trihm") {
  withinfo <- sortlistin %>% 
    left_join(info %>% rename("gene" = "olfrname"), by = "gene") %>% 
    mutate(sortearly = ifelse(sortrank <= max(sortrank)/2, "early", "late"),
           class_fct = as_factor(class),
           rtp_fct = as_factor(RTP))
  
  stats_oe <- withinfo %>% 
    filter(!is.na(oe_region)) %>% 
    select(sortearly, oe_region) %>% 
    group_by(sortearly) %>% 
    count(oe_region) %>% 
    arrange(desc(n)) %>%
    pivot_wider(names_from = sortearly, values_from = n) %>% 
    as.matrix() 
  stats_oe <- stats_oe[,-1]
  class(stats_oe) <- "numeric"
  test_oe <- fisher.test(stats_oe)
  pval_oe <- format(test_oe$p.value, digits = 5)
  odds_oe <- format(test_oe$estimate, digits = 3)
  
  stats_class <- withinfo %>% 
    filter(!is.na(class)) %>% 
    select(sortearly, class) %>% 
    group_by(sortearly) %>% 
    count(class) %>% 
    arrange(desc(n)) %>%
    pivot_wider(names_from = sortearly, values_from = n) %>% 
    as.matrix()
  stats_class <- stats_class[,-1]
  class(stats_class) <- "numeric"
  test_class <- fisher.test(stats_class)
  pval_class <- format(test_class$p.value, digits = 5)
  odds_class <- format(test_class$estimate, digits = 3)
  
  stats_fi <- withinfo %>% 
    filter(!is.na(fisurface)) %>% 
    select(sortearly, fisurface) %>% 
    group_by(sortearly) %>% 
    count(fisurface) %>% 
    arrange(desc(n)) %>%
    pivot_wider(names_from = sortearly, values_from = n) %>% 
    as.matrix()
  stats_fi <- stats_fi[,-1]
  class(stats_fi) <- "numeric"
  test_fi <- fisher.test(stats_fi)
  pval_fi <- format(test_fi$p.value, digits = 5)
  odds_fi <- format(test_fi$estimate, digits = 3)
  
  stats_tanvd <- withinfo %>% 
    filter(!is.na(tz_vd)) %>% 
    select(sortearly, tz_vd) %>% 
    group_by(sortearly) %>% 
    count(tz_vd) %>% 
    arrange(desc(n)) %>%
    pivot_wider(names_from = sortearly, values_from = n) %>% 
    as.matrix()
  stats_tzvd <- stats_tanvd[,-1]
  class(stats_tzvd) <- "numeric"
  test_tzvd <- fisher.test(stats_tzvd)
  pval_tzvd <- format(test_tzvd$p.value, digits = 5)
  odds_tzvd <- format(test_tzvd$estimate, digits = 3)
  
  stats_rtp <- withinfo %>% 
    filter(!is.na(RTP)) %>% 
    filter(RTP != "ns") %>%
    select(sortearly, RTP) %>% 
    group_by(sortearly) %>% 
    count(RTP) %>% 
    arrange(desc(n)) %>%
    pivot_wider(names_from = sortearly, values_from = n) %>% 
    as.matrix() 
  stats_rtp <- stats_rtp[,-1]
  class(stats_rtp) <- "numeric"
  test_rtp <- fisher.test(stats_rtp)
  pval_rtp <- format(test_rtp$p.value, digits = 5)
  odds_rtp <- format(test_rtp$estimate, digits = 3)
  
  sorttanzone <- withinfo %>% 
    select(gene, tz_val) %>% 
    filter(!is.na(tz_val)) %>% 
    as.matrix()
  rownames(sorttanzone) <- sorttanzone[,1]
  trimtanzone <- sorttanzone[,2]
  trimtanzone <- as.numeric(trimtanzone)
  
  df_df <- as.data.frame(trimtanzone)
  df_df$target = rownames(df_df)
  df_melt <- melt(df_df, id.vars ="target")
  names(df_melt)[2:3] <- c("Feature", "value")
  #make genelist into character vector
  df_melt$target <- as.character(df_melt$target)
  df_melt$target <- factor(df_melt$target, levels=unique(df_melt$target))
  
  tan <- ggplot(df_melt, aes(Feature, target, height = 3)) + 
    geom_tile(aes(fill = value), color = "grey") +
    scale_fill_viridis() +
    guides(fill = guide_legend(nrow = 1)) +
    theme(plot.margin = unit(c(1,1,1,1), "cm"),
          legend.title = element_text(size = 12, vjust = 1),
          legend.text = element_text(size = 6),
          legend.position = c(0.5, -0.1),
          plot.title = element_text(size = 20),
          axis.title.x = element_blank(),
          axis.title.y = element_text(size = 15, vjust = 0, face="bold"),
          axis.text.x = element_text(angle = 330, hjust = 0, vjust = 2, size = 0),
          axis.text.y = element_text(size = 0),
          axis.ticks = element_blank()) +
    labs(fill = "MiyIndex")
  
  tan_vd <- ggplot(withinfo %>% filter(!is.na(tz_vd)),
                  aes(sortrank, tz_vd, fill = tz_vd)) +
    geom_violin() +
    coord_flip() +
    ylab(paste("Tan Zone Dorsal Ventral", 
               paste0("p = ", pval_oe), 
               paste0("odds = ", odds_oe),
               sep = "\n")) +
    scale_fill_manual(values=c("red", "black")) +
    theme_cowplot() +
    theme(legend.position = "none",
          axis.line = element_blank(),
          axis.ticks.x = element_blank(),
          axis.ticks.y = element_blank(),
          axis.text.y = element_blank(),
          axis.title.y = element_blank())
  
  oe <- ggplot(withinfo %>% filter(!is.na(oe_region)), 
               aes(sortrank, oe_region, fill = oe_region)) +
    geom_violin() +
    coord_flip() +
    ylab(paste("Matsunami OE DiffE", 
               paste0("p = ", pval_tzvd), 
               paste0("odds = ", odds_tzvd),
               sep = "\n")) +
    scale_fill_manual(values=c("red", "black")) +
    theme_cowplot() +
    theme(legend.position = "none",
          axis.line = element_blank(),
          axis.ticks.x = element_blank(),
          axis.ticks.y = element_blank(),
          axis.text.y = element_blank(),
          axis.title.y = element_blank())
  
  class <- ggplot(withinfo %>% filter(!is.na(class_fct)), 
                  aes(sortrank, class_fct, fill = class_fct)) +
    geom_violin() +
    coord_flip() +
    ylab(paste("Class", 
               paste0("p = ", pval_class), 
               paste0("odds = ", odds_class),
               sep = "\n")) +
    scale_fill_manual(values=c("red", "black")) +
    theme_cowplot() +
    theme(legend.position = "none",
          axis.line = element_blank(),
          axis.ticks.x = element_blank(),
          axis.ticks.y = element_blank(),
          axis.text.y = element_blank(),
          axis.title.y = element_blank())
  
  FIsurface <- ggplot(withinfo %>% filter(!is.na(fisurface)),
                      aes(sortrank, fisurface, fill = fisurface)) +
    geom_violin() +
    coord_flip() +
    ylab(paste("FI surface enriched", 
               paste0("p = ", pval_fi), 
               paste0("odds = ", odds_fi),
               sep = "\n")) +
    scale_fill_manual(values=c("red", "black")) +
    theme_cowplot() +
    theme(legend.position = "none",
          axis.line = element_blank(),
          axis.ticks.x = element_blank(),
          axis.ticks.y = element_blank(),
          axis.text.y = element_blank(),
          axis.title.y = element_blank())
  
  #typically used as a quick control
  RTP <- ggplot(withinfo %>% 
                  filter(!is.na(RTP)) %>% 
                  filter(RTP != "ns"),
                      aes(sortrank, rtp_fct, fill = rtp_fct)) +
    geom_violin() +
    coord_flip() +
    ylab(paste("RTP u/oOR", 
               paste0("p = ", pval_rtp), 
               paste0("odds = ", odds_rtp),
               sep = "\n")) +
    scale_fill_manual(values=c("red", "black", "blue")) +
    theme_cowplot() +
    theme(legend.position = "none",
          axis.line = element_blank(),
          axis.ticks.x = element_blank(),
          axis.ticks.y = element_blank(),
          axis.text.y = element_blank(),
          axis.title.y = element_blank())
  
  if (featureOut == "trihm") {
    trihm <- tan + oe + class
    return(trihm)
  } else if (featureOut == "trivio") {
    trivio <- tan_vd + oe + class
    return(trivio)
  } else if (featureOut == "tan1") {
    return(tan)
  } else if (featureOut == "tan2") {
    return(tanvd)
  } else if (featureOut == "rtp") {
    return(RTP)
  } else if (featureOut == "oe") {
    return(oe)
  } else if (featureOut == "class") {
    return(class)
  } else if (featureOut == "fi") {
    return(FIsurface)
  } else if (featureOut == "quadbox") {
    quadbox <- (tan_vd + oe) / (class + FIsurface)
    return(quadbox)
  } #endif 
} #end PlotFeatures

#Plot two sets of random ORs as a mental control for other features
PlotRandom <- function(sortlistin, x = 200, y = 225, 
                       out = "plot") {
  withinfo <- sortlistin %>% 
    left_join(info %>% 
                rename("gene" = "olfrname"), by = "gene") %>% 
    mutate(sortearly = ifelse(sortrank <= max(sortrank)/2, "early", "late"))
  
  #get random OR sets and give them a feature
  random_names <- withinfo %>% pull(gene)
  getxy <- sample(random_names, size = x+y)
  xnames <- getxy[1:x]
  ynames <- getxy[(x+1):y]
  
  random_set <- withinfo %>%
    mutate(inXinY = ifelse(gene %in% xnames, "X", 
                           ifelse(gene %in% ynames, "Y", NA))) %>%
    filter(!is.na(inXinY))
  
  stats_random <- random_set %>% 
    select(sortearly, inXinY) %>% 
    group_by(sortearly) %>% 
    count(inXinY) %>% 
    arrange(desc(n)) %>%
    pivot_wider(names_from = sortearly, values_from = n) %>% 
    as.matrix()
  stats_random <- stats_random[,-1]
  class(stats_random) <- "numeric"
  test_random <- fisher.test(stats_random)
  pval_random <- format(test_random$p.value, digits = 5)
  odds_random <- format(test_random$estimate, digits = 3)
  random_stat <- c(pval_random, odds_random)
  
  random_violin <- ggplot(random_set, 
                          aes(inXinY, sortrank, fill = inXinY)) +
    geom_violin() +
    ylab(paste("Randomly picked ORs", 
               paste0("p = ", pval_random), 
               paste0("odds = ", odds_random),
               sep = "\n")) +
    scale_fill_manual(values=c("red", "black")) +
    theme_cowplot() +
    theme(legend.position = "none",
          axis.line = element_blank(),
          axis.ticks.x = element_blank(),
          axis.ticks.y = element_blank(),
          axis.text.y = element_blank(),
          axis.title.y = element_blank())
  
  if (str_detect(tolower(out), "stat"))
    return(random_stat)
  else {
    return(random_violin)
  } #endif
} #end PlotRandom

#check if expression of OR across sections is uniform in distribution using ks.test(x, "punif)$p.value > 0.05
UniformityFilter <- function(df, ranklist = NA, rounds = 10, out = "plot") {
  kstest_df <- matrix(nrow = nrow(df), ncol = 2*rounds)
  colnames(kstest_df) <- rep(c("Pval", "Dval"), each = rounds)
  
  for (r in 1:rounds) {
    #deal with ties for ks.test() input
    pseudocount <- runif(ncol(df), 0, 0.0001)
    
    for (i in 1:nrow(df)) {
      orexp <- df[i,]
      orexpPC <- orexp + pseudocount
      test <- ks.test(orexpPC, "punif")
      kstest_df[i,r] <- test$p.value
      kstest_df[i,rounds + r] <- test$statistic
    } #endfor
  } #endfor
  
  #distribution is uniform if kstest p value is not significant
  stat_all <- suppressMessages(as_tibble(kstest_df, .name_repair = "unique"))
  pmean <- stat_all %>% 
    select(contains("Pval")) %>%
    mutate(meanPval = rowMeans(.)) %>%
    pull(meanPval)
  
  dmean <- stat_all %>%
    select(contains("Dval")) %>%
    mutate(meanDval = rowMeans(.)) %>%
    pull(meanDval)

  stat_tib <- tibble(gene = rownames(df),
                   meanPval = pmean,
                   meanDval = dmean,
                   isUniform = ifelse(meanPval > 0.05, T, F))
  
  uniform_ORs <- stat_tib %>% filter(isUniform == T) %>% pull(gene)
  notuniform_ORs <- stat_tib %>% filter(isUniform == F) %>% pull(gene)
  
  uniform_len <- length(uniform_ORs)
  notuniform_len <- length(notuniform_ORs)
  
  uniform_df <- df[which(stat_tib$isUniform == T),]
  notuniform_df <- df[-which(stat_tib$isUniform == T),]
  
  if (str_detect(tolower(out), "plot")) {
    if (nrow(uniform_df) == 0) {
      return(stat_tib)
    }  #endif
    
    uniform_sort <- SortByList(uniform_df, ranklist)
    notuniform_sort <- SortByList(notuniform_df, ranklist)
    uniform_hm <- MakeBlackRedHeatmap(uniform_sort, 
                                      title = paste0(uniform_len, 
                                                     " Uniform ORs"))
    notuniform_hm <- MakeBlackRedHeatmap(notuniform_sort, 
                                         title = paste0(notuniform_len, 
                                                        " Non-Uniform ORs"))
    both_hm <- uniform_hm + notuniform_hm
    return(both_hm)
  } else if (str_detect(tolower(out), "stat")) {
    return(stat_tib)
  } else {
    return(uniform_ORs)
  } #endif
} #end UniformityFilter

#function for distance between 3d points
Dist3d <- function(df, row1, row2) {
  distance <- sqrt((df$AntPos[row1] - df$AntPos[row2])^2 +
                  (df$MedLat[row1] - df$MedLat[row2])^2 + 
                  (df$VenDor[row1] - df$VenDor[row2])^2)
  return(distance)
} #end Dist3d

#data for correlation matrix split into X many tree groups
PlotBranches <- function(branches) {
  dend <- hclust(dist(meanpos_cor, method = "euclidean"), method = "complete")
  
  #2d stuff
  branchheatmaps <- vector("list", length = branches)
  branchsize <- vector("numeric", length = branches)
  branchcor <- vector("numeric", length = branches)
  #3d stuff
  branch3d <- vector("list", length = branches)
  branchMdist <- vector("list", length = branches)
  branchLdist <- vector("list", length = branches)
  branchavgdist <- vector("numeric", length = branches)
  
  for (i in 1:branches) {
    if (i %% 50 == 0) {
      message(i)
    } #endif
    names <- names(which(cutree(dend, k = branches) == i))
    
    #remove ectopic
    ectopic <- c("Olfr32", "Olfr287")
    #checking for 1 at a time due to warning: the condition has length > 1 and only first element will be used
    if (max(ectopic %in% names) == 1) {
      names <- names[-which(names %in% ectopic)]
    } #endif
    
    branchsize[i] <- length(names)
    
    if(branchsize[i] > 1){
      title_ap <- paste("AP-branch", as.character(i), sep = "")
      title_ml <- paste("ML-branch", as.character(i), sep = "")
      title_vd <- paste("VD-branch", as.character(i), sep = "")
      
      ap_plot <- MakeBLRHtree4(m4_Vnorm, m6_Vnorm, m7_Vnorm, m10_Vnorm, 
                               names, title = title_ap)
      ml_plot <- MakeBLRHtree3(m8_Vnorm, m11_Vnorm, m16_Vnorm, 
                               names, title = title_ml)
      vd_plot <- MakeBLRHtree3(m9_Vnorm, m12_Vnorm, m14_Vnorm, 
                               names, title = title_vd)
      branchheatmaps[[i]] <- ap_plot + ml_plot + vd_plot
      
      #branchcor mean includes 1.0 values from ORx vs ORx
      branchcor[i] <- mean(meanpos_cor[which(rownames(meanpos_cor) %in% names),
                                       which(colnames(meanpos_cor) %in% names)])
      branch3d[[i]] <- ListML(names, "point")
      
      #calculate distance between points for ORs in group
      data4dist <- ListDorML(names, out = "data") %>% 
        filter(p50 == clustmaxp)
      dataMdist <- data4dist %>% 
        filter(side == "Medial") %>%
        unique()
      dataLdist <- data4dist %>%
        filter(side == "Lateral") %>%
        unique()
      permM <- permutations(n=nrow(dataMdist), r=2)
      permL <- permutations(n=nrow(dataLdist), r=2)
      someM <- vector("numeric", length = nrow(permM))
      someL <- vector("numeric", length = nrow(permL))
      
      for (m in 1:nrow(permM)) {
        someM[m] <- Dist3d(dataMdist, permM[m,1], permM[m,2])
      } #endfor
      
      for (l in 1:nrow(permL)) {
        someL[l] <- Dist3d(dataLdist, permL[l,1], permL[l,2])
      } #endfor
      
      branchMdist[[i]] <- someM
      branchLdist[[i]] <- someL
      branchavgdist[i] <- mean(c(someM, someL))
      
    } else {
      branchcor[i] <- NA
      branchheatmaps[[i]] <- NA
      branchMdist[[i]] <- NA
      branchLdist[[i]] <- NA
      branchavgdist[i] <- NA
      branch3d[[i]] <- DorsalML(olfr = names)
    } #endif
  } #endfor
  
  branchout <- list(branchsize, branchcor, branchheatmaps, 
                    branchMdist, branchLdist, branchavgdist, branch3d)
  names(branchout) <- c("size", "cor", "heatmaps", 
                        "med_dist", "lat_dist", "avg_dist", "3D")
  return(branchout)
} #end PlotBranches

SaveHeatmaps <- function(..., path = NA) {
  arguments <- list(...)
  argnames <- as.list(match.call())
  for (i in 1:length(arguments)) {
    ggsave(plot = arguments[[i]], 
           filename = paste0(path, "/", argnames[i + 1], ".png"),
           width = 6.41, height = 6.8, units = "in")
    print(paste0("heatmap saved: ", path, "/", argnames[i + 1], ".png"))
  }
}

SaveSquarePlots <- function(..., path = NA) {
  arguments <- list(...)
  argnames <- as.list(match.call())
  for (i in 1:length(arguments)) {
    ggsave(plot = arguments[[i]], 
           filename = paste0(path, "/", argnames[i + 1], ".png"),
           width = 6.6, height = 6.6, units = "in")
    print(paste0("plot saved: ", path, "/", argnames[i + 1], ".png"))
  }
}

#load data
allmice_tpm <- read_csv("~/Desktop/obmap/r_analysis/3dimOB/input/covarintactchemo_over_samples_200923.csv") %>% select(-X1)
allmice_covars <- read_csv("~/Desktop/obmap/r_analysis/3dimOB/input/allmice_covariates_trim_voxweights_v3.csv")
symline <- read_csv("~/Desktop/obmap/r_analysis/3dimOB/input/symline.csv")
good_3dpoint <- read_csv("~/Desktop/obmap/r_analysis/3dimOB/output/goodpointORs972_allchemo_listdorML_200820.csv")
info <- read_csv("~/Desktop/obmap/r_analysis/3dimOB/input/info_201129.csv") %>%
  rename("olfrname" = "gene")

set.seed(711)
```

# mean tpm across all sections

``` r
allmice_no16 <- allmice_tpm %>% filter(rep != 16) %>% arrange(rep)
allmice_ortaar <- allmice_no16[,1:1135]

#calculate average TPM of an OR across all sections
or_mean_all <- vector(length = ncol(allmice_ortaar), mode = "numeric")
for (i in 1:ncol(allmice_ortaar)) {
  if (i < 6) {
    or_mean_all[i] <- NA
  } else {
  or_mean_all[i] <- sum(allmice_ortaar[,i])/nrow(allmice_ortaar)
  }
}

#plot how many ORs below each TPM cutoff
ors_below <- vector(length = 100, mode = "numeric")
below_num <- vector(length = 100, mode = "numeric")
for (i in 1:100) {
  ors_below[i] <- length(which(or_mean_all < 2*i))
  below_num[i] <- 2*i
}
ors_below_num <- data.frame(below_num, ors_below)

ggplot(ors_below_num, aes(below_num, ors_below)) + 
  geom_point()
```

![](ob_heatmaps_files/figure-gfm/unnamed-chunk-1-1.png)<!-- -->

# filter ORs by: uniformity, zero in a replicate, mean of all sections &lt; 1tpm

``` r
#some ORs in some replicates have only 0 tpm
mean_byrep <- matrix(nrow = ncol(allmice_ortaar) - 5, 
                     ncol = length(unique(allmice_ortaar$rep)))
reps <- unique(allmice_ortaar$rep)
for (i in 1:length(reps)) {
  rep_tib <- allmice_ortaar %>% filter(rep == reps[i])
  for (j in 6:ncol(rep_tib)) {
    mean_byrep[j-5,i] <- sum(rep_tib[,j]/nrow(rep_tib))
  }
}
rownames(mean_byrep) <- colnames(allmice_ortaar)[6:ncol(allmice_ortaar)]
colnames(mean_byrep) <- reps

zero_inrep <- which(mean_byrep == 0, arr.ind = T)
zero_inrep_tb <- tibble(gene = rownames(zero_inrep),
                        rep = zero_inrep[,2])

#ectopic ORs
ectopic_ORs <- c("Olfr287", "Olfr32", "Olfr361")

#ORs with a mean below 1tpm in all sections
below1_ORs <- colnames(allmice_ortaar)[which(or_mean_all < 1)]

#ORs with only 0s in a replicate
zero_inrep_ORs <- unique(zero_inrep_tb$gene)

#ORs with uniform distribution in a replicate
uniform_ORs <- readRDS("~/Desktop/obmap/r_analysis/heatmaps/output/alldim_uniORs.RDS")

filter_ORs <- unique(c(ectopic_ORs, below1_ORs, zero_inrep_ORs, uniform_ORs))
length(filter_ORs)
```

    ## [1] 235

``` r
#saveRDS(filter_ORs, "~/Desktop/obmap/r_analysis/heatmaps/output/filterORs.RDS")

allmice_final <- allmice_ortaar %>% select(-all_of(filter_ORs))

#using 24 to start as it falls early on the curve and removes 350 ORs
#previous version used [,-which(or_mean_all < 10)] 
alldim_ap_tpm_cut <- allmice_final %>% filter(dim == "AntPos")
alldim_ml_tpm_cut <- allmice_final %>% filter(dim == "MedLat")
alldim_vd_tpm_cut <- allmice_final %>% filter(dim == "VenDor")
```

# ANTERIOR-POSTERIOR

``` r
#normalize values
m4_tpm <- alldim_ap_tpm_cut %>% filter(rep == 4) %>% select(-name, -rep, -dim, -dimrep)
m6_tpm <- alldim_ap_tpm_cut %>% filter(rep == 6) %>% select(-name, -rep, -dim, -dimrep)
m7_tpm <- alldim_ap_tpm_cut %>% filter(rep == 7) %>% select(-name, -rep, -dim, -dimrep)
m10_tpm <- alldim_ap_tpm_cut %>% filter(rep == 10) %>% select(-name, -rep, -dim, -dimrep)
m13_tpm <- alldim_ap_tpm_cut %>% filter(rep == 13) %>% select(-name, -rep, -dim, -dimrep)
m15_tpm <- alldim_ap_tpm_cut %>% filter(rep == 15) %>% select(-name, -rep, -dim, -dimrep)

m4_norm <- NormTPM(m4_tpm)
m6_norm <- NormTPM(m6_tpm)
m7_norm <- NormTPM(m7_tpm)
m10_norm <- NormTPM(m10_tpm)
m13_norm <- NormTPM(m13_tpm)
m15_norm <- NormTPM(m15_tpm)

#apply voxel factors
#make sort for voxel adjusted TPM
ap_voxweights <- allmice_covars %>% filter(rep == 10) %>% select(slice, weight)

m4_Vnorm <- Voxnormify(m4_tpm, ap_voxweights)
m6_Vnorm <- Voxnormify(m6_tpm, ap_voxweights)
m7_Vnorm <- Voxnormify(m7_tpm, ap_voxweights)
m10_Vnorm <- Voxnormify(m10_tpm, ap_voxweights)
m13_Vnorm <- Voxnormify(m13_tpm, ap_voxweights)
m15_Vnorm <- Voxnormify(m15_tpm, ap_voxweights)

#heatmaply_cor(m10_Vnorm, Rowv=T, Colv=F)

m46710_voxsort <- MakeSort4(m4_Vnorm, m6_Vnorm, m7_Vnorm, m10_Vnorm, method = "posmean")

m4_voxsort <- SortByList(m4_Vnorm, m46710_voxsort)
m6_voxsort <- SortByList(m6_Vnorm, m46710_voxsort)
m7_voxsort <- SortByList(m7_Vnorm, m46710_voxsort)
m10_voxsort <- SortByList(m10_Vnorm, m46710_voxsort)
m13_voxsort <- SortByList(m13_Vnorm, m46710_voxsort)
m15_voxsort <- SortByList(m15_Vnorm, m46710_voxsort)

#heatmaply_cor(m4_voxsort, Rowv = F, Colv = F)
#heatmaply_cor(m6_voxsort, Rowv = F, Colv = F)
#heatmaply_cor(m7_voxsort, Rowv = F, Colv = F)
#heatmaply_cor(m10_voxsort, Rowv = F, Colv = F)

ap4 <- MakeBlackRedHeatmap(m4_voxsort)
ap6 <- MakeBlackRedHeatmap(m6_voxsort)
ap7 <- MakeBlackRedHeatmap(m7_voxsort)
ap10 <- MakeBlackRedHeatmap(m10_voxsort)
ap13 <- MakeBlackRedHeatmap(m13_voxsort)
ap15 <- MakeBlackRedHeatmap(m15_voxsort)

#examine mergeDF used to make sort
m46710mergeV <- normNormDF4(m4_Vnorm, m6_Vnorm, m7_Vnorm, m10_Vnorm)
m46710mergeV_voxnormsort <- SortByList(m46710mergeV, m46710_voxsort)
apmerge <- MakeBlackRedHeatmap(m46710mergeV_voxnormsort)

apcombohm1 <- ap4 + ap6 + ap7
apcombohm2 <- ap10 + ap13 + ap15

#feature plots
#OR features
apfeats <- PlotFeatures(m46710_voxsort, featureOut = "quadbox")
#RTP control
aprtp <- PlotFeatures(m46710_voxsort, featureOut = "rtp")
#random control
aprandom <- PlotRandom(m46710_voxsort)
#Taar positions
aptaar <- m46710_voxsort %>% 
  mutate(type = ifelse(str_detect(gene, "Olfr"), "Olfr", "Taar")) %>%
  ggplot(aes(type, wavgs, fill = type)) +
  geom_violin() +
  theme_cowplot() +
  scale_fill_manual(values = c("red", "black")) +
    theme(axis.line = element_blank(),
          axis.ticks.x = element_blank(),
          axis.ticks.y = element_blank(),
          axis.text.y = element_blank(),
          axis.title.y = element_blank(),
          legend.position = "none")

apmerge + ggtitle("AP Merge")
```

![](ob_heatmaps_files/figure-gfm/unnamed-chunk-3-1.png)<!-- -->

``` r
apcombohm1
```

![](ob_heatmaps_files/figure-gfm/unnamed-chunk-3-2.png)<!-- -->

``` r
apcombohm2
```

![](ob_heatmaps_files/figure-gfm/unnamed-chunk-3-3.png)<!-- -->

``` r
apother <- aptaar + plot_spacer() + aprtp + aprandom
apfeats
```

![](ob_heatmaps_files/figure-gfm/unnamed-chunk-3-4.png)<!-- -->

``` r
apother
```

![](ob_heatmaps_files/figure-gfm/unnamed-chunk-3-5.png)<!-- -->

``` r
#test 200 random selections
random_iter <- 200
ap_random_stats <- matrix(nrow = random_iter, ncol = 2)
for (i in 1:random_iter) {
  i_out <- as.numeric(PlotRandom(m46710_voxsort, out = "stat"))
  ap_random_stats[i,1] <- i_out[1]
  ap_random_stats[i,2] <- i_out[2]
}

#mean pval is:
mean(ap_random_stats[,1])
```

    ## [1] 0.5634342

``` r
#mean odds ratio is:
mean(ap_random_stats[,2])
```

    ## [1] 1.219405

# DORSAL-VENTRAL

``` r
## VD normalize values
m9_tpm <- alldim_vd_tpm_cut %>% filter(rep == 9, slice <= 22) %>% select(-name, -rep, -dim, -dimrep)
m12_tpm <- alldim_vd_tpm_cut %>% filter(rep == 12, slice <= 22) %>% select(-name, -rep, -dim, -dimrep)
m14_tpm <- alldim_vd_tpm_cut %>% filter(rep == 14, slice <= 22) %>% select(-name, -rep, -dim, -dimrep)

m9_norm <- NormTPM(m9_tpm)
m12_norm <- NormTPM(m12_tpm)
m14_norm <- NormTPM(m14_tpm)

#hierarchical clustering of both axes
#heatmaply(m12_norm)
#heatmaply_cor(m12_norm)
#plots dont mean much without sorting, good to look at hclust of section axis however. ideally y should move sequentially, may be worth investigating sections that appear out of order for OR tpm total, # missing ORs, # high ORs, etc.

## VD apply voxel factors
# VD make sort for voxel adjusted TPM
vd_voxweights <- allmice_covars %>% filter(rep == 9) %>% select(slice, weight)

m9_Vnorm <- Voxnormify(m9_tpm, vd_voxweights)
m12_Vnorm <- Voxnormify(m12_tpm, vd_voxweights)
m14_Vnorm <- Voxnormify(m14_tpm, vd_voxweights)

#heatmaply_cor(m9_Vnorm, Rowv=T, Colv=F)

#sort and flip to make vd into dv
m91214_voxsort <- MakeSort3(m9_Vnorm, m12_Vnorm, m14_Vnorm, 
                            method = "posmean") %>%
  arrange(desc(sortrank)) %>%
  rename(vdrank = sortrank) %>%
  mutate(sortrank = 1:n()) %>%
  select(-vdrank)

m9_voxsort <- SortByList(m9_Vnorm, m91214_voxsort)
m12_voxsort <- SortByList(m12_Vnorm, m91214_voxsort)
m14_voxsort <- SortByList(m14_Vnorm, m91214_voxsort)

#heatmaply_cor(m4_voxsort, Rowv = F, Colv = F)
#heatmaply_cor(m6_voxsort, Rowv = F, Colv = F)
#heatmaply_cor(m7_voxsort, Rowv = F, Colv = F)
#heatmaply_cor(m10_voxsort, Rowv = F, Colv = F)

dv9 <- MakeBlackRedHeatmap(m9_voxsort)
dv12 <- MakeBlackRedHeatmap(m12_voxsort)
dv14 <- MakeBlackRedHeatmap(m14_voxsort)

#examine mergeDF used to make sort
m91214mergeV <- normNormDF3(m9_Vnorm, m12_Vnorm, m14_Vnorm)
m91214mergeV_voxnormsort <- SortByList(m91214mergeV, m91214_voxsort)
dvmerge <- MakeBlackRedHeatmap(m91214mergeV_voxnormsort)

dvcombohm <- dv9 + dv12 + dv14

dvfeats <- PlotFeatures(m91214_voxsort, featureOut = "quadbox")
dvrtp <- PlotFeatures(m91214_voxsort, featureOut = "rtp")
dvrandom <- PlotRandom(m91214_voxsort)

#Taar positions
dvtaar <- m91214_voxsort %>% 
  mutate(type = ifelse(str_detect(gene, "Olfr"), "Olfr", "Taar")) %>%
  ggplot(aes(type, wavgs, fill = type)) +
  geom_violin() +
  theme_cowplot() +
  scale_fill_manual(values = c("red", "black")) +
    theme(legend.position = "none",
          axis.line = element_blank(),
          axis.ticks.x = element_blank(),
          axis.ticks.y = element_blank(),
          axis.text.y = element_blank(),
          axis.title.y = element_blank())

dvmerge + ggtitle("DV Merge")
```

![](ob_heatmaps_files/figure-gfm/unnamed-chunk-4-1.png)<!-- -->

``` r
dvcombohm
```

![](ob_heatmaps_files/figure-gfm/unnamed-chunk-4-2.png)<!-- -->

``` r
dvother <- dvtaar + plot_spacer() + dvrtp + dvrandom
dvfeats
```

![](ob_heatmaps_files/figure-gfm/unnamed-chunk-4-3.png)<!-- -->

``` r
dvother
```

![](ob_heatmaps_files/figure-gfm/unnamed-chunk-4-4.png)<!-- -->

# MEDIAL-LATERAL

## ML tpm cutoff

## avoiding m16 because it is significantly different from the other replicates

``` r
## ML normalize values
m8_tpm <- alldim_ml_tpm_cut %>% filter(rep == 8, slice <= 19) %>% select(-name, -rep, -dim, -dimrep) #olfr361 is NA after normalization
m11_tpm <- alldim_ml_tpm_cut %>% filter(rep == 11, slice <= 19) %>% select(-name, -rep, -dim, -dimrep)
# m16_tpm <- alldim_ml_tpm_cut %>% filter(rep == 16, slice <= 19) %>% select(-name, -rep, -dim, -dimrep)

m8_norm <- NormTPM(m8_tpm)
m11_norm <- NormTPM(m11_tpm)
# m16_norm <- NormTPM(m16_tpm)

# m16_reorder <- c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11",
#                  "12", "13", "14", "15", "16", "17", "18", "19")
# m16_norm <- m16_norm[,m16_reorder]

#hierarchical clustering of both axes
#heatmaply(m16_norm)
#heatmaply_cor(m8_norm)
#plots dont mean much without sorting, good to look at hclust of section axis however. ideally y should move sequentially, may be worth investigating sections that appear out of order for OR tpm total, # missing ORs, # high ORs, etc.


## ML apply voxel factors
# ML make sort for voxel adjusted TPM
ml_voxweights <- allmice_covars %>% filter(rep == 8) %>% select(slice, weight)

m8_Vnorm <- Voxnormify(m8_tpm, ml_voxweights)
m11_Vnorm <- Voxnormify(m11_tpm, ml_voxweights)
# m16_Vnorm <- Voxnormify(m16_tpm, ml_voxweights) 

#heatmaply_cor(m9_Vnorm, Rowv=T, Colv=F)

m811_voxsort <- MakeSort4(m8_Vnorm, m8_Vnorm, m11_Vnorm, m11_Vnorm, method = "posmean")

m8_voxsort <- SortByList(m8_Vnorm, m811_voxsort)
m11_voxsort <- SortByList(m11_Vnorm, m811_voxsort)
#m16_voxsort <- SortByList(m16_Vnorm, m811_voxsort)

#heatmaply_cor(m8_voxsort, Rowv = F, Colv = F)

ml8 <- MakeBlackRedHeatmap(m8_voxsort)
ml11 <- MakeBlackRedHeatmap(m11_voxsort)
#ml16 <- MakeBlackRedHeatmap(m16_voxsort)

#examine mergeDF used to make sort
m811mergeV <- normNormDF4(m8_Vnorm, m8_Vnorm, m11_Vnorm, m11_Vnorm)
m811mergeV_voxnormsort <- SortByList(m811mergeV, m811_voxsort)
mlmerge <- MakeBlackRedHeatmap(m811mergeV_voxnormsort)

mlcombohm <- ml8 + ml11

mlfeats <- PlotFeatures(m811_voxsort, featureOut = "quadbox")
mlrtp <- PlotFeatures(m811_voxsort, featureOut = "rtp")
mlrandom <- PlotRandom(m811_voxsort)

#Taar positions
mltaar <- m811_voxsort %>% 
  mutate(type = ifelse(str_detect(gene, "Olfr"), "Olfr", "Taar")) %>%
  ggplot(aes(type, wavgs, fill = type)) +
  geom_violin() +
  theme_cowplot() +
  scale_fill_manual(values = c("red", "black")) +
    theme(legend.position = "none",
          axis.line = element_blank(),
          axis.ticks.x = element_blank(),
          axis.ticks.y = element_blank(),
          axis.text.y = element_blank(),
          axis.title.y = element_blank())

mlmerge + ggtitle("ML Merge")
```

![](ob_heatmaps_files/figure-gfm/unnamed-chunk-5-1.png)<!-- -->

``` r
mlcombohm
```

![](ob_heatmaps_files/figure-gfm/unnamed-chunk-5-2.png)<!-- -->

``` r
mlother <- mltaar + plot_spacer() + mlrtp + mlrandom
mlfeats
```

![](ob_heatmaps_files/figure-gfm/unnamed-chunk-5-3.png)<!-- -->

``` r
mlother
```

![](ob_heatmaps_files/figure-gfm/unnamed-chunk-5-4.png)<!-- -->

# Save plots

``` r
SaveHeatmaps(ap4, ap6, ap7, ap10, ap13, ap15,
             dv9, dv12, dv14, ml8, ml11,
             apmerge, dvmerge, mlmerge,
             path = "~/Desktop/obmap/r_analysis/heatmaps/output/below1uni/")

SaveSquarePlots(apfeats, apother, dvfeats, dvother, mlfeats, mlother, 
                path = "~/Desktop/obmap/r_analysis/heatmaps/output/below1uni/")
```

# matrix similarity tests

try sorting with kmeans clustering

``` r
voxsort_cor <- MultiSpearMatrix(m4_voxsort, m6_voxsort, m7_voxsort, m10_voxsort)
#heatmaply_cor(voxsort_cor, Rowv=F, Colv=F)

km <- kmeans(scale(m4_voxsort), 23, iter.max = 1000)
km_sort <- tibble("gene" = names(km$cluster), "sortrank" = km$cluster, wavgs = 2)
m4_kmsort <- SortByList(m4_voxsort, km_sort)
#heatmaply_cor(m4_kmsort, Rowv=F, Colv=F)
#need to add km_sort as well as voxsort and compare, perhaps possible to determine order of km clusters from AtoP
```

# Which ORs covary?

PCA, tsne, dimred plots For each OR, create a vector holding
meanposition for all replicates and all dimensions Calculate all
pairwise correlations for these vectors and determine a cutoff for
covariance, can also compare lowest distance by subtracting or something

``` r
alldim_m4 <- alldim_ap_tpm_cut %>% filter(rep == 4) %>% select(-name, -rep, -dim, -dimrep)
alldim_m6 <- alldim_ap_tpm_cut %>% filter(rep == 6) %>% select(-name, -rep, -dim, -dimrep)
alldim_m7 <- alldim_ap_tpm_cut %>% filter(rep == 7) %>% select(-name, -rep, -dim, -dimrep)
alldim_m10 <- alldim_ap_tpm_cut %>% filter(rep == 10) %>% select(-name, -rep, -dim, -dimrep)

alldim_m8 <- alldim_ml_tpm_cut %>% filter(rep == 8) %>% select(-name, -rep, -dim, -dimrep)
alldim_m11 <- alldim_ml_tpm_cut %>% filter(rep == 11) %>% select(-name, -rep, -dim, -dimrep)
# alldim_m16 <- alldim_ml_tpm_cut %>% filter(rep == 16) %>% select(-name, -rep, -dim, -dimrep)

alldim_m9 <- alldim_vd_tpm_cut %>% filter(rep == 9) %>% select(-name, -rep, -dim, -dimrep)
alldim_m12 <- alldim_vd_tpm_cut %>% filter(rep == 12) %>% select(-name, -rep, -dim, -dimrep)
alldim_m14 <- alldim_vd_tpm_cut %>% filter(rep == 14) %>% select(-name, -rep, -dim, -dimrep)

m4_vnorm_ad <- Voxnormify(alldim_m4, ap_voxweights)
m6_vnorm_ad <- Voxnormify(alldim_m6, ap_voxweights)
m7_vnorm_ad <- Voxnormify(alldim_m7, ap_voxweights)
m10_vnorm_ad <- Voxnormify(alldim_m10, ap_voxweights)

m8_vnorm_ad <- Voxnormify(alldim_m8, ml_voxweights)
m11_vnorm_ad <- Voxnormify(alldim_m11, ml_voxweights)
# m16_vnorm_ad <- Voxnormify(alldim_m16, ml_voxweights)

m9_vnorm_ad <- Voxnormify(alldim_m9, vd_voxweights)
m12_vnorm_ad <- Voxnormify(alldim_m12, vd_voxweights)
m14_vnorm_ad <- Voxnormify(alldim_m14, vd_voxweights)

#for each mouse, for each row
dfs <- list(m4_vnorm_ad, m6_vnorm_ad, m7_vnorm_ad, m8_vnorm_ad, 
            m9_vnorm_ad, m10_vnorm_ad, m11_vnorm_ad, m12_vnorm_ad, 
            m14_vnorm_ad) #, m16_vnorm_ad

#problem is dfs are different lengths, may need to use a different 
meanpos_mtx <- matrix(ncol = length(dfs), nrow = nrow(dfs[[1]]))
for (i in 1:length(dfs)) {
  df_norm <- dfs[[i]]
  secs <- c(1:23)
  wavgs <- vector(mode = "numeric", length = nrow(df_norm))
  for (j in 1:nrow(df_norm)) {
    vals <- df_norm[j,]
    sumvals <- sum(vals)
    prods <- vector(mode = "numeric", length = length(vals))
    for (k in 1:length(vals)) {
      prods[k] <- vals[k] * secs[k]
    }
    wavgs[j] <- sum(prods)/sumvals
  }
  meanpos_mtx[,i] <- wavgs
}

rownames(meanpos_mtx) <- rownames(m4_vnorm_ad)
colnames(meanpos_mtx) <- c("m4_AP", "m6_AP", "m7_AP", "m8_ML", 
                           "m9_VD", "m10_AP", "m11_ML", "m12_VD",
                           "m14_VD") #, "m16_ML"

meanpos_tb <- as_tibble(meanpos_mtx) %>% mutate(gene_id = rownames(meanpos_mtx))

#write_csv(meanpos_tb, "~/Desktop/obmap/r_analysis/heatmaps/output/meanpositions.csv")

heatmaply_cor(meanpos_mtx)
```

    ## Warning in heatmaply.heatmapr(hm, colors = colors, limits = limits,
    ## scale_fill_gradient_fun = scale_fill_gradient_fun, : Upper limit is not >=
    ## highest value in x, max of limits is set to the max of the range (otherwise,
    ## colors will be broken!)

![](ob_heatmaps_files/figure-gfm/unnamed-chunk-8-1.png)<!-- -->

``` r
meanpos_cor <- cor(t(meanpos_mtx))
meanpos_cov <- cov(t(meanpos_mtx))

lancet <- read_csv("~/Downloads/lancet_mouseORnames.csv")
```

    ## 
    ## ── Column specification ────────────────────────────────────────────────────────
    ## cols(
    ##   Symbol = col_character(),
    ##   Pseudo = col_double(),
    ##   Olfrname = col_character(),
    ##   alias = col_character()
    ## )

``` r
corinfo <- tibble(Olfrname = rownames(meanpos_mtx)) %>%
  left_join(lancet, by = "Olfrname") %>% 
  rename(olfrname = Olfrname) %>%
  left_join(info, by = "olfrname") %>%
  select(olfrname:class, fisurface, tan_zone) %>%
  mutate(symfam = as.numeric(str_extract(Symbol, "[0-9]+")),
         orsort = 1:n()) %>%
  arrange(class) %>%
  mutate(csort = 1:n()) %>%
  arrange(symfam) %>%
  mutate(sfsort = 1:n()) %>%
  arrange(orsort)

#make heatmaply_cor with class/family dendros and labels
meanpos_cor_sc <- cbind(meanpos_cor, corinfo$class)
```

    ## Warning in cbind(meanpos_cor, corinfo$class): number of rows of result is not a
    ## multiple of vector length (arg 2)

``` r
# heatmaply_cor(meanpos_cor_sc[,-863], 
#               Rowv = corinfo$class, Colv = "Rowv",
#               dendrogram = "none",
#               row_side_colors = meanpos_cor_sc[,863])

#if still sucks, plot ORs on umap (avoids another fucking heatmap in my paper) and color by class and fam
```

# 10 group covariance

``` r
heatmaply_cor(meanpos_cor, k_row = 10, dist_method = "euclidean", hclust_method = "complete")

dend <- hclust(dist(meanpos_cor, method = "euclidean"), method = "complete")
cutree(dend, k = 10)

tree1 <- names(which(cutree(dend, k = 10) == 1))
tree2 <- names(which(cutree(dend, k = 10) == 2))
tree3 <- names(which(cutree(dend, k = 10) == 3))
tree4 <- names(which(cutree(dend, k = 10) == 4))
tree5 <- names(which(cutree(dend, k = 10) == 5))
tree6 <- names(which(cutree(dend, k = 10) == 6))
tree7 <- names(which(cutree(dend, k = 10) == 7))
tree8 <- names(which(cutree(dend, k = 10) == 8))
tree9 <- names(which(cutree(dend, k = 10) == 9))
tree10 <- names(which(cutree(dend, k = 10) == 10))

#AP
ap1 <- MakeBLRHtree4(m4_Vnorm, m6_Vnorm, m7_Vnorm, m10_Vnorm, tree1, "AP-tree1")
ap2 <- MakeBLRHtree4(m4_Vnorm, m6_Vnorm, m7_Vnorm, m10_Vnorm, tree2, "AP-tree2")
ap3 <- MakeBLRHtree4(m4_Vnorm, m6_Vnorm, m7_Vnorm, m10_Vnorm, tree3, "AP-tree3")
ap4 <- MakeBLRHtree4(m4_Vnorm, m6_Vnorm, m7_Vnorm, m10_Vnorm, tree4, "AP-tree4")
ap5 <- MakeBLRHtree4(m4_Vnorm, m6_Vnorm, m7_Vnorm, m10_Vnorm, tree5, "AP-tree5")
ap6 <- MakeBLRHtree4(m4_Vnorm, m6_Vnorm, m7_Vnorm, m10_Vnorm, tree6, "AP-tree6")
ap7 <- MakeBLRHtree4(m4_Vnorm, m6_Vnorm, m7_Vnorm, m10_Vnorm, tree7, "AP-tree7")
ap8 <- MakeBLRHtree4(m4_Vnorm, m6_Vnorm, m7_Vnorm, m10_Vnorm, tree8, "AP-tree8")
ap9 <- MakeBLRHtree4(m4_Vnorm, m6_Vnorm, m7_Vnorm, m10_Vnorm, tree9, "AP-tree9")
ap10 <- MakeBLRHtree4(m4_Vnorm, m6_Vnorm, m7_Vnorm, m10_Vnorm, tree10, "AP-tree10")

#ML
ml1 <- MakeBLRHtree3(m8_Vnorm, m11_Vnorm, m16_Vnorm, tree1, "ML-tree1")
ml2 <- MakeBLRHtree3(m8_Vnorm, m11_Vnorm, m16_Vnorm, tree2, "ML-tree2")
ml3 <- MakeBLRHtree3(m8_Vnorm, m11_Vnorm, m16_Vnorm, tree3, "ML-tree3")
ml4 <- MakeBLRHtree3(m8_Vnorm, m11_Vnorm, m16_Vnorm, tree4, "ML-tree4")
ml5 <- MakeBLRHtree3(m8_Vnorm, m11_Vnorm, m16_Vnorm, tree5, "ML-tree5")
ml6 <- MakeBLRHtree3(m8_Vnorm, m11_Vnorm, m16_Vnorm, tree6, "ML-tree6")
ml7 <- MakeBLRHtree3(m8_Vnorm, m11_Vnorm, m16_Vnorm, tree7, "ML-tree7")
ml8 <- MakeBLRHtree3(m8_Vnorm, m11_Vnorm, m16_Vnorm, tree8, "ML-tree8")
ml9 <- MakeBLRHtree3(m8_Vnorm, m11_Vnorm, m16_Vnorm, tree9, "ML-tree9")
ml10 <- MakeBLRHtree3(m8_Vnorm, m11_Vnorm, m16_Vnorm, tree10, "ML-tree10")

#VD
vd1 <- MakeBLRHtree3(m9_Vnorm, m12_Vnorm, m14_Vnorm, tree1, "VD-tree1")
vd2 <- MakeBLRHtree3(m9_Vnorm, m12_Vnorm, m14_Vnorm, tree2, "VD-tree2")
vd3 <- MakeBLRHtree3(m9_Vnorm, m12_Vnorm, m14_Vnorm, tree3, "VD-tree3")
vd4 <- MakeBLRHtree3(m9_Vnorm, m12_Vnorm, m14_Vnorm, tree4, "VD-tree4")
vd5 <- MakeBLRHtree3(m9_Vnorm, m12_Vnorm, m14_Vnorm, tree5, "VD-tree5")
vd6 <- MakeBLRHtree3(m9_Vnorm, m12_Vnorm, m14_Vnorm, tree6, "VD-tree6")
vd7 <- MakeBLRHtree3(m9_Vnorm, m12_Vnorm, m14_Vnorm, tree7, "VD-tree7")
vd8 <- MakeBLRHtree3(m9_Vnorm, m12_Vnorm, m14_Vnorm, tree8, "VD-tree8")
vd9 <- MakeBLRHtree3(m9_Vnorm, m12_Vnorm, m14_Vnorm, tree9, "VD-tree9")
vd10 <- MakeBLRHtree3(m9_Vnorm, m12_Vnorm, m14_Vnorm, tree10, "VD-tree10")

ap1 + ml1 + vd1
ap2 + ml2 + vd2
ap3 + ml3 + vd3
ap4 + ml4 + vd4
ap5 + ml5 + vd5
ap6 + ml6 + vd6
ap7 + ml7 + vd7
ap8 + ml8 + vd8
ap9 + ml9 + vd9
ap10 + ml10 + vd10

length(tree1)
length(tree2)
length(tree3)
length(tree4)
length(tree5)
length(tree6)
length(tree7)
length(tree8)
length(tree9)
length(tree10)

#removing ectopics from the treesets
#tree7 <- tree7[-which(tree7 == "Olfr287")]
#tree7 <- tree7[-which(tree7 == "Olfr32")]
# 
#ListDorML(tree1, "point")
# Plot_listML(tree2, "point")
# Plot_listML(tree3, "point")
# Plot_listML(tree4, "point")
# Plot_listML(tree5, "point")
# Plot_listML(tree6, "point")
# Plot_listML(tree7, "point")
# Plot_listML(tree8, "point")
# Plot_listML(tree9, "point")
# Plot_listML(tree10, "point")
```

# examine covariance

``` r
meanpos_dist <- matrix(nrow = nrow(meanpos_mtx), ncol = nrow(meanpos_mtx))
for (i in 1:nrow(meanpos_mtx)) {
  for (j in 1:nrow(meanpos_mtx)) {
    meanpos_dist[i,j] <- sum(abs(meanpos_mtx[i,] - meanpos_mtx[j,]))/10
  }
}

#are high correlation pairs close in distance?
sd(meanpos_cor)
bestcor <- which(meanpos_cor > 0.6)
sd(meanpos_dist)
bestdistmean <- which(meanpos_dist < 2.5)

bestdistmean[which(bestdistmean < 783)]

bdm_rank <- matrix(nrow = nrow(meanpos_dist), ncol = ncol(meanpos_dist))
for (i in 1:nrow(meanpos_dist)) {
  bdm_rank[i,] <- min_rank(meanpos_dist[i,])
}
which(bdm_rank == 2)
rownames(bdm_rank) <- ornames
colnames(bdm_rank) <- ornames

#should we group by dim, determine avgpos of dim, and then compare distance/correlation?
ap_meanpos <- meanpos_tb %>% select(gene_id, m4_AP, m6_AP, m7_AP, m10_AP) %>% mutate(mean_AP = (m4_AP + m6_AP + m7_AP + m10_AP)/4)
vd_meanpos <- meanpos_tb %>% select(gene_id, m9_VD, m12_VD, m14_VD) %>% mutate(mean_VD = (m9_VD + m12_VD + m14_VD)/3)
#excluding m16 due to poor correlation with m8 and m11
#cor(ml_meanpos$m8_ML, ml_meanpos$m16_ML) = -0.2465, -.409 for m11vsm16
ml_meanpos <- meanpos_tb %>% select(gene_id, m8_ML, m11_ML) %>% mutate(mean811_ML = (m8_ML + m11_ML)/2)

dimmeanpos <- data.frame(ap_meanpos$mean_AP, vd_meanpos$mean_VD, ml_meanpos$mean811_ML) %>% data.matrix
dimmeanpos_tb <- ap_meanpos %>% select(gene_id, mean_AP) %>% mutate(mean_VD = vd_meanpos$mean_VD, mean_ML = ml_meanpos$mean811_ML)
#plot as 3D scatter
plot_ly(type = "scatter3d", mode = "markers") %>% 
  add_trace(data=dimmeanpos_tb, x= ~mean_AP, y= ~mean_ML, z= ~mean_VD,
            text = ~paste('Gene: ', gene_id,
                         '<br>AP:', mean_AP, 
                         '<br>ML:', mean_ML,
                         '<br>VD:', mean_VD),
            marker = list(size = 6, line = list(color = 'black', width = 0.5))) %>%
  layout(scene = list(xaxis = list(title = 'Anterior-Posterior'),
                      yaxis = list(title = 'Medial-Lateral'),
                      zaxis = list(title = 'Ventral-Dorsal')))

dim_dist <- matrix(nrow = nrow(dimmeanpos), ncol = nrow(dimmeanpos))
for (i in 1:nrow(dimmeanpos)) {
  for (j in 1:nrow(dimmeanpos)) {
    dim_dist[i,j] <- sum(abs(dimmeanpos[i,] - dimmeanpos[j,]))/3
  }
}

#using the distance between dimensional mean positions, define clusters
#for each OR, rank 1 = least distance
dd_rank <- matrix(nrow = nrow(dim_dist), ncol = nrow(dim_dist))
for (i in 1:nrow(dim_dist)) {
  dd_rank[i,] <- min_rank(dim_dist[i,])
}

#pick 2 to 6 since the diagonal is 0 which becomes rank 1
top5 <- which(between(x = dd_rank, left = 2, right = 6)) 
dd_bin <- matrix(nrow = nrow(dd_rank), ncol = nrow(dd_rank), data = 0)
dd_bin[top5] <- 1

#cluster matrix
cluster_list <- rep(NA, nrow(dd_bin))
for (k in 1:ncol(dd_bin)) {
  #skip points that are already in a cluster
  if (is.na(cluster_list[k])) {
    round <- 1
    #do a bunch of rounds, find a way to have it run until it stops finding new points
    while (round < 3) {
      neigh <- which(dd_bin[,k] == 1)
      #for each neighbor of previous round of neighbors, find new neighbors
      for (l in 1:length(neigh)) {
        neigh <- c(neigh, which(dd_bin[,neigh[l]] == 1))
        neigh <- neigh[-which(duplicated(neigh))]
        neigh <- sort(neigh)
      } #endforl
      round <- round + 1
    } #endwhile
    cluster_list[neigh] <- k
  } else {
    next
  } #endif
} #endfork
```

## Compare Ant and Post Peak Positions across mice

``` r
appp_vnorm <- MultiAPPP(m4_Vnorm, m6_Vnorm, m7_Vnorm, m10_Vnorm, "_m4", "_m6", "_m7", "_m10")

forcomboant <- appp_vnorm %>% select(antpeak_m4, antpeak_m6) %>%
  rename(combo4 = antpeak_m4, combo6 = antpeak_m6)
forcombopost <- appp_vnorm %>% select(postpeak_m4, postpeak_m6) %>%
  rename(combo4 = postpeak_m4, combo6 = postpeak_m6)
forcombo <- bind_rows(forcomboant, forcombopost)

#plot ant vs ant, post vs post
ggplot(appp_vnorm) +
  geom_jitter(aes(antpeak_m4, antpeak_m6), color = "red", alpha = 0.4) +
  geom_jitter(aes(postpeak_m4, postpeak_m6), color = "blue", alpha = 0.4) +
  geom_abline(slope = 1) + 
  geom_smooth(data = forcombo, aes(combo4, combo6), color = "purple", method = "lm", se=F)
```

    ## `geom_smooth()` using formula 'y ~ x'

![](ob_heatmaps_files/figure-gfm/unnamed-chunk-11-1.png)<!-- -->

``` r
# #above as boxes or violins
# ggplot(appp_vnorm %>% mutate(apm4 = antpeak_m4 + 100, ppm4 = postpeak_m4 + 100)) +
#   geom_boxplot(aes(as.character(apm4), antpeak_m6), color = "red", alpha = 0.3) +
#   geom_boxplot(aes(as.character(ppm4), postpeak_m6), color = "blue", alpha = 0.3)
# 
# #above as faceted for ant and post
# m46facet_a <- appp_vnorm %>% select(gene, antpeak_m4, antpeak_m6) %>% 
#   mutate(ap = "ant") %>% 
#   rename(peak_4 = antpeak_m4, peak_6 = antpeak_m6)
# m46facet_p <- appp_vnorm %>% select(gene, postpeak_m4, postpeak_m6) %>% 
#   mutate(ap = "post") %>% 
#   rename(peak_4 = postpeak_m4, peak_6 = postpeak_m6)
# m46facet <- bind_rows(m46facet_a, m46facet_p)
# 
# ggplot(m46facet) + 
#   geom_boxplot(aes(as.character(peak_4 + 100), peak_6)) +
#   facet_grid(ap ~ .)


# examine distances for ant and post peaks between mice
# calculate distances
appp_vnorm_dist <- appp_vnorm %>% mutate(antdist_46 = abs(antpeak_m4 - antpeak_m6),
                                         antdist_47 = abs(antpeak_m4 - antpeak_m7),
                                         antdist_410 = abs(antpeak_m4 - antpeak_m10),
                                         antdist_67 = abs(antpeak_m6 - antpeak_m7),
                                         antdist_610 = abs(antpeak_m6 - antpeak_m10),
                                         antdist_710 = abs(antpeak_m7 - antpeak_m10),
                                         antdistmean_46710 = (antdist_46 + antdist_47 + 
                                                              antdist_410 + antdist_67 +
                                                              antdist_610 + antdist_710)/6,
                                         postdist_46 = abs(postpeak_m4 - postpeak_m6),
                                         postdist_47 = abs(postpeak_m4 - postpeak_m7),
                                         postdist_410 = abs(postpeak_m4 - postpeak_m10),
                                         postdist_67 = abs(postpeak_m6 - postpeak_m7),
                                         postdist_610 = abs(postpeak_m6 - postpeak_m10),
                                         postdist_710 = abs(postpeak_m7 - postpeak_m10),
                                         postdistmean_46710 = (postdist_46 + postdist_47 + 
                                                              postdist_410 + postdist_67 +
                                                              postdist_610 + postdist_710)/6,
                                         aptotal_46 = antdist_46 + postdist_46,
                                         aptotal_47 = antdist_47 + postdist_47,
                                         aptotal_410 = antdist_410 + postdist_410,
                                         aptotal_67 = antdist_67 + postdist_67,
                                         aptotal_610 = antdist_610 + postdist_610,
                                         aptotal_710 = antdist_710 + postdist_710,
                                         aptotal_mean = (aptotal_46 + aptotal_47 +
                                                         aptotal_410 + aptotal_67 + 
                                                         aptotal_610 + aptotal_710)/6,
                                         apmean_46 = (antdist_46 + postdist_46)/2,
                                         apmean_47 = (antdist_47 + postdist_47)/2,
                                         apmean_410 = (antdist_410 + postdist_410)/2,
                                         apmean_67 = (antdist_67 + postdist_67)/2,
                                         apmean_610 = (antdist_610 + postdist_610)/2,
                                         apmean_710 = (antdist_710 + postdist_710)/2,
                                         apmean_mean = (apmean_46 + apmean_47 +
                                                        apmean_410 + apmean_67 +
                                                        apmean_610 + apmean_710)/6)

#test correlation
cor(appp_vnorm$antpeak_m4, appp_vnorm$antpeak_m6, method = "spearman")
```

    ## [1] 0.5440827

``` r
#how many ORs have identical distance between ant and post peak (could be diff directions thou)
length(which(appp_vnorm_dist$antdist_46 == appp_vnorm_dist$postdist_46))
```

    ## [1] 157

``` r
#plot distance between ant peaks for m4 and m6
#for these two mice, most ORs have anterior peaks that are within 2 sections
ggplot(appp_vnorm_dist) + 
  geom_histogram(aes(antdist_46), fill = "red")
```

    ## `stat_bin()` using `bins = 30`. Pick better value with `binwidth`.

![](ob_heatmaps_files/figure-gfm/unnamed-chunk-11-2.png)<!-- -->

``` r
#plot distance between post peaks for m4 and m6
#for these two mice, most ORs have posterior peaks that are within 2.5 sections
ggplot(appp_vnorm_dist) + 
  geom_histogram(aes(postdist_46), fill = "blue")
```

    ## `stat_bin()` using `bins = 30`. Pick better value with `binwidth`.

![](ob_heatmaps_files/figure-gfm/unnamed-chunk-11-3.png)<!-- -->

``` r
#plot mean distance between all mice for all ant and all post peaks for each OR
#across all mice, most ORs have anterior and posterior peaks within 5 sections 
ggplot(appp_vnorm_dist) + 
  geom_histogram(aes(apmean_mean))
```

    ## `stat_bin()` using `bins = 30`. Pick better value with `binwidth`.

![](ob_heatmaps_files/figure-gfm/unnamed-chunk-11-4.png)<!-- -->

``` r
#plot ant vs post distance for m4 vs m6
#for these two mice, most ORs have anterior and posterior peaks within 2 sections
#46 = 2
#47 = 4
#410 = 4
#67 = 3
#610 = 4
#710 = 4
ggplot(appp_vnorm_dist %>% group_by(antdist_46, postdist_46) %>% mutate(count = n()))  + 
  geom_point(aes(antdist_46, postdist_46, size = count), alpha = 0.25) + 
  geom_smooth(aes(antdist_46, postdist_46))
```

    ## `geom_smooth()` using method = 'loess' and formula 'y ~ x'

![](ob_heatmaps_files/figure-gfm/unnamed-chunk-11-5.png)<!-- -->

``` r
#plot ant vs post distance mean for all mice
ggplot(appp_vnorm_dist %>% group_by(antdistmean_46710, postdistmean_46710) %>% mutate(count = n())) +
  geom_point(aes(antdistmean_46710, postdistmean_46710, size = count), alpha = 0.25) + 
  geom_smooth(aes(antdistmean_46710, postdistmean_46710))
```

    ## `geom_smooth()` using method = 'loess' and formula 'y ~ x'

![](ob_heatmaps_files/figure-gfm/unnamed-chunk-11-6.png)<!-- -->

# Use heatmap information to determine OB symmetry line

Should be slanted across the OB moving more laterally as it moves more
posterior

``` r
#hmm im not really seeing what i expect
#basically says an ML plane ~13 is the symmetry line
#SymLiner(m46710mergeV_voxnormsort, m91214mergeV_voxnormsort, m81116mergeV_voxnormsort, out="plot")

#lets try a single replicate
dr1 <- SymLiner(m10_Vnorm, m9_Vnorm, m8_Vnorm)
dr2 <- SymLiner(m13_Vnorm, m12_Vnorm, m11_Vnorm)
#dr3 <- SymLiner(m15_Vnorm, m14_Vnorm, m16_Vnorm) #ml for m16 is flipped so actually lm

fit1 <- SymLiner(m10_Vnorm, m9_Vnorm, m8_Vnorm, out = "fit")
fit2 <- SymLiner(m13_Vnorm, m12_Vnorm, m11_Vnorm, out = "fit")
#fit3 <- SymLiner(m15_Vnorm, m14_Vnorm, m16_Vnorm, out = "fit")

SymLiner(m10_Vnorm, m9_Vnorm, m8_Vnorm, out = "plot")
```

![](ob_heatmaps_files/figure-gfm/unnamed-chunk-12-1.png)<!-- -->

``` r
SymLiner(m13_Vnorm, m12_Vnorm, m11_Vnorm, out = "plot")
```

![](ob_heatmaps_files/figure-gfm/unnamed-chunk-12-2.png)<!-- -->

``` r
#SymLiner(m15_Vnorm, m14_Vnorm, m16_Vnorm, out = "plot")

ggplot() +
  geom_blank() +
  geom_abline(aes(slope = 0.32655, intercept = 7.61119, color = "rep1"), size = 2) +
  geom_abline(aes(slope = 0.24922, intercept = 8.83890, color = "rep2"), size = 2) +
  xlim(0,23) +
  ylim(0,22) + 
  ggtitle("Glomeruli symmetry line as calculated by single dimension heatmap peaks") +
  xlab("Mean position of top 2 AP peaks") +
  ylab("Mean position of top 2 ML peaks") +
  theme_cowplot() +
  xlab("Anterior-Posterior Peak Intersections") +
  ylab("Medial-Lateral Peak Intersections")
```

![](ob_heatmaps_files/figure-gfm/unnamed-chunk-12-3.png)<!-- -->

``` r
apvals <- c(1:28)
mlvals <- vector(length = length(apvals), mode = "numeric")
for (i in 1:length(apvals)) {
  mlvals[i] <- round(8.225045 + 0.287885 * apvals[i], digits = 2)
}

#output symline positions for voxels
symline <- tibble(apvals, mlvals)
write_csv(symline, "~/Desktop/obmap/r_analysis/heatmaps/output/symline.csv")

#for dr1: m10 and m8, goes from 8.83 to 14.3 while AP moves 1->23 (slope = 0.24922)
#for dr2: m11 and m13, ML symmetry line goes from 7.61 to 14.8 while AP moves 1->23 (slope = 0.32655)
#i think the ml and ap ones are trimmed to fit multiple replicates, may need to use original to see if any differences in symline arise

#test other indiv replicates as well as all ap vs all ml
#export results for v20 model


#perhaps I should not be using mice averaged peaks and instead do AP1 vs ML1 etc.
#or i could be discounting ORs whose glomeruli lie very close the symmetry line in the ML dim since AntPeakPostPeak function has a minimum distance of separation
```

# Get Peaks for all AP and ML dims in order to make mock 3D projection

DV location comes from tanindex

``` r
ap_allreps <- MultiAPPP(m4_Vnorm, m6_Vnorm, m7_Vnorm, m10_Vnorm, 
                        "_ap4", "_ap6", "_ap7", "_ap10")
# ap_first should be more lateral
ap_first <- ap_allreps %>% select(gene, contains("ant")) %>%
  rename(AntPos4 = antpeak_ap4, AntPos6 = antpeak_ap6, 
         AntPos7 = antpeak_ap7, AntPos10 = antpeak_ap10,
         APval4 = antval_ap4, APval6 = antval_ap6,
         APval7 = antval_ap7, APval10 = antval_ap10) %>%
  rowwise() %>%
  mutate(side = "Lateral",
         AntPosAvg = mean(c(AntPos4, AntPos6, AntPos7, AntPos10))) %>%
  ungroup()
# ap_second should be more medial
ap_second <- ap_allreps %>% select(gene, contains("post"))%>%
  rename(AntPos4 = postpeak_ap4, AntPos6 = postpeak_ap6, 
         AntPos7 = postpeak_ap7, AntPos10 = postpeak_ap10,
         APval4 = postval_ap4, APval6 = postval_ap6,
         APval7 = postval_ap7, APval10 = postval_ap10) %>%
  mutate(side = "Medial",
         AntPosAvg = mean(c(AntPos4, AntPos6, AntPos7, AntPos10))) %>%
  ungroup()

ml_allreps <- MultiAPPP(m8_voxsort, m11_voxsort, m13_voxsort, m13_voxsort, 
                        "_ml8", "_ml11", "_ap13", "_ap14") %>% 
  select(-antpeak_ap14, -postpeak_ap14, -antval_ap14, -postval_ap14)
# ml_first is more medial
ml_first <- ml_allreps %>% 
  select(gene, antpeak_ml8, antpeak_ml11, postpeak_ap13, 
         antval_ml8, antval_ml11, postval_ap13) %>%
  rename(MedLat8 = antpeak_ml8, MedLat11 = antpeak_ml11, 
         AntPos13 = postpeak_ap13,
         MLval8 = antval_ml8, MLval11 = antval_ml11,
         APval13 = postval_ap13) %>%
  rowwise() %>%
  mutate(MedLatAvg = mean(c(MedLat8, MedLat11))) %>%
  ungroup()
# ml_second is more lateral
ml_second <- ml_allreps %>% 
  select(gene, postpeak_ml8, postpeak_ml11,  
         antpeak_ap13, postval_ml8, postval_ml11, antval_ap13) %>%
  rename(MedLat8 = postpeak_ml8, MedLat11 = postpeak_ml11,
         AntPos13 = antpeak_ap13, MLval8 = postval_ml8, 
         MLval11 = postval_ml11, APval13 = antval_ap13) %>%
  rowwise() %>%
  mutate(MedLatAvg = mean(c(MedLat8, MedLat11))) %>%
  ungroup()

vd_allreps <- MultiAPPP(m9_voxsort, m12_voxsort, m14_voxsort, m15_voxsort, 
                        "_vd9", "_vd12", "_vd14", "_ap15")
vd_first <- vd_allreps %>% select(gene, contains("ant")) %>%
  rename(VenDor9 = antpeak_vd9, VenDor12 = antpeak_vd12, 
         VenDor14 = antpeak_vd14, AntPos15 = antpeak_ap15,
         VDval9 = antval_vd9, VDval12 = antval_vd12,
         VDval14 = antval_vd14, APval15 = antval_ap15) %>%
  rowwise() %>%
  mutate(VenDorAvg = mean(c(VenDor9, VenDor12, VenDor14))) %>%
  ungroup()
vd_second <- vd_allreps %>% select(gene, contains("post")) %>%
  rename(VenDor9 = postpeak_vd9, VenDor12 = postpeak_vd12, 
         VenDor14 = postpeak_vd14, AntPos15 = postpeak_ap15,
         VDval9 = postval_vd9, VDval12 = postval_vd12,
         VDval14 = postval_vd14, APval15 = postval_ap15) %>%
  rowwise() %>%
  mutate(VenDorAvg = mean(c(VenDor9, VenDor12, VenDor14))) %>%
  ungroup() 

#using vd first for both
lateral_join <- left_join(ap_first, ml_second, by = "gene") %>% left_join(vd_first, by = "gene")
medial_join <- left_join(ap_second, ml_first, by = "gene") %>% left_join(vd_second, by = "gene")

heatmap_peaks <- bind_rows(lateral_join, medial_join) %>% 
  rename(olfrname = gene) %>%
  rowwise() %>%
  mutate(AntPosAvg = mean(c(AntPos4, AntPos6, AntPos7, AntPos10, AntPos13, AntPos15))) %>%
  ungroup() %>%
  left_join(info, by = "olfrname") %>%
  arrange(olfrname) %>%
  group_by(olfrname) %>%
  mutate(VDavg = mean(c(VenDorAvg)))

#write_csv(heatmap_peaks, "~/Desktop/obmap/r_analysis/heatmaps/output/heatmap_peaks.csv")
```

``` r
ap_vox <- m46710_voxsort %>% rename(AntPos = sortrank) %>% rename(ap_mp = wavgs, ap_rank = AntPos)
vd_vox <- m91214_voxsort %>% rename(VenDor = sortrank) %>% rename(vd_mp = wavgs, vd_rank = VenDor)
ml_vox <- m811_voxsort %>% rename(MedLat = sortrank) %>% rename(ml_mp = wavgs, ml_rank = MedLat)
ranklist <- left_join(ap_vox, vd_vox, by = "gene")
ranklist2 <- left_join(ranklist, ml_vox, by = "gene")
#write_csv(ranklist2, "~/Desktop/obmap/r_analysis/heatmaps/output/dim_ranklist_201129.csv")
```

## compare DV mean position vs tan/luis indexes, look at outliers, see if i can create something that relates to 3D DV

``` r
ls_idx <- read_csv("~/Desktop/obmap/r_analysis/3dimOB/input/LS_3Dindexes_real_pred.csv")

ap_meanpos <- m46710_voxsort %>%
  rename("olfrname" = "gene",
         "ap_wavgs" = "wavgs",
         "ap_sortrank" = "sortrank")

dv_compare <- m91214_voxsort %>% 
  rename("olfrname" = "gene") %>% 
  left_join(info, by = "olfrname") %>%
  left_join(ap_meanpos, by = "olfrname") %>%
  filter(!is.na(tz_val)) %>% 
  filter(!is.na(DPT_index)) %>%
  filter(!is.na(oe_region)) %>%
  mutate(tzflip = 6 - as.numeric(tan_zone), 
         wavgbin = round(wavgs), 
         isDor = ifelse(oe_region == "Dorsal", 1, 0),
         isTanDori = ifelse(tz_vd == "Dorsal", 1, 0)) %>% 
  filter(!is.na(tz_val)) %>% 
  group_by(wavgbin) %>%
  mutate(mean_dpt = mean(DPT_index), 
         sd_dpt = sd(DPT_index) * 1.645, 
         dptsdhi = mean_dpt + sd_dpt, 
         dptsdlo = mean_dpt - sd_dpt) %>%
  ungroup() %>%
  rowwise() %>%
  mutate(checkout = ifelse(DPT_index > dptsdhi, "ahi", 
                           ifelse(DPT_index < dptsdlo, "zlo", 
                                  "mid")))

ggplot(dv_compare, aes(wavgs, tzflip, color = fisurface)) + 
  geom_point() +
  geom_smooth() +
  ggtitle("Tan Index")
```

![](ob_heatmaps_files/figure-gfm/unnamed-chunk-15-1.png)<!-- -->

``` r
ggplot(dv_compare, aes(wavgs, DPT_index)) + 
  geom_point() +
  geom_smooth() +
  ggtitle("Luis Index")
```

![](ob_heatmaps_files/figure-gfm/unnamed-chunk-15-2.png)<!-- -->

``` r
ggplot(dv_compare, aes(wavgbin, DPT_index, color = checkout)) +
  geom_point() +
  geom_smooth() +
  xlab("Mean Position - Ventral to Dorsal") +
  ggtitle("Luis Index - 1.645SD outliers")
```

![](ob_heatmaps_files/figure-gfm/unnamed-chunk-15-3.png)<!-- -->

``` r
ggplot(dv_compare, aes(wavgs, tzflip, color = ap_wavgs)) +
  geom_point() +
  scale_color_viridis_c() +
  ggtitle("More dorsal ORs are also more anterior")
```

![](ob_heatmaps_files/figure-gfm/unnamed-chunk-15-4.png)<!-- -->

``` r
plot_ly(type = "scatter", mode = "markers") %>%
  add_trace(data = dv_compare, x = ~wavgbin, y = ~DPT_index, 
            color = ~checkout, size = 2,
            text = ~paste('Gene:', olfrname, 
                          '<br>wavgbin:', wavgbin, 
                          '<br>DPTindex:', DPT_index),
            marker = list(size = 6))
```

    ## Warning: `arrange_()` is deprecated as of dplyr 0.7.0.
    ## Please use `arrange()` instead.
    ## See vignette('programming') for more help
    ## This warning is displayed once every 8 hours.
    ## Call `lifecycle::last_warnings()` to see where this warning was generated.

![](ob_heatmaps_files/figure-gfm/unnamed-chunk-15-5.png)<!-- -->

``` r
#wavg = 3 is quite high in terms of DPTindex
#wavg = 20 is quite low in terms of DPTindex

dv_compare %>% group_by(wavgbin) %>% summarise(count = n(),
                                                 meanDPT = mean(DPT_index), 
                                                 meantz = mean(tzflip),
                                                 propDor = sum(isDor)/count,
                                                 propTan = sum(isTanDori)/count) %>%
  ggplot() +
  geom_line(aes(wavgbin, propTan), color = "red") +
  geom_line(aes(wavgbin, propDor), color = "blue") +
  geom_point(aes(wavgbin, propTan), size = 2, color = "red") +
  geom_point(aes(wavgbin, propDor), size = 2, color = "blue") +
  ggtitle("Proportion of ORs that are Tan/Matsu Dorsal in heatmap section",
          subtitle = "Tan = Red, Matsunami = Blue") +
  theme_cowplot()
```

![](ob_heatmaps_files/figure-gfm/unnamed-chunk-15-6.png)<!-- -->

``` r
dv_compare %>% group_by(wavgbin) %>% summarise(count = n(),
                                                 meanDPT = mean(DPT_index), 
                                                 meantz = mean(tzflip),
                                                 propDor = sum(isDor)/count) %>%
  ggplot() + 
  geom_line(aes(wavgbin, meanDPT), color = "red") + 
  geom_line(aes(wavgbin, meantz*10), color = "blue") +
  geom_line(aes(wavgbin, count), color = "green") +
  geom_line(aes(wavgbin, propDor * 100), color = "purple") +
  geom_point(aes(wavgbin, meanDPT), color = "red", size = 2) + 
  geom_point(aes(wavgbin, meantz*10), color = "blue", size = 2) +
  geom_point(aes(wavgbin, count), color = "green", size = 2) +
  geom_point(aes(wavgbin, propDor * 100), color = "purple", size = 2) +
  ggtitle("Tan (x10, blue), Luis (red), MatsuDE (purple), # ORs (green)")
```

![](ob_heatmaps_files/figure-gfm/unnamed-chunk-15-7.png)<!-- -->

``` r
#idea: define set of best predicted 3D DVreceptors (good relation to tan/luis/VDwavg), use to re-predict "bad" predictions?  or perhaps this just leads to more avg values?  would need to do for both medial and lateral halves, perhaps ignore tan and use DPT since more ORs have DPT index
#predictors would be ORs across samples matrix with training set having AP,ML,VD predictions
#do i need to predict one dim at a time?  also need to avoid non-scaffold spaces
#perhaps some sort of dim reduction to "unwrap" OB scaffold into 2d space, then transform tpm matrix into dimred space?
#would need to check if "good" predictions have multiple expression patterns
```

# Filter for uniformity and TPM

Run on non-filtered ORs

``` r
uni_heatmaps_ap4 <- UniformityFilter(m4_Vnorm, ranklist = m46710_voxsort)

uni_heatmaps_dv9 <- UniformityFilter(m9_Vnorm, ranklist = m91214_voxsort)

uni_heatmaps_ml8 <- UniformityFilter(m8_Vnorm, ranklist = m811_voxsort)

uni_heatmaps_ap4
uni_heatmaps_dv9
uni_heatmaps_ml8


#ap
m4_uni <- UniformityFilter(m4_Vnorm, out = "stat")
m6_uni <- UniformityFilter(m6_Vnorm, out = "stat")
m7_uni <- UniformityFilter(m7_Vnorm, out = "stat")
m10_uni <- UniformityFilter(m10_Vnorm, out = "stat")
ap_uni <- tibble(gene = m4_uni$gene,
                 uni4 = m4_uni$isUniform,
                 uni6 = m6_uni$isUniform,
                 uni7 = m7_uni$isUniform,
                 uni10 = m10_uni$isUniform)

#dv
m9_uni <- UniformityFilter(m9_Vnorm, out = "stat")
m12_uni <- UniformityFilter(m12_Vnorm, out = "stat")
m14_uni <- UniformityFilter(m14_Vnorm, out = "stat")
dv_uni <- tibble(gene = m9_uni$gene,
                 uni9 = m9_uni$isUniform,
                 uni12 = m12_uni$isUniform,
                 uni14 = m14_uni$isUniform)

#ml
m8_uni <- UniformityFilter(m8_Vnorm, out = "stat")
m11_uni <- UniformityFilter(m11_Vnorm, out = "stat")
ml_uni <- tibble(gene = m8_uni$gene,
                 uni8 = m8_uni$isUniform,
                 uni11 = m11_uni$isUniform)

all_uni <- left_join(ap_uni, dv_uni, by = "gene") %>% 
  left_join(ml_uni, by = "gene") %>%
  mutate(ap_score = uni4 + uni6 + uni7 + uni10,
         dv_score = uni9 + uni12 + uni14,
         ml_score = uni8 + uni11,
         all_score = ap_score + dv_score + ml_score)

#counts of uniform ORs
all_uni %>% 
  group_by(all_score) %>% 
  summarise(count = n()) %>% 
  arrange(desc(all_score))

all_uni %>% 
  group_by(ap_score) %>% 
  summarise(count = n()) %>% 
  arrange(desc(ap_score))

all_uni %>% 
  group_by(dv_score) %>% 
  summarise(count = n()) %>% 
  arrange(desc(dv_score))

all_uni %>% 
  group_by(ml_score) %>% 
  summarise(count = n()) %>% 
  arrange(desc(ml_score))

ggplot(all_uni, aes(all_score, ap_score)) +
  geom_point()

#how many ORs are uniform in all dims
ap_uniORs <- all_uni %>% filter(ap_score > 0) %>% pull(gene)
dv_uniORs <- all_uni %>% filter(dv_score > 0) %>% pull(gene)
ml_uniORs <- all_uni %>% filter(ml_score > 0) %>% pull(gene)

alldim_uniORs <- unique(c(ap_uniORs, dv_uniORs, ml_uniORs))
#96 ORs were uniform in a single dim

length(ap_uniORs)
length(dv_uniORs)
length(ml_uniORs)

which(dv_uniORs %in% ap_uniORs)
length(which(dv_uniORs %in% ap_uniORs)) 
#22 of 39 dv uni ORs are in the 71 ap uni ORs
apdv_uniORs <- dv_uniORs[which(dv_uniORs %in% ap_uniORs)]

which(ml_uniORs %in% apdv_uniORs)
length(which(ml_uniORs %in% apdv_uniORs))
#8 of 22 ml uni ORs are in the 22 apdv uni ORs
alldim_overlap_uniORs <- ml_uniORs[which(ml_uniORs %in% apdv_uniORs)]
```

### Legacy

# Load 3D model parameters and data

functions loaded from v21\_gen25.Rmd

``` r
# 3d interactive scatterplots with plotly
#given an OR, pick a number of high probability voxels based on signal to noise ratios
Scat_rank <- function(olfr, topX = 72, out = "plot", title = NA) {
  reranked <- ranked %>% 
    filter(olfrname == olfr) %>% 
    mutate(rankofrank = min_rank(voxRankSNR),
           rankcol = min_rank(desc(ifelse(rankofrank <= topX, rankofrank, NA))),
           isRanked = is.na(rankcol))
  
  besties <- reranked %>% filter(isRanked == 0) %>% arrange(voxRankSNR)
  worsties <- reranked %>% filter(isRanked == 1)
  all <- bind_rows(besties, worsties)
  
  if (out == "data") {
    return(all)
  } else {
    p <- plot_ly(type = "scatter3d", mode = "markers") %>% 
      add_trace(data=worsties, x=~AntPos, y=~MedLat, z=~VenDor, 
                color=~rankcol, opacity=0.15,
                text = ~paste('Gene:', olfrname, 
                              '<br>voxRankSNR:', voxRankSNR, 
                              '<br>voxSNRdim', voxSNRdim),
                marker = list(size = 6)) %>%
      add_trace(data=besties, x=~AntPos, y=~MedLat, z=~VenDor, color=~rankcol,
                text = ~paste('Gene:', olfrname, 
                              '<br>voxRankSNR:', voxRankSNR, 
                              '<br>voxSNRdim', voxSNRdim),
                marker = list(size = 6, 
                              line = list(color = 'black', width = 0.5))) %>%
      layout(title = title,
             scene = list(xaxis = list(title = 'Anterior-Posterior'),
                          yaxis = list(title = 'Medial-Lateral'),
                          zaxis = list(title = 'Ventral-Dorsal')))
    return(p)
  } #endif
} #endfunction

#plot p50 values for an olfrs top X voxels
p50plot <- function(olfr, rankx = 50, title = NA) {
  bestp50 <- ranked %>% 
    filter(olfrname == olfr) %>% 
    mutate(rankp1 = min_rank(p50), rankp = min_rank(desc(rankp1))) %>% 
    filter(rankp <= rankx)
  worstp50 <- ranked %>% 
    filter(olfrname == olfr) %>% mutate(rankp = min_rank(p50)) %>% 
    filter(rankp > rankx) %>% mutate(rankna = NA)
  plot_ly(type = "scatter3d", mode = "markers") %>% 
    add_trace(data=worstp50, x=~AntPos, y=~MedLat, z=~VenDor, 
              color='rgb(10,10,10)', opacity=0.1,
              text = ~paste('Gene:', olfr, 
                            '<br>voxRankSNR:', voxRankSNR, 
                            '<br>voxSNRdim', voxSNRdim),
              marker = list(size = 6)) %>%
    add_trace(data=bestp50, x=~AntPos, y=~MedLat, z=~VenDor, color=~rankp,
              text = ~paste('Gene:', olfr, 
                            '<br>voxRankSNR:', voxRankSNR, 
                            '<br>voxSNRdim', voxSNRdim),
              marker = list(size = 6, line = list(color = 'black', width = 0.5))) %>%
    layout(title = title,
           scene = list(xaxis = list(title = 'Anterior-Posterior'),
                        yaxis = list(title = 'Medial-Lateral'),
                        zaxis = list(title = 'Ventral-Dorsal')))
} #endfunction

#define clusters of best X p50 points using a pairwise matrix and ability to output plots and data
#given a number of top ranking positions for an OR, cluster the points based on spatial position
Cluster <- function (olfr, topX = 72, minClustSize = 5, out = "plot", title = NA) {
  df <- ranked %>% 
    filter(olfrname == olfr) %>% 
    mutate(rankofrank = min_rank(desc(min_rank(p50)))) %>% 
    filter(rankofrank <= topX) %>% 
    arrange(desc(p50))
  
  cut <- ranked %>% 
    filter(olfrname == olfr) %>% 
    mutate(rankofrank = min_rank(desc(min_rank(p50)))) %>% 
    filter(rankofrank > topX) %>% 
    arrange(desc(p50))
  
  #make pairwise of matrixes
  dmatrix <- matrix(data = 0, nrow = nrow(df), ncol = nrow(df))
  for (i in 1:nrow(df)) {
    distances <- vector("numeric", length = nrow(df))
    toprankvox <- df$voxel[i]
    ap <- df$AntPos[i]
    ml <- df$MedLat[i]
    vd <- df$VenDor[i]
    for (j in 1:nrow(df)) {
      if (i == j) {
        distances[j] <- 0
      } else {
      distances[j] <- sqrt((df$AntPos[i] - df$AntPos[j])^2 + 
                        (df$MedLat[i] - df$MedLat[j])^2 +
                        (df$VenDor[i] - df$VenDor[j])^2)
      } #endif
    } #endforj
    neighbors <- which(distances <= sqrt(3))
    dmatrix[neighbors, i] <- 1
  } #endfori
  
  #cluster matrix
  cluster_list <- rep(NA, nrow(df))
  for (k in 1:ncol(dmatrix)) {
    #skip points that are already in a cluster
    if (is.na(cluster_list[k])) {
      round <- 1
      #do a bunch of rounds, find a way to have it run until it stops finding new points
      while (round < topX/3) {
        neigh <- which(dmatrix[,k] == 1)
        #for each neighbor of previous round of neighbors, find new neighbors
        for (l in 1:length(neigh)) {
          neigh <- c(neigh, which(dmatrix[,neigh[l]] == 1))
          neigh <- neigh[-which(duplicated(neigh))]
          neigh <- sort(neigh)
        } #endforl
        round <- round + 1
      } #endwhile
      cluster_list[neigh] <- k
    } else {
      next
    } #endif
  } #endfork
  
  df_out <- df %>% mutate(rawclust = cluster_list) %>% 
    group_by(rawclust) %>% 
    mutate(clustmaxp = max(p50), 
           clustminp = min(p50), 
           clustmeanp = mean(p50)) %>%
    add_tally() %>%
    ungroup() %>%
    mutate(clustmaxprank = dense_rank(desc(clustmaxp)), 
           clustmeanprank = dense_rank(desc(clustmeanp)), 
           clustsizerank = dense_rank(desc(n))) %>% 
    arrange(clustsizerank) %>%
    mutate(isCluster = ifelse(n >= minClustSize, 1, 0), 
           clust_unique = dense_rank(desc(clustmaxp * isCluster))) %>%
    select(p2.5:ORrankpervox, rankofrank, rawclust:clustsizerank, clust_unique, isCluster)
  
  df_clustered <- df_out %>% filter(isCluster == 1)
  too_small <- df_out %>% filter(isCluster == 0.5)

  cut_out <- cut %>% mutate(rawclust = NA,
                            clustmaxp = NA,
                            clustminp = NA,
                            clustmeanp = NA,
                            n = NA,
                            clustmeanprank = NA,
                            clustmaxprank = NA,
                            clustsizerank = NA,
                            isCluster = 0,
                            clust_unique = NA) %>%
                      select(p2.5:ORrankpervox, 
                             rankofrank, 
                             rawclust:clustsizerank, 
                             clust_unique, 
                             isCluster) %>% 
                      bind_rows(too_small)
  
  all_out <- bind_rows(df_clustered, cut_out)
  
  #return various things 'rgb(10,10,10)'
  if (out == "data") {
    return(all_out)
  } else {
    p <-  plot_ly(type = "scatter3d", mode = "markers") %>% 
      add_trace(data=cut_out, x=~AntPos, y=~MedLat, z=~VenDor, 
                color="shell", opacity=0.15,
                text = ~paste('Gene:', olfrname, 
                              '<br>C_size_rank:', clustsizerank, 
                              '<br>C_mean_p50:', clustmeanp,
                              '<br>C_max_p50:', clustmaxp,
                              '<br>Cluster:', clust_unique),
                marker = list(size = 6)) %>%
      add_trace(data=df_out, x=~AntPos, y=~MedLat, z=~VenDor, color=~clust_unique,
                text = ~paste('Gene:', olfrname, 
                              '<br>C_size_rank:', clustsizerank, 
                              '<br>C_mean_p50:', clustmeanp,
                              '<br>C_max_p50:', clustmaxp,
                              '<br>Cluster:', clust_unique),
                marker = list(size = 6, line = list(color = 'black', width = 0.5))) %>%
      layout(title = title,
             scene = list(xaxis = list(title = 'Anterior-Posterior'),
                          yaxis = list(title = 'Medial-Lateral'),
                          zaxis = list(title = 'Ventral-Dorsal')))
    return(p)
  } #endif
} #endfunction

#run cluster given an incrementing number of voxels to cluster from
#needs to also output some sort of summary statistic to define the optimal number of initial voxels that returns the "best" clusters 
BestML <- function(olfr, topMin = 50, topMax = 200, topBy = 25, minSize = 5,
                        clustersPerHalfBulb = 1, out = "plot", title = NA) {
  cphb <- 0
  topStep <- topMin
  while (cphb < clustersPerHalfBulb) {
    df_in <- Cluster(olfr, topX = topStep, minClustSize = minSize, out = "data")
    df_ml <- df_in %>% 
      filter(isCluster == 1) %>% 
      group_by(clust_unique) %>% 
      mutate(meanML = mean(MedLat),
             meanAP = mean(AntPos),
             meanVD = mean(VenDor),
             side = ifelse(meanML >= symline$mlvals[which(symline$apvals == round(meanAP))],
                           "Lateral", "Medial")) %>% 
      ungroup()
    
    df_bestM <- df_ml %>% 
      filter(side == "Medial") %>%
      mutate(MedRank = dense_rank(clust_unique),
             LatRank = NA,
             TopStep = topStep,
             sideRank = MedRank) %>%
      filter(MedRank <= clustersPerHalfBulb)
    
    df_bestL <- df_ml %>%
      filter(side == "Lateral") %>%
      mutate(MedRank = NA,
             LatRank = dense_rank(clust_unique),
             TopStep = topStep,
             sideRank = LatRank) %>%
      filter(LatRank <= clustersPerHalfBulb)
    
    clust_bestM <- unique(df_bestM$clust_unique)
    clust_bestL <- unique(df_bestL$clust_unique)
    
    #checks and counters
    cphb <- min(length(clust_bestM), length(clust_bestL))
    topStep <- topStep + topBy
  } #endwhile
  
  clust_bestML <- c(clust_bestM, clust_bestL)
  which_bestML <- df_in[-which(df_ml$clust_unique %in% clust_bestML),]
  
  df_notbest <- df_in %>% 
    mutate(sideRank = NA) %>% 
    select(-clust_unique) %>% 
    mutate(clust_unique = NA)
  df_best <- bind_rows(df_bestM, df_bestL)
  df_all <- bind_rows(df_best, df_notbest)
  
  #output
  if (out == "data") {
    return(df_all)
  } else if (out == "best") {
    return(df_best)
  } else if (out == "notbest") {
    return(df_notbest)
    } else {
    p <- plot_ly(type = "scatter3d", mode = "markers") %>% 
      add_trace(data=df_notbest, x=~AntPos, y=~MedLat, z=~VenDor, 
                color="shell", opacity=0.15,
                text = ~paste('Gene:', olfrname, 
                              '<br>C_size_rank:', clustsizerank, 
                              '<br>C_mean_p50:', clustmeanp,
                              '<br>C_max_p50:', clustmaxp,
                              '<br>Cluster:', clust_unique),
                marker = list(size = 6)) %>%
      add_trace(data=df_best, x=~AntPos, y=~MedLat, z=~VenDor, color=~clust_unique,
                text = ~paste('Gene:', olfrname, 
                              '<br>C_size_rank:', clustsizerank, 
                              '<br>C_mean_ML:', meanML,
                              '<br>C_max_p50:', clustmaxp,
                              '<br>Cluster:', clust_unique,
                              '<br>SideRank:', sideRank),
                marker = list(size = 6, line = list(color = 'black', width = 0.5))) %>%
      layout(title = title,
             scene = list(xaxis = list(title = 'Anterior-Posterior'),
                          yaxis = list(title = 'Medial-Lateral'),
                          zaxis = list(title = 'Ventral-Dorsal')))
    return(p)
  } #endif
} #endfunction

#given a list of Olfr names, output 1 medial and 1 lateral cluster for each name
#perhaps add an if or arg for really big lists to run a reduced step
ListML <- function(x, out = "plot", title = NA) {
  #use a list to build a df of unknown size instead of bind_row each iteration
  list_out <- vector("list", length = length(x))
  for (i in 1:length(x)) {
    list_out[[i]] <- BestML(x[i], topMin = 50, topMax = 150, topBy = 50, 
                                 minSize = 5, clustersPerHalfBulb = 1, out = "best")
  } #endfor
  
  #always need notbest for the shape shell1
  notbest <- BestML(x[1], topMin = 100, topMax = 150, topBy = 50, minSize = 5, 
                         clustersPerHalfBulb = 1, out = "notbest")
  df_out <- bind_rows(list_out)
  
  #output
  if (out == "data") {
    return(df_out)
  } else if (out == "point") {
    df_point <- df_out %>% filter(p50 == clustmaxp)
    p <- plot_ly(type = "scatter3d", mode = "markers") %>% 
      add_trace(data=notbest, x=~AntPos, y=~MedLat, z=~VenDor, 
                color="shell", opacity=0.15,
                text = ~paste('Gene:', olfrname, 
                              '<br>C_size_rank:', clustsizerank, 
                              '<br>C_mean_p50:', clustmeanp,
                              '<br>C_max_p50:', clustmaxp,
                              '<br>Cluster:', clust_unique),
                marker = list(size = 6)) %>%
      add_trace(data=df_point, x=~AntPos, y=~MedLat, z=~VenDor, color=~olfrname,
                text = ~paste('Gene:', olfrname, 
                              '<br>C_size_rank:', clustsizerank, 
                              '<br>C_mean_ML:', meanML,
                              '<br>C_max_p50:', clustmaxp,
                              '<br>Cluster:', clust_unique,
                              '<br>SideRank:', sideRank),
                marker = list(size = 6, line = list(color = 'black', width = 0.5))) %>%
      layout(title = title,
             scene = list(xaxis = list(title = 'Anterior-Posterior'),
                          yaxis = list(title = 'Medial-Lateral'),
                          zaxis = list(title = 'Ventral-Dorsal')))
    return(p)
  } else {
    p <- plot_ly(type = "scatter3d", mode = "markers") %>% 
      add_trace(data=notbest, x=~AntPos, y=~MedLat, z=~VenDor, 
                color="shell", opacity=0.15,
                text = ~paste('Gene:', olfrname, 
                              '<br>C_size_rank:', clustsizerank, 
                              '<br>C_mean_p50:', clustmeanp,
                              '<br>C_max_p50:', clustmaxp,
                              '<br>Cluster:', clust_unique),
                marker = list(size = 6)) %>%
      add_trace(data=df_out, x=~AntPos, y=~MedLat, z=~VenDor, color=~olfrname,
                text = ~paste('Gene:', olfrname, 
                              '<br>C_size_rank:', clustsizerank, 
                              '<br>C_mean_ML:', meanML,
                              '<br>C_max_p50:', clustmaxp,
                              '<br>Cluster:', clust_unique,
                              '<br>SideRank:', sideRank),
                marker = list(size = 6, line = list(color = 'black', width = 0.5))) %>%
      layout(title = title,
             scene = list(xaxis = list(title = 'Anterior-Posterior'),
                          yaxis = list(title = 'Medial-Lateral'),
                          zaxis = list(title = 'Ventral-Dorsal')))
    return(p)
  } #endif
} #endfunction

#run Cluster and check if either top ranked M or L glom is dorsal and has known dorsal expression or is class 1. If True, find a dorsal glom for both M and L
DorsalML <- function(olfr, topMin = 100, topBy = 100, minSize = 3,
                     clustIn = 5, clustOut = 1, out = "plot") {
  print(olfr)
  clustFound <- 0
  topStep <- topMin
  while (clustFound != clustOut) {
    df_in <- Cluster(olfr, topX = topStep, minClustSize = minSize, out = "data")
    df_ml <- df_in %>% 
      filter(isCluster == 1) %>% 
      group_by(clust_unique) %>% 
      mutate(meanML = mean(MedLat),
             meanAP = mean(AntPos),
             meanVD = mean(VenDor)) %>%
      ungroup() %>%
      rowwise() %>%
      mutate(side = ifelse(meanML >= symline$mlvals[which(symline$apvals == round(meanAP))],
                           "Lateral", "Medial")) %>% 
      ungroup() %>%
      left_join(info, by = "olfrname") %>%
      rowwise() %>% 
      mutate(dorsalRating = sum(ifelse(class == 1, 2, 0), 
                                ifelse(tzsimple < 2, 1, 0),
                                ifelse(oe_region == "Dorsal", 1, 0),
                                na.rm = T)) %>%
      ungroup()
    
    df_bestM <- df_ml %>% 
      filter(side == "Medial") %>%
      mutate(MedRank = dense_rank(clust_unique),
             LatRank = NA,
             TopStep = topStep,
             minSize = minSize,
             clustIn = clustIn,
             sideRank = MedRank) %>%
      filter(MedRank <= clustIn) %>%
      arrange(MedRank)
    
    df_bestL <- df_ml %>%
      filter(side == "Lateral") %>%
      mutate(MedRank = NA,
             LatRank = dense_rank(clust_unique),
             TopStep = topStep,
             minSize = minSize,
             clustIn = clustIn,
             sideRank = LatRank) %>%
      filter(LatRank <= clustIn) %>%
      arrange(LatRank)
    
    max <- max(c(df_bestM$meanVD[1], df_bestL$meanVD[1]))
    min <- min(c(df_bestM$meanVD[1], df_bestL$meanVD[1]))
    
    #would be nice to remake using dplyr::case_when()
    if (is.na(max)) {
      print(paste("NA", as.character(topStep)))
      topStep <- topStep + topBy
      next
    } else {
      if (max >= 13) {
        if (min < 13) {
          if (df_bestM$dorsalRating[1] >= 2) {
            #OR has dorsal info, constrain both med/lat to have dorsal glom
            if (df_bestM$meanVD[1] == max) {
              #Medial was dorsal, find dorsal lateral glom and viceversa
              df_bestM <- df_bestM %>%
                filter(MedRank == 1) %>% mutate(test = "1m")
              df_bestL <- df_bestL %>%
                filter(meanVD >= 13) %>%
                filter(LatRank == min(LatRank)) %>% mutate(test = "1m")
              
            } else if (df_bestL$meanVD[1] == max) {
              df_bestL <- df_bestL %>%
                filter(LatRank == 1) %>% mutate(test = "1l")
              df_bestM <- df_bestM %>% 
                filter(meanVD >= 13) %>%
                filter(MedRank == min(MedRank)) %>% mutate(test = "1l")
            } #endif df_bestM$meanVD[1] == max
          } else if (df_bestM$dorsalRating[1] == 1) {
            #OR has mixed info, return best
            df_bestM <- df_bestM %>% filter(MedRank == 1)  %>% mutate(test = "2")
            df_bestL <- df_bestL %>% filter(LatRank == 1)  %>% mutate(test = "2")
          } else {
            #OR has ventral info, constrain both med/lat to have ventral glom
            if (df_bestM$meanVD[1] == min) {
              #Medial was ventral, find ventral ateral glom and viceversa
              df_bestM <- df_bestM %>%
                filter(MedRank == 1)   %>% mutate(test = "3m")
              df_bestL <- df_bestL %>%
                filter(meanVD < 13) %>%
                filter(LatRank == min(LatRank))  %>% mutate(test = "3m")
              
            } else if (df_bestL$meanVD[1] == min) {
              df_bestL <- df_bestL %>%
                filter(LatRank == 1)  %>% mutate(test = "3l")
              df_bestM <- df_bestM %>% 
                filter(meanVD < 13) %>%
                filter(MedRank == min(MedRank))  %>% mutate(test = "3l")
            } #endif (df_bestM$meanVD[1] == min)
          } #endif (df_bestM$dorsalRating[1] >= 2)
        } #endif (min < 13)
      } else {
        #both gloms are ventral, check if dorsal info
        if (df_bestM$dorsalRating[1] >= 2) {
          #both gloms are ventral, but OR has dorsal info, constrain to dorsal
          df_bestM <- df_bestM %>%
            filter(meanVD >= 13) %>%
            filter(MedRank == min(MedRank))  %>% mutate(test = "4")
          df_bestL <- df_bestL %>%
            filter(meanVD >= 13) %>%
            filter(LatRank == min(LatRank))  %>% mutate(test = "4")
        } else {
          df_bestM <- df_bestM %>% filter(MedRank == 1)   %>% mutate(test = "5")
          df_bestL <- df_bestL %>% filter(LatRank == 1)  %>% mutate(test = "5")
        } #endif (df_bestM$dorsalRating[1] >= 2)
      } #endif (max >= 13)
    } #endif (is.nax(max))
    
    df_bestM <- df_bestM %>% filter(MedRank == min(MedRank))  %>% mutate(test = "suck")
    df_bestL <- df_bestL %>% filter(LatRank == min(LatRank)) %>% mutate(test = "suck")
    
    clust_bestM <- unique(df_bestM$clust_unique)
    clust_bestL <- unique(df_bestL$clust_unique)
    
    #checks and counters
    clustFound <- min(c(length(clust_bestM), length(clust_bestL)))
    print(topStep)
    topStep <- topStep + topBy
    clustIn <- round(clustIn * 1.5)
  } #endwhile
  
  clust_bestML <- c(clust_bestM, clust_bestL)
  which_bestML <- df_in[-which(df_ml$clust_unique %in% clust_bestML),]
  
  df_best <- bind_rows(df_bestM, df_bestL) %>% unique()
  df_notbest <- df_in %>% 
    mutate(sideRank = NA) %>% 
    select(-clust_unique) %>% 
    mutate(clust_unique = NA) %>%
    unique()
  df_all <- bind_rows(df_best, df_notbest) %>% unique()
  
  #output
  if (out == "data") {
    return(df_all)
  } else if (out == "best") {
    return(df_best)
  } else if (out == "notbest") {
    return(df_notbest)
  } else {
    p <- plot_ly(type = "scatter3d", mode = "markers") %>%
      add_trace(data=df_notbest, x=~AntPos, y=~MedLat, z=~VenDor,
                color="shell", opacity=0.15,
                text = ~paste('Gene:', olfrname,
                              '<br>C_size_rank:', clustsizerank,
                              '<br>C_mean_p50:', clustmeanp,
                              '<br>C_max_p50:', clustmaxp,
                              '<br>Cluster:', clust_unique),
                marker = list(size = 6)) %>%
      add_trace(data=df_best, x=~AntPos, y=~MedLat, z=~VenDor, color=~clust_unique,
                text = ~paste('Gene:', olfrname,
                              '<br>C_size_rank:', clustsizerank,
                              '<br>C_mean_ML:', meanML,
                              '<br>C_max_p50:', clustmaxp,
                              '<br>Cluster:', clust_unique,
                              '<br>SideRank:', sideRank),
                marker = list(size = 6, line = list(color = 'black', width = 0.5))) %>%
      layout(title = title,
             scene = list(xaxis = list(title = 'Anterior-Posterior'),
                          yaxis = list(title = 'Medial-Lateral'),
                          zaxis = list(title = 'Ventral-Dorsal')))
    return(p)
  } #endif
} #endfunction

#given a list of Olfr names, output 1 medial and 1 lateral cluster for each name
ListDorML <- function(x, out = "plot", title = NA) {
  #use a list to build a df of unknown size instead of bind_row each iteration
  list_out <- vector("list", length = length(x))
  for (i in 1:length(x)) {
    list_out[[i]] <- DorsalML(x[i], topBy = 200, out = "best")
  } #endfor
  
  #always need notbest for the shape shell1
  notbest <- DorsalML(x[1], out = "notbest")
  df_out <- bind_rows(list_out)
  
  #output
  if (out == "data") {
    return(df_out)
  } else if (out == "point") {
    df_point <- df_out %>% filter(p50 == clustmaxp)
    p <- plot_ly(type = "scatter3d", mode = "markers") %>% 
      add_trace(data=notbest, x=~AntPos, y=~MedLat, z=~VenDor, 
                color="shell", opacity=0.15,
                text = ~paste('Gene:', olfrname, 
                              '<br>C_size_rank:', clustsizerank, 
                              '<br>C_mean_p50:', clustmeanp,
                              '<br>C_max_p50:', clustmaxp,
                              '<br>Cluster:', clust_unique),
                marker = list(size = 6)) %>%
      add_trace(data=df_point, x=~AntPos, y=~MedLat, z=~VenDor, color=~olfrname,
                text = ~paste('Gene:', olfrname, 
                              '<br>C_size_rank:', clustsizerank, 
                              '<br>C_mean_ML:', meanML,
                              '<br>C_max_p50:', clustmaxp,
                              '<br>Cluster:', clust_unique,
                              '<br>SideRank:', sideRank),
                marker = list(size = 6, line = list(color = 'black', width = 0.5))) %>%
      layout(title = title,
             scene = list(xaxis = list(title = 'Anterior-Posterior'),
                          yaxis = list(title = 'Medial-Lateral'),
                          zaxis = list(title = 'Ventral-Dorsal')))
    return(p)
  } else {
    p <- plot_ly(type = "scatter3d", mode = "markers") %>% 
      add_trace(data=notbest, x=~AntPos, y=~MedLat, z=~VenDor, 
                color="shell", opacity=0.15,
                text = ~paste('Gene:', olfrname, 
                              '<br>C_size_rank:', clustsizerank, 
                              '<br>C_mean_p50:', clustmeanp,
                              '<br>C_max_p50:', clustmaxp,
                              '<br>Cluster:', clust_unique),
                marker = list(size = 6)) %>%
      add_trace(data=df_out, x=~AntPos, y=~MedLat, z=~VenDor, color=~olfrname,
                text = ~paste('Gene:', olfrname, 
                              '<br>C_size_rank:', clustsizerank, 
                              '<br>C_mean_ML:', meanML,
                              '<br>C_max_p50:', clustmaxp,
                              '<br>Cluster:', clust_unique,
                              '<br>SideRank:', sideRank),
                marker = list(size = 6, line = list(color = 'black', width = 0.5))) %>%
      layout(title = title,
             scene = list(xaxis = list(title = 'Anterior-Posterior'),
                          yaxis = list(title = 'Medial-Lateral'),
                          zaxis = list(title = 'Ventral-Dorsal')))
    return(p)
  } #endif
} #endfunction


kzY <- allmice_tpm %>% filter(dimrep == 1 | dimrep == 2)
ectopic <- c("Olfr287","Olfr32")
kzYgood <- dplyr::select(kzY, -ectopic) %>% dplyr::select(-name, -rep, -slice, -dim, -dimrep)
kzmY <- as.matrix(kzYgood)
rownames(kzmY) <- paste0("sample", 1:dim(kzY)[1])

#rdirichlet requires vector not dataframe so matrix the tibble
kzX <- read_csv("~/Desktop/obmap/r_analysis/3dimOB/input/allmice_covariates_trim_voxweights_v3.csv") %>% filter(name %in% kzY$name)
N <- dim(kzmY)[1] #number of sections
D <-  dim(kzmY)[2] #number of genes

# Voxel System
d1 <- 23 #antpost
d2 <- 22 #medlat
d3 <- 23 #vendor

#load voxels from bulb_in_cube 
# voxalls_test <- read_csv("~/Desktop/obmap/r_analysis/3dimOB/input/blankOBcoords200922.csv") %>% mutate(type = "latfill") #incorporating fill for 1 voxel thick layer on most lateral surface
# voxalls_0 <- readRDS("~/Desktop/obmap/r_analysis/3dimOB/input/191010_outer_inner_coords_ml23nohole.RDS")
# voxalls_0 <- bind_rows(voxalls_0, voxalls_test)
# voxalls_1 <- voxalls_0[-which(duplicated(voxalls_0)),]
# vox_ap <- voxalls_1[which(voxalls_1$AntPos <= d1),]
# vox_ml <- vox_ap[which(vox_ap$MedLat <= d2),]
# vox_vd <- vox_ml[which(vox_ml$VenDor <= d3),]
# voxalls <- vox_vd %>% select(-type) %>% unique()
voxalls <- read_csv("~/Desktop/obmap/r_analysis/3dimOB/input/voxalls_200924.csv")

tidy_result <- readRDS("~/Desktop/obmap/r_analysis/3dimOB/tidyresults/200i_allchemo_dimrep12_allvox.RDS")

#figure out number of rows in tidy_result corresponding to each Olfr
#repeat olfrname for each voxel, total rows/1100 (#ofORs) returns voxel #
ranked <- tidy_result %>%
  arrange(gene) %>%
  mutate(olfrname = rep(colnames(kzmY), each=max(tidy_result$voxel))) %>%
  group_by(olfrname) %>%
  mutate(voxrankperOR = rank(desc(p50))) %>%
  ungroup() %>%
  group_by(voxel) %>%
  mutate(ORrankpervox = rank(desc(p50))) %>%
  ungroup() %>% group_by(gene) %>% 
  mutate(allVoxAvgp50 = mean(p50)) %>% 
  ungroup() %>% group_by(gene, AntPos) %>% 
  mutate(AP_Avgp50 = mean(p50), AP_SNR = p50/AP_Avgp50) %>%
  ungroup() %>% group_by(gene, MedLat) %>% 
  mutate(ML_Avgp50 = mean(p50), ML_SNR = p50/ML_Avgp50) %>% 
  ungroup() %>% group_by(gene, VenDor) %>% 
  mutate(VD_Avgp50 = mean(p50), VD_SNR = p50/VD_Avgp50) %>% 
  ungroup() %>% 
  mutate(voxSNRdim = (AP_SNR + ML_SNR + VD_SNR)/3) %>% 
  group_by(gene) %>% 
  mutate(geneRankSNR = rank(desc(voxSNRdim))) %>% 
  ungroup() %>% group_by(voxel) %>% 
  mutate(voxRankSNR = rank(desc(voxSNRdim))) %>% 
  ungroup() %>% 
  select(p2.5, p50, p97.5, AntPos:ORrankpervox, geneRankSNR, 
         voxRankSNR, AP_Avgp50:voxSNRdim, voxel)
```

# Split into n groups

``` r
#branch400 <- PlotBranches(400)
#saveRDS(branch400, "~/Desktop/obmap/r_analysis/heatmaps/output/branch400.RDS")
branch400 <- readRDS("~/Desktop/obmap/r_analysis/heatmaps/output/branch400.RDS")

#400 branches, only 127 have >1 OR
#300 brnaches, only 122 have > 1 OR
#200 branches, only 118 have > 1 OR

goodcor400 <- which(branch400$cor > 0.95)
gooddist400 <- which(branch400$avg_dist < 12)
goodset <- goodcor400[which(goodcor400 %in% gooddist400)]

bestcor400 <- which(branch400$cor == max(branch400$cor, na.rm=T))
branch400$heatmaps[bestcor400]

bestdist400 <- which(branch400$avg_dist == min(branch400$avg_dist, na.rm=T))
branch400$`3D`[107]
#find something to define the branch with the highest ranked cor that is also the highest ranked 3d
```
