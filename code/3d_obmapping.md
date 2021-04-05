OBmap 3D
================
kanazian
2021-04-05

# Goal:

Create a 3D reconstruction of single dimension targeted transcriptomics
data in order to calculate the position (not necessarily a specific
voxel but moreso a domain-like area) of an ORs two glomeruli.

## Prior Steps:

Input TPM data is gencode v25 aligned, RSEM quantified and merged using
make\_obmtx.Rmd file. Input voxels are based on MRI model which is
transformed to 3D coordinates using mri\_OB\_cubed.Rmd file. Weights are
calculated based on the number of voxels in that section for each
dimension. Line of glomeruli mirror symmetry was calculated from
two-peak AP and ML positions using ob\_heatmaps.Rmd file.

## Plan:

Use best 2 replicates (1 and 2) to generate the statistical model with
200 iterations (more iter will require running on server or adding more
swap space). Model involves ILR transform of dirichlet sampled TPM data
which is then used to calculate Bayesian posterior probabilities for
each OR in each voxel. p50 quantile value will be used to calculate most
likely positions for each OR. To avoid picking a point that is an
intersection of high expression sections in each dimension, a clustering
function is used to determine high expression areas. A dorsal-ventral
constraint function based on prior knowledge of an ORs position in the
olfactory epithelium is combined with predicted OB positions in order to
ensure 2 matching DV predictions if at least one predicted position
matches with olfactory epithelium position.

# Packages and functions

``` r
#knitr::opts_chunk$set(warning=F)
library(tictoc)
library(MCMCpack)
library(plotly) 
library(tidyverse)
library(driver) # devtools::install_github("jsilve24/driver")
library(Rcpp)
library(ggsci)
library(cowplot)
#detach("package:vegan", unload=TRUE)

# functions -----------------------------------------------------------
# model prep
SampleData <- function(Y, alpha, iter=1){
  d <- array(0, dim=c(N, D, iter))
  for (i in 1:nrow(Y)){
    d[i,,] <- t(rdirichlet(iter, Y[i,] + alpha[i,]))
  } #endfor
  return(d)
} #end SampleData

# calculate quantiles
sourceCpp("~/Desktop/obmap/r_analysis/data/3dimOB/model/kzcolQuant.cpp")

#function for distance between 3d points
Dist3d <- function(df, row1, row2) {
  thing <- sqrt((df$AntPos[row1] - df$AntPos[row2])^2 +
                  (df$MedLat[row1] - df$MedLat[row2])^2 + 
                  (df$VenDor[row1] - df$VenDor[row2])^2)
  return(thing)
} #end Dist3d

# 3d interactive scatterplot
Scat3d <- function(df, ml, ap, vd, color="voxrankperOR") {
  plot_ly(df, 
          x = ml, 
          y = ap, 
          z = vd, 
          color = color,
          text = ~paste('Gene:', olfrname, 
                        '<br>Rank:', voxrankperOR,
                        '<br>p50:', p50),
          marker = list(size = 6,
                        line = list(color = 'black',
                                    width = 0.5)),
          type = 'scatter3d',
          mode = 'markers') %>%
    layout(scene = list(xaxis = list(title = 'Medial-Lateral'),
                        yaxis = list(title = 'Anterior-Posterior'),
                        zaxis = list(title = 'Ventral-Dorsal')))
} #end Scat3d

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
                marker = list(size = 5, color = "grey")) %>%
      add_trace(data=besties, x=~AntPos, y=~MedLat, z=~VenDor, color=~rankcol,
                text = ~paste('Gene:', olfrname, 
                              '<br>voxRankSNR:', voxRankSNR, 
                              '<br>voxSNRdim', voxSNRdim),
                marker = list(size = 6, line = list(color = 'black', width = 0.5))) %>%
      layout(title = title,
             scene = list(xaxis = list(title = 'Anterior-Posterior'),
                          yaxis = list(title = 'Medial-Lateral'),
                          zaxis = list(title = 'Ventral-Dorsal')))
    return(p)
  } #endif
} #end Scat_rank

#plot p50 values for an olfrs top X voxels
p50plot <- function(olfr, rankx = 50, title = NA) {
  bestp50 <- ranked %>% 
    filter(olfrname == olfr) %>% 
    mutate(rankp = min_rank(desc(p50))) %>% filter(rankp <= rankx)
  worstp50 <- ranked %>% 
    filter(olfrname == olfr) %>% mutate(rankp = min_rank(p50)) %>% 
    filter(rankp > rankx) %>% mutate(rankna = NA)
  plot_ly(type = "scatter3d", mode = "markers") %>% 
    add_trace(data=worstp50, x=~AntPos, y=~MedLat, z=~VenDor, 
              color="shell", opacity=0.15,
              text = ~paste('Gene:', olfr, 
                            '<br>voxRankSNR:', voxRankSNR, 
                            '<br>voxSNRdim', voxSNRdim),
              marker = list(size = 5, color = "grey")) %>%
    add_trace(data=bestp50, x=~AntPos, y=~MedLat, z=~VenDor, color=~rankp,
              text = ~paste('Gene:', olfr, 
                            '<br>voxRankSNR:', voxRankSNR, 
                            '<br>voxSNRdim', voxSNRdim),
              marker = list(size = 6, line = list(color = 'black', width = 0.5))) %>%
    layout(title = title,
           scene = list(xaxis = list(title = 'Anterior-Posterior'),
                        yaxis = list(title = 'Medial-Lateral'),
                        zaxis = list(title = 'Ventral-Dorsal')))
} #end p50plot

#define clusters of best X p50 points using a pairwise matrix and ability to output plots and data
#given a number of top ranking positions for an OR, cluster the points based on spatial position
#why the hell am i ranking just to rerank with desc???
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
  
  #make pairwise matrixes
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
  
  all_out <- bind_rows(df_clustered, cut_out) %>% unique()
  
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
                marker = list(size = 5, color = "grey")) %>%
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
} #end Cluster

#run cluster given an incrementing number of voxels to cluster from
#needs to also output some sort of summary statistic to define the optimal number of initial voxels that returns the "best" clusters 
BestML <- function(olfr, topMin = 50, topMax = 200, topBy = 25, minSize = 2,
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
             clustrank = min_rank(desc(p50))) %>% 
      ungroup() %>%
      rowwise() %>%
      mutate(fadenotbest = ifelse(clustrank == 1, 1, 0),
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
                marker = list(size = 5, color = "grey")) %>%
      add_trace(data=df_best %>% filter(fadenotbest == 1), x=~AntPos, y=~MedLat, z=~VenDor, 
                color="Best Probability",
                text = ~paste('Gene:', olfrname, 
                              '<br>C_size_rank:', clustsizerank, 
                              '<br>C_mean_ML:', meanML,
                              '<br>C_max_p50:', clustmaxp,
                              '<br>Cluster:', clust_unique,
                              '<br>SideRank:', sideRank),
                marker = list(size = 6, color = "red", line = list(color = 'black', width = 0.5))) %>%
      add_trace(data=df_best %>% filter(fadenotbest == 0), x=~AntPos, y=~MedLat, z=~VenDor, 
                color="High Probability", opacity = 0.4,
                text = ~paste('Gene:', olfrname, 
                              '<br>C_size_rank:', clustsizerank, 
                              '<br>C_mean_ML:', meanML,
                              '<br>C_max_p50:', clustmaxp,
                              '<br>Cluster:', clust_unique,
                              '<br>SideRank:', sideRank),
                marker = list(size = 6, color = "orange", line = list(color = 'black', width = 0.5))) %>%
      layout(title = title,
             scene = list(xaxis = list(title = 'Anterior-Posterior'),
                          yaxis = list(title = 'Medial-Lateral'),
                          zaxis = list(title = 'Ventral-Dorsal')))
    return(p)
  } #endif
} #end BestML

#given a list of Olfr names, output 1 medial and 1 lateral cluster for each name
#perhaps add an if or arg for really big lists to run a reduced step
ListML <- function(x, out = "plot", title = NA) {
  #use a list to build a df of unknown size instead of bind_row each iteration
  list_out <- vector("list", length = length(x))
  for (i in 1:length(x)) {
    list_out[[i]] <- BestML(x[i], topMin = 50, topMax = 150, topBy = 50, 
                            minSize = 2, clustersPerHalfBulb = 1, out = "best")
  } #endfor
  
  #always need notbest for the shape shell1
  notbest <- BestML(x[1], topMin = 100, topMax = 150, topBy = 50, minSize = 2, 
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
                marker = list(size = 5, color = "grey")) %>%
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
} #end ListML

#run Cluster and check if either top ranked M or L glom is dorsal and has known dorsal expression or is class 1. If True, find a dorsal glom for both M and L
DorsalML <- function(olfr, topMin = 50, topBy = 50, minSize = 3,
                     clustIn = 5, clustOut = 1, out = "plot",
                     title = NA) {
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
             meanVD = mean(VenDor),
             clustrank = min_rank(desc(p50))) %>%
      ungroup() %>%
      rowwise() %>%
      mutate(fadenotbest = ifelse(clustrank == 1, 1, 0),
             side = ifelse(meanML >= symline$mlvals[which(symline$apvals == round(meanAP))],
                           "Lateral", "Medial")) %>% 
      ungroup() %>%
      left_join(info, by = "olfrname") %>%
      rowwise() %>% 
      mutate(dorsalRating = sum(ifelse(class == 1, 2, 0), 
                                ifelse(tz_val < 2, 1, 0),
                                ifelse(oe_region == "Dorsal", 1, 0),
                                ifelse(str_detect(olfrname, "Olfr"), 0, 2),
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
      add_trace(data=blankdata, x=~AntPos, y=~MedLat, z=~VenDor,
                color="shell", opacity=0.15,
                marker = list(size = 5, color = "grey")) %>%
      add_trace(data=df_best %>% filter(fadenotbest == 1), x=~AntPos, y=~MedLat, z=~VenDor, 
                color="Best Probability",
                text = ~paste('Gene:', olfrname,
                              '<br>C_size_rank:', clustsizerank,
                              '<br>C_mean_ML:', meanML,
                              '<br>C_max_p50:', clustmaxp,
                              '<br>Cluster:', clust_unique,
                              '<br>SideRank:', sideRank),
                marker = list(size = 6, color = "red", line = list(color = 'black', width = 0.5))) %>%
      add_trace(data=df_best %>% filter(fadenotbest == 0), x=~AntPos, y=~MedLat, z=~VenDor, 
                color="High Probability",  opacity=0.6,
                text = ~paste('Gene:', olfrname,
                              '<br>C_size_rank:', clustsizerank,
                              '<br>C_mean_ML:', meanML,
                              '<br>C_max_p50:', clustmaxp,
                              '<br>Cluster:', clust_unique,
                              '<br>SideRank:', sideRank),
                marker = list(size = 6, color = "orange", line = list(color = 'black', width = 0.5))) %>%
      layout(title = title,
             scene = list(xaxis = list(title = 'Anterior-Posterior'),
                          yaxis = list(title = 'Medial-Lateral'),
                          zaxis = list(title = 'Ventral-Dorsal')))
    return(p)
  } #endif
} #end DorsalML

#given a list of Olfr names, output 1 medial and 1 lateral cluster for each name
ListDorML <- function(x, by, out = "plot", title = NA) {
  #use a list to build a df of unknown size instead of bind_row each iteration
  list_out <- vector("list", length = length(x))
  for (i in 1:length(x)) {
    list_out[[i]] <- DorsalML(x[i], topBy = by, out = "best")
  } #endfori
  
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
                marker = list(size = 5, color = "grey")) %>%
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
                marker = list(size = 5, color = "grey")) %>%
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
} #end ListDorML
```

# Load info and features

Both sexes represented in each dimensional group. Line of glomeruli
symmetry determined from 1D AP,ML heatmap mean of top 2 peaks

``` r
info <- read_csv("~/Desktop/obmap/r_analysis/data/3dimOB/model/info_210127.csv", 
                 col_names = TRUE) %>% 
  rename("olfrname" = "gene")
```

    ## 
    ## ── Column specification ────────────────────────────────────────────────────────
    ## cols(
    ##   gene = col_character(),
    ##   class = col_double(),
    ##   oe_region = col_character(),
    ##   fisurface = col_logical(),
    ##   tan_zone = col_character(),
    ##   tz_val = col_double(),
    ##   tz_vd = col_character(),
    ##   RTP = col_character(),
    ##   known = col_logical(),
    ##   DPT_index = col_double(),
    ##   Miy_index = col_double(),
    ##   Zolfr_Momb = col_double(),
    ##   ap_mp = col_double(),
    ##   ap_rank = col_double(),
    ##   vd_mp = col_double(),
    ##   vd_rank = col_double(),
    ##   ml_mp = col_double(),
    ##   ml_rank = col_double(),
    ##   protlength = col_double()
    ## )

``` r
symline <- read_csv("~/Desktop/obmap/r_analysis/data/3dimOB/model/symline_3rep.csv")
```

    ## 
    ## ── Column specification ────────────────────────────────────────────────────────
    ## cols(
    ##   apvals = col_double(),
    ##   mlvals = col_double()
    ## )

``` r
blankdata <- readRDS("~/Desktop/obmap/r_analysis/data/mri_to_R/210404_voxV4b_strip25_straightMdoubled_outeronly_rescaled.RDS") %>% unique()

filterORs <- readRDS("~/Desktop/obmap/r_analysis/data/heatmaps/filterORs.RDS")

kzY <- read_csv("~/Desktop/obmap/r_analysis/data/3dimOB/model/covarintactchemo2_over_samples_210227.csv", col_names = TRUE) %>% 
  select(-contains("Vmn"))
```

    ## 
    ## ── Column specification ────────────────────────────────────────────────────────
    ## cols(
    ##   .default = col_double(),
    ##   name = col_character(),
    ##   dim = col_character()
    ## )
    ## ℹ Use `spec()` for the full column specifications.

``` r
#TRY MAKING A MOMBOR
#make a testOR - Anterior, Dorsal, Medial and Middle, Ventral, Middle
testORs <- kzY %>% 
  select(name:dimrep) %>% 
  mutate(Olfr3_3_18_na = case_when(dim == "AntPos" & slice == 3 ~ 500,
                                dim == "AntPos" & slice == 2 ~ 250,
                                dim == "AntPos" & slice == 4 ~ 250,
                                dim == "MedLat" & slice == 3 ~ 500,
                                dim == "MedLat" & slice == 2 ~ 250,
                                dim == "MedLat" & slice == 4 ~ 250,
                                dim == "VenDor" & slice == 18 ~ 500,
                                dim == "VenDor" & slice == 17 ~ 250,
                                dim == "VenDor" & slice == 19 ~ 250,
                                dim == "AntPos" & slice == 2 ~ 500,
                                dim == "AntPos" & slice == 1 ~ 250,
                                dim == "AntPos" & slice == 3 ~ 250,
                                dim == "MedLat" & slice == 9 ~ 500,
                                dim == "MedLat" & slice == 8 ~ 250,
                                dim == "MedLat" & slice == 10 ~ 250,
                                dim == "VenDor" & slice == 17 ~ 500,
                                dim == "VenDor" & slice == 18 ~ 250,
                                dim == "VenDor" & slice == 16 ~ 250),
         Olfr3_3_18 = ifelse(is.na(Olfr3_3_18_na), 0, Olfr3_3_18_na)) %>%
  select(name, Olfr3_3_18)

kzYprep <- kzY %>%
  left_join(testORs, by = "name") %>% #adding testORs
  filter(rep != 14) %>% #remove newCap VD
  filter(rep != 13) %>% #remove newCap AP Dimrep 2
  filter(rep != 15) %>% #remove newCap AP Dimrep 3
  filter(rep != 16) #remove ML Dimrep3
  
kzYgood <- kzYprep %>% select(-name, -rep, -slice, -dim, -dimrep)
kzmY <- as.matrix(kzYgood)
rownames(kzmY) <- paste0("sample", 1:dim(kzYgood)[1])

N <- dim(kzmY)[1] #number of sections
D <-  dim(kzmY)[2] #number of genes

# Voxel System
d1 <- 23 #antpost
d2 <- 22 #medlat
d3 <- 23 #vendor
```

# calculate weights

``` r
#rdirichlet requires vector not dataframe so matrix the tibble
glomcounts <- read_csv("~/Desktop/obmap/r_analysis/data/3dimOB/model/allmice_covariates_trim_voxweights_v4.csv") %>% filter(name %in% kzYprep$name) %>% select(-rep, -slice, -dim) %>%
  rename("count.prop" = "weight")
```

    ## 
    ## ── Column specification ────────────────────────────────────────────────────────
    ## cols(
    ##   name = col_character(),
    ##   rep = col_double(),
    ##   slice = col_double(),
    ##   dim = col_character(),
    ##   count = col_double(),
    ##   weight = col_double(),
    ##   glomarea.apOMP.mlvdseb = col_double(),
    ##   totalarea.seb = col_double(),
    ##   count.clear = col_double(),
    ##   count.unclear = col_double(),
    ##   count.total = col_double(),
    ##   glom.seb = col_double(),
    ##   glom.kz = col_double()
    ## )

``` r
#make covariates from 3d model
cntAP <- blankdata %>% group_by(AntPos) %>% summarise(voxels = n(), vox_norms_ratio2total = voxels/nrow(blankdata)) %>% mutate(dim = "AntPos") %>% rename("slice" = "AntPos")
cntML <- blankdata %>% group_by(MedLat) %>% summarise(voxels = n(), vox_norms_ratio2total = voxels/nrow(blankdata)) %>% mutate(dim = "MedLat") %>% rename("slice" = "MedLat")
cntVD <- blankdata %>% group_by(VenDor) %>% summarise(voxels = n(), vox_norms_ratio2total = voxels/nrow(blankdata)) %>% mutate(dim = "VenDor") %>% rename("slice" = "VenDor")
cntall <- bind_rows(cntAP, cntML) %>% bind_rows(cntVD)
voxweights <- kzYprep %>% 
  left_join(cntall, by = c("dim", "slice")) %>% 
  select(name:dimrep, voxels, vox_norms_ratio2total) %>%
  group_by(dim) %>%
  mutate(voxels_dim_mean = mean(voxels)) %>%
  ungroup() %>%
  mutate(vox_norms_ratio2mean = voxels / voxels_dim_mean) %>%
  rename("vox_r2t" = "vox_norms_ratio2total",
         "vox_r2m" = "vox_norms_ratio2mean") %>%
  select(-voxels_dim_mean)

totaltpm_norms <- readRDS("~/Desktop/obmap/r_analysis/data/make_TPMmtx/totalORtpm_normalizations_justname_210404.RDS") %>%
  rename("rep_r2m" = "rep_norms_ratio2mean",
         "rep_r2t" = "rep_norms_ratio2total",
         "dim_r2m" = "dim_norms_ratio2mean",
         "dim_r2t" = "dim_norms_ratio2total")

#dim_norms and rep_norms is the total OR tpm (without Olfr287,32) of that section divided by the mean or total OR tpm for all sections of that dim or rep
allfeats <- voxweights %>% 
  left_join(totaltpm_norms, by = "name") %>% 
  left_join(glomcounts, by = "name")

tpmvoxnorms <- allfeats %>% 
  select(name:dim_r2t) %>%
  mutate(rt_vt = rep_r2t / vox_r2t,
         dt_vt = dim_r2t / vox_r2t,
         voxflip = 1/vox_r2t) %>%
  select(name:dimrep, rt_vt:voxflip)
#rep_ratio2mean/vox_ratio2mean is equal to rep_ratio2total/vox_ratio2total

tpmvoxnorms %>% 
  gather(key = "norm", value = "factor", rt_vt, dt_vt) %>%
  filter(dimrep == 1) %>%
  ggplot() +
  geom_line(aes(slice, factor, color = norm)) +
  facet_wrap(~ dim, nrow = 1, ncol = 3)
```

![](3d_obmapping_files/figure-gfm/unnamed-chunk-2-1.png)<!-- -->

``` r
kzX <- tpmvoxnorms %>%
  mutate(weight = dt_vt)

saveRDS(kzX, "~/Desktop/obmap/r_analysis/data/3dimOB/model/210404_vox4B_wt_dtvt_OConly.RDS")
```

# Build model

Input: TPM (gene length and sequence depth normalized expression unit)
for 1088 intact ORs across 269 samples across 6 Anterior-Posterior, 3
Medial-Lteral, and 3 Ventral-Dorsal 100um vibratome sliced single OBs.

``` r
tic()
voxalls <- blankdata

# pick prior and Monte Carlo Sample 
alpha <- matrix(0.65, N, D)
iter <- 500
pi_post <- SampleData(kzmY, alpha, iter) 

check <- ifelse(t(voxalls[,kzX$dim]) == kzX$slice, T, F)
weights <- ifelse(check, kzX$weight, 0)
weights <- miniclo(t(weights)) 
colnames(weights) <- rownames(kzmY) 

# We are going to use the ILR basis for computation
contrast.matrix <- create_default_ilr_base(D)
eta_post <- ilr_array(pi_post, contrast.matrix, 2)

# now calculate weighted composition
q_ilr <- array(0, dim=c(nrow(voxalls), D-1, iter))
for (i in 1:iter) q_ilr[,,i] <- weights %*% eta_post[,,i] #matrix multiplication operator

#clear RAM
remove(eta_post)
gc()

#proportional composition rather than in ILR, using driver 18-08-16 new ilrInv_array with transposed q_ilr and coord == 1
ilr_prop <- aperm(q_ilr, c(2,1,3))

remove(q_ilr)
gc() 
q_proportions <- ilrInv_array(ilr_prop, contrast.matrix, 1)

remove(ilr_prop)
gc()

q_proptidy <- aperm(q_proportions, c(3,2,1))
d <- dim(q_proptidy)
q_proptidy <- matrix(q_proptidy, d[1], prod(d[-1]))
quantiles <- c(0.025, .25, .5, .75, .975)
q_prop_quantiles <- t(kzcolQuant(q_proptidy, quantiles))
colnames(q_prop_quantiles) <- paste0("p", 100*quantiles)
voxel=rep(1:d[2], times=D)
coord.ap = voxalls[voxel, "AntPos"] 
coord.ml = voxalls[voxel, "MedLat"]
coord.vd = voxalls[voxel, "VenDor"]
tidy_result <- data.frame(voxel=rep(1:d[2], times=D), 
                          gene=rep(1:d[3], each=nrow(voxalls)), 
                          q_prop_quantiles) %>% 
  bind_cols(coord.ap, coord.ml, coord.vd)
toc()

remove(q_proptidy)
gc()

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

# load save model output, no longer saving tidyresults

``` r
saveRDS(ranked, "~/Desktop/obmap/r_analysis/data/3dimOB/model/210404_vox4B_rk_dtvt_OConly.RDS")
```

# Find gloms for all ORs and get brief statistics

``` r
allORs <- ranked %>% select(olfrname) %>% unique() %>% pull()

allDors <- ListDorML(allORs, by = 50, out = "data")
saveRDS(allDors, "~/Desktop/obmap/r_analysis/data/3dimOB/model/210404_vox4B_ad50_dtvt_OConly.RDS")
```

``` r
ranked <- readRDS("~/Desktop/obmap/r_analysis/data/3dimOB/model/210404_vox4B_rk_dtvt_OConly.RDS")
allDors <- readRDS("~/Desktop/obmap/r_analysis/data/3dimOB/model/210404_vox4B_ad50_dtvt_OConly.RDS")
allDorpoint <- allDors %>% 
  filter(p50 == clustmaxp) %>%
  filter(olfrname != "Olfr3_3_18")

#how many uniquely positioned voxels? duplicated 
allDorpoint %>% 
  select(AntPos:VenDor) %>% 
  unique() %>%
  nrow()
```

    ## [1] 1198

``` r
#how many in ML1
allDorpoint %>% 
  filter(MedLat == 1) %>%
  nrow()
```

    ## [1] 169

``` r
#how many in ML1 are unique?
allDorpoint %>% 
  filter(MedLat == 1) %>%
  select(AntPos:VenDor) %>%
  unique() %>%
  nrow()
```

    ## [1] 122

``` r
#how many in VD1
allDorpoint %>% 
  filter(VenDor == 1) %>%
  nrow()
```

    ## [1] 64

``` r
#1/voxProp: 1151,509,213,65
#voxProp/dimTotal: 1115,77,59,77
#dimTotal/voxProp: 1162,483,220,86
#V4Adtvt: 1111,537,234,88
#V4Bdtvt: 1198,169,122,64 looks v good

#how many predictions per slice per dim, also for unique
ad_ap <- allDorpoint %>% count(AntPos) %>% rename("slice" = "AntPos") %>% mutate(dim = "AntPos", type = "all")
ad_apu <- allDorpoint %>% select(AntPos:VenDor) %>% unique() %>% count(AntPos) %>% rename("slice" = "AntPos") %>% mutate(dim = "AntPos", type = "unique")

ad_ml <- allDorpoint %>% count(MedLat) %>% rename("slice" = "MedLat") %>% mutate(dim = "MedLat", type = "all")
ad_mlu <- allDorpoint %>% select(AntPos:VenDor) %>% unique() %>% count(MedLat) %>% rename("slice" = "MedLat") %>% mutate(dim = "MedLat", type = "unique")

ad_vd <- allDorpoint %>% count(VenDor) %>% rename("slice" = "VenDor") %>% mutate(dim = "VenDor", type = "all")
ad_vdu <- allDorpoint %>% select(AntPos:VenDor) %>% unique() %>% count(VenDor) %>% rename("slice" = "VenDor") %>% mutate(dim = "VenDor", type = "unique")

ad_dimwise <- bind_rows(ad_ap, ad_apu) %>%
  bind_rows(ad_ml) %>%
  bind_rows(ad_mlu) %>%
  bind_rows(ad_vd) %>%
  bind_rows(ad_vdu)

ad_dimwise %>%
  ggplot() +
  geom_line(aes(slice, n, color = type)) + 
  facet_wrap(~ dim) +
  ggtitle("DorsalML prediction positions by dim") +
  theme(legend.position = "bottom")
```

![](3d_obmapping_files/figure-gfm/unnamed-chunk-6-1.png)<!-- -->

``` r
#plot all ORs, both gloms
shellpoints <- DorsalML("Olfr299", out = "notbest")
```

    ## [1] "Olfr299"
    ## [1] 50

``` r
plot_ly(type = "scatter3d", mode = "markers") %>% 
      add_trace(data=shellpoints, x=~AntPos, y=~MedLat, z=~VenDor, 
                color="shell", opacity=0.15,
                text = ~paste('Gene:', olfrname, 
                              '<br>C_size_rank:', clustsizerank, 
                              '<br>C_mean_p50:', clustmeanp,
                              '<br>C_max_p50:', clustmaxp,
                              '<br>Cluster:', clust_unique),
                marker = list(size = 5, color = "grey")) %>%
      add_trace(data=allDorpoint, x=~AntPos, y=~MedLat, z=~VenDor, color=~olfrname,
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
```

    ## Warning in RColorBrewer::brewer.pal(N, "Set2"): n too large, allowed maximum for palette Set2 is 8
    ## Returning the palette you asked for with that many colors

    ## Warning in RColorBrewer::brewer.pal(N, "Set2"): n too large, allowed maximum for palette Set2 is 8
    ## Returning the palette you asked for with that many colors

    ## Fontconfig warning: "/etc/fonts/conf.avail/53-monospace-lcd-filter.conf", line 10: Having multiple values in <test> isn't supported and may not work as expected
    ## Fontconfig warning: "/etc/fonts/conf.avail/53-monospace-lcd-filter.conf", line 10: Having multiple values in <test> isn't supported and may not work as expected

![](3d_obmapping_files/figure-gfm/unnamed-chunk-7-1.png)<!-- -->

## Dorsal constraint function

``` r
DorsalML("Olfr881", topMin = 50, topBy = 50, title = "Olfr881 - DorsalML50")
```

    ## [1] "Olfr881"
    ## [1] "NA 50"
    ## [1] "NA 100"
    ## [1] "NA 150"
    ## [1] 200

    ## Fontconfig warning: "/etc/fonts/conf.avail/53-monospace-lcd-filter.conf", line 10: Having multiple values in <test> isn't supported and may not work as expected
    ## Fontconfig warning: "/etc/fonts/conf.avail/53-monospace-lcd-filter.conf", line 10: Having multiple values in <test> isn't supported and may not work as expected

![](3d_obmapping_files/figure-gfm/unnamed-chunk-8-1.png)<!-- -->

``` r
DorsalML("Olfr1377", topMin = 50, topBy = 50, title = "Olfr1377 - DorsalML50")
```

    ## [1] "Olfr1377"
    ## [1] "NA 50"
    ## [1] 100

    ## Fontconfig warning: "/etc/fonts/conf.avail/53-monospace-lcd-filter.conf", line 10: Having multiple values in <test> isn't supported and may not work as expected
    ## Fontconfig warning: "/etc/fonts/conf.avail/53-monospace-lcd-filter.conf", line 10: Having multiple values in <test> isn't supported and may not work as expected

![](3d_obmapping_files/figure-gfm/unnamed-chunk-8-2.png)<!-- -->

``` r
DorsalML("Olfr16", topMin = 50, topBy = 50, title = "Olfr16/MOR23 - DorsalML50")
```

    ## [1] "Olfr16"
    ## [1] 50

    ## Fontconfig warning: "/etc/fonts/conf.avail/53-monospace-lcd-filter.conf", line 10: Having multiple values in <test> isn't supported and may not work as expected
    ## Fontconfig warning: "/etc/fonts/conf.avail/53-monospace-lcd-filter.conf", line 10: Having multiple values in <test> isn't supported and may not work as expected

![](3d_obmapping_files/figure-gfm/unnamed-chunk-8-3.png)<!-- -->

``` r
DorsalML("Olfr17", topMin = 50, topBy = 50, title = "Olfr17/P2 - DorsalML50")
```

    ## [1] "Olfr17"
    ## [1] 50

    ## Fontconfig warning: "/etc/fonts/conf.avail/53-monospace-lcd-filter.conf", line 10: Having multiple values in <test> isn't supported and may not work as expected
    ## Fontconfig warning: "/etc/fonts/conf.avail/53-monospace-lcd-filter.conf", line 10: Having multiple values in <test> isn't supported and may not work as expected

![](3d_obmapping_files/figure-gfm/unnamed-chunk-8-4.png)<!-- -->

``` r
DorsalML("Olfr15", topMin = 50, topBy = 50, title = "Olfr15/MOR256-17 - DorsalML50")
```

    ## [1] "Olfr15"
    ## [1] 50

    ## Fontconfig warning: "/etc/fonts/conf.avail/53-monospace-lcd-filter.conf", line 10: Having multiple values in <test> isn't supported and may not work as expected
    ## Fontconfig warning: "/etc/fonts/conf.avail/53-monospace-lcd-filter.conf", line 10: Having multiple values in <test> isn't supported and may not work as expected

![](3d_obmapping_files/figure-gfm/unnamed-chunk-8-5.png)<!-- -->

``` r
DorsalML("Olfr160", topMin = 50, topBy = 50, title = "Olfr160/M72 - DorsalML50")
```

    ## [1] "Olfr160"
    ## [1] "NA 50"
    ## [1] "NA 100"
    ## [1] "NA 150"
    ## [1] 200

    ## Fontconfig warning: "/etc/fonts/conf.avail/53-monospace-lcd-filter.conf", line 10: Having multiple values in <test> isn't supported and may not work as expected
    ## Fontconfig warning: "/etc/fonts/conf.avail/53-monospace-lcd-filter.conf", line 10: Having multiple values in <test> isn't supported and may not work as expected

![](3d_obmapping_files/figure-gfm/unnamed-chunk-8-6.png)<!-- -->

``` r
DorsalML("Olfr155", topMin = 50, topBy = 50, title = "Olfr155/MOR37a - DorsalML50")
```

    ## [1] "Olfr155"
    ## [1] 50

    ## Fontconfig warning: "/etc/fonts/conf.avail/53-monospace-lcd-filter.conf", line 10: Having multiple values in <test> isn't supported and may not work as expected
    ## Fontconfig warning: "/etc/fonts/conf.avail/53-monospace-lcd-filter.conf", line 10: Having multiple values in <test> isn't supported and may not work as expected

![](3d_obmapping_files/figure-gfm/unnamed-chunk-8-7.png)<!-- -->

``` r
DorsalML("Olfr1507", topMin = 50, topBy = 50, title = "Olfr1507/MOR28 - DorsalML50")
```

    ## [1] "Olfr1507"
    ## [1] 50

    ## Fontconfig warning: "/etc/fonts/conf.avail/53-monospace-lcd-filter.conf", line 10: Having multiple values in <test> isn't supported and may not work as expected
    ## Fontconfig warning: "/etc/fonts/conf.avail/53-monospace-lcd-filter.conf", line 10: Having multiple values in <test> isn't supported and may not work as expected

![](3d_obmapping_files/figure-gfm/unnamed-chunk-8-8.png)<!-- -->

``` r
DorsalML("Olfr3_3_18", topMin = 50, topBy = 50, title = "testOR") 
```

    ## [1] "Olfr3_3_18"
    ## [1] "NA 50"
    ## [1] "NA 100"
    ## [1] 150

    ## Fontconfig warning: "/etc/fonts/conf.avail/53-monospace-lcd-filter.conf", line 10: Having multiple values in <test> isn't supported and may not work as expected
    ## Fontconfig warning: "/etc/fonts/conf.avail/53-monospace-lcd-filter.conf", line 10: Having multiple values in <test> isn't supported and may not work as expected

![](3d_obmapping_files/figure-gfm/unnamed-chunk-8-9.png)<!-- -->

``` r
# testOR med glom at 3,3,18 lat glom at 2,9,17
```

# 3D plots of OR probability as selected and constrained using my own functions

## Plot X best voxels for an OR as ranked by SNR

Each voxel holds a probability for every OR, as determined by the
expression of that ORs in the dimensional slices that intersect at that
voxel. The total of the probability for all ORs in a voxel sums to 1.
SNR is signal-to-noise ratio, essentially how much does an individual
voxel stand out from the average probability of its 3 dimensional slices
(all voxels at that single dimension position). Idea was to eliminate
ectopics whose high expression across the whole OB dominated the voxel
probability for nearly all voxels. 72 is \~2.5% of total voxel positions
and presents a robust visualization, meaning that increasing this number
typically results in adding neighboring voxels.

``` r
Scat_rank("Olfr881", 72, title = "Olfr881 - Top 72 SnR voxels")
Scat_rank("Olfr1377", 72, title = "Olfr1377 - Top 72 SnR voxels")
zapORs <- c("Olfr15", "Olfr16", "Olfr17", "Olfr155", "Olfr160", "Olfr1507")
#other zapiec & mombaerts labeled ORs
# Scat_rank("Olfr16",15)
# Scat_rank("Olfr15", 15)
# Scat_rank("Olfr17", 50)
# Scat_rank("Olfr160", 15)
# Scat_rank("Olfr155", 15)
# Scat_rank("Olfr1507", 30)
```

## Plot best X voxels based on raw p50 (probability)

Using raw probability here, its pretty much the same as SNR after
removing ectopics.

``` r
#ranked %>% filter(olfrname == "Olfr160") %>% filter(MedLat < 3) %>% arrange(desc(p50)) %>% select(voxrankperOR, AntPos:VenDor)

p50plot("Olfr881", 72, title = "Olfr881 - Top 72 probability voxels")
p50plot("Olfr1377", 72, title = "Olfr1377 - Top 72 probability voxels")
# p50plot("Olfr16", 60)
# p50plot("Olfr17", 70)
# p50plot("Olfr15", 36)
# p50plot("Olfr160", 36)
# p50plot("Olfr155", 36)
# p50plot("Olfr1507", 36)
```

## Clustering top ranked points

Given a number of top probability ranked points, cluster spatial
neighbors.

``` r
Cluster("Olfr881", 1000, minClustSize = 2, out = "plot", title = "Olfr881 - Clustered top 100 probability voxels")
Cluster("Olfr1377", 125, minClustSize = 2, out = "plot", title = "Olfr1377 - Clustered top 100 probability voxels")
# Cluster("Olfr16", out = 100, "plot")
# Cluster("Olfr17", out = 100, "plot")
# Cluster("Olfr15", out = 100, "plot")
# Cluster("Olfr160", out = 100, "plot")
# Cluster("Olfr155", out = 100, "plot")
# Cluster("Olfr1507", out = 100, "plot")
```

## Pick the best Medial and Lateral cluster

Plot best Medial and best Lateral cluster from the set of clusters
across a symmetry line that is calculated using single dimension heatmap
data For clarity, it is not picking the 2 best clusters but the best
cluster from the medial half and the best cluster from the lateral half.

``` r
BestML("Olfr881", clustersPerHalfBulb = 1, title = "Olfr881 - Highest Probability Medial and Lateral Clusters")
BestML("Olfr1377", clustersPerHalfBulb = 1, title = "Olfr1377 - Highest Probability Medial and Lateral Cluster")
# BestML("Olfr16", clustersPerHalfBulb = 1)
# BestML("Olfr17", clustersPerHalfBulb = 1)
# BestML("Olfr15", clustersPerHalfBulb = 1)
# BestML("Olfr160", clustersPerHalfBulb = 1)
# BestML("Olfr155", clustersPerHalfBulb = 1)
# BestML("Olfr1507", clustersPerHalfBulb = 1)
```

# features from intro, moved here for simplicity during making of many models

``` r
#symline visualization
ggplot() +
  geom_blank() +
  geom_abline(aes(slope = 0.32655, intercept = 7.61119, color = "rep1")) +
  geom_abline(aes(slope = 0.24922, intercept = 8.83890, color = "rep2")) +
  geom_abline(aes(slope = 0.287885, intercept = 8.225045, color = "mean")) +
  xlim(0,23) +
  ylim(0,22) + 
  ggtitle("Glomeruli symmetry line as calculated by single dimension heatmap peaks") +
  xlab("Mean position of top 2 AP peaks") +
  ylab("Mean position of top 2 ML peaks")
```

![](3d_obmapping_files/figure-gfm/unnamed-chunk-13-1.png)<!-- -->

``` r
#blankOB colored by DV position on a black background, see title slide of 2020 presentations
dvml_blank <- blankdata %>%
  mutate(dimdv = ifelse(VenDor >= 13, 
                        "Dorsal", 
                        "Ventral"),
         dimml = ifelse(MedLat > 17, 
                        "Lateral",
                        ifelse(MedLat < 5,
                               "Medial",
                               NA)),
         dim = ifelse(is.na(dimml), dimdv, dimml)) %>%
  rowwise() %>%
  mutate(symdim = ifelse(MedLat == round(symline$mlvals[which(symline$apvals == AntPos)]),"SymLine", dim)) %>%
  ungroup()

#blankOB
plot_ly(dvml_blank, x = ~MedLat, y = ~AntPos, z = ~VenDor, color = ~VenDor, 
        text = ~paste('Dim: ', dim,
                      '<br>AP:', AntPos, 
                      '<br>ML:', MedLat,
                      '<br>VD:', VenDor),
        marker = list(size = 8, symbol = "circle"),
        type = 'scatter3d',
        mode = 'markers',
        hoverinfo = "none",
        hovertext="none") %>%
  layout(scene = list(xaxis = list(title = 'Medial-Lateral',
                                   backgroundcolor="rgb(200, 200, 230",
                                   gridcolor="rgb(255,255,255)",
                                   showbackground=TRUE,
                                   zerolinecolor="rgb(255,255,255"),
                      yaxis = list(title = 'Anterior-Posterior',
                                   backgroundcolor="rgb(200, 200, 230",
                                   gridcolor="rgb(255,255,255)",
                                   showbackground=TRUE,
                                   zerolinecolor="rgb(255,255,255)"),
                      zaxis = list(title = 'Ventral-Dorsal',
                                   backgroundcolor="rgb(200, 200, 230",
                                   gridcolor="rgb(255,255,255)",
                                   showbackground=TRUE,
                                   zerolinecolor="rgb(255,255,255)")))
```

    ## Fontconfig warning: "/etc/fonts/conf.avail/53-monospace-lcd-filter.conf", line 10: Having multiple values in <test> isn't supported and may not work as expected
    ## Fontconfig warning: "/etc/fonts/conf.avail/53-monospace-lcd-filter.conf", line 10: Having multiple values in <test> isn't supported and may not work as expected

![](3d_obmapping_files/figure-gfm/unnamed-chunk-13-2.png)<!-- -->

``` r
#rankedolfr is a tidy_result ranked tibble using all intact olfr tpm
#rankedchemo is a tidy_result ranked tibble using all intact chemosensory gene tpm
#orbyorcor <- tibble(olfrname = ranked_olfr2$olfrname, orp50 = ranked_olfr2$p50, chemp50 = ranked_chemORs$p50) %>% group_by(olfrname) %>% summarise(cor = cor(orp50, chemp50)) %>% left_join(info, by = "olfrname")
#orXor_cor <- read_csv("~/Desktop/obmap/r_analysis/3dimOB/output/chemoVolfr_p50.csv")
#ggplot(orXor_cor) + geom_violin(aes(oe_region, cor)) + 
#  ggtitle("150iter, using all Chemo does not significantly alter voxel probability rankings")
```
