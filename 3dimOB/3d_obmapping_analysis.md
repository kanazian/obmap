Analysis of 3D predictions
================
kanazian
Oct 25, 2020

# Goal:

Examine OBmap 3D model predicted positions for concordance with OE DV,
FI surface, and other features as well as known OR positions.

# Plan:

1.  Better definition of dorsal-ventral zone based on OE data related to
    predicted position.
2.  Define area of functional imaging using wachowiak functional imaging
    area enriched ORs.
3.  Annotate positions of Mombaerts ORs onto model and compare predicted
    positions and distances.

<!-- end list -->

``` r
knitr::opts_chunk$set(warning=F)
library(plotly) 
library(tidyverse)
library(cowplot)
library(patchwork)

# functions -----------------------------------------------------------
#https://stackoverflow.com/questions/49215193/r-error-cant-join-on-because-of-incompatible-types
MatchColClasses <- function(df1, df2) {

  sharedColNames <- names(df1)[names(df1) %in% names(df2)]
  sharedColTypes <- sapply(df1[,sharedColNames], class)

  for (n in sharedColNames) {
     class(df2[, n]) <- sharedColTypes[n]
  }

  return(df2)
}

#define clusters of best X p50 points using a pairwise matrix and ability to output plots and data
#given a number of top ranking positions for an OR, cluster the points based on spatial position
Cluster <- function (olfr, topX = 72, minClustSize = 5, chooseOut = "plot", title = NA) {
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
  
  all_out <- bind_rows(df_clustered, cut_out) %>% unique()
  
  #return various things 'rgb(10,10,10)'
  if (chooseOut == "data") {
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
BestML <- function(olfr, topMin = 50, topMax = 200, topBy = 25, minSize = 2,
                   clustersPerHalfBulb = 1, chooseOut = "plot", title = NA) {
  cphb <- 0
  topStep <- topMin
  while (cphb < clustersPerHalfBulb) {
    df_in <- Cluster(olfr, topX = topStep, minClustSize = minSize, chooseOut = "data")
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
  if (chooseOut == "data") {
    return(df_all)
  } else if (chooseOut == "best") {
    return(df_best)
  } else if (chooseOut == "notbest") {
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
ListML <- function(x, chooseOut = "plot", title = NA) {
  #use a list to build a df of unknown size instead of bind_row each iteration
  list_out <- vector("list", length = length(x))
  for (i in 1:length(x)) {
    list_out[[i]] <- BestML(x[i], topMin = 50, topMax = 150, topBy = 50, 
                            minSize = 2, clustersPerHalfBulb = 1, chooseOut = "best")
  } #endfor
  
  #always need notbest for the shape shell1
  notbest <- BestML(x[1], topMin = 100, topMax = 150, topBy = 50, minSize = 2, 
                    clustersPerHalfBulb = 1, chooseOut = "notbest")
  df_out <- bind_rows(list_out)
  
  #output
  if (chooseOut == "data") {
    return(df_out)
  } else if (chooseOut == "point") {
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
DorsalML <- function(olfr, topMin = 100, topBy = 100, minSize = 2,
                     clustIn = 5, clustOut = 1, chooseOut = "plot",
                     title = NA) {
  print(olfr)
  clustFound <- 0
  topStep <- topMin
  while (clustFound != clustOut) {
    df_in <- Cluster(olfr, topX = topStep, minClustSize = minSize, chooseOut = "data")
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
  if (chooseOut == "data") {
    return(df_all)
  } else if (chooseOut == "best") {
    return(df_best)
  } else if (chooseOut == "notbest") {
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
ListDorML <- function(x, chooseOut = "plot", title = NA) {
  #use a list to build a df of unknown size instead of bind_row each iteration
  list_out <- vector("list", length = length(x))
  for (i in 1:length(x)) {
    list_out[[i]] <- DorsalML(x[i], topBy = 200, chooseOut = "best")
  } #endfor
  
  #always need notbest for the shape shell1
  notbest <- DorsalML(x[1], chooseOut = "notbest")
  df_out <- bind_rows(list_out)
  
  #output
  if (chooseOut == "data") {
    return(df_out)
  } else if (chooseOut == "point") {
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

#using heatmap peak data, create a 3d raw point or find nearest blank OB shell point to that raw point, options to use tan zone to assign a VD position or use average of VD peaks
Heat3D <- function(olfr, dimrep = 1, dv = "oe", raw = F, chooseOut = "plot") {
  df_in <- heatmap_peaks %>% filter(olfrname %in% olfr) %>% arrange(side)
  
  if (nrow(df_in) == 0) {
    print(paste("No", olfr, "Found"))
    next
  } #endif
  
  if (dv == "oe") {
    df_dv <- df_in %>%
      mutate(VD_selected = VDtz,
             VD_value = "Tan Based")
  } else {
    df_dv <- df_in %>%
      mutate(VD_selected = VDavg,
             VD_value = "Avg Based")
  } #endif
  
  if (dimrep == 1) {
    df_out <- df_dv %>% 
      select(olfrname, AntPos10, APval10, MedLat8, MLval8, VD_selected, VD_value, side) %>%
      rename(AP_selected = AntPos10, AP_value = APval10,
             ML_selected = MedLat8, ML_value = MLval8)
  } else if (dimrep == 2) {
    df_out <- df_dv %>% 
      select(olfrname, AntPos13, APval13, MedLat11, MLval11, VD_selected, VD_value, side) %>%
      rename(AP_selected = AntPos13, AP_value = APval13,
             ML_selected = MedLat11, ML_value = MLval11)
  } else if (dimrep == 3) {
    df_out <- df_dv %>% 
      select(olfrname, AntPos15, APval15, MedLat16, MLval16, VD_selected, VD_value, side) %>%
      rename(AP_selected = AntPos15, AP_value = APval15,
             ML_selected = MedLat16, ML_value = MLval16)
  } else if (dimrep == "123") {
    df_out <- df_dv %>% rowwise() %>%
      mutate(AP_selected = mean(c(AntPos10, AntPos13, AntPos15)),
             AP_value = mean(c(APval10, APval13, APval15)),
             ML_selected = mean(c(MedLat8, MedLat11, MedLat16)),
             ML_value = mean(c(MLval8, MLval11, MLval16))) %>%
      ungroup()
  } else if (dimrep == "12") {
    df_out <- df_dv %>% rowwise() %>%
      mutate(AP_selected = mean(c(AntPos10, AntPos13)),
             AP_value = mean(c(APval10, APval13)),
             ML_selected = mean(c(MedLat8, MedLat11)),
             ML_value = mean(c(MLval8, MLval11))) %>%
      ungroup()
  } #endif
  
  if (raw == T) {
    if (chooseOut == "plot") {
      plot <- plot_ly(type = "scatter3d", mode = "markers") %>%
        add_trace(data = blankdata,
                  x = ~AntPos, y=~MedLat, z=~VenDor,
                  color="shell", opacity=0.15,
                  text = ~paste('AntPos:', AntPos,
                                '<br>MedLat:', MedLat,
                                '<br>VenDor:', VenDor),
                  marker = list(size = 6)) %>%
        add_trace(data = df_out, 
                  x = ~AP_selected, 
                  y = ~ML_selected, 
                  z = ~VD_selected, 
                  color = ~olfrname,
                  text = ~paste('Gene:', olfrname, 
                                '<br>AntPos:', AP_selected,
                                '<br>APvalue', AP_value,
                                '<br>MedLat:', ML_selected,
                                '<br>MLvalue:', ML_value,
                                '<br>VenDor:', VD_selected,
                                '<br>VDvalue:', VD_value,
                                '<br>Side:', side,
                                '<br>Point:', "Raw"),
                  marker = list(size = 6, 
                                line = list(color='black', width = 0.5))) %>%
        layout(scene = list(xaxis = list(title = 'Anterior-Posterior'),
                            yaxis = list(title = 'Medial-Lateral'),
                            zaxis = list(title = 'Ventral-Dorsal')))
      return(plot)
    } else if (chooseOut == "data") {
      return(df_out)
    }
  } else if (raw == F) { 
    df_close <- tibble("AntPos" = numeric(), "MedLat" = numeric(), "VenDor" = numeric())
    for (i in 1:nrow(df_out)) {
      distances <- vector(mode = "numeric", length = nrow(blankdata))
      sidevec <- c("Lateral", "Medial")
      for (j in 1:nrow(blankdata)) {
        distances[j] <- sqrt((df_out$AP_selected[i] - blankdata$AntPos[j])^2 +
                               (df_out$ML_selected[i] - blankdata$MedLat[j])^2 + 
                               (df_out$VD_selected[i] - blankdata$VenDor[j])^2)
        df_close[i,] <- blankdata[which(rank(distances, ties.method = "first") == 1),]
      }
    }
    
    df_close_out <- df_close %>% 
      mutate(side = sidevec) %>%
      left_join(df_out, by = "side") %>%
      select(olfrname, AntPos, VenDor, MedLat, side, everything())
    
    if (chooseOut == "plot") {
      plot <- plot_ly(type = "scatter3d", mode = "markers") %>%
        add_trace(data = blankdata,
                  x = ~AntPos, y=~MedLat, z=~VenDor,
                  color="shell", opacity=0.15,
                  text = ~paste('AntPos:', AntPos,
                                '<br>MedLat:', MedLat,
                                '<br>VenDor:', VenDor),
                  marker = list(size = 6)) %>%
        add_trace(data = df_close_out, 
                  x = ~AntPos, 
                  y = ~MedLat, 
                  z = ~VenDor, 
                  color = ~olfrname,
                  text = ~paste('Gene:', olfrname, 
                                '<br>AntPos:', AP_selected,
                                '<br>APvalue', AP_value,
                                '<br>MedLat:', ML_selected,
                                '<br>MLvalue:', ML_value,
                                '<br>VenDor:', VD_selected,
                                '<br>VDvalue:', VD_value,
                                '<br>Side:', side,
                                '<br>Point:', "Nearest Shell"),
                  marker = list(size = 6, 
                                line = list(color='black', width = 0.5))) %>%
        layout(scene = list(xaxis = list(title = 'Anterior-Posterior'),
                            yaxis = list(title = 'Medial-Lateral'),
                            zaxis = list(title = 'Ventral-Dorsal')))
      return(plot)
    } else if (chooseOut == "data") {
      return(df_close_out)
    } #endif
  } #endif
} #endfunction

#compute distances between heatmap peak point and good point data featuring eurovision movie memes
DistHeat3D <- function(olfr_list, heat_dimrep = 1, heat_dv = "oe", heat_raw = F) {
  doubletrouble <- tibble(olfrname = character(), side = character(), 
                          distance = numeric(), ap_ml_vd_3d = character(), 
                          ap_ml_vd_heat = character(), .rows = 0)
  
  for (i in 1:length(olfr_list)) {
    heatpoint <- Heat3D(olfr_list[i], dimrep = heat_dimrep, 
                        dv = heat_dv, raw = heat_raw, chooseOut = "data") %>%
      select(olfrname, AntPos, MedLat, VenDor, side) %>% 
      mutate(origin = "heat",
             ap_ml_vd_heat = paste(AntPos, MedLat, VenDor,sep = "_"))
    
    goodpoint <- good_point %>% 
      filter(olfrname == olfr_list[i]) %>% 
      select(olfrname, AntPos, MedLat, VenDor, side) %>% 
      mutate(origin = "3d",
             ap_ml_vd_3d = paste(AntPos, MedLat, VenDor,sep = "_")) 
    
    doubletrouble <- bind_rows(doubletrouble, 
                                  left_join(heatpoint, goodpoint, 
                                            by = c("olfrname", "side"), 
                                            suffix = c("_h", "_3")) %>%
                                    mutate(distance = sqrt((AntPos_h - AntPos_3)^2 +
                                                             (MedLat_h - MedLat_3)^2 +
                                                             (VenDor_h - VenDor_3)^2)) %>%
                                    select(olfrname, side, distance, 
                                           ap_ml_vd_3d, ap_ml_vd_heat))
  } #endfor
  doubletrouble <- doubletrouble %>% 
    mutate(dimrep = heat_dimrep,
           vd_assignment = heat_dv,
           use_raw = heat_raw)
  return(doubletrouble)
} #endfunction


Plot_predictions <- function(gene_vector, varcolor=~olfrname, chooseOut = "point", title = NA) {
  clust_df <- filter_preds %>% filter(olfrname %in% gene_vector)
  if (chooseOut == "point") {
    point_df <- clust_df %>% filter(p50 == clustmaxp)
    point <- plot_ly(type = "scatter3d", mode = "markers") %>% 
      add_trace(data=blankdata, x=~AntPos, y=~MedLat, z=~VenDor, 
                color="shell", opacity=0.15,
                marker = list(size = 6)) %>%
      add_trace(data=point_df, x=~AntPos, y=~MedLat, z=~VenDor, color=varcolor,
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
    return(point)
  } else {
    clust <- plot_ly(type = "scatter3d", mode = "markers") %>% 
      add_trace(data=blankdata, x=~AntPos, y=~MedLat, z=~VenDor, 
                color="shell", opacity=0.15,
                marker = list(size = 6)) %>%
      add_trace(data=clust_df, x=~AntPos, y=~MedLat, z=~VenDor, color=~olfrname,
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
    return(clust)
  } #endif chooseOut
} #endfunction
```

# Load data

``` r
kzY <- read_csv("~/Desktop/obmap/r_analysis/3dimOB/input/covarintactchemo_over_samples_200923.csv", col_names = TRUE) %>% select(-X1) %>% filter(dimrep == 1 | dimrep == 2)
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
ectopic <- c("Olfr287","Olfr32")
kzYgood <- dplyr::select(kzY, -ectopic) %>% dplyr::select(-name, -rep, -slice, -dim, -dimrep)
```

    ## Note: Using an external vector in selections is ambiguous.
    ## ℹ Use `all_of(ectopic)` instead of `ectopic` to silence this message.
    ## ℹ See <https://tidyselect.r-lib.org/reference/faq-external-vector.html>.
    ## This message is displayed once per session.

``` r
kzmY <- as.matrix(kzYgood)

#feature info
info <- read_csv("~/Desktop/obmap/r_analysis/3dimOB/input/knowntanwavgFI.csv", col_names = TRUE) %>% 
  rename("olfrname" = "gene") %>%
  select(olfrname:RTP, known, lowTPM)
```

    ## 
    ## ── Column specification ────────────────────────────────────────────────────────
    ## cols(
    ##   .default = col_double(),
    ##   gene = col_character(),
    ##   tan_zone = col_character(),
    ##   oe_region = col_character(),
    ##   RTP = col_character(),
    ##   fisurface = col_logical(),
    ##   ish_id = col_logical(),
    ##   lacz_xg = col_logical(),
    ##   fp_xg = col_logical(),
    ##   med_glom = col_character(),
    ##   lat_glom = col_character(),
    ##   citation = col_character(),
    ##   known = col_logical(),
    ##   lowTPM = col_logical()
    ## )
    ## ℹ Use `spec()` for the full column specifications.

``` r
#line of bulb symmetry as found using single-dimension heatmap data
#the average of the two calculated lines is used to call whether a predicted glomeruli position is medial or lateral
symline <- read_csv("~/Desktop/obmap/r_analysis/3dimOB/input/symline.csv")
```

    ## 
    ## ── Column specification ────────────────────────────────────────────────────────
    ## cols(
    ##   apvals = col_double(),
    ##   mlvals = col_double()
    ## )

``` r
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

![](3d_obmapping_analysis_files/figure-gfm/unnamed-chunk-1-1.png)<!-- -->

``` r
# Example 3D left OB for future orientation of dimensions
#blankdata <- Scat_rank("Olfr881", 1, "data") %>% select(AntPos:VenDor) %>% arrange(AntPos, MedLat, VenDor)
#write_csv(blankdata, "~/Desktop/rproj/obmap/allmice/v21_gen25/base_files/blankOBcoords.csv")
#blankOBcoords200922
blankdata <- read_csv("~/Desktop/obmap/r_analysis/3dimOB/input/blankOBcoords200922.csv")
```

    ## 
    ## ── Column specification ────────────────────────────────────────────────────────
    ## cols(
    ##   AntPos = col_double(),
    ##   MedLat = col_double(),
    ##   VenDor = col_double()
    ## )

``` r
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

![](3d_obmapping_analysis_files/figure-gfm/unnamed-chunk-1-2.png)<!-- -->

``` r
#load 3d model predictions (clusters) from ListDorML function with non-OLfr genes expected to be Dorsal
all_predictions <- read_csv("~/Desktop/obmap/r_analysis/3dimOB/output/allchemo_dornonor_ldML_201012.csv")
```

    ## 
    ## ── Column specification ────────────────────────────────────────────────────────
    ## cols(
    ##   .default = col_double(),
    ##   olfrname = col_character(),
    ##   side = col_character(),
    ##   tan_zone = col_character(),
    ##   oe_region = col_character(),
    ##   RTP = col_character(),
    ##   known = col_logical(),
    ##   lowTPM = col_logical(),
    ##   test = col_character()
    ## )
    ## ℹ Use `spec()` for the full column specifications.

``` r
#allORs2 <- read_csv("~/Desktop/obmap/r_analysis/3dimOB/output/allchemo_ldML_200925.csv") #%>% filter(p50 == clustmaxp)
#aORcor <- tibble(olfrname = allORs$olfrname, oneAP = allORs$AntPos, oneML = allORs$MedLat, oneVD = allORs$VenDor, twoAP = allORs2$AntPos, twoML = allORs2$MedLat, twoVD = allORs2$VenDor) %>% filter(!str_detect(olfrname, "Olfr")) %>% summarise(corAP = cor(oneAP, twoAP), corML = cor(oneML, twoML), corVD = cor(oneVD, twoVD))
# main change is corVD, due to addition of dorsal rating for VMNrs and TAARs
```

# Filter

``` r
#get summary statistics for each OR
kzmY[1:5,1:5]
```

    ##      Olfr299 Olfr109 Olfr281 Olfr1015 Olfr1347
    ## [1,]    0.00    0.00    0.33     5.77      0.0
    ## [2,]  277.35    0.00  602.80   228.19      0.3
    ## [3,]    0.00   99.66 2451.54     0.11      0.0
    ## [4,]    0.00    0.00 1882.25   340.79      0.0
    ## [5,]    0.00    0.00 2612.72     0.00      0.0

``` r
ornames <- colnames(kzmY)
max_tpm <- vector(mode = "numeric", length = ncol(kzmY))
mean_tpm <- vector(mode = "numeric", length = ncol(kzmY))
sd_tpm <- vector(mode = "numeric", length = ncol(kzmY))
for (i in 1:ncol(kzmY)) {
  max_tpm[i] <- max(kzmY[,i])
  mean_tpm[i] <- mean(kzmY[,i])
  sd_tpm[i] <- sd(kzmY[,i])
}

# remove ORs with low max/mean values and low max tpm values, can probably do better
metrics <- tibble(ornames, max_tpm, mean_tpm, sd_tpm) %>%
  rowwise() %>%
  mutate(max2mean = max_tpm/mean_tpm,
         max2mean_lo = ifelse(max2mean < 5, T, F),
         maxSDmean = (max_tpm - mean_tpm)/sd_tpm,
         maxSDmean_lo = ifelse(maxSDmean < 3, T, F),
         coefvar = sd_tpm/mean_tpm,
         coefvar_lo = ifelse(coefvar < 1, T, F),
         maxtpm_lo = ifelse(max_tpm < 10, T, F))

goodORs <- metrics %>%
  filter(max2mean_lo == F) %>%
  filter(maxSDmean_lo == F) %>%
  filter(coefvar_lo == F) %>%
  filter(maxtpm_lo == F) %>%
  select(ornames) %>%
  as_vector() 

filter_preds <- all_predictions
filter_olfrs <- filter_preds %>% filter(str_detect(olfrname, "Olfr")) %>% select(olfrname) %>% unique() %>% as_vector()
filter_vmnrs <- filter_preds %>% filter(str_detect(olfrname, "Vmn")) %>% select(olfrname) %>% unique() %>% as_vector()
filter_taars <- filter_preds %>% filter(str_detect(olfrname, "Taar")) %>% select(olfrname) %>% unique() %>% as_vector()
filter_all <- filter_preds %>% select(olfrname) %>% unique() %>% as_vector()
```

# Class 1 vs Class 2 positions

Expect Class 1 OR positions to be primarily dorsal-anterior to
dorsal-central. Filtered out some ORs based on max to mean TPM ratio
lower than 10 (968 ORs remaining). These could be ORs that are poorly
enriched/dropout or have non-traditional expression across the OB
(requires closer
investigation).

``` r
#lets look at class 1 vs class 2 ORs, note that is it possible for a voxel to hold multiple OR cluster points
filter_preds <- filter_preds %>% 
  mutate(class_fct = as_factor(class))
Plot_predictions(filter_olfrs, varcolor = ~class_fct)
```

![](3d_obmapping_analysis_files/figure-gfm/unnamed-chunk-3-1.png)<!-- -->

``` r
#a look at proportions to deal with point density
class_props <- filter_preds %>% 
  filter(!is.na(class)) %>%
  filter(p50 == clustmaxp) %>%
  group_by(olfrname) %>%
  mutate(VDmeanpos = round(mean(VenDor))) %>%
  ungroup() %>%
  mutate(propC1 = ifelse(class == 1, T, F)) %>%
  group_by(VDmeanpos) %>%
  summarise(count = n(),
            sum_class1 = sum(propC1),
            proportion_class1 = sum_class1/count)

ggplot(class_props) + 
  geom_bar(aes(VDmeanpos, proportion_class1), stat = "identity") + 
  ggtitle("Proportion of Class 1 ORs across Ventral-Dorsal Axis") + 
  xlab("Ventral  <<<    100um sections   >>>  Dorsal")
```

![](3d_obmapping_analysis_files/figure-gfm/unnamed-chunk-3-2.png)<!-- -->

``` r
ggplot(class_props) + 
  geom_bar(aes(VDmeanpos, sum_class1), stat = "identity") + 
  ggtitle("Number of Class 1 ORs across Ventral-Dorsal Axis") + 
  xlab("Ventral  <<<    100um sections   >>>  Dorsal")
```

![](3d_obmapping_analysis_files/figure-gfm/unnamed-chunk-3-3.png)<!-- -->

# OE region positions (zonal expression of OR as determined by Matsunami Lab DiffE)

3 samples of dorsal OE vs 3 samples of ventral OE Could also examine
relation to more discrete tan et al. zone indices but this is more
readable.

``` r
Plot_predictions(filter_olfrs, varcolor = ~oe_region, title = "OE Zone of cluster points for 1115 ORs")
```

![](3d_obmapping_analysis_files/figure-gfm/unnamed-chunk-4-1.png)<!-- -->

``` r
#a look at proportions to deal with point density
oe_props <- filter_preds %>% 
  filter(str_detect(olfrname, "Olfr")) %>%
  filter(!is.na(oe_region)) %>%
  group_by(olfrname) %>%
  mutate(VDmeanpos = round(mean(VenDor))) %>%
  ungroup() %>%
  mutate(isdor = ifelse(oe_region == "Dorsal", T, F)) %>%
  group_by(VDmeanpos) %>%
  summarise(count = n(),
            sum_dorsal = sum(isdor),
            proportion_dorsal = sum_dorsal/count,
            proportion_ventral = 1-proportion_dorsal)

ggplot(oe_props) + 
  geom_bar(aes(VDmeanpos, proportion_dorsal), stat = "identity") + 
  ggtitle("Proportion of Dorsal OE Zone ORs across Ventral-Dorsal Axis") + 
  xlab("Ventral  <<<    100um sections   >>>  Dorsal")
```

![](3d_obmapping_analysis_files/figure-gfm/unnamed-chunk-4-2.png)<!-- -->

``` r
ggplot(oe_props) + 
  geom_bar(aes(VDmeanpos, sum_dorsal), stat = "identity") + 
  ggtitle("Number of Dorsal OE Zone ORs across Ventral-Dorsal Axis") + 
  xlab("Ventral  <<<    100um sections   >>>  Dorsal")
```

![](3d_obmapping_analysis_files/figure-gfm/unnamed-chunk-4-3.png)<!-- -->

# Tan zone indexes

``` r
filter_preds <- filter_preds %>% 
  mutate(tzsimplest = ifelse(tzsimple <= 5, 6-tzsimple, NA),
         tzbins = round(tzsimplest/0.5)) %>%
  filter(!is.na(tzsimplest))
Plot_predictions(filter_olfrs, varcolor = ~tzsimplest)
```

![](3d_obmapping_analysis_files/figure-gfm/unnamed-chunk-5-1.png)<!-- -->

``` r
Plot_predictions(filter_olfrs, varcolor = ~tzbins)
```

![](3d_obmapping_analysis_files/figure-gfm/unnamed-chunk-5-2.png)<!-- -->

``` r
#predictions for medial vs lateral are probably something good to discuss in paper
filter_preds %>% 
  filter(p50 == clustmaxp) %>%
  filter(side == "Medial") %>%
  ggplot() +
  geom_jitter(aes(tzbins, VenDor)) +
  geom_smooth(aes(tzbins, VenDor)) +
  ggtitle("Medial points - tan and prediction position")
```

    ## `geom_smooth()` using method = 'loess' and formula 'y ~ x'

![](3d_obmapping_analysis_files/figure-gfm/unnamed-chunk-6-1.png)<!-- -->

``` r
filter_preds %>% 
  filter(p50 == clustmaxp) %>%
  filter(side == "Lateral") %>%
  ggplot() +
  geom_jitter(aes(tzbins, VenDor)) +
  geom_smooth(aes(tzbins, VenDor)) +
  ggtitle("Lateral points - tan and prediction position")
```

    ## `geom_smooth()` using method = 'loess' and formula 'y ~ x'

![](3d_obmapping_analysis_files/figure-gfm/unnamed-chunk-6-2.png)<!-- -->

``` r
filter_preds %>% 
  filter(p50 == clustmaxp) %>%
  filter(side == "Lateral") %>%
  ggplot() +
  geom_boxplot(aes(as_factor(tzbins), VenDor)) +
  ggtitle("Lateral points - tan and prediction position")
```

![](3d_obmapping_analysis_files/figure-gfm/unnamed-chunk-6-3.png)<!-- -->

``` r
#find outliers in terms of DV voxel and DV index (tan, luis, etc.) relationship (outside of 1.645 SDs from mean or ~90%)
#index is a colname, using {{colname}} to pass colname as function arg
Find_DVoutliers <- function(index, sd_from_mean = 1.645) {
  for (i in min(filter_preds$VenDor):max(filter_preds$VenDor)) {
    filtered_lat <- filter_preds %>%
      filter(p50 == clustmaxp) %>%
      filter(side == "Lateral") %>%
      filter(VenDor == i) %>%
      mutate(temp_idx = {{index}})
    lat_sum <- filtered_lat %>% summarise(avg_idx = mean(temp_idx), sd_idx = sd(temp_idx))
    lat_midsdlo <- lat_sum$avg_idx[1] - sd_from_mean * lat_sum$sd_idx[1]
    lat_midsdhi <- lat_sum$avg_idx[1] + sd_from_mean *lat_sum$sd_idx[1]
    filtered_lat <- filtered_lat %>%
      mutate(checkhi = ifelse({{index}} > lat_midsdhi, 500, 0),
             checklo = ifelse({{index}} < lat_midsdlo, 50, 0),
             checkmid = ifelse(between({{index}}, lat_midsdlo, lat_midsdhi), 5, 0),
             barhi = lat_midsdhi,
             barlo = lat_midsdlo,
             barmid = lat_sum$avg_idx[1]) %>%
      rowwise() %>%
      mutate(pack = checkhi + checklo + checkmid) %>%
      ungroup() %>%
      select(olfrname, side, VenDor, p50, {{index}}, checkhi, checklo, checkmid, barhi, barlo, barmid, pack) %>%
      filter(pack > 0)
    
    filtered_med <- filter_preds %>%
      filter(p50 == clustmaxp) %>%
      filter(side == "Medial") %>%
      filter(VenDor == i) %>%
      mutate(temp_idx = {{index}})
    med_sum <- filtered_med %>% summarise(avg_idx = mean(temp_idx), sd_idx = sd(temp_idx))
    med_midsdlo <- med_sum$avg_idx[1] - sd_from_mean * med_sum$sd_idx[1]
    med_midsdhi <- med_sum$avg_idx[1] + sd_from_mean * med_sum$sd_idx[1]
    filtered_med <- filtered_med %>%
      mutate(checkhi = ifelse({{index}} > med_midsdhi, 500, 0),
             checklo = ifelse({{index}} < med_midsdlo, 50, 0),
             checkmid = ifelse(between({{index}}, med_midsdlo, med_midsdhi), 5, 0),
             barhi = med_midsdhi,
             barlo = med_midsdlo,
             barmid = med_sum$avg_idx[1]) %>%
      rowwise() %>%
      mutate(pack = checkhi + checklo + checkmid) %>%
      ungroup() %>%
      select(olfrname, side, VenDor, p50, {{index}}, checkhi, checklo, checkmid, barhi, barlo, barmid, pack) %>%
      filter(pack > 0)
    
    both_latmed <- bind_rows(filtered_lat, filtered_med)
    
    #combine
    if (i == max(filter_preds$VenDor)) {
      filtered_all <- bind_rows(filtered_all, both_latmed) %>%
        mutate(type = ifelse(pack == 500, "ahi",
                             ifelse(pack == 50, "zlo",
                                    ifelse(pack == 5, "mid", "woops"))))
      return(filtered_all)
    } else if (i == min(filter_preds$VenDor)) {
      filtered_all <- both_latmed
    } else {
      filtered_all <- bind_rows(filtered_all, both_latmed)
    }#endif
  } #endfor
} #endfunc

#outliers of VD by bin is now outliers_tzbin_vd until close

outliers <- Find_DVoutliers(tzsimplest, 1)

Analyze_DVoutliers <- function(df, index, chooseOut = "Lat") {
  if(str_detect(tolower(chooseOut), "lat")) {
    lat_indiv <- df %>%
      filter(side == "Lateral") %>% 
      ggplot() + 
      geom_point(aes(VenDor, {{index}}, color = type), alpha = 0.5) +
      geom_point(aes(VenDor, barhi), color = "black", shape = 4) +
      geom_point(aes(VenDor, barlo), color = "black", shape = 4) +
      geom_point(aes(VenDor, barmid), color = "black", shape = 4) + 
      ggtitle("Lateral predictions") +
      theme(legend.position = "none")
    
    lat_groups <- df %>%
      filter(side == "Lateral") %>% 
      group_by(VenDor, type) %>% 
      summarise(avgp50 = mean(p50),
                group_size = n()) %>% 
      ggplot() + 
      geom_point(aes(type, avgp50, size = group_size)) +
      facet_wrap(~ VenDor) +
      theme(legend.position = "none")
    
    lat_plots <- lat_indiv + lat_groups
    return(lat_plots)
  } else if (str_detect(tolower(chooseOut), "med")) {
    med_indiv <- df %>%
      filter(side == "Medial") %>% 
      ggplot() + 
      geom_point(aes(VenDor, {{index}}, color = type), alpha = 0.5) +
      geom_point(aes(VenDor, barhi), color = "black", shape = 4) +
      geom_point(aes(VenDor, barlo), color = "black", shape = 4) +
      geom_point(aes(VenDor, barmid), color = "black", shape = 4) + 
      ggtitle("Medial predictions") +
      theme(legend.position = "none")
    
    med_groups <- df %>%
      filter(side == "Medial") %>% 
      group_by(VenDor, type) %>% 
      summarise(avgp50 = mean(p50),
                group_size = n()) %>% 
      ggplot() + 
      geom_point(aes(type, avgp50, size = group_size)) +
      facet_wrap(~ VenDor) +
      theme(legend.position = "none")
    
    med_plots <- med_indiv + med_groups
    return(med_plots)
  } else {
    df_out <- df %>%
      group_by(VenDor, side, type) %>%
      mutate(avgp50 = mean(p50),
             group_size = n()) %>%
      ungroup() %>%
      mutate()
    return(df_out)
  }
}

Analyze_DVoutliers(outliers, tzsimplest, "med")
```

    ## `summarise()` regrouping output by 'VenDor' (override with `.groups` argument)

![](3d_obmapping_analysis_files/figure-gfm/unnamed-chunk-7-1.png)<!-- -->

``` r
Analyze_DVoutliers(outliers, tzsimplest, "lat")
```

    ## `summarise()` regrouping output by 'VenDor' (override with `.groups` argument)

![](3d_obmapping_analysis_files/figure-gfm/unnamed-chunk-7-2.png)<!-- -->

# DV indexes from Luis 3D OE project

``` r
#olfr1204 had a dpt of 0, changed to 1
ls_idx <- read_csv("~/Desktop/obmap/r_analysis/3dimOB/input/LS_3Dindexes_real_pred.csv") %>%
  mutate(logDPT = log2(DPT_index + 0.1),
         rankDPT = min_rank(DPT_index))
```

    ## 
    ## ── Column specification ────────────────────────────────────────────────────────
    ## cols(
    ##   olfrname = col_character(),
    ##   Rfpred_index = col_double(),
    ##   DPT_index = col_double(),
    ##   MiyReal_index = col_double(),
    ##   Tan_index = col_character(),
    ##   Miy_index = col_double(),
    ##   Zolfr_real = col_character(),
    ##   Zolfr_Noreal = col_double(),
    ##   Zolfr_Momb = col_double()
    ## )

``` r
filter_preds <- filter_preds %>%
  left_join(ls_idx, by = "olfrname")

Plot_predictions(filter_olfrs, varcolor = ~logDPT)
```

![](3d_obmapping_analysis_files/figure-gfm/unnamed-chunk-8-1.png)<!-- -->

``` r
test <- filter_preds %>%
  filter(p50 == clustmaxp) %>%
  filter(side == "Medial") %>%
  filter(tzbins == 2) %>%
  mutate(vdrank = min_rank(desc(VenDor)),
         pack = ifelse(vdrank <= 10, "dorsal", 
                       ifelse(vdrank > (max(vdrank) - 10), "ventral", 
                              ifelse(between(vdrank, (mean(vdrank) - 5), (mean(vdrank) + 5)),
                                     "mid", NA)))) %>%
  select(olfrname, side, VenDor, tzsimple, tzbins, p50, pack) %>%
  filter(!is.na(pack))
```

# Topics from Luis 3D OE project

``` r
#Mayra/Antonio/Luis topics
topics <- read_csv("~/Desktop/obmap/r_analysis/3dimOB/input/LS_DegreesOfBelonging_201015.csv") 
```

    ## 
    ## ── Column specification ────────────────────────────────────────────────────────
    ## cols(
    ##   olfrname = col_character(),
    ##   T1 = col_double(),
    ##   T2 = col_double(),
    ##   T3 = col_double(),
    ##   T4 = col_double(),
    ##   T5 = col_double(),
    ##   max_topic = col_character()
    ## )

``` r
# topic_max <- topics_raw %>%
#   rowwise() %>%
#   arrange(desc(T5)) %>%
#   mutate(maxval = pmax(T1, T2, T3, T4, T5)) %>%
#   pivot_longer(T1:T5, names_to = "topic", values_to = "prob") %>%
#   rowwise() %>%
#   mutate(max_topic = ifelse(maxval == prob, topic, NA)) %>%
#   select(olfrname, max_topic) %>%
#   filter(!is.na(max_topic)) %>%
#   ungroup() 
# 
# topics_join <- left_join(topics, topic_max, by = "olfrname")
# write_csv(topics_join, "~/Desktop/obmap/r_analysis/3dimOB/input/LS_DegreesOfBelonging_201015.csv")

filter_preds <- filter_preds %>% 
  left_join(topics, by = "olfrname") %>%
  filter(!is.na(max_topic))

topic_ors <- filter_preds$olfrname

Plot_predictions(topic_ors, varcolor = ~max_topic)
```

![](3d_obmapping_analysis_files/figure-gfm/unnamed-chunk-9-1.png)<!-- -->

``` r
topic_X <- filter_preds %>% filter(max_topic == "T2") %>%
  select(olfrname) %>% as_vector()
Plot_predictions(topic_X, varcolor = ~max_topic)
```

![](3d_obmapping_analysis_files/figure-gfm/unnamed-chunk-9-2.png)<!-- -->

``` r
#1 is all over
#2 is dorsal, class 1 enriched (67 class 1 ORs, 99 class 2 ORs)
#3 is ventral
#4 is only 1 Olfr338, dorsal
#5 is posterior
```

# Find positions for 50 ORs enriched in Wachowiak Functional Imaging surface samples

Shawn Burton provided samples from 8 OBs from 4 mice, each OB cut into 2
pieces with 1 piece representing the functional imaging surface. Using
the above algorithm for picking the best medial and lateral cluster for
a given OR. Need to update list for newest alignment updates and to
include all FI ORs now that code was improved in terms of
speed.

``` r
olfr_result <- read_csv("~/Desktop/rproj/obmap_inactive/wach_diffe/starrsem_aligned/out/wach_v16model_top25FIenriched.csv") %>% select(-Used) %>% mutate(FIsur = ifelse(logFC > 1, ifelse(FDR < 0.05, 1, 0), 0))

func_sig <- olfr_result %>% 
  filter(FIsur == 1) %>% 
  filter(!str_detect(Gene_name, "-ps")) %>% 
  arrange(FDR)
func_sig_olfr <- func_sig$Gene_name
func_sig_olfr
```

    ##  [1] "Olfr376"  "Olfr629"  "Olfr57"   "Olfr1046" "Olfr1122" "Olfr937" 
    ##  [7] "Olfr994"  "Olfr228"  "Olfr570"  "Olfr558"  "Olfr969"  "Olfr402" 
    ## [13] "Olfr683"  "Olfr597"  "Olfr19"   "Olfr578"  "Olfr550"  "Olfr566" 
    ## [19] "Olfr1496" "Olfr974"  "Olfr478"  "Olfr677"  "Olfr147"  "Olfr561" 
    ## [25] "Olfr957"  "Olfr971"  "Olfr633"  "Olfr231"  "Olfr1032" "Olfr922" 
    ## [31] "Olfr1023" "Olfr635"  "Olfr1031" "Olfr690"  "Olfr1010" "Olfr1019"
    ## [37] "Olfr609"  "Olfr874"  "Olfr510"  "Olfr1134" "Olfr160"  "Olfr197" 
    ## [43] "Olfr467"  "Olfr150"  "Olfr935"  "Olfr64"   "Olfr5"    "Olfr506" 
    ## [49] "Olfr1377" "Olfr338"  "Olfr1086" "Olfr1020" "Olfr202"  "Olfr1339"
    ## [55] "Olfr225"  "Olfr691"  "Olfr895"  "Olfr20"   "Olfr982"  "Olfr146" 
    ## [61] "Olfr1328" "Olfr1129" "Olfr432"  "Olfr152"  "Olfr1449" "Olfr490" 
    ## [67] "Olfr557"  "Olfr51"   "Olfr488"  "Olfr1448" "Olfr1154" "Olfr1128"
    ## [73] "Olfr1511" "Olfr881"  "Olfr133"

``` r
ggplot(olfr_result) + 
  geom_point(aes(logFC,-log10(FDR), alpha = 0.25, color = as.factor(FIsur), size = 1.3)) +
  geom_vline(xintercept = 0) + 
  geom_hline(yintercept = -log10(0.05)) + 
  theme_cowplot() + 
  theme(legend.position = "none") + 
  xlab("nonFIsurface  <<<  log2FoldChange  >>>  FIsurface")
```

![](3d_obmapping_analysis_files/figure-gfm/unnamed-chunk-10-1.png)<!-- -->

# Plot significantly enriched FI surface ORs

Best in this case refers to highest FDR (aka how consistency enriched in
functional imaging surface). Perhaps color by an adjusted FDR?

``` r
#plot only highest p50 voxel of each cluster
Plot_predictions(func_sig_olfr, varcolor=~p50, chooseOut = "point", title = "Medial/Lateral cluster points for the top 30 FI surface enriched ORs")
```

![](3d_obmapping_analysis_files/figure-gfm/unnamed-chunk-11-1.png)<!-- -->

# plot heatmap peaks and calc dist to DorML

``` r
heatmap_peaks <- read_csv("~/Desktop/obmap/r_analysis/3dimOB/input/heatmap_peaks.csv")

Heat3D("Olfr881", dimrep = 1, dv = "oe", raw = F, chooseOut = "plot")

all_heat_in_good <- heatmap_peaks$olfrname[which(heatmap_peaks$olfrname %in% filter_olfrs)]
#ahig_dist <- DistHeat3D(all_heat_in_good)
#write_csv(ahig_dist, "~/Desktop/rproj/obmap/allmice/v21_gen25/heatmapORs_distance_to3D.csv")
ahig_dist <- read_csv("~/Desktop/obmap/r_analysis/3dimOB/output/heatmapORs_distance_to3D.csv")
ahig_dist %>% ggplot(aes(side, distance)) + geom_violin()
```

## Check if new Dorsal constraint function improves output

New function checks if OR has evidence of dorsal OE expression (class 1
or miyamichi zone \< 2, matsunami diffE = dorsal). If so, check if
either the Medial or Lateral predicted position is dorsal. If OR is
likely dorsal and either Medial or Lateral glomerulus prediction is
dorsal but other halfbulb glom is not dorsal, find a dorsal glom for
that halfbulb In independent images for the more lateral glomerulus,
1377 seems slightly more anterior lateral than 881