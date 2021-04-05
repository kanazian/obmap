Extracting OB shape from MRI model
================
kanazian
2021-04-04

# Goal:

Generate 3D voxel coordinates representing the outer surface (roughly
the glomerular layer) for incorporation as a scaffold for the
statistical reconstruction of 3D OR glomeruli positions from single
dimension targeted transcriptomics data.

## Prior Steps:

1.  Downloaded MRI model (.nii) from Duke CIVM voxport, Microscopic DTI
    of mouse brain (PMC3085633)
2.  Converted .nii to .stl (common 3d print format) following this
    guide: <https://layerfmri.com/2018/07/25/3d-printing-nii-data/>
3.  Loaded .stl into Blender
4.  Created an array of cubes matching the \# of sections taken for
    targeted transcriptomics of each dimension, 24 Anterior-Posterior,
    22 Medial-Lateral, 23 Ventral-Dorsal. This cube array was stretched
    in each dimension to encompass the entire OB due to difference in
    age between MRI study (p66-78) and our (p21) mice. OB was tilted
    slightly in order to better represent the less erect and more
    elongated nature of the OB after the OE is removed in our dissection
5.  Output OB and cubegrid as .ply files

## Plan:

3D models are comprised of polygons. The corner points of these polygons
can be found as 3D coordinates in the .ply files. By finding the
coordinates of the corners of each cube in the grid, I can determine
which cubes contain OB polygon corner coordinates. These cubes are the
ones containing surface of the OB where OR glomeruli are located.

## Potential problems:

The model includes both OBs with no gap between them so there is
technically no medial surface to model, can correct using the most
dorsal and ventral medial coordinates and filling in.

There may be overly large polygons. MRI scans are taken in sections
while the OB is a curved surface so even the best model is kind of like
building a ball with legos. This may lead to overly large rectangles
whose edges may span multiple cubes leading to corners not being present
in certain cubes which do contain OB surface. Can correct by adding
neighbors to initial set of surface containing cubes in order to fill in
gaps where surface should be.

# Notes:

As this analysis involves a lot of 3D models, I use a lot of plotly in
order to have interactive 3D plots which does not display in the github
document (and it doesn’t seem there is an easy way to do so). Until that
time, if you are interested in seeing the 3D outputs, please download
the .html files where you can freely rotate and zoom the 3D models.

### load packages

``` r
knitr::opts_chunk$set(warning=F)
library(geomorph)
library(plotly)
library(tidyverse)
library(knitr)

OBplot <- function(df) {
  obplot <- plot_ly(df, x=~MedLat, y=~AntPos, z=~VenDor, 
                    color=~type, 
                    marker=list(size = 6, 
                                line = list(color = 'black', width = 0.5)), 
                    text=~paste("AP:", AntPos,
                                "<br>ML:", MedLat,
                                "<br>VD:", VenDor,
                                "<br>Type:", type), 
                    type='scatter3d', mode='markers') 
  return(obplot)
} #end OBplot


#Find_Misses is a function that will find voxels that don't contain OB 3D model vertices but contain surface. Can be used to progressively fill in misses.
#hits is a dataframe with cube_number, AP, ML, DV coordinates, and additional info pertaining to the presence of vertices within a voxel. Coordinate names are AntPos, MedLat, VenDor
#all voxels is a expand.grid voxel framework for all possible AP, ML, DV coordinates. Coordinate names are: cube_ap, cube_ml, cube_dv
#output is a single dataframe with AntPos, MedLat, VenDor coordinates as well as a new variable: type, which indicates whether the point was originally a hit or miss

FindMisses <- function(hits, all_voxels) {
  #remove coordinates with hits
  nohit_which <- vector("double", length = dim(all_voxels)[1])
  for (cube in 1:dim(all_voxels)[1]) {
  nohit_which[cube] <- length(which(hits$AntPos == all_voxels$cube_ap[cube] & 
                                    hits$MedLat == all_voxels$cube_ml[cube] &
                                    hits$VenDor == all_voxels$cube_dv[cube]))
  } #endfor
  nohits <- all_voxels[which(nohit_which == 0),] %>% as_tibble()
  
  #vectors to hold neighboring hits
  neighs_ap <- vector("double", length = dim(nohits)[1])
  neighs_ml <- vector("double", length = dim(nohits)[1])
  neighs_vd <- vector("double", length = dim(nohits)[1])
  
  #find non-hit voxels with two neighbors in a single cardinal direction
  for (cube in 1:dim(nohits)[1]) {
  neighs_ap[cube] <- length(which(between(hits$AntPos, nohits$cube_ap[cube]-1, nohits$cube_ap[cube]+1) &
        hits$AntPos != nohits$cube_ap[cube] &
        hits$MedLat == nohits$cube_ml[cube] &
        hits$VenDor == nohits$cube_dv[cube]))
  
  neighs_ml[cube] <- length(which(between(hits$MedLat, nohits$cube_ml[cube]-1, nohits$cube_ml[cube]+1) &
        hits$MedLat != nohits$cube_ml[cube] &
        hits$AntPos == nohits$cube_ap[cube] &
        hits$VenDor == nohits$cube_dv[cube]))
    
  neighs_vd[cube] <- length(which(between(hits$VenDor, nohits$cube_dv[cube]-1, nohits$cube_dv[cube]+1) &
        hits$VenDor != nohits$cube_dv[cube] &
        hits$MedLat == nohits$cube_ml[cube] &
        hits$AntPos == nohits$cube_ap[cube]))
  } #endfor
  
  neigh_grid <- as_tibble(cbind(nohits, neighs_ap, neighs_ml, neighs_vd))
  
  two_neighs <- filter(neigh_grid, neighs_ap == 2 | neighs_ml == 2 | neighs_vd == 2)
  
  simpleHits <- hits %>% select(AntPos:VenDor) %>% mutate(type = "Hit")
  simpleMiss <- two_neighs %>% rename("AntPos" = cube_ap, "MedLat" = cube_ml, "VenDor" = cube_dv) %>% select(AntPos:VenDor) %>% mutate(type = "Missed")
  
  hitandmiss <- bind_rows(simpleHits, simpleMiss)
} #end FindMisses

#df is a tibble containing at least three columns named AntPos, MedLat, and VenDor whose values exceed desired limits
#ap, ml, vd are integers indicating the desired length for that dimension
Rescale <- function(df, ap, ml, vd) {
  apscale <- ap/max(df$AntPos)
  mlscale <- ml/max(df$MedLat)
  vdscale <- vd/max(df$VenDor)
  rescaled <- df %>% 
    select(AntPos, MedLat, VenDor) %>% 
    mutate(newAP = AntPos * apscale,
           newML = MedLat * mlscale,
           newVD = VenDor * vdscale,
           APint = round(newAP),
           MLint = round(newML),
           VDint = round(newVD),
           APd2i = round(newAP) - newAP,
           MLd2i = round(newML) - newML,
           VDd2i = round(newVD) - newVD)
  return(rescaled)
} #end Rescale
```

## Blender view of OB in cube grid

Trimmed OB in half due to crashes

``` r
include_graphics(".3dbraincubes.png") 
```

<img src=".3dbraincubes.png" width="75%" height="75%" />

``` r
#image needs to have same path as .Rmd file
```

### Using geomorph to read .ply 3D data files

``` r
#,showSpecimen=T if you want to view in rgl
cubes <- read.ply("~/Desktop/obmap/r_analysis/data/mri_to_R/v1_cube.ply") 
brain <- read.ply("~/Desktop/obmap/r_analysis/data/mri_to_R/v2_partial_brain.ply")

cubeverts <- as_tibble(t(cubes$vb)) %>% mutate(position = 1:n())
brainverts <- as_tibble(t(brain$vb)) %>% mutate(position = 1:n())

saveRDS(cubeverts, "~/Desktop/obmap/r_analysis/data/mri_to_R/cubeverts.RDS")
saveRDS(brainverts, "~/Desktop/obmap/r_analysis/data/mri_to_R/brainverts.RDS")
```

### load vertices (polygon corner coordinates)

``` r
cubeverts <- readRDS("~/Desktop/obmap/r_analysis/data/mri_to_R/cubeverts.RDS")
cubeverts %>% ggplot(aes(xpts, ypts)) + geom_point()
```

![](extracting_mri_OB_files/figure-gfm/unnamed-chunk-3-1.png)<!-- -->

``` r
brainverts <- readRDS("~/Desktop/obmap/r_analysis/data/mri_to_R/brainverts.RDS")
brainverts %>% ggplot(aes(ypts, zpts)) + geom_point()
```

![](extracting_mri_OB_files/figure-gfm/unnamed-chunk-3-2.png)<!-- -->

# Initial investigation and dealing with problems

Examining the vertices of the cube file (simpler than brain) indicates
that each corner of a cube has multiple vertices. Since this file breaks
objects down to triangles, a single cube could have as few as 3 or as
many as 6 vertices. Not quite sure how closely neighboring cubes in the
array are positioned, this file indicates small differences. Goal of
this chunk is to determine which vertices are “real” corners by finding
vertices with relatively large and equally spaced distances between them
in each plane. I also know the number of cubes alongside each dimension
and the number of corner vertices should be equal to that number. Also
note that the original cube is CENTERED at 0,0,0 and not at a vertex

``` r
#X, the medial-lateral dimension
cv_x1 <- cubeverts %>% group_by(xpts) %>% count()

cv_x2 <- vector(mode = "double", length = dim(cv_x1)[1])
for (i in 1:dim(cv_x1)[1]) {
  if (i != 1) {
    cv_x2[i] <- abs(cv_x1$xpts[i] - cv_x1$xpts[i-1]) 
  } else {
    cv_x2[i] <- 0
  } #endif
} #endfor

cv_x3 <- unlist(cv_x2)
cv_x <- add_column(cv_x1, cv_x3) %>% as_tibble() %>% mutate(axisorder = 1:n())


#Y, the anterior-posterior dimension
cv_y1 <- cubeverts %>% group_by(ypts) %>% count()

cv_y2 <- vector(mode = "double", length = dim(cv_y1)[1])
for (i in 1:dim(cv_y1)[1]) {
  if (i != 1) {
    cv_y2[i] <- abs(cv_y1$ypts[i] - cv_y1$ypts[i-1]) 
  } else {
    cv_y2[i] <- 0
  } #endif
} #endfor

cv_y3 <- unlist(cv_y2)
cv_y <- add_column(cv_y1, cv_y3) %>% as_tibble() %>% mutate(axisorder = 1:n())


#Z, the ventral-dorsal dimension
cv_z1 <- cubeverts %>% group_by(zpts) %>% count()

cv_z2 <- vector(mode = "double", length = dim(cv_z1)[1])
for (i in 1:dim(cv_z1)[1]) {
  if (i != 1) {
    cv_z2[i] <- abs(cv_z1$zpts[i] - cv_z1$zpts[i-1]) 
  } else {
    cv_z2[i] <- 0
  } #endif
} #endfor

cv_z3 <- unlist(cv_z2)
cv_z <- add_column(cv_z1, cv_z3) %>% as_tibble() %>% mutate(axisorder = 1:n())

#distance between cube edges should be similar to:
(max(cv_y$ypts)-min(cv_y$ypts))/27
```

    ## [1] 0.04412704

``` r
#seems that in cv_y, the 1st and 8th coordinates, 9th and 15th, 16th and 22nd, etc. are vertices that have a similar distance to the cube size for that dimension as well as follow a pattern of aabcbaa|aabcbaa|aabcb....cbaa
#for vertice coordinate arranged from min to max, seq(1,length(x),7) and seq(9,length(x),7) will produce list positions that are "true" corners
trueCornersInOrder <- c(seq(1, 190, 7), seq(9,190,7))
```

# Make df to hold cube grid

Colnames: cube.coord.ml, cube.coord.ap, cube.coord.dv, minML, maxML,
minAP, maxAP, minDV, maxDV

``` r
#Medial Lateral is X dim in 3D file, low values are more lateral, higher more medial. In OBmap low values are more medial, higher more lateral. Hence need to flip.
mutML <- 24
mutAP <- 27
mutDV <- 23
true_ml <- cv_x %>% filter(axisorder %in% trueCornersInOrder) %>% mutate(cube_ml = rep(seq(mutML,1), each = 2), lim_ml = rep(c("minML", "maxML"),mutML)) %>% select(xpts, cube_ml, lim_ml) %>% spread(lim_ml, xpts)

#Anterior Posterior is Y dim in 3D file, low Y is more posterior, higher more anterior. In OBmap low values are anterior, higher more posterior. Hence need to flip.
true_ap <- cv_y %>% filter(axisorder %in% trueCornersInOrder) %>% mutate(cube_ap = rep(seq(mutAP,1), each = 2), lim_ap = rep(c("minAP", "maxAP"),mutAP)) %>% select(ypts, cube_ap, lim_ap) %>% spread(lim_ap, ypts)

#Dorsal Ventral is Z dim in 3D file, low Z is more ventral, higher more dorsal. In OBmap low values are ventral, higher more dorsal. Hence no need to flip.
true_dv <- cv_z %>% filter(axisorder %in% trueCornersInOrder) %>% mutate(cube_dv = rep(seq(mutDV), each = 2), lim_dv = rep(c("minDV", "maxDV"),mutDV)) %>% select(zpts, cube_dv, lim_dv) %>% spread(lim_dv, zpts)


#make a coordinate grid like in OBMap, however note that the the numbering will be off since there are more cubes in this model (27x24x23) compared to OBMap (24,23,22)
cube_grid <- expand.grid(1:mutAP, 1:mutML, 1:mutDV)
colnames(cube_grid) <- c("cube_ap", "cube_ml", "cube_dv")

grid_coords <- as_tibble(cube_grid) %>% left_join(true_ap, by = "cube_ap") %>% left_join(true_ml, by = "cube_ml") %>% left_join(true_dv, by = "cube_dv") %>% mutate(cube = 1:n()) %>% select(cube, everything()) %>% dplyr::rename("AntPos" = cube_ap, "MedLat" = cube_ml, "VenDor" = cube_dv)
```

### Using vertices from the brain file, check which vertices are within the global cube grid limits. I will use only these vertices for assignments into the cube grid

``` r
xmax <- max(cv_x1$xpts)
xmin <- min(cv_x1$xpts)
ymax <- max(cv_y1$ypts)
ymin <- min(cv_y1$ypts)
zmax <- max(cv_z1$zpts)
zmin <- min(cv_z1$zpts)

betweenX <- which(between(brainverts$xpts, xmin, xmax)) #1,227,597
betweenY <- which(between(brainverts$ypts, ymin, ymax)) #1,113,114
betweenZ <- which(between(brainverts$zpts, zmin, zmax)) #1,159,631
betweenXY <- intersect(betweenX, betweenY) #541,372
betweenXYZ <- intersect(betweenXY, betweenZ) #298,953

brainV_in_cubes <- brainverts[betweenXYZ,]

brainV_in_cubes %>% ggplot(aes(xpts, ypts)) + geom_point()
```

![](extracting_mri_OB_files/figure-gfm/unnamed-chunk-6-1.png)<!-- -->

## Found \~300k vertices that fall within the grid of cubes

Now lets find exactly which cubes have vertices within them.

``` r
hits <- vector("double", length = dim(grid_coords)[1])
for (cube in 1:dim(grid_coords)[1]) {
  verts_in_cube <- length(which(between(brainV_in_cubes$xpts, 
                                        grid_coords$minML[cube], 
                                        grid_coords$maxML[cube]) &
                                  between(brainV_in_cubes$ypts, 
                                          grid_coords$minAP[cube], 
                                          grid_coords$maxAP[cube]) &
                                  between(brainV_in_cubes$zpts, 
                                          grid_coords$minDV[cube], 
                                          grid_coords$maxDV[cube])))
  
  hits[cube] <- verts_in_cube
}
```

## Plot the cubes with vertices

``` r
grid_hits <- grid_coords %>% mutate(hitcount = hits, hitlog = ifelse(hits > 0, TRUE, FALSE))
only_hits <- filter(grid_hits, hitlog == TRUE)

saveRDS(only_hits, "~/Desktop/obmap/r_analysis/data/mri_to_R/hits_within_cube.RDS")

#a little bit of optical nerve seems to have entered the cube grid
weirdbit <- only_hits %>% filter(AntPos > 23 & MedLat > 19 & VenDor < 3) #14 points
getWeird <- only_hits %>% mutate(isWeird = ifelse(AntPos > 23 & MedLat > 19 & VenDor < 3, TRUE, FALSE)) %>% filter(isWeird == FALSE) %>% select(-isWeird) #now missing 14 points

plot_hits <- plot_ly(getWeird, x=~MedLat, y=~AntPos, z=~VenDor, marker=list(size = 6, line = list(color = 'black', width = 0.5)), text=~paste("Hits:", hitcount), type='scatter3d', mode='markers')
plot_hits
```

    ## Fontconfig warning: "/etc/fonts/conf.avail/53-monospace-lcd-filter.conf", line 10: Having multiple values in <test> isn't supported and may not work as expected
    ## Fontconfig warning: "/etc/fonts/conf.avail/53-monospace-lcd-filter.conf", line 10: Having multiple values in <test> isn't supported and may not work as expected

![](extracting_mri_OB_files/figure-gfm/unnamed-chunk-8-1.png)<!-- -->

## Remove cubes with very few hits prior to filling missing spaces

``` r
thresh <- 25
stripper <- getWeird %>%
  mutate(type = case_when(AntPos < 22 & hitcount < thresh ~ "Low",
                          AntPos > 24 & hitcount < thresh ~ "Low",
                          MedLat < 9 & hitcount < thresh ~ "Low",
                          MedLat > 13 & hitcount < thresh ~ "Low",
                          VenDor < 20 & hitcount < thresh ~ "Low",
                          TRUE ~ "High")) #for strip thresh25

OBplot(stripper)
```

    ## Fontconfig warning: "/etc/fonts/conf.avail/53-monospace-lcd-filter.conf", line 10: Having multiple values in <test> isn't supported and may not work as expected
    ## Fontconfig warning: "/etc/fonts/conf.avail/53-monospace-lcd-filter.conf", line 10: Having multiple values in <test> isn't supported and may not work as expected

![](extracting_mri_OB_files/figure-gfm/unnamed-chunk-9-1.png)<!-- -->

``` r
stripped <- stripper %>%
  filter(type == "High") %>%
  select(-type)
```

## Some cubes appear to be missing

This is likely due to relatively flat faces on the object leading to
larger triangles leading to more dispersed vertices. Will try to correct
this by examining neighbors in the cardinal directions and filling spots
that have 2 cardinal neighbors. Manual addition will be done after if
needed.

``` r
findMissOut_1 <- FindMisses(stripped, cube_grid)

#plot round 1  
OBplot(findMissOut_1)
```

    ## Fontconfig warning: "/etc/fonts/conf.avail/53-monospace-lcd-filter.conf", line 10: Having multiple values in <test> isn't supported and may not work as expected
    ## Fontconfig warning: "/etc/fonts/conf.avail/53-monospace-lcd-filter.conf", line 10: Having multiple values in <test> isn't supported and may not work as expected

![](extracting_mri_OB_files/figure-gfm/unnamed-chunk-10-1.png)<!-- -->

### Find Misses Round 2: The Sequel

``` r
findMissOut_2 <- FindMisses(findMissOut_1, cube_grid)
OBplot(findMissOut_2)
```

    ## Fontconfig warning: "/etc/fonts/conf.avail/53-monospace-lcd-filter.conf", line 10: Having multiple values in <test> isn't supported and may not work as expected
    ## Fontconfig warning: "/etc/fonts/conf.avail/53-monospace-lcd-filter.conf", line 10: Having multiple values in <test> isn't supported and may not work as expected

![](extracting_mri_OB_files/figure-gfm/unnamed-chunk-11-1.png)<!-- -->

``` r
#adding points that DONT have 2 cardinal neighbors aka wont be found in Find_Misses
findMissOut_3 <- findMissOut_2 %>% 
  add_row(AntPos = 22, MedLat = 3, VenDor = 21, type="add") %>% #for no stripping
  add_row(AntPos = 1, MedLat = 6, VenDor = 16, type="add") %>% #for strip thresh25
  add_row(AntPos = 1, MedLat = 7, VenDor = 16, type="add") %>% #for strip thresh25
  add_row(AntPos = 1, MedLat = 6, VenDor = 15, type="add") %>% #for strip thresh25
  add_row(AntPos = 1, MedLat = 7, VenDor = 15, type="add") %>% #for strip thresh25
  add_row(AntPos = 19, MedLat = 22, VenDor = 11, type="add") %>% #for strip thresh25
  add_row(AntPos = 20, MedLat = 22, VenDor = 11, type="add") %>% #for strip thresh25
  add_row(AntPos = 18, MedLat = 4, VenDor = 1, type="add") %>% #for strip thresh25
  add_row(AntPos = 18, MedLat = 3, VenDor = 1, type="add") %>% #for strip thresh25
  add_row(AntPos = 18, MedLat = 2, VenDor = 1, type="add") %>% #for strip thresh25
  add_row(AntPos = 18, MedLat = 1, VenDor = 1, type="add") %>% #for strip thresh25
  add_row(AntPos = 19, MedLat = 1, VenDor = 1, type="add")#for strip thresh25
OBplot(findMissOut_3)
```

    ## Fontconfig warning: "/etc/fonts/conf.avail/53-monospace-lcd-filter.conf", line 10: Having multiple values in <test> isn't supported and may not work as expected
    ## Fontconfig warning: "/etc/fonts/conf.avail/53-monospace-lcd-filter.conf", line 10: Having multiple values in <test> isn't supported and may not work as expected

![](extracting_mri_OB_files/figure-gfm/unnamed-chunk-11-2.png)<!-- -->

# Declare a national emergency and build a border wall

Originally written in 2018

``` r
most_medial <- findMissOut_3 %>% filter(MedLat == 1)
#for each AP coord, if more than 2 points in AP coord, add a bunch of vd points to fill therange
border_wall <- expand.grid(0,0,0)

for (i in seq(range(most_medial$AntPos)[1],range(most_medial$AntPos)[2])) {
  med_ap_count <- most_medial %>% filter(AntPos == i)
  
  if (dim(med_ap_count)[1] >= 2) {
    min_vd <- min(med_ap_count$VenDor)
    max_vd <- max(med_ap_count$VenDor)
    build_the_wall <- expand.grid(med_ap_count$AntPos[1],med_ap_count$MedLat[1],min_vd:max_vd)
  }
  
  border_wall <- bind_rows(border_wall, build_the_wall)
}

border_wall <- as_tibble(border_wall) %>% rename("AntPos" = Var1, "MedLat" = Var2, "VenDor" = Var3) %>% mutate(type = "Wall") %>% filter(AntPos > 0)

walled <- bind_rows(findMissOut_3, border_wall)
OBplot(walled)
```

    ## Fontconfig warning: "/etc/fonts/conf.avail/53-monospace-lcd-filter.conf", line 10: Having multiple values in <test> isn't supported and may not work as expected
    ## Fontconfig warning: "/etc/fonts/conf.avail/53-monospace-lcd-filter.conf", line 10: Having multiple values in <test> isn't supported and may not work as expected

![](extracting_mri_OB_files/figure-gfm/unnamed-chunk-12-1.png)<!-- -->

# Add an interior layer

``` r
#define a central coordinate in the middle of the OB, for each hit, note direction toward central coordinate and add a point
#try to avoid multi dimensions, lets instead move along the anterior posterior 
inner_layer <- tibble(AntPos = 0, MedLat = 0, VenDor = 0, type = "init")

for (cube in 1:nrow(walled)) {
  if (walled$AntPos[cube] == 1) {
    #if front wall, add a cube directly posterior
    inner_layer <- inner_layer %>% add_row(AntPos = 2, MedLat = walled$MedLat[cube], VenDor = walled$VenDor[cube], type = "F")
  } else if (walled$MedLat[cube] >= 10 && walled$VenDor[cube] >= 13) {
    #if more Lateral and more Dorsal, add a cube more ventral and more medial
    inner_layer <- inner_layer %>% add_row(AntPos = walled$AntPos[cube], MedLat = walled$MedLat[cube], VenDor = walled$VenDor[cube]-1, type = "A") %>% add_row(AntPos = walled$AntPos[cube], MedLat = walled$MedLat[cube]-1, VenDor = walled$VenDor[cube], type = "A")
  } else if (walled$MedLat[cube] >= 10  &&  walled$VenDor[cube] < 13) {
    #if more Lateral and more Ventral, add ...
    inner_layer <- inner_layer %>% add_row(AntPos = walled$AntPos[cube], MedLat = walled$MedLat[cube], VenDor = walled$VenDor[cube]+1, type = "B") %>% add_row(AntPos = walled$AntPos[cube], MedLat = walled$MedLat[cube]-1, VenDor = walled$VenDor[cube], type = "B")
  } else if (walled$MedLat[cube] < 10  &&  walled$VenDor[cube] >= 13) {
    #if more Medial and more Dorsal, add ...
    inner_layer <- inner_layer %>% add_row(AntPos = walled$AntPos[cube], MedLat = walled$MedLat[cube], VenDor = walled$VenDor[cube]-1, type = "C") %>% add_row(AntPos = walled$AntPos[cube], MedLat = walled$MedLat[cube]+1, VenDor = walled$VenDor[cube], type = "C")
  } else if (walled$MedLat[cube] < 10  &&  walled$VenDor[cube] < 13) {
    #if more Medial and more Ventral, add ...
    inner_layer <- inner_layer %>% add_row(AntPos = walled$AntPos[cube], MedLat = walled$MedLat[cube], VenDor = walled$VenDor[cube]+1, type = "D") %>% add_row(AntPos = walled$AntPos[cube], MedLat = walled$MedLat[cube]+1, VenDor = walled$VenDor[cube], type = "D")
  } else {
    #else make a note
    inner_layer <- inner_layer %>% add_row(AntPos = walled$AntPos[cube], MedLat = walled$MedLat[cube], VenDor = walled$VenDor[cube], type = "E")
  }
}
inner_layer %>% group_by(type) %>% count()
```

    ## # A tibble: 6 x 2
    ## # Groups:   type [6]
    ##   type      n
    ##   <chr> <int>
    ## 1 A      1202
    ## 2 B       964
    ## 3 C      1046
    ## 4 D       890
    ## 5 F        26
    ## 6 init      1

``` r
inners_dups <- inner_layer %>% filter(type != "init")
inner <- inners_dups[-which(duplicated(inners_dups)==TRUE),]

outer <- walled %>% select(-type) %>% mutate(type = "outer") 
like_an_onion <- bind_rows(outer, inner)

OBplot(like_an_onion)
```

    ## Fontconfig warning: "/etc/fonts/conf.avail/53-monospace-lcd-filter.conf", line 10: Having multiple values in <test> isn't supported and may not work as expected
    ## Fontconfig warning: "/etc/fonts/conf.avail/53-monospace-lcd-filter.conf", line 10: Having multiple values in <test> isn't supported and may not work as expected

![](extracting_mri_OB_files/figure-gfm/unnamed-chunk-13-1.png)<!-- -->

``` r
#saveRDS(like_an_onion, "~/Desktop/obmap/r_analysis/mri_to_R/output/190218_outer_inner_coords.RDS")
```

# adjusting for alt models

Scale from full dim down to model dim Find way to add gloms to medial
wall - add concave curve to medial surface

``` r
#take anterior of AP26, remove odd looking strip from Dorsal-Medial edge
trim_to25 <- like_an_onion %>% filter(AntPos <= 25) %>% filter(type == "outer") %>% mutate(oddstrip = ifelse(MedLat == 1, ifelse(VenDor == 21, T, F), F)) %>% filter(oddstrip == F) %>% select(-oddstrip) %>% mutate(origin = "oniontrim25") %>% mutate(type2 = ifelse(VenDor == 1, "outer", type)) %>% select(-type) %>% rename(type = type2) %>% filter(!(AntPos == 20 & MedLat == 2 & VenDor == 22)) %>% filter(!(AntPos == 24 & MedLat == 8 & VenDor == 23))
OBplot(trim_to25)
```

    ## Fontconfig warning: "/etc/fonts/conf.avail/53-monospace-lcd-filter.conf", line 10: Having multiple values in <test> isn't supported and may not work as expected
    ## Fontconfig warning: "/etc/fonts/conf.avail/53-monospace-lcd-filter.conf", line 10: Having multiple values in <test> isn't supported and may not work as expected

![](extracting_mri_OB_files/figure-gfm/unnamed-chunk-14-1.png)<!-- -->

``` r
# actually looks a little too thick 
# check14 <- trim_to25 %>% filter(AntPos == 14)
# plot_ly(data=check14, x=~MedLat, y=~AntPos, z=~VenDor, color=~type, marker=list(size = 6, line = list(color = 'black', width = 0.5)), text=~paste("Type:", type), type='scatter3d', mode='markers')

#fill in posterior bits missing between AP20-25
backbit <- trim_to25 %>% filter(AntPos > 20)
fillpost <- tibble(AntPos = 0, MedLat = 0, VenDor = 0, type = "init")
for (i in 1:nrow(backbit)) {
  positions <- c(backbit$AntPos[i], backbit$MedLat[i], backbit$VenDor[i])
  typei <- backbit$type[i]
  if (positions[3] == 23) {
    next
  } else if (positions[2] == 24) {
    next
  }
  if (backbit %>% filter(AntPos == positions[1] + 1 & 
                         MedLat == positions[2] & 
                         VenDor == positions[3]) %>% nrow() == 0) {
    for (j in positions[1]:25) {
      fillpost <- fillpost %>% add_row(AntPos = j, MedLat = positions[2], VenDor = positions[3], type = typei)
    } #endfor
  } #endif
} #endfor

fillpostout <- fillpost %>% mutate(origin = "fillpost") %>% filter(type != "init")

trimpost <- bind_rows(trim_to25, fillpostout)
OBplot(trimpost)
```

    ## Fontconfig warning: "/etc/fonts/conf.avail/53-monospace-lcd-filter.conf", line 10: Having multiple values in <test> isn't supported and may not work as expected
    ## Fontconfig warning: "/etc/fonts/conf.avail/53-monospace-lcd-filter.conf", line 10: Having multiple values in <test> isn't supported and may not work as expected

![](extracting_mri_OB_files/figure-gfm/unnamed-chunk-14-2.png)<!-- -->

``` r
#AP26 will be projected posterior and interior to AP25
just25 <- trimpost %>% filter(AntPos == 25)

#create an inside ring based on direction to mean ML/VD point of AP section
#more specifically, the medial face will move straight back but the lateral hemisphere will begin curving inward
ring26 <- just25 %>% select(-AntPos) %>% mutate(AntPos = 26) %>% select(AntPos, everything())
meanVD26 <- mean(ring26$VenDor)
meanML26 <- mean(ring26$MedLat)

inside26 <- tibble(AntPos = 0, MedLat = 0, VenDor = 0, type = "init")
#set VD1 ML bounds for AP26
ventral_minML26 <- ring26 %>% filter(VenDor == 2) %>% pull(MedLat) %>% min()
ventral_maxML26 <- ring26 %>% filter(VenDor == 2) %>% pull(MedLat) %>% max()

for (i in 1:nrow(ring26)) {
  positions <- c(ring26$AntPos[i], ring26$MedLat[i], ring26$VenDor[i])
  sign_newVD <- sign(positions[3] - meanVD26) #if positive, new point goes ventral
  sign_newML <- sign(positions[2] - meanML26) #if positive, new point goes medial
  pos_newVD <- ifelse(sign_newVD == 1, positions[3]-1, positions[3]+1)
  #push medial side straight back, but curve lateral side
  pos_newML <- ifelse(positions[2] < 10, 
                      positions[2], ifelse(sign_newML == 1, positions[2]-1, positions[2]+1))
  
  if (positions[3] > 2) {
    inside26 <- inside26 %>% add_row(AntPos = 26, MedLat = pos_newML, VenDor = pos_newVD, type = "outer")
  } else if (positions[3] == 2) {
    inside26 <- inside26 %>% add_row(AntPos = 26, MedLat = pos_newML, VenDor = pos_newVD, type = "outer")
    inside26 <- inside26 %>% add_row(AntPos = 26, MedLat = pos_newML, VenDor = 2, type = "outer")
  }
}
ven26floor <- tibble(AntPos = 26, MedLat = (ventral_minML26 + 1):(ventral_maxML26 - 1), 
                   VenDor = 1, type = "outer")

inside26out <- bind_rows(inside26, ven26floor) %>% mutate(origin = "add26") %>% filter(type != "init")

#do the same for AP27
ring27 <- inside26out %>% select(-AntPos) %>% mutate(AntPos = 27) %>% select(AntPos, everything())
meanVD27 <- mean(ring27$VenDor)
meanML27 <- mean(ring27$MedLat)

inside27 <- tibble(AntPos = 0, MedLat = 0, VenDor = 0, type = "init")
#set VD1 ML bounds for AP27
ventral_minML27 <- ring27 %>% filter(VenDor == 1) %>% pull(MedLat) %>% min()
ventral_maxML27 <- ring27 %>% filter(VenDor == 1) %>% pull(MedLat) %>% max()
for (i in 1:nrow(ring27)) {
  positions <- c(ring27$AntPos[i], ring27$MedLat[i], ring27$VenDor[i])
  sign_newVD <- sign(positions[3] - meanVD27) #if positive, new point goes ventral
  sign_newML <- sign(positions[2] - meanML27) #if positive, new point goes medial
  pos_newVD <- ifelse(sign_newVD == 1, positions[3]-1, positions[3]+1)
  pos_newML <- ifelse(positions[2] < 10, 
                      positions[2], ifelse(sign_newML == 1, positions[2]-1, positions[2]+1))
  if (positions[3] > 2) {
    inside27 <- inside27 %>% add_row(AntPos = 27, MedLat = pos_newML, VenDor = pos_newVD, type = "outer")
  } else if (positions[3] == 2) {
    inside27 <- inside27 %>% add_row(AntPos = 27, MedLat = pos_newML, VenDor = pos_newVD, type = "outer")
    inside27 <- inside27 %>% add_row(AntPos = 27, MedLat = pos_newML, VenDor = 2, type = "outer")
  }
}
ven27floor <- tibble(AntPos = 27, MedLat = (ventral_minML27 + 1):(ventral_maxML27 - 1), 
                   VenDor = 1, type = "outer")

inside27out <- bind_rows(inside27, ven27floor) %>% mutate(origin = "add27") %>% filter(type != "init")

check_expand <- bind_rows(just25, inside26out) %>% bind_rows(inside27out)
OBplot(check_expand)
```

    ## Fontconfig warning: "/etc/fonts/conf.avail/53-monospace-lcd-filter.conf", line 10: Having multiple values in <test> isn't supported and may not work as expected
    ## Fontconfig warning: "/etc/fonts/conf.avail/53-monospace-lcd-filter.conf", line 10: Having multiple values in <test> isn't supported and may not work as expected

![](extracting_mri_OB_files/figure-gfm/unnamed-chunk-14-3.png)<!-- -->

``` r
trim_expand <- bind_rows(trimpost, check_expand)
OBplot(trim_expand)
```

    ## Fontconfig warning: "/etc/fonts/conf.avail/53-monospace-lcd-filter.conf", line 10: Having multiple values in <test> isn't supported and may not work as expected
    ## Fontconfig warning: "/etc/fonts/conf.avail/53-monospace-lcd-filter.conf", line 10: Having multiple values in <test> isn't supported and may not work as expected

![](extracting_mri_OB_files/figure-gfm/unnamed-chunk-14-4.png)<!-- -->

``` r
trim_expand_add <- trim_expand %>%
  add_row(AntPos = 26, MedLat = 10, VenDor = 21, origin = "add", type = "outer") %>%
  add_row(AntPos = 26, MedLat = 1, VenDor = 1, origin = "add", type = "outer") %>%
  add_row(AntPos = 27, MedLat = 2, VenDor = 1, origin = "add", type = "outer") %>%
  add_row(AntPos = 27, MedLat = 1, VenDor = 1, origin = "add", type = "outer")

OBplot(trim_expand_add)
```

    ## Fontconfig warning: "/etc/fonts/conf.avail/53-monospace-lcd-filter.conf", line 10: Having multiple values in <test> isn't supported and may not work as expected
    ## Fontconfig warning: "/etc/fonts/conf.avail/53-monospace-lcd-filter.conf", line 10: Having multiple values in <test> isn't supported and may not work as expected

![](extracting_mri_OB_files/figure-gfm/unnamed-chunk-14-5.png)<!-- -->

``` r
#saveRDS(trim_expand_add, "~/Desktop/obmap/r_analysis/data/mri_to_R/temp_strip25_straightM_outeronly.RDS")

#double the medial face
ml2 <- trim_expand_add %>% filter(MedLat == 1) %>% select(-MedLat) %>% mutate(MedLat = 2)
thicc <- bind_rows(trim_expand_add, ml2)
OBplot(thicc)
```

    ## Fontconfig warning: "/etc/fonts/conf.avail/53-monospace-lcd-filter.conf", line 10: Having multiple values in <test> isn't supported and may not work as expected
    ## Fontconfig warning: "/etc/fonts/conf.avail/53-monospace-lcd-filter.conf", line 10: Having multiple values in <test> isn't supported and may not work as expected

![](extracting_mri_OB_files/figure-gfm/unnamed-chunk-14-6.png)<!-- -->

\#make concave medial wall COULD ALSO REDUCE MEDIAL WALL SIZE BY SLOPING
THE MORE MEDIAL SECTIONS causing the medial wall to begin at a lower
dorsal ceiling and higher ventral floor, would probably remove \~40
voxels or \~12%

``` r
medwallout <- trim_expand %>% filter(MedLat == 1)
medwallin <- trim_expand %>% filter(MedLat == 2)
medwallin_outtype <- medwallin %>% filter(type == "outer")

medwalloutcc <- medwallout %>% mutate(cc = ifelse(between(VenDor, 13, 15), 
                                                  ifelse(between(AntPos, 3, 26),
                                                         T, F), F))
medwallcc <- medwalloutcc %>% filter(cc == T) %>% select(-MedLat) %>% mutate(MedLat = 2)

newmedwallout <- medwalloutcc %>% filter(cc == F)

newmedwall <- bind_rows(newmedwallout, medwallcc)

trimpost_nomed <- trim_expand %>% filter(MedLat >= 3)
trimpost_newmed <- bind_rows(trimpost_nomed, newmedwall) %>% bind_rows(medwallin_outtype)
OBplot(trimpost_newmed)
```

    ## Fontconfig warning: "/etc/fonts/conf.avail/53-monospace-lcd-filter.conf", line 10: Having multiple values in <test> isn't supported and may not work as expected
    ## Fontconfig warning: "/etc/fonts/conf.avail/53-monospace-lcd-filter.conf", line 10: Having multiple values in <test> isn't supported and may not work as expected

![](extracting_mri_OB_files/figure-gfm/unnamed-chunk-15-1.png)<!-- -->

# scale to AP23, ML22, VD23

``` r
#rescale OB with concave medial wall
rescaledcc <- Rescale(trimpost_newmed, 23, 22, 23)
rescaledcc_out <- tibble(AntPos = rescaledcc$APint, 
                         MedLat = rescaledcc$MLint, 
                         VenDor = rescaledcc$VDint, type = "concave_medwall") %>%
  unique()
OBplot(rescaledcc_out)
```

    ## Fontconfig warning: "/etc/fonts/conf.avail/53-monospace-lcd-filter.conf", line 10: Having multiple values in <test> isn't supported and may not work as expected
    ## Fontconfig warning: "/etc/fonts/conf.avail/53-monospace-lcd-filter.conf", line 10: Having multiple values in <test> isn't supported and may not work as expected

![](extracting_mri_OB_files/figure-gfm/unnamed-chunk-16-1.png)<!-- -->

``` r
#saveRDS(rescaledcc_out, "~/Desktop/obmap/r_analysis/mri_to_R/output/210310_cubegrid_rescaled_concaveMedial_outeronly.RDS")

#rescale OB with straight medial wall
rescaled <- Rescale(thicc, 23, 22, 23)
rescaled_out <- tibble(AntPos = rescaled$APint, 
                         MedLat = rescaled$MLint, 
                         VenDor = rescaled$VDint, type = "straight_medwall") %>%
  unique()
OBplot(rescaled_out)
```

    ## Fontconfig warning: "/etc/fonts/conf.avail/53-monospace-lcd-filter.conf", line 10: Having multiple values in <test> isn't supported and may not work as expected
    ## Fontconfig warning: "/etc/fonts/conf.avail/53-monospace-lcd-filter.conf", line 10: Having multiple values in <test> isn't supported and may not work as expected

![](extracting_mri_OB_files/figure-gfm/unnamed-chunk-16-2.png)<!-- -->

``` r
#saveRDS(rescaled_out, "~/Desktop/obmap/r_analysis/data/mri_to_R/210404_voxV4_strip25_straightMdoubled_outeronly_rescaled.RDS")
```