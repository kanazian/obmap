---
title: "OBmap single cell style analysis"
author: "kanazian"
date: "Sep 24, 2020"
output: 
  html_document: 
    keep_md: true
---
# Goal: 
Perform scRNAseq cluster DiffE style analysis on replicates of single dimension OBmap data to examine gene differential expression between positional clusters and dimensional reduction for positional samples. Following: https://bioconductor.org/packages/release/bioc/vignettes/SC3/inst/doc/SC3.html#de-genes

## Prep tpm data for sc3



#Prep single cell experiment

```r
allchemo_dim <- allchemo_tpm %>% filter(dim == "AntPos")
sample_names <- allchemo_dim$name
covars <- allchemo_dim %>% select(name:dimrep)

allchemo_trim <- allchemo_dim %>% select(-name, -rep, -slice, -dim, -dimrep) %>% t()
colnames(allchemo_trim) <- sample_names
allchemo_trim[1:5,1:5]
```

```
##          M10S01 M10S02  M10S03  M10S04  M10S05
## Olfr299    0.00 277.35    0.00    0.00    0.00
## Olfr109    0.00   0.00   99.66    0.00    0.00
## Olfr281    0.33 602.80 2451.54 1882.25 2612.72
## Olfr1015   5.77 228.19    0.11  340.79    0.00
## Olfr1347   0.00   0.30    0.00    0.00    0.00
```

```r
allchemo_coldata <- covars %>% select(name, slice)
coldata_ac <- data.frame("cell_type1" = allchemo_coldata$slice)
rownames(coldata_ac) <- allchemo_coldata$name

sce <- SingleCellExperiment(
    assays = list(
        counts = as.matrix(allchemo_trim),
        logcounts = log2(as.matrix(allchemo_trim) + 1)
    ), 
    colData = coldata_ac
)

rowData(sce)$feature_symbol <- rownames(sce)
```


#visualize dimreduction of the ~24 AP slices

```r
sce_pca <- runPCA(sce)
plotPCA(sce_pca, colour_by = "cell_type1")
```

![](obmap_scstyle_files/figure-html/unnamed-chunk-2-1.png)<!-- -->

```r
sce_tsne <- runTSNE(sce)
plotTSNE(sce_tsne, colour_by = "cell_type1")
```

![](obmap_scstyle_files/figure-html/unnamed-chunk-2-2.png)<!-- -->

```r
sce_umap <- runUMAP(sce)
plotUMAP(sce_umap, colour_by = "cell_type1")
```

![](obmap_scstyle_files/figure-html/unnamed-chunk-2-3.png)<!-- -->

```r
umap_rd <- reducedDim(sce_umap) %>% as_tibble()
```

```
## Warning: The `x` argument of `as_tibble.matrix()` must have unique column names if `.name_repair` is omitted as of tibble 2.0.0.
## Using compatibility `.name_repair`.
## This warning is displayed once every 8 hours.
## Call `lifecycle::last_warnings()` to see where this warning was generated.
```

```r
umap_tibble <- tibble(name = rownames(reducedDim(sce_umap)), 
                      UMAP_1 = umap_rd$V1, UMAP_2 = umap_rd$V2) %>%
  left_join(covars, by = "name") %>%
  mutate(slice_fct = as_factor(slice))

umapplot <- ggplot(umap_tibble) + 
  geom_point(aes(UMAP_1, UMAP_2, fill = slice),
             size = 6, pch = 21, color = "black") +
  scale_fill_viridis_c() +
  theme_cowplot() +
  ggtitle("UMAP of OR Gene Expression", 
          subtitle = paste0(as.character(dim(allchemo_trim)[2]), 
                            " Anterior-Posterior sections from 6 mice"))
ggsave(umapplot, filename = "~/Desktop/obmap/r_analysis/sc_style/output/ap_umap.png", device = "png")
```

```
## Saving 7 x 5 in image
```

```r
umapplot
```

![](obmap_scstyle_files/figure-html/unnamed-chunk-2-4.png)<!-- -->

```r
sce_umapcomp <- runUMAP(sce, ncomponents = 5)
plotUMAP(sce_umapcomp, ncomponents = 5, colour_by = "cell_type1")
```

![](obmap_scstyle_files/figure-html/unnamed-chunk-2-5.png)<!-- -->


# run sc3
Interactive is very good, but keeping plot code for output

```r
sce_sc3 <- sc3(sce, ks = 3:10, biology = TRUE)
```

```
## Setting SC3 parameters...
```

```
## Warning: 'isSpike' is deprecated.
## See help("Deprecated")
```

```
## Calculating distances between the cells...
```

```
## Performing transformations and calculating eigenvectors...
```

```
## Performing k-means clustering...
```

```
## Calculating consensus matrix...
```

```
## Calculating biology...
```

```r
#sc3_interactive(sce_sc3)
sc3_export_results_xls(sce_sc3, filename = "~/Desktop/obmap/r_analysis/sc_style/output/ap_sc3.xls")

plotPCA(
    sce_sc3, 
    colour_by = "sc3_5_clusters", 
    size_by = "sc3_5_log2_outlier_score")
```

```
## Warning: call 'runPCA' explicitly to compute results
```

![](obmap_scstyle_files/figure-html/unnamed-chunk-3-1.png)<!-- -->

```r
sc3_plot_consensus(
    sce_sc3, k = 5, 
    show_pdata = c(
        "cell_type1", 
        "log10_total_features",
        "sc3_5_clusters", 
        "sc3_5_log2_outlier_score"))
```

```
## Provided columns 'log10_total_features' do not exist in the phenoData table!
```

![](obmap_scstyle_files/figure-html/unnamed-chunk-3-2.png)<!-- -->

```r
sc3_plot_silhouette(sce_sc3, k = 5)
```

![](obmap_scstyle_files/figure-html/unnamed-chunk-3-3.png)<!-- -->

```r
sc3_plot_expression(sce_sc3, k = 5)
```

![](obmap_scstyle_files/figure-html/unnamed-chunk-3-4.png)<!-- -->

```r
sc3_plot_cluster_stability(sce_sc3, k = 5)
```

```
## Warning: Use of `d$Cluster` is discouraged. Use `Cluster` instead.
```

```
## Warning: Use of `d$Stability` is discouraged. Use `Stability` instead.
```

![](obmap_scstyle_files/figure-html/unnamed-chunk-3-5.png)<!-- -->

```r
sc3_plot_de_genes(
    sce_sc3, k = 5, 
    show_pdata = c(
        "cell_type1", 
        "log10_total_features",
        "sc3_5_clusters", 
        "sc3_5_log2_outlier_score"))
```

```
## Provided columns 'log10_total_features' do not exist in the phenoData table!
```

![](obmap_scstyle_files/figure-html/unnamed-chunk-3-6.png)<!-- -->

```r
sc3_plot_markers(
    sce_sc3, k = 5, 
    show_pdata = c(
        "cell_type1", 
        "log10_total_features",
        "sc3_3_clusters", 
        "sc3_3_log2_outlier_score"))
```

```
## Provided columns 'log10_total_features' do not exist in the phenoData table!
```

![](obmap_scstyle_files/figure-html/unnamed-chunk-3-7.png)<!-- -->
