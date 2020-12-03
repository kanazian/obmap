wachR2\_gen25
================
kanazian
Nov 23 2020

Goal: Re-examine the differential expression analysis of Functional
Imaging OB Surface vs non-FIsurface samples from Wachowiak lab (second
set).

``` r
library(edgeR)
library(DESeq)
library(merTools)
library(plotly)
library(tidyverse)
library(ggrepel)
library(cowplot)

#Set working directory here
#setwd("~/Desktop/obmap/r_analysis/func_surface/")

#Differential Expression function
DiffExp <- function (targets, countsTable) {
  Treat <- factor(targets$Treatment);Subject <- factor(targets$Subject);design <- model.matrix(~Subject+Treat)
  cds <- newCountDataSet(countsTable,Treat);cds <- estimateSizeFactors(cds);cds <- estimateDispersions(cds);d <- nbinomTest(cds,"0","1")
  e.litter <- DGEList(counts=countsTable)
  e.litter <- estimateGLMCommonDisp(e.litter,design)
  e.litter <- estimateGLMTrendedDisp(e.litter,design)
  e.litter <- estimateGLMTagwiseDisp(e.litter,design)
  fit <- glmFit(e.litter, design);lrt <- glmLRT(fit);diff <- topTags(lrt,n=dim(lrt)[1])$table
  result <- merge(merge(diff,countsTable,by=0,sort=F),d, by.x="Row.names", by.y="id",sort=F)
  colnames(result)[1] <- "id"
  return(result)
}
```

## Load files and run DiffE

``` r
#group1 - nonFI surface/ventral OB
g1a <- read.table("~/Desktop/obmap/r_analysis/func_surface/sr_gen25_out/1LV_wR2.genes.results", header=T)
g1b <- read.table("~/Desktop/obmap/r_analysis/func_surface/sr_gen25_out/1RV_wR2.genes.results", header=T)
g1c <- read.table("~/Desktop/obmap/r_analysis/func_surface/sr_gen25_out/2LV_wR2.genes.results", header=T)
g1d <- read.table("~/Desktop/obmap/r_analysis/func_surface/sr_gen25_out/2RV_wR2.genes.results", header=T)
g1e <- read.table("~/Desktop/obmap/r_analysis/func_surface/sr_gen25_out/3LV_wR2.genes.results", header=T)
g1f <- read.table("~/Desktop/obmap/r_analysis/func_surface/sr_gen25_out/3RV_wR2.genes.results", header=T)
g1g <- read.table("~/Desktop/obmap/r_analysis/func_surface/sr_gen25_out/4LV_wR2.genes.results", header=T)
g1h <- read.table("~/Desktop/obmap/r_analysis/func_surface/sr_gen25_out/4RV_wR2.genes.results", header=T)

#group2 - FI surface/dorsal OB
g2a <- read.table("~/Desktop/obmap/r_analysis/func_surface/sr_gen25_out/1LD_wR2.genes.results", header=T)
g2b <- read.table("~/Desktop/obmap/r_analysis/func_surface/sr_gen25_out/1RD_wR2.genes.results", header=T)
g2c <- read.table("~/Desktop/obmap/r_analysis/func_surface/sr_gen25_out/2LD_wR2.genes.results", header=T)
g2d <- read.table("~/Desktop/obmap/r_analysis/func_surface/sr_gen25_out/2RD_wR2.genes.results", header=T)
g2e <- read.table("~/Desktop/obmap/r_analysis/func_surface/sr_gen25_out/3LD_wR2.genes.results", header=T)
g2f <- read.table("~/Desktop/obmap/r_analysis/func_surface/sr_gen25_out/3RD_wR2.genes.results", header=T)
g2g <- read.table("~/Desktop/obmap/r_analysis/func_surface/sr_gen25_out/4LD_wR2.genes.results", header=T)
g2h <- read.table("~/Desktop/obmap/r_analysis/func_surface/sr_gen25_out/4RD_wR2.genes.results", header=T)

#make table, takes a couple mins
countsTable <- round(data.frame(one1=g1a$expected_count, 
                                one2=g1b$expected_count, 
                               one3=g1c$expected_count, 
                               one4=g1d$expected_count, 
                               one5=g1e$expected_count, 
                               one6=g1f$expected_count, 
                               one7=g1g$expected_count, 
                               one8=g1h$expected_count, 
                               two1=g2a$expected_count, 
                               two2=g2b$expected_count, 
                               two3=g2c$expected_count,  
                               two4=g2d$expected_count, 
                               two5=g2e$expected_count, 
                               two6=g2f$expected_count, 
                               two7=g2g$expected_count, 
                               two8=g2h$expected_count, 
                               row.names=g1a$gene_id))

#Use this line if you want to run DiffE using only Olfr genes
#countsTable2 <- countsTable[which(str_detect(rownames(countsTable), "Olfr")),]

#make the labels, needs to match countsTable
diffE_labels <- data.frame(FileName=c("one1","one2","one3","one4","one5","one6","one7","one8", 
                                      "two1","two2","two3","two4","two5","two6","two7","two8"),
                           Subject=c("A","B","C","D","E","F","G","H",
                                     "A","B","C","D","E","F","G","H"),
                           Treatment=c(0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1))

#run differential expression with all genes, check dimension of dataframe and view head
diffE_result <- DiffExp(diffE_labels, countsTable) %>% rename("gene_id" = "id")
```

    ## Warning in parametricDispersionFit(means, disps): Dispersion fit did not
    ## converge.

``` r
#dim(diffE_result)
diffE_result[1:5,1:19]
```

    ##                gene_id     logFC   logCPM       LR        PValue           FDR
    ## 1 ENSMUSG00000111590.2 -11.30811 4.175700 545.7242 1.071886e-120 5.947576e-116
    ## 2 ENSMUSG00000049334.3 -11.00785 3.980512 354.3144  4.871073e-79  1.351406e-74
    ## 3 ENSMUSG00000063137.5 -12.37683 4.225478 332.8950  2.250454e-74  4.162365e-70
    ## 4 ENSMUSG00000049573.3 -10.40797 3.557266 303.4543  5.823703e-68  8.078496e-64
    ## 5 ENSMUSG00000107748.4 -10.96619 2.797568 300.5796  2.463220e-67  2.733533e-63
    ##   one1 one2 one3 one4 one5 one6 one7 one8 two1 two2 two3 two4 two5
    ## 1  502  628 1058  292  337  683  755 1214    0    0    0    0    0
    ## 2    0 1988    0  852    0 1192    0  966    0    1    0    0    0
    ## 3    0    0    4 1902  108   32 1957 1486    0    0    0    0    0
    ## 4    0   39    0  896   42  248    0 2374    0    0    0    0    0
    ## 5    0  102  466  267  858    0    0  398    0    0    0    0    0

``` r
#id = gene name
#left side is edgeR output
#logFC = log fold change (Between group1 and group2)
#logCPM = log Counts Per Million
#Pvalue = raw significance
#FDR = multiple comparison corrected significance - use this, typically with 0.05 to start
#one1:two3 = counts per sample
#right side is DEseq output
#baseMean/A/B = Mean of counts
#padj = DEseq MCcorrected significance
```

## Prep results and plot all genes in volcano

``` r
#tibbles are nicer dataframes
der <-as_tibble(diffE_result)

###I hate dealing with attributes so will reload the file as tibble and then left_join names
write_csv(der, "~/Downloads/temp_results.csv")
der <- read_csv("~/Downloads/temp_results.csv") %>%
  rowwise() %>%
  mutate(trim_id = str_split(gene_id, "[.]")[[1]][1]) %>%
  select(-gene_id) %>%
  select(trim_id, everything()) %>%
  rename(gene_id = trim_id) %>%
  ungroup()

#change ensembl_geneIDs to gene names
#since this conversion file also have ensembl transcript IDs, remove duplicates first
ens_names <- read_csv("~/Desktop/obmap/r_analysis/inputs/convert_ensID_geneName.csv") %>% 
  select(ens_id, gene_name) %>% 
  rename("gene_id" = "ens_id") %>%
  unique()

der_named <- left_join(der, ens_names, by = "gene_id") %>% 
  select(gene_name, logFC:padj) %>% 
  rowwise() %>%
  mutate(isOlfr = ifelse(str_detect(gene_name, "Olfr"), T, F),
         isTaar = ifelse(str_detect(gene_name, "Taar"), T, F),
         isVmn = ifelse(str_detect(gene_name, "Vmn"), T, F),
         isPS = ifelse(str_detect(gene_name, "-ps"), T, F),
         isChemo = ifelse(sum(isOlfr, isTaar, isVmn) > 0, T, F),
         mkate2OR = ifelse(gene_name == "Olfr881", "Olfr881", 
                              ifelse(gene_name == "Olfr1377", 
                                     "Olfr1377", "Olfr"))) %>%
  ungroup()

write_csv(der_named, "~/Desktop/obmap/r_analysis/func_surface/allgenes_results.csv")

#quick volcano plot using edgeR
ggplot(der_named) + 
  geom_point(aes(logFC,-log10(FDR), alpha = 0.25, color = isOlfr)) + 
  geom_vline(xintercept = 0)  + 
  geom_hline(yintercept = -log10(0.05), linetype = 2) + 
  scale_color_manual(values = c("#000000","#FF0000"))
```

    ## Warning: Removed 19 rows containing missing values (geom_point).

![](wachR2_gen25_files/figure-gfm/unnamed-chunk-2-1.png)<!-- -->

``` r
##interactively to ask where do ORs fall among all genes
plot_ly(der_named, x = ~logFC, y = ~-log10(FDR), 
        color = ~isOlfr, 
        text = ~paste("ID: ", gene_name, 
                      "<br>logFC: ", logFC,                                                             "<br>FDR", FDR))
```

    ## Warning: `arrange_()` is deprecated as of dplyr 0.7.0.
    ## Please use `arrange()` instead.
    ## See vignette('programming') for more help
    ## This warning is displayed once every 8 hours.
    ## Call `lifecycle::last_warnings()` to see where this warning was generated.

    ## Warning in RColorBrewer::brewer.pal(N, "Set2"): minimal value for n is 3, returning requested palette with 3 different levels
    
    ## Warning in RColorBrewer::brewer.pal(N, "Set2"): minimal value for n is 3, returning requested palette with 3 different levels

![](wachR2_gen25_files/figure-gfm/unnamed-chunk-2-2.png)<!-- -->

## plot only chemoreceptors and calculate new FDR for subsets of all genes

``` r
#itnact chemoreceptors
only_chemo <- der_named %>% 
  filter(isChemo == T) %>%
  filter(isPS == F)

only_or <- only_chemo %>% 
  filter(isOlfr == T)
only_taar <- only_chemo %>% 
  filter(isTaar == T)
only_vmn <- only_chemo %>% 
  filter(isVmn == T)

#recalculate FDR
only_or$FDRor <- p.adjust(only_or$PValue, method = "fdr")
only_taar$FDRtaar <- p.adjust(only_taar$PValue, method = "fdr")
only_vmn$FDRvmn <- p.adjust(only_vmn$PValue, method = "fdr")

makers <- only_or %>% mutate(fisurface = ifelse(logFC > 0, ifelse(FDRor < 0.05, T, F), F))
fisigor <- only_or$gene_name[which(only_or$logFC > 0 & only_or$FDRor < 0.05)]
only_taar$gene_name[which(only_taar$logFC > 0 & only_taar$FDRtaar < 0.05)]
```

    ## [1] "Taar8a" "Taar5"  "Taar8b" "Taar7f" "Taar3"  "Taar9"  "Taar7d" "Taar7e"
    ## [9] "Taar2"

``` r
only_vmn$gene_name[which(only_vmn$logFC > 0 & only_vmn$FDRvmn < 0.05)]
```

    ##  [1] "Vmn1r201" "Vmn1r84"  "Vmn1r232" "Vmn1r24"  "Vmn2r23"  "Vmn1r57" 
    ##  [7] "Vmn1r18"  "Vmn1r198" "Vmn2r59"  "Vmn2r21"  "Vmn1r40"  "Vmn2r8"  
    ## [13] "Vmn2r33"  "Vmn1r214" "Vmn1r168" "Vmn2r3"   "Vmn1r87"  "Vmn2r75" 
    ## [19] "Vmn1r142" "Vmn2r112" "Vmn2r82"  "Vmn2r68"  "Vmn1r160" "Vmn1r73"

``` r
write_csv(only_or, 
          "~/Desktop/obmap/r_analysis/func_surface/olfr_fdradj_201123.csv")
write_csv(only_taar, 
          "~/Desktop/obmap/r_analysis/func_surface/taar_fdradj_201123.csv")
write_csv(only_vmn, 
          "~/Desktop/obmap/r_analysis/func_surface/vmn_fdradj_201123.csv")

only_or %>% 
  filter(isOlfr == T) %>%
  ggplot(aes(logFC,-log10(FDRor), alpha = 0.25)) +
  geom_point() + 
  geom_vline(xintercept = 0)  + 
  geom_hline(yintercept = -log10(0.05), linetype = 2) +
  ggtitle("Olfrs from All Genes DE", subtitle = "FDR re-calculated") +
  xlab("Not FI surface  <<<   logFC   >>>     FI surface")
```

![](wachR2_gen25_files/figure-gfm/unnamed-chunk-3-1.png)<!-- -->

``` r
only_taar %>% 
  filter(isTaar == T) %>%
  ggplot(aes(logFC, -log10(FDRtaar), alpha = 0.25)) + 
  geom_point() + 
  geom_vline(xintercept = 0) + 
  geom_hline(yintercept = -log10(0.05)) +
  ggtitle("Taars from All Genes DE", subtitle = "FDR re-calculated") +
  xlab("Not FI surface  <<<   logFC   >>>     FI surface")
```

![](wachR2_gen25_files/figure-gfm/unnamed-chunk-3-2.png)<!-- -->

``` r
only_vmn %>% 
  filter(isVmn == T) %>%
  ggplot(aes(logFC, -log10(FDRvmn), alpha = 0.25)) + 
  geom_point() + 
  geom_vline(xintercept = 0) + 
  geom_hline(yintercept = -log10(0.05)) +
  ggtitle("Vmnrs from All Genes DE", subtitle = "FDR re-calculated") +
  xlab("Not FI surface  <<<   logFC   >>>     FI surface")
```

![](wachR2_gen25_files/figure-gfm/unnamed-chunk-3-3.png)<!-- -->

``` r
plot_ly(only_chemo, x = ~logFC, y = ~-log10(FDR), 
        color = ~isTaar,
        text = ~paste("ID: ", gene_name,
                      "<br>logFC: ", logFC,
                      "<br>FDR", FDR))
```

    ## No trace type specified:
    ##   Based on info supplied, a 'scatter' trace seems appropriate.
    ##   Read more about this trace type -> https://plot.ly/r/reference/#scatter

    ## No scatter mode specifed:
    ##   Setting the mode to markers
    ##   Read more about this attribute -> https://plot.ly/r/reference/#scatter-mode

    ## Warning in RColorBrewer::brewer.pal(N, "Set2"): minimal value for n is 3, returning requested palette with 3 different levels
    
    ## Warning in RColorBrewer::brewer.pal(N, "Set2"): minimal value for n is 3, returning requested palette with 3 different levels

![](wachR2_gen25_files/figure-gfm/unnamed-chunk-3-4.png)<!-- -->

``` r
info <- read_csv("~/Desktop/obmap/r_analysis/3dimOB/input/knowntanwavgFI.csv") %>%
  rename("gene_name" = "gene") %>%
  select(gene_name:RTP)
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
feats <- only_or %>% 
  select(gene_name:FDR, FDRor, isOlfr:mkate2OR) %>% 
  left_join(info, by = "gene_name") %>%
  filter(isPS == F)

# how many from each class
feats %>% 
  filter(FDRor < 0.05) %>%
  filter(logFC > 0) %>%
  group_by(class) %>%
  summarise(count = n())
```

    ## `summarise()` ungrouping output (override with `.groups` argument)

    ## # A tibble: 2 x 2
    ##   class count
    ##   <dbl> <int>
    ## 1     1    27
    ## 2     2    94

``` r
# how many from each OE region
feats %>% 
  filter(FDRor < 0.05) %>%
  filter(logFC > 0) %>%
  group_by(oe_region) %>%
  summarise(count = n())
```

    ## `summarise()` ungrouping output (override with `.groups` argument)

    ## # A tibble: 3 x 2
    ##   oe_region count
    ##   <chr>     <int>
    ## 1 Dorsal       99
    ## 2 Ventral      21
    ## 3 <NA>          1

``` r
# how many TAARs enriched in FI surface
only_taar %>% 
  filter(FDRtaar < 0.05) %>%
  filter(logFC > 0) %>% 
  select(gene_name, PValue, FDR, FDRtaar)
```

    ## # A tibble: 9 x 4
    ##   gene_name   PValue      FDR  FDRtaar
    ##   <chr>        <dbl>    <dbl>    <dbl>
    ## 1 Taar8a    1.44e-28 1.82e-25 2.16e-27
    ## 2 Taar5     5.09e- 6 6.87e- 4 3.82e- 5
    ## 3 Taar8b    1.37e- 5 1.68e- 3 6.83e- 5
    ## 4 Taar7f    1.07e- 4 1.07e- 2 4.01e- 4
    ## 5 Taar3     1.73e- 4 1.64e- 2 5.19e- 4
    ## 6 Taar9     4.61e- 4 3.88e- 2 1.13e- 3
    ## 7 Taar7d    5.28e- 4 4.38e- 2 1.13e- 3
    ## 8 Taar7e    1.84e- 3 1.25e- 1 3.46e- 3
    ## 9 Taar2     8.14e- 3 4.10e- 1 1.36e- 2

# Subset count table to just ORs then run DE

``` r
gene_list <- as.character(g1a$gene_id) %>% 
  as_tibble() %>% 
  rename("gene_id" = "value") %>% 
  rowwise() %>% 
  mutate(trim_id = str_split(gene_id, "[.]")[[1]][1]) %>% 
  rename("dot_id" = "gene_id", "gene_id" = "trim_id") %>% 
  left_join(ens_names, by = "gene_id") %>% 
  mutate(isOlfr = str_detect(gene_name, "Olfr")) %>% 
  filter(isOlfr == T)

orct <- countsTable[which(rownames(countsTable) %in% gene_list$dot_id),]

#run differential expression, check dimension of dataframe and view head
olfr_result <- DiffExp(diffE_labels, orct) %>% 
  as_tibble() %>% 
  rename("gene_id" = "id") %>% 
  mutate(isSignif = ifelse(FDR < 0.05, "FDR<0.05", "Notsig")) %>%
  rowwise() %>% 
  mutate(trim_id = str_split(gene_id, "[.]")[[1]][1]) %>% 
  rename("dot_id" = "gene_id", "gene_id" = "trim_id") %>% 
  left_join(ens_names, by = "gene_id") %>%
  select(gene_name, logFC:isSignif) %>%
  mutate(FI_sig = ifelse(FDR < 0.05 & logFC > 0, T, F),
         RM_sig = ifelse(FDR < 0.05 & logFC < 0, T, F),
         mkate2OR = ifelse(gene_name == "Olfr881", "Olfr881", 
                              ifelse(gene_name == "Olfr1377", "Olfr1377", 
                                     "Olfr")))

olfr_result %>% filter(FI_sig == T)

ggplot(olfr_result) + 
  geom_point(aes(logFC,-log10(FDR), alpha = 0.25, color = mkate2OR)) + 
  geom_vline(xintercept = 0) + 
  geom_hline(yintercept = -log10(0.05))

write_csv(olfr_result, "~/Desktop/obmap/r_analysis/func_surface/orOnlyDiffE_info.csv")
```
