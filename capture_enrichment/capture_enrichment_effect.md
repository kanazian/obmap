obmap\_enrichment\_effect
================
kanazian
Nov 23, 2020

\#Goal Determine how consistent enrichment is between biological and
technical replicates of mouse OBs. Determine the range of enrichment
fold for all ORs.

# setup and load data

``` r
library(tidyverse)
library(heatmaply)

rawdat <- read_csv("~/Desktop/obmap/r_analysis/capture_enrichment/obmap_enrichmenteffect_gen25.csv") %>%
  filter(!is.na(gene_name)) %>%
  mutate(isOR = str_detect(gene_name, "Olfr"),
         isPS = str_detect(gene_name, "-ps"))

ectopic <- c("Olfr287","Olfr32")
olfrs_nops <- rawdat %>% filter(isOR == T) %>% filter(isPS == F)
olfrs <- olfrs_nops %>% filter(!(gene_name %in% ectopic)) %>% select(-ens_raw, -ens_id)
```

# check enrichment effect

``` r
olfrs %>% ggplot() +
  geom_point(aes(M6nocapOB, M6OB))
```

![](capture_enrichment_effect_files/figure-gfm/unnamed-chunk-1-1.png)<!-- -->

``` r
print("M6 stats")
```

    ## [1] "M6 stats"

``` r
#number of ORs with reads
sum(length(which(olfrs$M6nocapOB > 0)))
```

    ## [1] 412

``` r
sum(length(which(olfrs$M6OB > 0)))
```

    ## [1] 844

``` r
mean(olfrs$M6nocapOB)
```

    ## [1] 0.0637276

``` r
mean(olfrs$M6OB)
```

    ## [1] 360.2719

``` r
median(olfrs$M6nocapOB)
```

    ## [1] 0

``` r
median(olfrs$M6OB)
```

    ## [1] 106.735

``` r
olfrs %>% ggplot() +
  geom_point(aes(M4OBnocap, M4OBcap))
```

![](capture_enrichment_effect_files/figure-gfm/unnamed-chunk-1-2.png)<!-- -->

``` r
#number of ORs with reads
print("M4 stats")
```

    ## [1] "M4 stats"

``` r
sum(length(which(olfrs$M4OBnocap > 0)))
```

    ## [1] 496

``` r
sum(length(which(olfrs$M4OBcap > 0)))
```

    ## [1] 992

``` r
mean(olfrs$M4OBnocap)
```

    ## [1] 0.1499821

``` r
mean(olfrs$M4OBcap)
```

    ## [1] 61.52858

``` r
median(olfrs$M4OBnocap)
```

    ## [1] 0

``` r
median(olfrs$M4OBcap)
```

    ## [1] 19.665

``` r
m6compare <- tibble(Olfrname = olfrs$gene_name,
                    m6nocap = olfrs$M6nocapOB,
                    m6cap = olfrs$M6OB) %>%
  mutate(enr = (m6cap + 1)/(m6nocap + 1))

# more abundant ORs were enriched more
ggplot(m6compare, aes(m6cap, enr)) +
  geom_point()
```

![](capture_enrichment_effect_files/figure-gfm/unnamed-chunk-1-3.png)<!-- -->

# biological replicates

``` r
bioreps <- olfrs %>% select(contains("OB")) %>% select(-M6nocapOB, -M4OBnocap)
bioreps <- olfrs %>% 
  select(M4OB, M14XOB, M13XOB, M12OB, 
         M16XOB, M9OBEC, M8OB, M10ECOB) 
# %>%
#   rowwise() %>% 
#   mutate(tots = sum(M4OB, M14XOB, M13XOB, M12OB, 
#                     M16XOB, M9OBEC, M8OB, M10ECOB)) %>%
#   filter(tots != 0) %>%
#   select(-tots)

#pairwise spearman correlation matrix
bio_pmtx <- matrix(data=0, nrow=ncol(bioreps), ncol=ncol(bioreps))
colnames(bio_pmtx) <- colnames(bioreps)
rownames(bio_pmtx) <- colnames(bioreps)
for (i in 1:ncol(bioreps)) {
  correlation <- vector("numeric", length = ncol(bioreps))
  for (j in 1:ncol(bioreps)) {
    if (i == j) {
      bio_pmtx[i,j] <- 1
    } else {
    bio_pmtx[i,j] <- cor(bioreps[i], bioreps[j], method = "s")
    }
  }
}

#calc mean cor for pairwise 
r1 <- sum(bio_pmtx[1,2:8])
r2 <- sum(bio_pmtx[1,3:8])
r3 <- sum(bio_pmtx[1,4:8])
r4 <- sum(bio_pmtx[1,5:8])
r5 <- sum(bio_pmtx[1,6:8])
r6 <- sum(bio_pmtx[1,7:8])
r7 <- sum(bio_pmtx[1,8])
sum(r1,r2,r3,r4,r5,r6,r7)/sum(7+6+5+4+3+2+1)
```

    ## [1] 0.7581894

``` r
heatmaply_cor(bio_pmtx)
```

![](capture_enrichment_effect_files/figure-gfm/unnamed-chunk-2-1.png)<!-- -->

``` r
ggplot(bioreps, aes(M4OB, M14XOB)) +
  geom_point(alpha = 0.25) +
  geom_smooth(method = "loess")
```

    ## `geom_smooth()` using formula 'y ~ x'

![](capture_enrichment_effect_files/figure-gfm/unnamed-chunk-2-2.png)<!-- -->

# technical replicates of capture enrichment (2 different arrays on same cDNA, same library)

``` r
techreps <- olfrs %>% select(NC_10Ar2, NC_10Br2, NC_4Ar2, NC_4Br2, OC_10Ar2, OC_10Br2, OC_4Ar2, OC_4Br2)

one <- cor(techreps$NC_10Ar2, techreps$OC_10Ar2, method = "spearman")
two <- cor(techreps$NC_10Br2, techreps$OC_10Br2, method = "spearman")
three <- cor(techreps$NC_4Ar2, techreps$OC_4Ar2, method = "spearman")
four <- cor(techreps$NC_4Br2, techreps$OC_4Br2, method = "spearman")
mean(c(one, two, three, four))
```

    ## [1] 0.9508835

``` r
#pairwise spearman correlation matrix
tech_pmtx <- matrix(data=0, nrow=ncol(techreps), ncol=ncol(techreps))
colnames(tech_pmtx) <- colnames(techreps)
rownames(tech_pmtx) <- colnames(techreps)
for (i in 1:ncol(techreps)) {
  correlation <- vector("numeric", length = ncol(techreps))
  for (j in 1:ncol(techreps)) {
    if (i == j) {
      tech_pmtx[i,j] <- 1
    } else {
    tech_pmtx[i,j] <- cor(techreps[i], techreps[j], method = "spearman")
    }
  }
}

heatmaply_cor(tech_pmtx)
```

![](capture_enrichment_effect_files/figure-gfm/unnamed-chunk-3-1.png)<!-- -->

``` r
ggplot(olfrs, aes(log2(M4OB + 1), log2(M4OBcap + 1))) + 
  geom_point(alpha = 0.25) + 
  geom_smooth(method = "lm")
```

    ## `geom_smooth()` using formula 'y ~ x'

![](capture_enrichment_effect_files/figure-gfm/unnamed-chunk-3-2.png)<!-- -->

``` r
ggplot(olfrs, aes(log2(NC_10Ar2 + 1), log2(OC_10Ar2 + 1))) + 
  geom_point(alpha = 0.25) + 
  geom_smooth(method = "lm")
```

    ## `geom_smooth()` using formula 'y ~ x'

![](capture_enrichment_effect_files/figure-gfm/unnamed-chunk-3-3.png)<!-- -->

``` r
ggplot(olfrs, aes(NC_10Ar2, OC_10Ar2)) + 
  geom_point(alpha = 0.25) + 
  geom_smooth(method = "lm")
```

    ## `geom_smooth()` using formula 'y ~ x'

![](capture_enrichment_effect_files/figure-gfm/unnamed-chunk-3-4.png)<!-- -->

``` r
#library(edgeR)
#maPlot(olfrs$NC_10Ar2, olfrs$OC_10Ar2)
#maPlot(rawdat$NC_10Ar2, rawdat$OC_10Ar2)
```

# ma plots

``` r
library(edgeR)

#best technical replicate - all genes
maPlot(rawdat$NC_10Ar2, rawdat$OC_10Ar2)
```

![](capture_enrichment_effect_files/figure-gfm/unnamed-chunk-4-1.png)<!-- -->

``` r
#best technical replicate - just ORs
maPlot(olfrs$NC_10Ar2, olfrs$OC_10Ar2)
```

![](capture_enrichment_effect_files/figure-gfm/unnamed-chunk-4-2.png)<!-- -->

``` r
#best biological replicate - just ORs
maPlot(olfrs$M4OB, olfrs$M4OBcap)
```

![](capture_enrichment_effect_files/figure-gfm/unnamed-chunk-4-3.png)<!-- -->

# fold enrichment using genomic DNA controls

``` r
b6dna <- olfrs %>% select(gene_id, b6cap100, b6nocap7, b6nocap9, m4cap, m4uncap, M6OB, M6nocapOB) %>% mutate(foldenr7 = b6cap100/(b6nocap7 + 0.1), foldenrm4 = m4cap/(m4uncap + 0.1), foldenrm6 = M6OB/(M6nocapOB + 0.1)) %>% arrange(desc(foldenrm6)) %>% mutate(order = factor(gene_id, gene_id))

ggplot(b6dna, aes(b6nocap7, b6cap100)) +
  geom_point()

ggplot(b6dna, aes(m4uncap, m4cap)) +
  geom_point()

ggplot(b6dna, aes(M6nocapOB, M6OB)) +
  geom_point()

ggplot(b6dna, aes(order, foldenrm6)) +
  geom_bar(stat="identity") +
  ggtitle("M6OB RNA fold enrichment")

ggplot(b6dna, aes(foldenrm4, foldenrm6)) + 
  geom_point() +
  geom_smooth(method = "lm")

ggplot(b6dna, aes(m4cap, M6OB)) + 
  geom_point() +
  geom_smooth(method = "lm")

ggplot(b6dna, aes(m4uncap, M6nocapOB)) + 
  geom_point() +
  geom_smooth(method = "lm")
```
