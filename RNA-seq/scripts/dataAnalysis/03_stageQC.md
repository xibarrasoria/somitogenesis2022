---
title: "Stage QC"
date: '08 May, 2020'
output:
  html_document:
    keep_md: true
    fig_width: 5
    fig_height: 5
    fig_caption: yes
    code_folding: hide
    toc: true
    toc_float: 
      collapsed: false
---



Hox genes are critical in establishing somite identity. In mouse, there are 13 different genes (1 to 13), each with up to four paralogs (A to D). As such, there is substantial redundancy in their expression patterns and function. Hox genes from each paralogous group are clustered and ordered from 1 to 13. The first Hox genes are expressed the earliest and subsequent genes are expressed later in development. That is, Hox1 genes are the first to be expressed and Hox13 are the last.

Each somite has a specific combination of Hox gene expression, with earlier somites having predominantly expression from early Hox genes, whilst these are downregulated in later somites that instead express later Hox genes. Such combinatorial code gives a specific identity to each somite and determines the type of structures it will develop into, i.e. cervical, thoracic, lumbar, sacral or caudal.

This combinatorial Hox code is therefore directly correlated with somite stage, and we can use it to make sure we have the correct stages reported in the experimental design.

### Data

We have the complete transcriptome of individual mouse somites spanning early (8 somites) to late (35 somites) development. At each stage, we have sequenced the last three somites that were segmented (with somiteI being the most recent, and somiteIII the one segmented two cycles before). 


```r
# metadata
meta <- read.table(paste0(dir, "RNA-seq/data/metadata_RNAseq.tsv"), stringsAsFactors = FALSE, header = TRUE)
meta.all <- meta
meta <- meta[meta$QC ==1,] ## use only good-qual samples

# normalised, batch-corrected data
dataNorm <- read.table(paste0(dir, "RNA-seq/data/geneCounts.NORM_batchCorrected_14PCs.tsv"), check.names = FALSE)
stopifnot(identical(colnames(dataNorm)[-1], meta$sample))

table(meta$somite, stage=meta$stage)
```

```
##       stage
##        8 18 21 25 27 35
##   SI   4  5  4  4  4  6
##   SII  4  5  3  3  4  4
##   SIII 4  5  4  4  3  6
```

Given the stages we have sampled we should have a good representation of the Hox code, with expression of most, if not all, Hox genes.

### Hox gene expression

Simply plotting the relative expression of all Hox genes across the samples we see good progression of expression from early to late Hox genes, as development proceeds. Stages are indicated at the top:

- green = cervical (stage8)
- light blue = thoracic (stag18)
- blue = thoracic (stage21)
- dark blue = thoracic (stage25)
- orange = lumbar (stage27)
- red = sacral (stage35)

And the paralogous Hox groups are coloured from yellow to orange to brown, from 1 to 13.


```r
hox <- dataNorm[grep("^Hox",dataNorm$gene),]
hox <- hox[grep("s", hox$gene, invert=TRUE),] # remove antisense and non-coding genes
row.names(hox) <- hox$gene
hox$gene <- NULL

hox.z <- t(apply(hox, 1, function(x) (x-mean(x))/sd(x) )) # visualise with z-scores

hox.num <- as.numeric(substr(row.names(hox.z), 5,6))
names(hox.num) <- row.names(hox.z)

cols.hox <- colorRampPalette(colors = c("yellow", "orange", "brown"))(length(unique(hox.num)))[hox.num]
cols.stage <- as.character(factor(meta$stage, labels = c("darkolivegreen4", "skyblue", "cornflowerblue", "dodgerblue4", "darkorange", "indianred1")))

heatmap.2(as.matrix(hox.z[order(hox.num), order(meta$stage)]), trace="none", col=brewer.pal(n=8, name="Blues"), ColSideColors = cols.stage[order(meta$stage)], RowSideColors = cols.hox[order(hox.num)], Colv = FALSE, dendrogram = "none", labCol = rep("",ncol(hox.z)), Rowv = FALSE, key.title = "", key.xlab = "z-score", key.ylab = "", density.info="none")
```

![](03_stageQC_files/figure-html/hoxGenes-1.png)<!-- -->

Cervical somites express only Hox genes 1-5 and these genes decrease in expression as development proceeds. Thoracic somites express genes 6-9. Lumbar somites switch to express genes 10 and 11. Finally, sacral somites express genes 10 to 13, with 12 and 13 being specific to them.

And given that the expression of the Hox genes is so consistent with stage, a PCA on only these 39 genes recapitulates the developmental progression of the samples. The first two PCs separate the cervical and sacral somites from the trunk somites, whereas the third PC separates the latest thoracic (stage 25) from the earlier ones (18 and 21).


```r
pca <- prcomp(t(hox))
eigs <- pca$sdev^2
prop.var <- eigs / sum(eigs) * 100

df <- as.data.frame(pca$x)
plots <- list()
plots[[1]] <- ggplot(df, aes(PC1, PC2, colour=as.factor(meta$stage), shape=meta$somite)) + geom_point(cex=2) + xlab(paste("PC1 -", round(prop.var[1],2), "%")) + ylab(paste("PC2 -", round(prop.var[2],2), "%")) + labs(colour="stage", shape="somite") + scale_color_manual(values = c("darkolivegreen4", "skyblue", "cornflowerblue", "dodgerblue4", "darkorange", "indianred1")) + th
plots[[2]] <- ggplot(df, aes(PC1, PC3, colour=as.factor(meta$stage), shape=meta$somite)) + geom_point(cex=2) + xlab(paste("PC1 -", round(prop.var[1],2), "%")) + ylab(paste("PC3 -", round(prop.var[3],2), "%")) + labs(colour="stage", shape="somite") + scale_color_manual(values = c("darkolivegreen4", "skyblue", "cornflowerblue", "dodgerblue4", "darkorange", "indianred1")) + th
ggarrange(plotlist = plots, ncol = 2, common.legend = TRUE, legend = "right")
```

![](03_stageQC_files/figure-html/hox_pca-1.png)<!-- -->

Two samples from stage 27 (both from the same embryo) separate substantially from the rest and associate more closely to the stage35 samples. These two show expression of late Hox genes like *Hoxd12* and *Hoxd13* that resemble more the expression of stage35 samples, compared to those from stage27.


```r
plots <- list()

df <- data.frame(expr=as.numeric(dataNorm[dataNorm$gene=="Hoxd12",-1]), stage=as.factor(meta$stage))
plots[[1]] <- ggplot(df, aes(stage, expr, colour=stage)) + geom_boxplot() + ggtitle("Hoxd12") + th + scale_color_manual(values = c("darkolivegreen4", "skyblue", "cornflowerblue", "dodgerblue4", "darkorange", "indianred1"))
df <- data.frame(expr=as.numeric(dataNorm[dataNorm$gene=="Hoxd13",-1]), stage=as.factor(meta$stage))
plots[[2]] <- ggplot(df, aes(stage, expr, colour=stage)) + geom_boxplot() + ggtitle("Hoxd13") + scale_color_manual(values = c("darkolivegreen4", "skyblue", "cornflowerblue", "dodgerblue4", "darkorange", "indianred1")) + th

ggarrange(plotlist = plots, ncol = 2, legend = "none")
```

![](03_stageQC_files/figure-html/outlier27-1.png)<!-- -->

This suggests that this embryo might have been misannotated and it is of a more advanced stage than 27 somites, although possibly not as far as 35 somites. Therefore, we remove this embryo from downstream analyses.


```r
## retrieve the two stage27 samples with high Hoxd12-13 expression
idx <- which(colnames(dataNorm) %in% meta[meta$stage==27,]$sample)
wrongStage <- colnames(dataNorm[,idx])[which(dataNorm[dataNorm$gene=="Hoxd12", idx]>2)]
stopifnot(identical(wrongStage, colnames(dataNorm[,idx])[which(dataNorm[dataNorm$gene=="Hoxd13", idx]>2)]))

meta.all$wrongStage <- ifelse(meta.all$sample %in% wrongStage, 1, 0)
write.table(meta.all, paste0(dir, "RNA-seq/data/metadata_RNAseq.tsv"), quote = FALSE, sep="\t", row.names = FALSE)
```



```r
sessionInfo()
```

```
## R version 3.6.1 (2019-07-05)
## Platform: x86_64-apple-darwin15.6.0 (64-bit)
## Running under: macOS High Sierra 10.13.6
## 
## Matrix products: default
## BLAS:   /Library/Frameworks/R.framework/Versions/3.6/Resources/lib/libRblas.0.dylib
## LAPACK: /Library/Frameworks/R.framework/Versions/3.6/Resources/lib/libRlapack.dylib
## 
## locale:
## [1] en_GB.UTF-8/en_GB.UTF-8/en_GB.UTF-8/C/en_GB.UTF-8/en_GB.UTF-8
## 
## attached base packages:
## [1] stats     graphics  grDevices utils     datasets  methods   base     
## 
## other attached packages:
## [1] gplots_3.0.1.1     RColorBrewer_1.1-2 ggpubr_0.2.5       magrittr_1.5      
## [5] ggplot2_3.3.0      edgeR_3.28.0       limma_3.42.0      
## 
## loaded via a namespace (and not attached):
##  [1] Rcpp_1.0.3         pillar_1.4.2       compiler_3.6.1     bitops_1.0-6      
##  [5] tools_3.6.1        digest_0.6.23      evaluate_0.14      lifecycle_0.1.0   
##  [9] tibble_2.1.3       gtable_0.3.0       lattice_0.20-38    pkgconfig_2.0.3   
## [13] rlang_0.4.2        rstudioapi_0.10    yaml_2.2.0         xfun_0.11         
## [17] gridExtra_2.3      withr_2.1.2        stringr_1.4.0      dplyr_0.8.3       
## [21] knitr_1.26         caTools_1.17.1.3   gtools_3.8.1       cowplot_1.0.0     
## [25] locfit_1.5-9.1     grid_3.6.1         tidyselect_0.2.5   glue_1.3.1        
## [29] R6_2.4.1           rmarkdown_1.18     gdata_2.18.0       farver_2.0.1      
## [33] purrr_0.3.3        scales_1.1.0       htmltools_0.4.0    assertthat_0.2.1  
## [37] colorspace_1.4-1   ggsignif_0.6.0     labeling_0.3       KernSmooth_2.23-16
## [41] stringi_1.4.3      munsell_0.5.0      crayon_1.3.4
```
