---
title: "Removal of technical variation from RNA-seq data"
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



We have processed and QCed RNA-seq data from mouse somites. All but one sample are of good enough quality for downstream analyses.


```r
meta <- read.table(paste0(dir, "RNA-seq/data/metadata_RNAseq.tsv"), stringsAsFactors = FALSE, header = TRUE)

data <- read.table(paste0(dir, "RNA-seq/data/geneCounts.RAW.tsv"), check.names = FALSE)
dataNorm <- read.table(paste0(dir, "RNA-seq/data/geneCounts.NORM_logCPM.tsv"), check.names = FALSE)

y <- readRDS(paste0(dir, "RNA-seq/results/01_edgeRobject.Rds"))
```

We observed strong technical effects relating to the date of sample collection. Since the experimental design is partially confounded with collection date, we cannot simply regress it out on the differential expression analysis. Thus, we need to identify and remove these technical variation.

One approach is to use PCA on the residuals after fitting the model including the biological effects of interest, and remove the amount of variation that is above that expected by chance. This will only preserve what is explicitly modelled in the design; in this case, the *stage* and *somite* of each sample. Any other effects, technical *and* biological, will be removed.

We fit a model including the design considering stage and somite, and apply PCA on the residuals of the fit. We then use the `parallelPCA` function from the `PCAtools` package to estimate the number of PCs to keep. This test is based on permuting the data to establish how much variation is explained under a random model; PCs that capture more variation than the random model are kept.


```r
## design matrix
meta$group <- as.factor(meta$group)
design <- model.matrix(~0+group, meta[meta$QC==1,]) ## only use samples that pass QC
colnames(design) <- paste0("stage", levels(meta$group))

## fit modoel
fit <- lmFit(dataNorm[,-1], design)
## obtain residuals
res <- residuals(fit, dataNorm[,-1])

## use a permutation test to estimate how many PCs explain more variation than expected by chance
n.pc <- parallelPCA(dataNorm[,-1])

## compute PCA on the residuals to identify major sources of variation unrelated to design
pcs <- prcomp(t(res))
```

For this dataset we retain the first 14 PCs.

Next, we regress out these PCs and run a PCA on the corrected data to check the batch effect. Samples separate by stage and somite better than before and, importantly, there is no segregation by collection date anymore, suggesting this procedure has identified and removed the variation associated wit the experiments' dates successfully.


```r
## regress out batch effect
dataNorm.corr <- removeBatchEffect(dataNorm[,-1], design = model.matrix(~0+group, meta[meta$QC==1,]), covariates = pcs$x[,1:n.pc$n])

## variance stablisation to obtain meaningful variable genes for the PCA
tmp <- data[row.names(dataNorm.corr), colnames(dataNorm.corr)]
data.vst <- vst(as.matrix(tmp))
data.vst <- removeBatchEffect(data.vst, design = model.matrix(~0+group, meta[meta$QC==1,]), covariates = pcs$x[,1:n.pc$n])

## select most variable genes
vars <- rowVars(data.vst)
names(vars) <- row.names(data.vst)
pca <- prcomp(t(dataNorm.corr[names(vars[order(vars, decreasing = TRUE)])[1:1000],]))

eigs <- pca$sdev^2
prop.var <- eigs / sum(eigs)
# plot(prop.var, ylab="proportion of variance explained", xlab="PC", bty="l")

df <- as.data.frame(pca$x)
plots <- list()
plots[[1]] <- ggplot(df, aes(PC1, PC2)) + geom_point(aes(col=as.factor(meta[meta$QC==1,]$stage))) + labs(col="stage") + ggtitle("stage") + xlab(paste("PC1 - ", round(prop.var[1]*100,2), "% var explained")) + ylab(paste("PC2 - ", round(prop.var[2]*100,2), "% var explained"))
plots[[2]] <- ggplot(df, aes(PC1, PC2)) + geom_point(aes(col=as.factor(meta[meta$QC==1,]$somite))) + labs(col="somite") + ggtitle("somite") + xlab(paste("PC1 - ", round(prop.var[1]*100,2), "% var explained")) + ylab(paste("PC2 - ", round(prop.var[2]*100,2), "% var explained"))
plots[[3]] <- ggplot(df, aes(PC1, PC2)) + geom_point(aes(col=as.factor(meta[meta$QC==1,]$date))) + labs(col="batch") + ggtitle("collection date") + xlab(paste("PC1 - ", round(prop.var[1]*100,2), "% var explained")) + ylab(paste("PC2 - ", round(prop.var[2]*100,2), "% var explained"))
ggarrange(plotlist = plots, ncol=2, nrow=2)
```

![](02_removeTechnicalVar_files/figure-html/plotPCA-1.png)<!-- -->

We can use these PCs as covariates to regress out when performing differential expression, to ensure that the technical effects are properly controlled for. 


```r
write.table(pcs$x, paste0(dir, "RNA-seq/results/02_pcs_residuals.tab"), quote = FALSE, sep="\t")

dataNorm.corr <- cbind(gene=dataNorm[,1], as.data.frame(dataNorm.corr))
write.table(dataNorm.corr, paste0(dir, "RNA-seq/data/geneCounts.NORM_batchCorrected_14PCs.tsv"), quote = FALSE, sep="\t")
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
## [1] parallel  stats4    stats     graphics  grDevices utils     datasets 
## [8] methods   base     
## 
## other attached packages:
##  [1] PCAtools_1.2.0              cowplot_1.0.0              
##  [3] lattice_0.20-38             reshape2_1.4.3             
##  [5] ggrepel_0.8.1               ggpubr_0.2.5               
##  [7] magrittr_1.5                ggplot2_3.3.0              
##  [9] RColorBrewer_1.1-2          DESeq2_1.26.0              
## [11] SummarizedExperiment_1.16.0 DelayedArray_0.12.0        
## [13] BiocParallel_1.20.0         matrixStats_0.55.0         
## [15] Biobase_2.46.0              GenomicRanges_1.38.0       
## [17] GenomeInfoDb_1.22.0         IRanges_2.20.1             
## [19] S4Vectors_0.24.1            BiocGenerics_0.32.0        
## [21] edgeR_3.28.0                limma_3.42.0               
## 
## loaded via a namespace (and not attached):
##  [1] bitops_1.0-6             bit64_0.9-7              tools_3.6.1             
##  [4] backports_1.1.5          R6_2.4.1                 irlba_2.3.3             
##  [7] rpart_4.1-15             Hmisc_4.3-0              DBI_1.0.0               
## [10] colorspace_1.4-1         nnet_7.3-12              withr_2.1.2             
## [13] tidyselect_0.2.5         gridExtra_2.3            bit_1.1-14              
## [16] compiler_3.6.1           htmlTable_1.13.3         labeling_0.3            
## [19] scales_1.1.0             checkmate_1.9.4          genefilter_1.68.0       
## [22] stringr_1.4.0            digest_0.6.23            foreign_0.8-72          
## [25] rmarkdown_1.18           XVector_0.26.0           base64enc_0.1-3         
## [28] pkgconfig_2.0.3          htmltools_0.4.0          htmlwidgets_1.5.1       
## [31] rlang_0.4.2              rstudioapi_0.10          RSQLite_2.1.3           
## [34] DelayedMatrixStats_1.8.0 farver_2.0.1             acepack_1.4.1           
## [37] dplyr_0.8.3              RCurl_1.95-4.12          BiocSingular_1.2.0      
## [40] GenomeInfoDbData_1.2.2   Formula_1.2-3            Matrix_1.2-18           
## [43] Rcpp_1.0.3               munsell_0.5.0            lifecycle_0.1.0         
## [46] stringi_1.4.3            yaml_2.2.0               zlibbioc_1.32.0         
## [49] plyr_1.8.4               grid_3.6.1               blob_1.2.0              
## [52] dqrng_0.2.1              crayon_1.3.4             splines_3.6.1           
## [55] annotate_1.64.0          locfit_1.5-9.1           zeallot_0.1.0           
## [58] knitr_1.26               pillar_1.4.2             ggsignif_0.6.0          
## [61] geneplotter_1.64.0       XML_3.99-0.3             glue_1.3.1              
## [64] evaluate_0.14            latticeExtra_0.6-28      data.table_1.12.6       
## [67] vctrs_0.2.0              gtable_0.3.0             purrr_0.3.3             
## [70] assertthat_0.2.1         xfun_0.11                rsvd_1.0.2              
## [73] xtable_1.8-4             survival_3.1-8           tibble_2.1.3            
## [76] AnnotationDbi_1.48.0     memoise_1.1.0            cluster_2.1.0
```

