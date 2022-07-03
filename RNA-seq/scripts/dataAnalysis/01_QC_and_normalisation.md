---
title: "RNA-seq data QC and normalisation"
date: '02 May, 2020'
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



We have sequenced 77 RNA-seq libraries from individual mouse somites, representing the three most recently segmented somites per embryo, across six different stages. Stages are named based on the total number of somites in the embryo. We have 3 to 6 biological replicates per stage.


```r
meta <- read.table(paste0(dir, "RNA-seq/data/metadata_RNAseq.tsv"), stringsAsFactors = FALSE, header = TRUE)
table(meta$somite, meta$stage)
```

```
##       
##        8 18 21 25 27 35
##   SI   4  5  4  4  4  6
##   SII  4  5  3  3  4  5
##   SIII 4  5  4  4  3  6
```

We have mapped the data to the mouse reference genome (`mm10`) and quantified expression by counting the number of aligned fragments per gene (using Ensembl annotation, version 96). 

### Quality control

Samples were sequenced at a median depth of nearly 17.5 million fragments; 75% of the data have library sizes between 14.9 and 23.2 million fragments.


```r
## data
data <- read.table(paste0(dir, "RNA-seq/data/geneCounts.RAW.tsv"), stringsAsFactors = FALSE, check.names = FALSE)
stopifnot(identical(colnames(data)[-1], meta$sample))

## mapping statistics (from STAR logs)
mapping.stats <- read.table(paste0(dir, "RNA-seq/data/mappingStatistics.tsv"), stringsAsFactors = FALSE, header = TRUE)
row.names(mapping.stats) <- meta[match(mapping.stats$sample, meta$library),]$sample
stopifnot(identical(colnames(data)[-1], row.names(mapping.stats)))

## counting statistics (from STAR gene counts)
counting.stats <- read.table(paste0(dir, "RNA-seq/data/countingStatistics.tsv"), stringsAsFactors = FALSE)
stopifnot(identical(colnames(data)[-1], row.names(counting.stats)))
counting.stats$uniquelyMapped <- colSums(data[,-1])

ggplot(mapping.stats, aes(1, total/1e6)) + geom_violin() + geom_dotplot(binaxis='y', stackdir='center', dotsize=0.5) + ylab("library size (million fragments)") + xlab("") + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
```

![](01_QC_and_normalisation_files/figure-html/depth-1.png)<!-- -->

Nearly 15 million of these fragments are uniquely mapped, 12.8 million of which are within annotated exons (interquartile range = 10.5 - 16.4 million).


```r
ggplot(counting.stats, aes(1, uniquelyMapped/1e6)) + geom_violin() + geom_dotplot(binaxis='y', stackdir='center', dotsize=0.5) + ylab("total fragments uniquely mapped to exons (millions)") + xlab("") + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
```

![](01_QC_and_normalisation_files/figure-html/inExons-1.png)<!-- -->

One sample, e17_SII-2, has a very small library with only 88342 total fragments and will be removed from downstream analyses. The next smallest library has 5.3 million fragments which might or might not be enough for a representative transcriptome.

In general, most samples show good mapping statistics with low proportion of unmapped and multimapped fragments, and of fragments mapped to non-exonic regions. There are a few outliers that show larger numbers, but since these have large libraries (over 20 million), the proportion of uniquely mapped fragments is still large.


```r
mapping.stats$unmapped <- mapping.stats$total-(mapping.stats$unique+mapping.stats$multimapped)
plots <- list()
plots[[1]] <- ggplot(mapping.stats, aes(1, unmapped/total*100)) + geom_violin(trim=FALSE) + geom_jitter(aes(1,unmapped/total*100, col=total/1e6), width = 0.05) + ylab("% unmapped") + scale_colour_gradientn(colours = (brewer.pal(n = 11, name = "PiYG"))) + labs(colour="libsize") + xlab("") + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())

plots[[2]] <- ggplot(mapping.stats, aes(1, multimapped/total*100)) + geom_violin(trim=FALSE) + geom_jitter(aes(1,multimapped/total*100, col=total/1e6), width = 0.05) + ylab("% multimapped") + scale_colour_gradientn(colours = (brewer.pal(n = 11, name = "PiYG"))) + labs(colour="libsize") + xlab("") + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())

plots[[3]] <- ggplot(counting.stats, aes(1, N_noFeature/mapping.stats$total*100)) + geom_violin(trim=FALSE) + geom_jitter(aes(1,N_noFeature/mapping.stats$total*100, col=mapping.stats$total/1e6), width = 0.05) + ylab("% outside exons") + scale_colour_gradientn(colours = (brewer.pal(n = 11, name = "PiYG"))) + labs(colour="libsize") + xlab("") + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())

plots[[4]] <- ggplot(counting.stats, aes(1, N_ambiguous/mapping.stats$total*100)) + geom_violin(trim=FALSE) + geom_jitter(aes(1,N_ambiguous/mapping.stats$total*100, col=mapping.stats$total/1e6), width = 0.05) + ylab("% ambiguous") + scale_colour_gradientn(colours = (brewer.pal(n = 11, name = "PiYG"))) + labs(colour="libsize") + xlab("") + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())

ggarrange(plotlist = plots, ncol = 2, nrow = 2)
```

![](01_QC_and_normalisation_files/figure-html/qc-1.png)<!-- -->

The two samples with the smallest libraries have good stats, with only ~7% of fragments unmapped, multimapped, outside exons and ambiguous. That leaves around 70% of fragments uniquely mapped and thus suggests the libraries are of good quality, but sequenced to lower depth.

Samples generally have around 22 thousand genes expressed, with a tight IQR between 21.1 and 23 thousand. The sample with the very small library is a clear outlier and needs to be removed. The second smallest also deviates a bit from the rest of the dataset, but has 18,565 detected genes, compared to 19,549 of the next sample.


```r
counting.stats$nGenes <- apply(data[,-1], 2, function(x) sum(x>0))
# summary(counting.stats$nGenes)

ggplot(counting.stats, aes(1, nGenes/1e3)) + geom_violin(trim=FALSE) + geom_jitter(aes(1,nGenes/1e3, col=mapping.stats$total/1e6), width = 0.02) + ylab("total genes detected (thousands)") + scale_colour_gradientn(colours = (brewer.pal(n = 11, name = "PiYG"))) + labs(colour="libsize") + xlab("") + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
```

![](01_QC_and_normalisation_files/figure-html/numGenes-1.png)<!-- -->

Given that the second smallest library shows a similar number of detected genes and good mapping statistics (reflecting good sample quality), it is likely that its data is useful, specially for genes expressed at moderate to high levels. Therefore, we will only remove the sample that almost failed sequencing and proceed with all others.


```r
bad.qual <- row.names(mapping.stats[which.min(mapping.stats$total),])
data <- data[,-which(colnames(data) == bad.qual)]
meta$QC <- ifelse(meta$sample == bad.qual, 0, 1)
```

### Normalisation

The first thing to account for is library size. Before estimating size factors, we filter very lowly expressed genes.

We retain around 40% of the annotated genes (46% of detected genes).


```r
## our factor of interest is the combination of stage and somite. 
# create a variable for it
meta$group <- paste(meta$stage, meta$somite, sep=".")

## create an edgeR object
y <- DGEList(counts=data[,-1], samples = meta[meta$QC==1,], genes = data[,1], group = meta[meta$QC==1,]$group)

## filter low abundance genes
means <- aveLogCPM(y)
keep <- filterByExpr(y)

plot(density(means[keep]), lwd=2, xlab="average log counts-per-million", main="", bty="l", ylim=c(0,1), xlim=c(-3,15))
lines(density(means[-keep]), lty=2, lwd=2)
legend("topright", legend = c("kept", "filtered"), lty=c(1,2), lwd=2)
```

![](01_QC_and_normalisation_files/figure-html/filter-1.png)<!-- -->

```r
summary(keep)
```

```
##    Mode   FALSE    TRUE 
## logical   35474   20062
```

```r
y <- y[keep, , keep.lib.sizes=FALSE]
```

With this reduced dataset of 20062 genes we calculate size factors to normalise the data. Size factors normalise for both library size and composition biases.

The method successfully standardises the libraries. The sample with second smallest library size and slightly fewer detected genes is highlighted in red; the normalisation seems to successfully bring it to the overall dataset distribution, so it should be ok to keep.


```r
## normalisation
y <- calcNormFactors(y)

dataNorm <- cpm(y, log=TRUE, prior.count = 1)

col <- c(rep("black",18), "red", rep("black",69)) # highlight the sample with small library size
par(mfrow=c(1,2))
boxplot(log2(data[keep,-1]+1), main="RAW", ylab=expression('log'[2]*' counts + 1'), axes=FALSE, border=col); box(bty="l"); axis(2, las=2)
boxplot(dataNorm, main="NORMALISED", ylab=expression('log'[2]*' counts-per-million + 1'), axes=FALSE, border=col); box(bty="l"); axis(2, las=2)
```

![](01_QC_and_normalisation_files/figure-html/norm-1.png)<!-- -->

### Exploratory analysis

With the normalised data, we can explore the general features of the dataset. Firstly, we compute a PCA on the thousand most variable genes to check if we can detect the biological effect of interest; we also check for any potential confounders.


The first PC is dominated by collection date, whereas the embryo stage gets pushed to PC2. This suggests that we have strong batch effects from sample processing date that need to be accounted for. There is no obvious grouping by somite, suggesting differences are much more subtle.


```r
## variance-stabilise data
tmp <- data[row.names(dataNorm), colnames(dataNorm)]
data.vst <- vst(as.matrix(tmp))

vars <- rowVars(data.vst)
names(vars) <- row.names(dataNorm)
pca <- prcomp(t(dataNorm[names(vars[order(vars, decreasing = TRUE)])[1:1000],]))

eigs <- pca$sdev^2
prop.var <- eigs / sum(eigs)
# plot(prop.var, ylab="proportion of variance explained", xlab="PC", bty="l")

df <- as.data.frame(pca$x)
plots <- list()
plots[[1]] <- ggplot(df, aes(PC1, PC2)) + geom_point(aes(col=as.factor(meta[meta$QC==1,]$stage))) + labs(col="stage") + ggtitle("stage") + xlab(paste("PC1 - ", round(prop.var[1]*100,2), "% var explained")) + ylab(paste("PC2 - ", round(prop.var[2]*100,2), "% var explained"))
plots[[2]] <- ggplot(df, aes(PC1, PC2)) + geom_point(aes(col=as.factor(meta[meta$QC==1,]$somite))) + labs(col="somite") + ggtitle("somite") + xlab(paste("PC1 - ", round(prop.var[1]*100,2), "% var explained")) + ylab(paste("PC2 - ", round(prop.var[2]*100,2), "% var explained"))
plots[[3]] <- ggplot(df, aes(PC1, PC2)) + geom_point(aes(col=as.factor(meta[meta$QC==1,]$date))) + labs(col="batch") + ggtitle("collection date") + xlab(paste("PC1 - ", round(prop.var[1]*100,2), "% var explained")) + ylab(paste("PC2 - ", round(prop.var[2]*100,2), "% var explained"))
ggarrange(plotlist = plots, ncol=3)
```

![](01_QC_and_normalisation_files/figure-html/pca-1.png)<!-- -->



```r
## metadata with the additional info of which samples are used in downstream analyses
write.table(meta, paste0(dir, "RNA-seq/data/metadata_RNAseq.tsv"), quote = FALSE, sep="\t", row.names = FALSE)

## save the normalised counts, but first add gene names
dataNorm <- cbind(data[match(row.names(dataNorm), row.names(data)),1], as.data.frame(dataNorm))
colnames(dataNorm)[1] <- "gene"
write.table(dataNorm, paste0(dir, "RNA-seq/data/geneCounts.NORM_logCPM.tsv"), quote = FALSE, sep="\t")

## also save edgeR object
saveRDS(y, paste0(dir, "RNA-seq/results/01_edgeRobject.Rds"))
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
##  [1] DESeq2_1.26.0               SummarizedExperiment_1.16.0
##  [3] DelayedArray_0.12.0         BiocParallel_1.20.0        
##  [5] matrixStats_0.55.0          Biobase_2.46.0             
##  [7] GenomicRanges_1.38.0        GenomeInfoDb_1.22.0        
##  [9] IRanges_2.20.1              S4Vectors_0.24.1           
## [11] BiocGenerics_0.32.0         ggpubr_0.2.5               
## [13] magrittr_1.5                ggplot2_3.3.0              
## [15] RColorBrewer_1.1-2          edgeR_3.28.0               
## [17] limma_3.42.0               
## 
## loaded via a namespace (and not attached):
##  [1] bit64_0.9-7            splines_3.6.1          Formula_1.2-3         
##  [4] assertthat_0.2.1       latticeExtra_0.6-28    blob_1.2.0            
##  [7] GenomeInfoDbData_1.2.2 yaml_2.2.0             RSQLite_2.1.3         
## [10] pillar_1.4.2           backports_1.1.5        lattice_0.20-38       
## [13] glue_1.3.1             digest_0.6.23          ggsignif_0.6.0        
## [16] XVector_0.26.0         checkmate_1.9.4        colorspace_1.4-1      
## [19] cowplot_1.0.0          htmltools_0.4.0        Matrix_1.2-18         
## [22] XML_3.99-0.3           pkgconfig_2.0.3        genefilter_1.68.0     
## [25] zlibbioc_1.32.0        purrr_0.3.3            xtable_1.8-4          
## [28] scales_1.1.0           tibble_2.1.3           htmlTable_1.13.3      
## [31] annotate_1.64.0        farver_2.0.1           withr_2.1.2           
## [34] nnet_7.3-12            survival_3.1-8         crayon_1.3.4          
## [37] memoise_1.1.0          evaluate_0.14          foreign_0.8-72        
## [40] tools_3.6.1            data.table_1.12.6      lifecycle_0.1.0       
## [43] stringr_1.4.0          munsell_0.5.0          locfit_1.5-9.1        
## [46] cluster_2.1.0          AnnotationDbi_1.48.0   compiler_3.6.1        
## [49] rlang_0.4.2            grid_3.6.1             RCurl_1.95-4.12       
## [52] rstudioapi_0.10        htmlwidgets_1.5.1      labeling_0.3          
## [55] bitops_1.0-6           base64enc_0.1-3        rmarkdown_1.18        
## [58] gtable_0.3.0           DBI_1.0.0              R6_2.4.1              
## [61] gridExtra_2.3          knitr_1.26             dplyr_0.8.3           
## [64] zeallot_0.1.0          bit_1.1-14             Hmisc_4.3-0           
## [67] stringi_1.4.3          Rcpp_1.0.3             geneplotter_1.64.0    
## [70] vctrs_0.2.0            rpart_4.1-15           acepack_1.4.1         
## [73] tidyselect_0.2.5       xfun_0.11
```
