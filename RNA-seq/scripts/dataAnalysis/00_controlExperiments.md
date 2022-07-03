---
title: "Control experiment to ensure left and right somites are equivalent"
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



Before creating matched RNA and ATAC-seq datasets using the two somites from each pair, it was important to determine that the two somites from a pair are truly equivalent.

To test this, we collected the first two pairs of somites from two different embryos from the same litter. All somites were used to produce RNA-seq libraries. Additionally, the RNA from a somite pair was used to produce libraries again, but using only a quarter of the recommended reagents ('bis' samples); this was to assess whether using less reagents resulted in the same quality data.


```r
## metadata
meta <- read.table(paste0(dir, "RNA-seq/data/metadata_RNAseq_CONTROLlibraries.tsv"), stringsAsFactors = FALSE, header = TRUE)
meta
```

```
##              sample library embryo somite side stage     date size     conc
## 1      ctrl_e1_SI-1 do16433     e1     SI    1 20-25 03.08.17  434 12.51092
## 2      ctrl_e1_SI-2 do16434     e1     SI    2 20-25 03.08.17  654 11.37400
## 3     ctrl_e1_SII-1 do16435     e1    SII    1 20-25 03.08.17  509 12.28483
## 4     ctrl_e1_SII-2 do16436     e1    SII    2 20-25 03.08.17  420 12.92770
## 5      ctrl_e2_SI-1 do16437     e2     SI    1 20-25 03.08.17  432 17.52016
## 6      ctrl_e2_SI-2 do16438     e2     SI    2 20-25 03.08.17  393 13.98273
## 7     ctrl_e2_SII-1 do16439     e2    SII    1 20-25 03.08.17  461 12.27829
## 8     ctrl_e2_SII-2 do16440     e2    SII    2 20-25 03.08.17  570 11.43267
## 9  ctrl_e1_SII-1bis do16441     e1    SII    1 20-25 03.08.17  682 34.50733
## 10 ctrl_e1_SII-2bis do16442     e1    SII    2 20-25 03.08.17  630 16.40458
```

### Quality control 

The ten samples were sequenced across one lane of an Illumina HiSeq 4000, producing single-end 50bp reads. All samples were sequenced successfully, producing around 37+-7.1 million reads.


```r
## counts
data <- read.table(paste0(dir, "RNA-seq/data/geneCounts_CONTROLlibraries.RAW.tsv"))
stopifnot(identical(meta$library, colnames(data)[-1]))

## mapping statistics (from STAR logs)
mapping.stats <- read.table(paste0(dir, "RNA-seq/data/mappingStatistics_CONTROLlibraries.tsv"), stringsAsFactors = FALSE, header = TRUE)
stopifnot(identical(colnames(data)[-1], mapping.stats$sample))

## counting statistics (from STAR gene counts)
counting.stats <- read.table(paste0(dir, "RNA-seq/data/countingStatistics_CONTROLlibraries.tsv"), stringsAsFactors = FALSE)
stopifnot(identical(colnames(data)[-1], row.names(counting.stats)))
counting.stats$uniquelyMapped <- colSums(data[,-1])

ggplot(mapping.stats, aes(1, total/1e6)) + geom_violin() + geom_dotplot(binaxis='y', stackdir='center', dotsize=0.5) + ylab("library size (million fragments)") + xlab("") + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
```

![](00_controlExperiments_files/figure-html/data-1.png)<!-- -->

From these, a median of 75.37% were uniquely mapped, plus 23.56% multimapped, which is expected for 50bp single-end reads. Most of the uniquely mapped reads could be unambiguously assigned to annotated exons (85.44%).


```r
# summary(mapping.stats$unique/mapping.stats$total*100)
# summary(mapping.stats$multimapped/mapping.stats$total*100)

ggplot(counting.stats, aes(1, uniquelyMapped/1e6)) + geom_violin() + geom_dotplot(binaxis='y', stackdir='center', dotsize=0.5) + ylab("total fragments uniquely mapped to exons (millions)") + xlab("") + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
```

![](00_controlExperiments_files/figure-html/inExons-1.png)<!-- -->

```r
# summary(counting.stats$uniquelyMapped/mapping.stats$unique*100)
```

All samples but one had very uniform mapping statistics. The outlier (`do16440`) has higher proportion of multimapped reads (30.2%) which results in a lower proportion of uniquely mapped reads (67.5%). But from those, the proportion within exons is equivalent to other samples.

Similarly, all samples have a very uniform number of genes detected, around 22 thousand,


```r
counting.stats$nGenes <- apply(data[,-1], 2, function(x) sum(x>0))
# summary(counting.stats$nGenes)

ggplot(counting.stats, aes(1, nGenes/1e3)) + geom_violin(trim=FALSE) + geom_jitter(aes(1,nGenes/1e3, col=mapping.stats$total/1e6), width = 0.02) + ylab("total genes detected (thousands)") + scale_colour_gradientn(colours = (brewer.pal(n = 11, name = "PiYG"))) + labs(colour="libsize") + xlab("") + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
```

![](00_controlExperiments_files/figure-html/numGenes-1.png)<!-- -->

Based on this, all samples are deemed of good quality.

### Normalisation

We normalise to account for differences in sequencing depth, using `edgeR`'s method.
We filter out lowly expressed genes:


```r
## create an edgeR object
meta$group <- paste(meta$embryo, meta$somite, sep=".")
y <- DGEList(counts=data[,-1], samples = meta, genes = data[,1], group = meta$group)

## filter low abundance genes
means <- aveLogCPM(y)
keep <- filterByExpr(y)

plot(density(means[keep]), lwd=2, xlab="average log counts-per-million", main="", bty="l", ylim=c(0,1), xlim=c(-5,15))
lines(density(means[-keep]), lty=2, lwd=2)
legend("topright", legend = c("kept", "filtered"), lty=c(1,2), lwd=2)
```

![](00_controlExperiments_files/figure-html/filter-1.png)<!-- -->

```r
summary(keep)
```

```
##    Mode   FALSE    TRUE 
## logical   38692   16844
```

```r
y <- y[keep, , keep.lib.sizes=FALSE]
```

And estimate normalisation factors, which successfully unify the samples.


```r
## normalisation
y <- calcNormFactors(y)

dataNorm <- cpm(y, log=TRUE, prior.count = 1)
## add gene names
dataNorm <- cbind(data[row.names(dataNorm),1], as.data.frame(dataNorm))
colnames(dataNorm)[1] <- "gene"
write.table(dataNorm, paste0(dir, "RNA-seq/data/geneCounts_CONTROLlibraries.NORM_logCPM.tsv"), quote = FALSE, sep="\t")

par(mfrow=c(1,2))
boxplot(log2(data[keep,-1]+1), main="RAW", ylab=expression('log'[2]*' counts + 1'), axes=FALSE); box(bty="l"); axis(2, las=2)
boxplot(dataNorm[,-1], main="NORMALISED", ylab=expression('log'[2]*' counts-per-million + 1'), axes=FALSE); box(bty="l"); axis(2, las=2)
```

![](00_controlExperiments_files/figure-html/norm-1.png)<!-- -->

A PCA on the normalised expression of the 1000 most variable genes reveals that the sample that was an outlier on mapping statistics also is an outlier in the PCA. Embryo of origin is the next largest source of variation.


```r
data.vst <- vst(as.matrix(data[keep,-1]))

vars <- rowVars(data.vst)
names(vars) <- row.names(data.vst)
pca <- prcomp(t(dataNorm[names(vars[order(vars, decreasing = TRUE)])[1:1000],-1]))

eigs <- pca$sdev^2
prop.var <- eigs / sum(eigs)
# plot(prop.var, ylab="proportion of variance explained", xlab="PC", bty="l")

df <- as.data.frame(pca$x)
ggplot(df, aes(PC1, PC2)) + geom_point(aes(col=as.factor(meta$embryo), shape=meta$somite)) + labs(col="embryo", shape="somite") + xlab(paste("PC1 - ", round(prop.var[1]*100,2), "% var explained")) + ylab(paste("PC2 - ", round(prop.var[2]*100,2), "% var explained"))
```

![](00_controlExperiments_files/figure-html/PCA-1.png)<!-- -->

Based on this, we remove this sample from downstream analyses.


```r
## remove aoutlier
outlier <- mapping.stats$sample[which.min(mapping.stats$unique/mapping.stats$total)]
meta <- meta[-which(meta$library==outlier),]
dataNorm <- dataNorm[,-which(colnames(dataNorm)==outlier)]

pca <- prcomp(t(dataNorm[names(vars[order(vars, decreasing = TRUE)])[1:1000],-1]))

eigs <- pca$sdev^2
prop.var <- eigs / sum(eigs)

df <- as.data.frame(pca$x)
ggplot(df, aes(PC1, PC2)) + geom_point(aes(col=as.factor(meta$embryo), shape=meta$somite)) + labs(col="embryo", shape="somite") + xlab(paste("PC1 - ", round(prop.var[1]*100,2), "% var explained")) + ylab(paste("PC2 - ", round(prop.var[2]*100,2), "% var explained"))
```

![](00_controlExperiments_files/figure-html/PCA_noOut-1.png)<!-- -->

### Differential expression analysis

We can use differential expression analysis to test whether the two somites from each pair have significant differences in gene expression. We use `edgeR` to estimate the dispersion.


```r
## design
design <- model.matrix(~0+somite+side, meta)

## remove outlier
y$samples <- y$samples[-which(y$samples$library == outlier),]
y$counts <- y$counts[,-which(colnames(y$counts)==outlier)]

## dispersion
y <- estimateDisp(y, design, robust = TRUE)
par(mfrow=c(1,2))
plotBCV(y)

## fit model
fit <- glmQLFit(y, design, robust=TRUE)
plotQLDisp(fit)
```

![](00_controlExperiments_files/figure-html/edger-1.png)<!-- -->

```r
# summary(fit$df.prior)
```

And we test for differences between side 1 and 2, while controlling for which somite was used. No significant genes are identified.


```r
## test
side <- glmQLFTest(fit, ncol(design))

## results
side <- as.data.frame(topTags(side, n=Inf))
# sum(side$FDR < 0.05)  # 0
write.table(side, paste0(dir, "RNA-seq/results/00_DEresults_side_allSomites.tsv"), quote = FALSE, sep="\t")

plot(side$logCPM, side$logFC, pch=16, cex=0.7, xlab="log2 mean expression", ylab="log2 fold-change", col=ifelse(side$FDR < 0.05, "red", "black"))
abline(h=0, col="red", lwd=2)
legend("topright", legend = c("DE", "non-DE"), pch=16, col=c("red", "black"))
```

![](00_controlExperiments_files/figure-html/somiteTrios_perStage-1.png)<!-- -->

Alternatively, we can test separately the somiteI and somiteII pairs, using the two embryos as biological replicates (for somite II we only have one replicate for one of the sides). Again, no significant genes are identified.


```r
## somite I only
y.SI <- DGEList(counts=data[,c(2:3,6:7)], genes=data[,1], samples=meta[c(1:2,5:6),], group=meta$group[c(1:2,5:6)])
y.SI <- y.SI[keep, , keep.lib.sizes=FALSE]
y.SI <- calcNormFactors(y.SI)

design.SI <- model.matrix(~side, meta[c(1:2,5:6),])
y.SI <- estimateDisp(y.SI, design.SI, robust = TRUE)
fit.SI <- glmQLFit(y.SI, design.SI, robust=TRUE)
side.SI <- glmQLFTest(fit.SI, ncol(design.SI))
side.SI <- as.data.frame(topTags(side.SI, n=Inf))

par(mfrow=c(1,2))
plot(side.SI$logCPM, side.SI$logFC, pch=16, cex=0.7, xlab="log2 mean expression", ylab="log2 fold-change", col=ifelse(side.SI$FDR < 0.05, "red", "black"), main="somite I")
abline(h=0, col="red", lwd=2)
legend("topright", legend = c("DE", "non-DE"), pch=16, col=c("red", "black"))

## somite II only
y.SII <- DGEList(counts=data[,c(4:5,8)], genes=data[,1], samples=meta[c(3:4,7),], group=meta$group[c(3:4,7)])
y.SII <- y.SII[keep, , keep.lib.sizes=FALSE]
y.SII <- calcNormFactors(y.SII)

design.SII <- model.matrix(~side, meta[c(3:4,7),])
y.SII <- estimateDisp(y.SII, design.SII, robust = TRUE)
fit.SII <- glmQLFit(y.SII, design.SII, robust=TRUE)
side.SII <- glmQLFTest(fit.SII, ncol(design.SII))
side.SII <- as.data.frame(topTags(side.SII, n=Inf))

plot(side.SII$logCPM, side.SII$logFC, pch=16, cex=0.7, xlab="log2 mean expression", ylab="log2 fold-change", col=ifelse(side.SII$FDR < 0.05, "red", "black"), main="somite II")
abline(h=0, col="red", lwd=2)
legend("topright", legend = c("DE", "non-DE"), pch=16, col=c("red", "black"))
```

![](00_controlExperiments_files/figure-html/perSomite-1.png)<!-- -->

Even though this experiment is not very well powered, any systematic clear differences between the two somites should be picked up. Thus, we can use the two somites as equivalent samples to perform matched RNA and ATAC experiments.

### Differences between reagent concentration

For one somite pair, the library for sequencing was prepared twice, using a quarter of the recommended reagents. Since there are no significant differences between the left and right somites from the same pair, we can use them as replicates to compare the reagent concentrations. It's not ideal but good enough for a sanity check. If reducing the reagents led to systematic differences in a group of genes (e.g., those expressed at low levels, or with high/low GC content) we should see them.

We don't. No genes are significantly DE between the full and quartered reagents libraries.


```r
y.reagents <- DGEList(counts=data[,c(4:5,9:10)], genes=data[,1], samples=meta[c(3:4,8:9),], group=meta$group[c(3:4,8:9)])
y.reagents <- y.reagents[keep, , keep.lib.sizes=FALSE]
y.reagents <- calcNormFactors(y.reagents)

design.reagents <- model.matrix(~c(0,0,1,1))
y.reagents <- estimateDisp(y.reagents, design.reagents, robust = TRUE)
fit.reagents <- glmQLFit(y.reagents, design.reagents, robust=TRUE)
reagents <- glmQLFTest(fit.reagents, ncol(design.reagents))
reagents <- as.data.frame(topTags(reagents, n=Inf))

plot(reagents$logCPM, reagents$logFC, pch=16, cex=0.7, xlab="log2 mean expression", ylab="log2 fold-change", col=ifelse(side.SI$FDR < 0.05, "red", "black"))
abline(h=0, col="red", lwd=2)
legend("topright", legend = c("DE", "non-DE"), pch=16, col=c("red", "black"))
```

![](00_controlExperiments_files/figure-html/reagents-1.png)<!-- -->

And indeed the expression estimates are pretty well matched, as expected for technical replicates.


```r
par(mfrow=c(1,2))
plot(dataNorm[,'do16435'], dataNorm[,'do16441'], pch=16, cex=0.5, main="e1-SII-side1", xlab="full", ylab="quarter", bty="l")
abline(0,1,col="red", lwd=2)
plot(dataNorm[,'do16436'], dataNorm[,'do16442'], pch=16, cex=0.5, main="e1-SII-side2", xlab="full", ylab="quarter", bty="l")
abline(0,1,col="red", lwd=2)
```

![](00_controlExperiments_files/figure-html/scatter-1.png)<!-- -->

### Conclusions

These control experiments show that it is valid to consider the left and right somites from the same pair as equivalent samples, with the same transcriptome. Also, that using a lower amount of reagents in the library prep step doesn't affect the library produced. 

This is illustrated below, with the correlation coefficients between the transcriptomes of all the samples. The technical replicates -that only differ in the amount of reagents used- have the highest coefficients, very close to 1; and these are significantly higher than the correlations between any other two samples (wilcoxon rank sum test p-value = 0.003175). 

The samples coming from the same somite-pair, that only differ on whether they are from the left or right side are the next highest correlations. However, these correlation can be as *low* as those between adjacent somites. Overall, the correlation coefficients decrease as the number of variables that are different between samples increases, as expected.


```r
corr <- cor(dataNorm[,-1], method="spearman")
# corr <- cor(dataNorm[names(vars[order(vars, decreasing = TRUE)])[1:1000],-1])

## group the different correlations
df <- data.frame(corr=unique(as.vector(corr))[-1])
## by the similarity of the two samples involved
df$group <- c("side","somite","somite","embryo","embryo-side","embryo-somite","somite","somite","somite","somite","embryo-side","embryo","embryo-somite","somite","somite","side","embryo-somite","embryo-somite","embryo","replicate","side","embryo-somite","embryo-somite","embryo","side","replicate","side","somite","embryo-somite","embryo-somite","somite","embryo-somite","embryo-somite","embryo","embryo-side","side")
## technical 'replicates' are the most similar
## followed by samples that only differ by 'side'
## then those that differ by the 'somite' taken (both sides considered equal)
## then those from a different 'embryo' but same somite and side
## versus 'embryo-side' where embryo and side are different
##lastly, 'embryo-somite' where both embryo and somite are different
df$group <- factor(as.character(df$group), levels=c("replicate", "side", "somite", "embryo", "embryo-side", "embryo-somite"))
ggplot(df, aes(group, corr)) + geom_violin() + geom_dotplot(binaxis='y', stackdir='center', dotsize=0.5)
```

![](00_controlExperiments_files/figure-html/unnamed-chunk-1-1.png)<!-- -->

```r
# wilcox.test(df[df$group=="replicate",]$corr, df[df$group!="replicate",]$corr) # 0.003175

# replicate corrs:  0.9856984 0.9805222 pearson
#                   0.9900903 0.9878620 spearman
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
##  [4] assertthat_0.2.1       statmod_1.4.32         latticeExtra_0.6-28   
##  [7] blob_1.2.0             GenomeInfoDbData_1.2.2 yaml_2.2.0            
## [10] RSQLite_2.1.3          pillar_1.4.2           backports_1.1.5       
## [13] lattice_0.20-38        glue_1.3.1             digest_0.6.23         
## [16] ggsignif_0.6.0         XVector_0.26.0         checkmate_1.9.4       
## [19] colorspace_1.4-1       htmltools_0.4.0        Matrix_1.2-18         
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

