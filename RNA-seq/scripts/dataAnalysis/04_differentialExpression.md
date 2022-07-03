---
title: "Differential expression analysis of mouse somites"
date: '26 September, 2020'
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





```r
meta <- read.table(paste0(dir, "RNA-seq/data/metadata_RNAseq.tsv"), stringsAsFactors = FALSE, header = TRUE)
meta.all <- meta # save complete metadata just in case

data <- read.table(paste0(dir, "RNA-seq/data/geneCounts.RAW.tsv"), check.names = FALSE)
stopifnot(identical(colnames(data)[-1], meta$sample))
y <- readRDS(paste0(dir, "RNA-seq/results/01_edgeRobject.Rds"))

## remove bad-quality and wrong-stage samples
meta <- meta[meta$QC == 1 & meta$wrongStage == 0,]

y$samples <- y$samples[which(y$samples$sample %in% meta$sample),]
y$counts <- y$counts[,which(colnames(y$counts) %in% meta$sample)]
stopifnot(identical(row.names(y$samples), colnames(y$counts)))
stopifnot(identical(colnames(y$counts), meta$sample))
```

To define the transcriptional changes that accompany somite maturation along development, we use differential expression analysis. To control for the batch effects observed in the data, we use the PCs computed in `02_removeTechnicalVar.Rmd` as covariates in the linear model.

We use `edgeR` to perform the differential expression (DE) testing. We first estimate the dispersion and fit the design containing the interaction of somite and stage.


```r
## design
meta$group <- as.factor(meta$group)
design <- model.matrix(~0+group, meta)
colnames(design) <- paste0("stage",levels(meta$group))

# add covariates
pcs <- read.table(paste0(dir, "RNA-seq/results/02_pcs_residuals.tab"))
# remove wrong stage samples
pcs <- pcs[row.names(pcs) %in% meta$sample,]
# add first 14 PCs to design
design <- cbind(design, pcs[,1:14])
colnames(design)[-c(1:length(levels(meta$group)))] <- paste0("PC", 1:14)

## dispersion
y <- estimateDisp(y, design, robust = TRUE)
par(mfrow=c(1,2))
plotBCV(y)

## fit model
fit <- glmQLFit(y, design, robust=TRUE)
plotQLDisp(fit)
```

![](04_differentialExpression_files/figure-html/edger-1.png)<!-- -->

```r
# summary(fit$df.prior)

saveRDS(y, file=paste0(dir, "RNA-seq/results/04_edgeRobject_pca.Rds"))
```

### Differences across somite trios

To test for changes in transcription as somites mature, we compare somites I, II and III. There are at least 3 replicates for each stage-somite combination.


```r
table(meta$stage, meta$somite)
```

```
##     
##      SI SII SIII
##   8   4   4    4
##   18  5   5    5
##   21  4   3    4
##   25  4   3    4
##   27  3   3    3
##   35  6   4    6
```

The changes could be conserved across development, or stage specific.

First, we test on a per-stage basis, testing all three pairwise comparisons from each stage at once with an anova-like approach. Therefore, an adjusted p-value lower than 0.05 indicates that the gene is significantly DE between at least a pair of somites. We further require an absolute fold-change of 1.5 or larger to consider a gene significantly DE. Several hundred genes are identified for each stage.


```r
my.contrasts <- makeContrasts(stage8.IIvsI    = stage8.SII-stage8.SI,
                              stage8.IIIvsI   = stage8.SIII-stage8.SI, 
                              stage8.IIIvsII  = stage8.SIII-stage8.SII, 
                              stage18.IIvsI   = stage18.SII-stage18.SI,
                              stage18.IIIvsI  = stage18.SIII-stage18.SI, 
                              stage18.IIIvsII = stage18.SIII-stage18.SII, 
                              stage21.IIvsI   = stage21.SII-stage21.SI,
                              stage21.IIIvsI  = stage21.SIII-stage21.SI, 
                              stage21.IIIvsII = stage21.SIII-stage21.SII, 
                              stage25.IIvsI   = stage25.SII-stage25.SI,
                              stage25.IIIvsI  = stage25.SIII-stage25.SI, 
                              stage25.IIIvsII = stage25.SIII-stage25.SII, 
                              stage27.IIvsI   = stage27.SII-stage27.SI,
                              stage27.IIIvsI  = stage27.SIII-stage27.SI, 
                              stage27.IIIvsII = stage27.SIII-stage27.SII, 
                              stage35.IIvsI   = stage35.SII-stage35.SI,
                              stage35.IIIvsI  = stage35.SIII-stage35.SI, 
                              stage35.IIIvsII = stage35.SIII-stage35.SII,
                              levels=design)

## test
somite.stage8 <- glmQLFTest(fit, contrast = my.contrasts[,1:3])
somite.stage18 <- glmQLFTest(fit, contrast = my.contrasts[,4:6])
somite.stage21 <- glmQLFTest(fit, contrast = my.contrasts[,7:9])
somite.stage25 <- glmQLFTest(fit, contrast = my.contrasts[,10:12])
somite.stage27 <- glmQLFTest(fit, contrast = my.contrasts[,13:15])
somite.stage35 <- glmQLFTest(fit, contrast = my.contrasts[,16:18])

## results
somite.stage8 <- as.data.frame(topTags(somite.stage8, n=Inf))
somite.stage8$maxFC <- apply(somite.stage8[,2:4], 1, function(x) x[which.max(abs(x))])
somite.stage18 <- as.data.frame(topTags(somite.stage18, n=Inf))
somite.stage18$maxFC <- apply(somite.stage18[,2:4], 1, function(x) x[which.max(abs(x))])
somite.stage21 <- as.data.frame(topTags(somite.stage21, n=Inf))
somite.stage21$maxFC <- apply(somite.stage21[,2:4], 1, function(x) x[which.max(abs(x))])
somite.stage25 <- as.data.frame(topTags(somite.stage25, n=Inf))
somite.stage25$maxFC <- apply(somite.stage25[,2:4], 1, function(x) x[which.max(abs(x))])
somite.stage27 <- as.data.frame(topTags(somite.stage27, n=Inf))
somite.stage27$maxFC <- apply(somite.stage27[,2:4], 1, function(x) x[which.max(abs(x))])
somite.stage35 <- as.data.frame(topTags(somite.stage35, n=Inf))
somite.stage35$maxFC <- apply(somite.stage35[,2:4], 1, function(x) x[which.max(abs(x))])

## define DE genes
fdr_thr <- 0.05
fc_thr <- 1.5

somite.stage8$DE = ifelse(somite.stage8$FDR < fdr_thr & abs(somite.stage8$maxFC) > log2(fc_thr), 1, 0)
somite.stage18$DE = ifelse(somite.stage18$FDR < fdr_thr & abs(somite.stage18$maxFC) > log2(fc_thr), 1, 0)
somite.stage21$DE = ifelse(somite.stage21$FDR < fdr_thr & abs(somite.stage21$maxFC) > log2(fc_thr), 1, 0)
somite.stage25$DE = ifelse(somite.stage25$FDR < fdr_thr & abs(somite.stage25$maxFC) > log2(fc_thr), 1, 0)
somite.stage27$DE = ifelse(somite.stage27$FDR < fdr_thr & abs(somite.stage27$maxFC) > log2(fc_thr), 1, 0)
somite.stage35$DE = ifelse(somite.stage35$FDR < fdr_thr & abs(somite.stage35$maxFC) > log2(fc_thr), 1, 0)

nDE <- t(data.frame(stage8 = sum(somite.stage8$DE),
                    stage18 = sum(somite.stage18$DE),
                    stage21 = sum(somite.stage21$DE),
                    stage25 = sum(somite.stage25$DE),
                    stage27 = sum(somite.stage27$DE),
                    stage35 = sum(somite.stage35$DE)))
colnames(nDE) <- "numberDE"
nDE
```

```
##         numberDE
## stage8       439
## stage18      767
## stage21      467
## stage25      152
## stage27      420
## stage35      862
```

The vast majority of genes are deemed significant in only one stage.


```r
somite.perStage <- data.frame(stage8 = ifelse(somite.stage8$DE == 1, 1, 0), 
                        stage18 = ifelse(somite.stage18[row.names(somite.stage8),]$DE == 1, 1, 0),
                        stage21 = ifelse(somite.stage21[row.names(somite.stage8),]$DE == 1, 1, 0),
                        stage25 = ifelse(somite.stage25[row.names(somite.stage8),]$DE == 1, 1, 0),
                        stage27 = ifelse(somite.stage27[row.names(somite.stage8),]$DE == 1, 1, 0),
                        stage35 = ifelse(somite.stage35[row.names(somite.stage8),]$DE == 1, 1, 0))
row.names(somite.perStage) <- row.names(somite.stage8)
upset(somite.perStage, nsets=6, nintersects = 100, sets.x.label = "number DE genes per stage")
```

![](04_differentialExpression_files/figure-html/overlapPerStage-1.png)<!-- -->

This could indicate that many expression changes are dependent on the developmental stage of the embryo. Alternatively, it could reflect a lack of power, with genes reaching significance in only some comparisons. 

We also test for differences that are conserved across stages, using somites from different stages but the same maturation level as replicates; this should increase power by having many more replicates. This test would identify changes that are conserved across development.

As expected, somites I vs III show the largest number of DE genes, as they should be the most different.


```r
my.contrasts <- makeContrasts(
  somiteIIvsI = (stage8.SII + stage18.SII + stage21.SII + stage25.SII + stage27.SII + stage35.SII)/6 - (stage8.SI + stage18.SI + stage21.SI + stage25.SI + stage27.SI + stage35.SI)/6,
  somiteIIIvsI = (stage8.SIII + stage18.SIII + stage21.SIII + stage25.SIII + stage27.SIII + stage35.SIII)/6 - (stage8.SI + stage18.SI + stage21.SI + stage25.SI + stage27.SI + stage35.SI)/6,
  somiteIIIvsII = (stage8.SIII + stage18.SIII + stage21.SIII + stage25.SIII + stage27.SIII + stage35.SIII)/6 - (stage8.SII + stage18.SII + stage21.SII + stage25.SII + stage27.SII + stage35.SII)/6,
  levels=design)

## test
somiteIIvsI.all <- glmQLFTest(fit, contrast = my.contrasts[,"somiteIIvsI"])
somiteIIIvsI.all <- glmQLFTest(fit, contrast = my.contrasts[,"somiteIIIvsI"])
somiteIIIvsII.all <- glmQLFTest(fit, contrast = my.contrasts[,"somiteIIIvsII"])

## results
somiteIIvsI.all <- as.data.frame(topTags(somiteIIvsI.all, n=Inf))
somiteIIIvsI.all <- as.data.frame(topTags(somiteIIIvsI.all, n=Inf))
somiteIIIvsII.all <- as.data.frame(topTags(somiteIIIvsII.all, n=Inf))

## define DE genes
somiteIIvsI.all$DE = ifelse(somiteIIvsI.all$FDR < fdr_thr & somiteIIvsI.all$logFC > log2(fc_thr), "up", 
                            ifelse(somiteIIvsI.all$FDR < fdr_thr & somiteIIvsI.all$logFC < -log2(fc_thr), "down", "NS"))
somiteIIIvsI.all$DE = ifelse(somiteIIIvsI.all$FDR < fdr_thr & somiteIIIvsI.all$logFC > log2(fc_thr), "up", 
                            ifelse(somiteIIIvsI.all$FDR < fdr_thr & somiteIIIvsI.all$logFC < -log2(fc_thr), "down", "NS"))
somiteIIIvsII.all$DE = ifelse(somiteIIIvsII.all$FDR < fdr_thr & somiteIIIvsII.all$logFC > log2(fc_thr), "up", 
                            ifelse(somiteIIIvsII.all$FDR < fdr_thr & somiteIIIvsII.all$logFC < -log2(fc_thr), "down", "NS"))

nDE <- t(data.frame(IvsII = sum(somiteIIvsI.all$DE != "NS"),
                    IIvsIII = sum(somiteIIIvsII.all$DE != "NS"),
                    IvsIII = sum(somiteIIIvsI.all$DE != "NS")))
colnames(nDE) <- "numberDE"
nDE
```

```
##         numberDE
## IvsII        275
## IIvsIII      404
## IvsIII       989
```


```r
plots <- list()

tmp <- somiteIIvsI.all[order(somiteIIvsI.all$DE),]
tmp <- tmp[order(tmp$FDR),]
tmp <- rbind(tmp[tmp$DE == "NS",], tmp[tmp$DE != "NS",])
plots[[1]] <- ggplot(tmp, aes(logCPM, logFC, col=as.factor(DE), label=genes)) + geom_point(size=1) +
  geom_label_repel(data = tmp[tmp$DE != "NS",][1:20,], size = 3, show.legend = FALSE) + 
  ylim(c(-4,4)) + labs(col="") + 
  scale_color_manual(values = c("dodgerblue4", "grey", "firebrick3")) + 
  geom_hline(yintercept = seq(-4,4,2), lty = 2, colour="lightgrey") + geom_hline(yintercept = 0) + 
  xlab(expression('log'[2]*' CPM')) + ylab(expression('log'[2]*' fold-change')) + 
  ggtitle("somite II vs I") + 
  th + guides(colour = guide_legend(reverse = TRUE)) + theme(legend.position = "bottom")

tmp <- somiteIIIvsII.all[order(somiteIIIvsII.all$DE),]
tmp <- tmp[order(tmp$FDR),]
tmp <- rbind(tmp[tmp$DE == "NS",], tmp[tmp$DE != "NS",])
plots[[2]] <- ggplot(tmp, aes(logCPM, logFC, col=as.factor(DE), label=genes)) + geom_point(size=1) +
  geom_label_repel(data = tmp[tmp$DE != "NS",][1:20,], size = 3, show.legend = FALSE) + 
  ylim(c(-4,4)) + labs(col="") + 
  scale_color_manual(values = c("dodgerblue4", "grey", "firebrick3")) + 
  geom_hline(yintercept = seq(-4,4,2), lty = 2, colour="lightgrey") + geom_hline(yintercept = 0) + 
  xlab(expression('log'[2]*' CPM')) + ylab(expression('log'[2]*' fold-change')) + 
  ggtitle("somite III vs II") + 
  th + guides(colour = guide_legend(reverse = TRUE)) + theme(legend.position = "bottom")

tmp <- somiteIIIvsI.all[order(somiteIIIvsI.all$DE),]
tmp <- tmp[order(tmp$FDR),]
tmp <- rbind(tmp[tmp$DE == "NS",], tmp[tmp$DE != "NS",])
plots[[3]] <- ggplot(tmp, aes(logCPM, logFC, col=as.factor(DE), label=genes)) + geom_point(size=1) +
  geom_label_repel(data = tmp[tmp$DE != "NS",][1:20,], size = 3, show.legend = FALSE) + 
  ylim(c(-4,4)) + labs(col="") + 
  scale_color_manual(values = c("dodgerblue4", "grey", "firebrick3")) + 
  geom_hline(yintercept = seq(-4,4,2), lty = 2, colour="lightgrey") + geom_hline(yintercept = 0) + 
  xlab(expression('log'[2]*' CPM')) + ylab(expression('log'[2]*' fold-change')) + 
  ggtitle("somite III vs I") + 
  th + guides(colour = guide_legend(reverse = TRUE)) + theme(legend.position = "bottom")

ggarrange(plotlist = plots, ncol = 3, common.legend = TRUE, legend = "bottom")
```

![](04_differentialExpression_files/figure-html/somiteMAplots-1.png)<!-- -->

The vast majority of the genes significantly different between somites I and II are also significant between somites I and III, suggesting monotonic changes that increase in magnitude as somites mature. This is also true for most DE genes between somites II and III. However, a large number of genes only reach significance when comparing somites I and III, suggesting that many of the changes are small; this is supported by the MA plots above.


```r
somite.average <- data.frame(IvsII   = ifelse(somiteIIvsI.all$DE != "NS", 1, 0), 
                             IIvsIII = ifelse(somiteIIIvsII.all[row.names(somiteIIvsI.all),]$DE != "NS", 1, 0),
                             IvsIII  = ifelse(somiteIIIvsI.all[row.names(somiteIIvsI.all),]$DE != "NS", 1, 0))
row.names(somite.average) <- row.names(somiteIIvsI.all)
upset(somite.average, sets.x.label = "number DE genes per stage")
```

![](04_differentialExpression_files/figure-html/overlapSomiteAve-1.png)<!-- -->

Finally, to combine all results, we define changes as stage-specific if they were identified in one of the per-stage changes, but not in the averaged test. Overall there are a couple thousand significantly different genes between the somite trios. 


```r
somiteTrios <- cbind(somite.average, somite.perStage[row.names(somite.average),])
somiteTrios$ave <- ifelse(rowSums(somiteTrios[,1:3])>0, 1, 0)
somiteTrios$stageSpecific <- ifelse(somiteTrios$ave == 0 & rowSums(somiteTrios[,-c(1:3,ncol(somiteTrios))])>0, 1, 0)

colSums(somiteTrios[,10:11])
```

```
##           ave stageSpecific 
##          1235          1742
```

```r
colSums(somiteTrios[somiteTrios$ave == 0 & somiteTrios$stageSpecific == 1,4:9])
```

```
##  stage8 stage18 stage21 stage25 stage27 stage35 
##     259     513     295      56     358     554
```

And most have small fold-changes.


```r
## capture the largest FC for each DE gene
somiteTrios.fc <- somiteTrios
somiteTrios.fc.de <- somiteTrios.fc[somiteTrios.fc$ave==1 | somiteTrios.fc$stageSpecific==1,]
somiteTrios.fc.nonde <- somiteTrios.fc[setdiff(row.names(somiteTrios.fc), row.names(somiteTrios.fc.de)),]

## DE genes
somiteTrios.fc.de$IvsII <- somiteIIvsI.all[row.names(somiteTrios.fc.de),]$logFC
somiteTrios.fc.de$IIvsIII <- somiteIIIvsII.all[row.names(somiteTrios.fc.de),]$logFC
somiteTrios.fc.de$IvsIII <- somiteIIIvsI.all[row.names(somiteTrios.fc.de),]$logFC
somiteTrios.fc.de$stage8 <- somite.stage8[row.names(somiteTrios.fc.de),]$maxFC
somiteTrios.fc.de$stage18 <- somite.stage18[row.names(somiteTrios.fc.de),]$maxFC
somiteTrios.fc.de$stage21 <- somite.stage21[row.names(somiteTrios.fc.de),]$maxFC
somiteTrios.fc.de$stage25 <- somite.stage25[row.names(somiteTrios.fc.de),]$maxFC
somiteTrios.fc.de$stage27 <- somite.stage27[row.names(somiteTrios.fc.de),]$maxFC
somiteTrios.fc.de$stage35 <- somite.stage35[row.names(somiteTrios.fc.de),]$maxFC

## non-DE genes
somiteTrios.fc.nonde$IvsII <- somiteIIvsI.all[row.names(somiteTrios.fc.nonde),]$logFC
somiteTrios.fc.nonde$IIvsIII <- somiteIIIvsII.all[row.names(somiteTrios.fc.nonde),]$logFC
somiteTrios.fc.nonde$IvsIII <- somiteIIIvsI.all[row.names(somiteTrios.fc.nonde),]$logFC
somiteTrios.fc.nonde$stage8 <- somite.stage8[row.names(somiteTrios.fc.nonde),]$maxFC
somiteTrios.fc.nonde$stage18 <- somite.stage18[row.names(somiteTrios.fc.nonde),]$maxFC
somiteTrios.fc.nonde$stage21 <- somite.stage21[row.names(somiteTrios.fc.nonde),]$maxFC
somiteTrios.fc.nonde$stage25 <- somite.stage25[row.names(somiteTrios.fc.nonde),]$maxFC
somiteTrios.fc.nonde$stage27 <- somite.stage27[row.names(somiteTrios.fc.nonde),]$maxFC
somiteTrios.fc.nonde$stage35 <- somite.stage35[row.names(somiteTrios.fc.nonde),]$maxFC

## use same rules to define average and stage-specific
somiteTrios.fc.de[somiteTrios.fc.de$ave==1,]$ave <- apply(somiteTrios.fc.de[somiteTrios.fc.de$ave==1,1:3], 1,
                                                          function(x) x[which.max(abs(x))])
somiteTrios.fc.de[somiteTrios.fc.de$stageSpecific==1,]$stageSpecific <- apply(
  somiteTrios.fc.de[somiteTrios.fc.de$stageSpecific==1,4:9], 1, function(x) x[which.max(abs(x))])

somiteTrios.fc.de$logFC <- somiteTrios.fc.de$ave + somiteTrios.fc.de$stageSpecific

## for non-DE genes, use the median fold-change for all tests
somiteTrios.fc.nonde$logFC <- apply(somiteTrios.fc.nonde[,1:9], 1, median)

## MA plot
# for DE genes, cap FC at |7|
df <- data.frame(m = somiteIIIvsI.all[row.names(somiteTrios.fc.de),]$logCPM, a = somiteTrios.fc.de$logFC)
df$pch <- 16
df[df$a > 7,]$pch <- 2
df[df$a > 7,]$a <- 7
df[df$a < -7,]$pch <- 6
df[df$a < -7,]$a <- -7

plot(somiteIIIvsI.all[row.names(somiteTrios.fc.nonde),]$logCPM, 
     somiteTrios.fc.nonde$logFC, 
     pch=16, col="grey",
     xlab=expression('log'[2]*' mean expression'),
     ylab=expression('log'[2]*' fold-change'),
     bty="l", xlim=range(df$m), ylim=c(-7,7))
points(df$m, df$a, pch=df$pch, col="steelblue")
abline(h=0, col=rgb(1,0,0,0.5), lwd=2)
abline(h=c(-log2(1.5), log2(1.5)), lty=2, col="grey20")
```

![](04_differentialExpression_files/figure-html/fold-change-1.png)<!-- -->

```r
# sum(abs(somiteTrios.fc.de$logFC)<1)/nrow(somiteTrios.fc.de)*100
# 49.6473
```


#### Stage-specific changes {.tabset}

When examining more closely the expression patterns of the genes that are significant in only one of the stage-specific tests, most genes show a consistent expression pattern across stages. However, the differences tend to be small, and thus only reach significance occasionally. This suggests that there are few genes that change in expression between maturing somites in only some stages. Generally, DE genes between somite trios show the same pattern of expression across development.


```r
dataNorm <- read.table(paste0(dir, "RNA-seq/data/geneCounts.NORM_batchCorrected_14PCs.tsv"), check.names = FALSE)
dataNorm <- dataNorm[,which(colnames(dataNorm) %in% c("gene",meta$sample)),]

retrieveGeneExpr <- function(dataNorm, meta, columns=c("stage", "somite", "date"), gene){
  d <- t(as.data.frame(dataNorm[dataNorm$gene==gene,-1]))
  stopifnot(identical(row.names(d), meta$sample))
  d <- cbind(d, meta[,columns])
  colnames(d)[1] <- "count"
  return(d)
}
```


##### Stage8


```r
## select stage-specific genes
tmp <- somite.stage8[row.names(somiteTrios[somiteTrios$stage8 == 1 & rowSums(somiteTrios[,4:9]) == 1 & somiteTrios$ave == 0,]),]

plots <- list()
for(i in 1:6){
  d <- retrieveGeneExpr(dataNorm, meta, gene=tmp[i,1])
  plots[[i]] <- ggplot(d, aes(x=as.factor(stage), y=count)) + 
    geom_boxplot(aes(fill=somite)) +
    ylab("log CPM") + xlab("stage") + ggtitle(tmp[i,1]) + 
    scale_fill_brewer(palette = "Purples") + 
    th
}
ggarrange(plotlist = plots, ncol = 3, nrow = 2, common.legend = TRUE, legend = "bottom")
```

![](04_differentialExpression_files/figure-html/s8.sp, -1.png)<!-- -->

##### Stage18


```r
## select stage-specific genes
tmp <- somite.stage18[row.names(somiteTrios[somiteTrios$stage18 == 1 & rowSums(somiteTrios[,4:9]) == 1 & somiteTrios$ave == 0,]),]

plots <- list()
for(i in 1:6){
  d <- retrieveGeneExpr(dataNorm, meta, gene=tmp[i,1])
  plots[[i]] <- ggplot(d, aes(x=as.factor(stage), y=count)) + 
    geom_boxplot(aes(fill=somite)) +
    ylab("log CPM") + xlab("stage") + ggtitle(tmp[i,1]) + 
    scale_fill_brewer(palette = "Purples") + 
    th
}
ggarrange(plotlist = plots, ncol = 3, nrow = 2, common.legend = TRUE, legend = "bottom")
```

![](04_differentialExpression_files/figure-html/s18.sp, -1.png)<!-- -->

##### Stage21


```r
## select stage-specific genes
tmp <- somite.stage21[row.names(somiteTrios[somiteTrios$stage21 == 1 & rowSums(somiteTrios[,4:9]) == 1 & somiteTrios$ave == 0,]),]

plots <- list()
for(i in 1:6){
  d <- retrieveGeneExpr(dataNorm, meta, gene=tmp[i,1])
  plots[[i]] <- ggplot(d, aes(x=as.factor(stage), y=count)) + 
    geom_boxplot(aes(fill=somite)) +
    ylab("log CPM") + xlab("stage") + ggtitle(tmp[i,1]) + 
    scale_fill_brewer(palette = "Purples") + 
    th
}
ggarrange(plotlist = plots, ncol = 3, nrow = 2, common.legend = TRUE, legend = "bottom")
```

![](04_differentialExpression_files/figure-html/s21.sp, -1.png)<!-- -->

##### Stage25


```r
## select stage-specific genes
tmp <- somite.stage25[row.names(somiteTrios[somiteTrios$stage25 == 1 & rowSums(somiteTrios[,4:9]) == 1 & somiteTrios$ave == 0,]),]

plots <- list()
for(i in 1:6){
  d <- retrieveGeneExpr(dataNorm, meta, gene=tmp[i,1])
  plots[[i]] <- ggplot(d, aes(x=as.factor(stage), y=count)) + 
    geom_boxplot(aes(fill=somite)) +
    ylab("log CPM") + xlab("stage") + ggtitle(tmp[i,1]) + 
    scale_fill_brewer(palette = "Purples") + 
    th
}
ggarrange(plotlist = plots, ncol = 3, nrow = 2, common.legend = TRUE, legend = "bottom")
```

![](04_differentialExpression_files/figure-html/s25.sp, -1.png)<!-- -->

##### Stage27


```r
## select stage-specific genes
tmp <- somite.stage27[row.names(somiteTrios[somiteTrios$stage27 == 1 & rowSums(somiteTrios[,4:9]) == 1 & somiteTrios$ave == 0,]),]

plots <- list()
for(i in 1:6){
  d <- retrieveGeneExpr(dataNorm, meta, gene=tmp[i,1])
  plots[[i]] <- ggplot(d, aes(x=as.factor(stage), y=count)) + 
    geom_boxplot(aes(fill=somite)) +
    ylab("log CPM") + xlab("stage") + ggtitle(tmp[i,1]) + 
    scale_fill_brewer(palette = "Purples") + 
    th
}
ggarrange(plotlist = plots, ncol = 3, nrow = 2, common.legend = TRUE, legend = "bottom")
```

![](04_differentialExpression_files/figure-html/s27.sp, -1.png)<!-- -->

##### Stage35


```r
## select stage-specific genes
tmp <- somite.stage35[row.names(somiteTrios[somiteTrios$stage35 == 1 & rowSums(somiteTrios[,4:9]) == 1 & somiteTrios$ave == 0,]),]

plots <- list()
for(i in 1:6){
  d <- retrieveGeneExpr(dataNorm, meta, gene=tmp[i,1])
  plots[[i]] <- ggplot(d, aes(x=as.factor(stage), y=count)) + 
    geom_boxplot(aes(fill=somite)) +
    ylab("log CPM") + xlab("stage") + ggtitle(tmp[i,1]) + 
    scale_fill_brewer(palette = "Purples") + 
    th
}
ggarrange(plotlist = plots, ncol = 3, nrow = 2, common.legend = TRUE, legend = "bottom")
```

![](04_differentialExpression_files/figure-html/s35.sp, -1.png)<!-- -->


### Differences across developmental stages

Now we test for expression changes as development proceeds. First we test for somite-specific changes. All pairwise comparisons between the six stages are tested at once, to avoid over testing. As before, a gene with an FDR of 0.05 or lower indicates that the gene is DE between at least a pair of stages; and we require a minimum fold-change of 1.5 to consider the gene significantly DE. Several thousand genes are significant.


```r
my.contrasts <- makeContrasts(stage8vs18 = stage18.SI - stage8.SI, stage8vs21 = stage21.SI - stage8.SI,
                              stage8vs25 = stage25.SI - stage8.SI, stage8vs27 = stage27.SI - stage8.SI, 
                              stage8vs35 = stage35.SI - stage8.SI, 
                              stage18vs21 = stage21.SI - stage18.SI, stage18vs25 = stage25.SI - stage18.SI, 
                              stage18vs27 = stage27.SI - stage18.SI, stage18vs35 = stage35.SI - stage18.SI, 
                              stage21vs25 = stage25.SI - stage21.SI, stage21vs27 = stage27.SI - stage21.SI, 
                              stage21vs35 = stage35.SI - stage21.SI, 
                              stage25vs27 = stage27.SI - stage25.SI, stage25vs35 = stage35.SI - stage25.SI, 
                              stage27vs35 = stage35.SI - stage27.SI, levels=design)
test <- glmQLFTest(fit, contrast = my.contrasts)
stage.somiteI <- as.data.frame(topTags(test, n=Inf))
stage.somiteI$maxFC <- apply(stage.somiteI[,2:16], 1, function(x) x[which.max(abs(x))])

my.contrasts <- makeContrasts(stage8vs18 = stage18.SII - stage8.SII, stage8vs21 = stage21.SII - stage8.SII,
                              stage8vs25 = stage25.SII - stage8.SII, stage8vs27 = stage27.SII - stage8.SII, 
                              stage8vs35 = stage35.SII - stage8.SII, 
                              stage18vs21 = stage21.SII - stage18.SII, stage18vs25 = stage25.SII - stage18.SII, 
                              stage18vs27 = stage27.SII - stage18.SII, stage18vs35 = stage35.SII - stage18.SII, 
                              stage21vs25 = stage25.SII - stage21.SII, stage21vs27 = stage27.SII - stage21.SII, 
                              stage21vs35 = stage35.SII - stage21.SII, 
                              stage25vs27 = stage27.SII - stage25.SII, stage25vs35 = stage35.SII - stage25.SII, 
                              stage27vs35 = stage35.SII - stage27.SII, levels=design)
test <- glmQLFTest(fit, contrast = my.contrasts)
stage.somiteII <- as.data.frame(topTags(test, n=Inf))
stage.somiteII$maxFC <- apply(stage.somiteII[,2:16], 1, function(x) x[which.max(abs(x))])

my.contrasts <- makeContrasts(stage8vs18 = stage18.SIII - stage8.SIII, stage8vs21 = stage21.SIII - stage8.SIII,
                              stage8vs25 = stage25.SIII - stage8.SIII, stage8vs27 = stage27.SIII - stage8.SIII, 
                              stage8vs35 = stage35.SIII - stage8.SIII, 
                              stage18vs21 = stage21.SIII - stage18.SIII, stage18vs25 = stage25.SIII - stage18.SIII, 
                              stage18vs27 = stage27.SIII - stage18.SIII, stage18vs35 = stage35.SIII - stage18.SIII, 
                              stage21vs25 = stage25.SIII - stage21.SIII, stage21vs27 = stage27.SIII - stage21.SIII, 
                              stage21vs35 = stage35.SIII - stage21.SIII, 
                              stage25vs27 = stage27.SIII - stage25.SIII, stage25vs35 = stage35.SIII - stage25.SIII, 
                              stage27vs35 = stage35.SIII - stage27.SIII, levels=design)
test <- glmQLFTest(fit, contrast = my.contrasts)
stage.somiteIII <- as.data.frame(topTags(test, n=Inf))
stage.somiteIII$maxFC <- apply(stage.somiteIII[,2:16], 1, function(x) x[which.max(abs(x))])

## define DE genes
stage.somiteI$DE = ifelse(stage.somiteI$FDR < fdr_thr & abs(stage.somiteI$maxFC) > log2(fc_thr), 1, 0)
stage.somiteII$DE = ifelse(stage.somiteII$FDR < fdr_thr & abs(stage.somiteII$maxFC) > log2(fc_thr), 1, 0)
stage.somiteIII$DE = ifelse(stage.somiteIII$FDR < fdr_thr & abs(stage.somiteIII$maxFC) > log2(fc_thr), 1, 0)

nDE <- t(data.frame(somiteI = sum(stage.somiteI$DE),
                    somiteII = sum(stage.somiteII$DE),
                    somiteIII = sum(stage.somiteIII$DE)))
colnames(nDE) <- "numberDE"
nDE
```

```
##           numberDE
## somiteI       4959
## somiteII      5119
## somiteIII     6951
```

The majority of significant genes are identified in at least two of the three somites, suggesting changes are well replicated.


```r
stage.perSomite <- data.frame(somiteI = ifelse(stage.somiteI$DE == 1, 1, 0), 
                              somiteII = ifelse(stage.somiteII[row.names(stage.somiteI),]$DE == 1, 1, 0),
                              somiteIII = ifelse(stage.somiteIII[row.names(stage.somiteI),]$DE == 1, 1, 0))
row.names(stage.perSomite) <- row.names(stage.somiteI)
upset(stage.perSomite, nsets=6, nintersects = 100, sets.x.label = "number DE genes per stage")
```

![](04_differentialExpression_files/figure-html/overlapPerSomite-1.png)<!-- -->

Next, we repeat the test but using the somite trios as replicates. This increases power and therefore the number of significant genes.


```r
my.contrasts <- makeContrasts(stage8vs18 = (stage18.SIII+stage18.SII+stage18.SI)/3 - (stage8.SIII+stage8.SII+stage8.SI)/3,
                              stage8vs21 = (stage21.SIII+stage21.SII+stage21.SI)/3 - (stage8.SIII+stage8.SII+stage8.SI)/3,
                              stage8vs25 = (stage25.SIII+stage25.SII+stage25.SI)/3 - (stage8.SIII+stage8.SII+stage8.SI)/3,
                              stage8vs27 = (stage27.SIII+stage27.SII+stage27.SI)/3 - (stage8.SIII+stage8.SII+stage8.SI)/3,
                              stage8vs35 = (stage35.SIII+stage35.SII+stage35.SI)/3 - (stage8.SIII+stage8.SII+stage8.SI)/3,
                              stage18vs21 = (stage21.SIII+stage21.SII+stage21.SI)/3 - (stage18.SIII+stage18.SII+stage18.SI)/3,
                              stage18vs25 = (stage25.SIII+stage25.SII+stage25.SI)/3 - (stage18.SIII+stage18.SII+stage18.SI)/3,
                              stage18vs27 = (stage27.SIII+stage27.SII+stage27.SI)/3 - (stage18.SIII+stage18.SII+stage18.SI)/3,
                              stage18vs35 = (stage35.SIII+stage35.SII+stage35.SI)/3 - (stage18.SIII+stage18.SII+stage18.SI)/3,
                              stage21vs25 = (stage25.SIII+stage25.SII+stage25.SI)/3 - (stage21.SIII+stage21.SII+stage21.SI)/3,
                              stage21vs27 = (stage27.SIII+stage27.SII+stage27.SI)/3 - (stage21.SIII+stage21.SII+stage21.SI)/3,
                              stage21vs35 = (stage35.SIII+stage35.SII+stage35.SI)/3 - (stage21.SIII+stage21.SII+stage21.SI)/3,
                              stage25vs27 = (stage27.SIII+stage27.SII+stage27.SI)/3 - (stage25.SIII+stage25.SII+stage25.SI)/3,
                              stage25vs35 = (stage35.SIII+stage35.SII+stage35.SI)/3 - (stage25.SIII+stage25.SII+stage25.SI)/3,
                              stage27vs35 = (stage35.SIII+stage35.SII+stage35.SI)/3 - (stage27.SIII+stage27.SII+stage27.SI)/3,
                              levels=design)

test <- glmQLFTest(fit, contrast = my.contrasts)
stage.all <- as.data.frame(topTags(test, n=Inf))
stage.all$maxFC <- apply(stage.all[,2:16], 1, function(x) x[which.max(abs(x))])
stage.all$DE <- ifelse(stage.all$FDR < fdr_thr & abs(stage.all$maxFC) > log2(fc_thr), 1, 0)
sum(stage.all$DE)
```

```
## [1] 7714
```

The vast majority of the genes identified on the per-somite tests are also significant when the somite trios are treated as replicates.


```r
stages <- cbind(average=stage.all$DE, stage.perSomite[row.names(stage.all),])
stages$somiteSpecific <- ifelse(stages$average == 0 & rowSums(stages[,-1])>0, 1, 0)

colSums(stages[,c(1,5)])
```

```
##        average somiteSpecific 
##           7714           2977
```

```r
colSums(stages[stages$average == 0 & stages$somiteSpecific == 1,2:4])
```

```
##   somiteI  somiteII somiteIII 
##       821      1069      1808
```

Overall, the somite-specific changes seem to be replicated in the other somites, but with much smaller differences, so they don't reach significance.

#### Somite-specific changes {.tabset}

As before, the genes that are only significant in one per-somite test tend to show consistent expression patterns in the other somites, reflecting small changes difficult to pick up, rather than true somite-specific changes.

##### Somite I


```r
## select somite-specific genes
tmp <- stage.somiteI[row.names(stages[stages$somiteI==1 & stages$average==0 & rowSums(stages[2:4])==1,]),]

plots <- list()
for(i in 1:6){
  d <- retrieveGeneExpr(dataNorm, meta, gene=tmp[i,1])
  plots[[i]] <- ggplot(d, aes(x=as.factor(somite), y=count)) + 
    geom_boxplot(aes(fill=as.factor(stage))) +
    ylab("log CPM") + xlab("somite") + ggtitle(tmp[i,1]) + 
    scale_fill_brewer(palette = "Oranges") + labs(fill="stage") + 
    th
}
ggarrange(plotlist = plots, ncol = 3, nrow = 2, common.legend = TRUE, legend = "bottom")
```

![](04_differentialExpression_files/figure-html/sI.sp, -1.png)<!-- -->

##### Somite II


```r
## select somite-specific genes
tmp <- stage.somiteII[row.names(stages[stages$somiteII==1 & stages$average==0 & rowSums(stages[2:4])==1,]),]

plots <- list()
for(i in 1:6){
  d <- retrieveGeneExpr(dataNorm, meta, gene=tmp[i,1])
  plots[[i]] <- ggplot(d, aes(x=as.factor(somite), y=count)) + 
    geom_boxplot(aes(fill=as.factor(stage))) +
    ylab("log CPM") + xlab("somite") + ggtitle(tmp[i,1]) + 
    scale_fill_brewer(palette = "Oranges") + labs(fill="stage") + 
    th
}
ggarrange(plotlist = plots, ncol = 3, nrow = 2, common.legend = TRUE, legend = "bottom")
```

![](04_differentialExpression_files/figure-html/sII.sp, -1.png)<!-- -->

##### Somite III


```r
## select somite-specific genes
tmp <- stage.somiteIII[row.names(stages[stages$somiteIII==1 & stages$average==0 & rowSums(stages[2:4])==1,]),]

plots <- list()
for(i in 1:6){
  d <- retrieveGeneExpr(dataNorm, meta, gene=tmp[i,1])
  plots[[i]] <- ggplot(d, aes(x=as.factor(somite), y=count)) + 
    geom_boxplot(aes(fill=as.factor(stage))) +
    ylab("log CPM") + xlab("somite") + ggtitle(tmp[i,1]) + 
    scale_fill_brewer(palette = "Oranges") + labs(fill="stage") + 
    th
}
ggarrange(plotlist = plots, ncol = 3, nrow = 2, common.legend = TRUE, legend = "bottom")
```

![](04_differentialExpression_files/figure-html/sIII.sp, -1.png)<!-- -->


####



```r
# somite trios, per stage
write.table(somite.stage8, paste0(dir, "RNA-seq/results/04_DEresults_somiteTrios_stage8.tsv"), quote = FALSE, sep="\t")
write.table(somite.stage18, paste0(dir, "RNA-seq/results/04_DEresults_somiteTrios_stage18.tsv"), quote = FALSE, sep="\t")
write.table(somite.stage21, paste0(dir, "RNA-seq/results/04_DEresults_somiteTrios_stage21.tsv"), quote = FALSE, sep="\t")
write.table(somite.stage25, paste0(dir, "RNA-seq/results/04_DEresults_somiteTrios_stage25.tsv"), quote = FALSE, sep="\t")
write.table(somite.stage27, paste0(dir, "RNA-seq/results/04_DEresults_somiteTrios_stage27.tsv"), quote = FALSE, sep="\t")
write.table(somite.stage35, paste0(dir, "RNA-seq/results/04_DEresults_somiteTrios_stage35.tsv"), quote = FALSE, sep="\t")

# somite trios, average
write.table(somiteIIvsI.all, paste0(dir, "RNA-seq/results/04_DEresults_somiteTrios_somiteIvsII.tsv"), quote = FALSE, sep="\t")
write.table(somiteIIIvsII.all, paste0(dir, "RNA-seq/results/04_DEresults_somiteTrios_somiteIIvsIII.tsv"), quote = FALSE, sep="\t")
write.table(somiteIIIvsI.all, paste0(dir, "RNA-seq/results/04_DEresults_somiteTrios_somiteIvsIII.tsv"), quote = FALSE, sep="\t")

# summary
write.table(somiteTrios, paste0(dir, "RNA-seq/results/04_DEresults_summary_somiteTrios.tsv"), quote = FALSE, sep="\t")

# stage, per somite
write.table(stage.somiteI, paste0(dir, "RNA-seq/results/04_DEresults_stage_somiteI.tsv"), quote = FALSE, sep="\t")
write.table(stage.somiteII, paste0(dir, "RNA-seq/results/04_DEresults_stage_somiteII.tsv"), quote = FALSE, sep="\t")
write.table(stage.somiteIII, paste0(dir, "RNA-seq/results/04_DEresults_stage_somiteIII.tsv"), quote = FALSE, sep="\t")

# stage, aveerage
write.table(stage.all, paste0(dir, "RNA-seq/results/04_DEresults_stage_average.tsv"), quote = FALSE, sep="\t")

# summary
write.table(stages, paste0(dir, "RNA-seq/results/04_DEresults_summary_stage.tsv"), quote = FALSE, sep="\t")
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
## [1] UpSetR_1.4.0       RColorBrewer_1.1-2 ggpubr_0.4.0       ggrepel_0.8.2     
## [5] ggplot2_3.3.2      edgeR_3.28.1       limma_3.42.2      
## 
## loaded via a namespace (and not attached):
##  [1] statmod_1.4.34    tidyselect_1.1.0  locfit_1.5-9.4    xfun_0.16        
##  [5] purrr_0.3.4       splines_3.6.1     haven_2.3.1       lattice_0.20-41  
##  [9] carData_3.0-4     colorspace_1.4-1  vctrs_0.3.2       generics_0.0.2   
## [13] htmltools_0.5.0   yaml_2.2.1        rlang_0.4.7       pillar_1.4.6     
## [17] foreign_0.8-72    glue_1.4.1        withr_2.2.0       readxl_1.3.1     
## [21] plyr_1.8.6        lifecycle_0.2.0   stringr_1.4.0     munsell_0.5.0    
## [25] ggsignif_0.6.0    gtable_0.3.0      cellranger_1.1.0  zip_2.0.4        
## [29] evaluate_0.14     labeling_0.3      knitr_1.29        rio_0.5.16       
## [33] forcats_0.5.0     curl_4.3          broom_0.7.0       Rcpp_1.0.5       
## [37] scales_1.1.1      backports_1.1.8   abind_1.4-5       farver_2.0.3     
## [41] gridExtra_2.3     hms_0.5.3         digest_0.6.25     stringi_1.4.6    
## [45] openxlsx_4.1.5    rstatix_0.6.0     dplyr_1.0.1       cowplot_1.0.0    
## [49] grid_3.6.1        tools_3.6.1       magrittr_1.5      tibble_3.0.3     
## [53] crayon_1.3.4      tidyr_1.1.1       car_3.0-8         pkgconfig_2.0.3  
## [57] ellipsis_0.3.1    data.table_1.12.8 rmarkdown_2.3     R6_2.4.1         
## [61] compiler_3.6.1
```
