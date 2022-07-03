---
title: "Differential accessibility analysis"
date: '11 August, 2020'
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



We have ATAC-seq data for 50 good-quality samples, comprising the last three generated somites of mouse embryos at six different developmental stages.


```r
## metadata
meta <- read.table(paste0(dir, "ATAC-seq/data/metadata_ATACseq_GQ.tsv"), stringsAsFactors = FALSE, header = TRUE)
meta <- meta[meta$QCpass==1,]

## data
filtered.data <- readRDS(paste0(dir, "ATAC-seq/results/03_windowCounts_filteredWindows.Rds"))
```

We have defined a set of peaks representing regions of open chromatin, and computed the sequencing data counts for 150bp windows covering such loci. We have used a trended-normalisation approach to remove systematic biases. We remove further unwanted technical variation by using PCA on the residuals of the model containing the somite-stage information; the first few PCs are correlated to FRiP and regressing out the first 18 components cleans up the data nicely.

### Differential accessibility analysis

Now we can perform differential analysis to obtain differences between the somite trios, and between stages, using the normalised data while regressing out the PCs that capture the remaining technical variation in the data.

We use `edgeR` to perform the differential testing. We first estimate the dispersion and fit the design containing the interaction of somite and stage, plus the 18 PCs.


```r
## PCA on the fit residuals
pcs <- read.table(paste0(dir, "ATAC-seq/results/03_pcs_residuals.tab"))

## edgeR object
y <- asDGEList(filtered.data)

## design
meta$group <- factor(paste(meta$stage, meta$somite, sep="."))
design <- model.matrix(~0+group, meta)
colnames(design) <- paste0("stage",levels(meta$group))

design <- cbind(design, pcs[,1:18])
colnames(design)[-c(1:18)] <- paste0("PC",1:18)

## dispersion
y <- estimateDisp(y, design)
par(mfrow=c(1,2))
plotBCV(y)

## fit model
fit <- glmQLFit(y, design, robust=TRUE)
plotQLDisp(fit)
```

![](04_differentialAccessibility_files/figure-html/dispersionFit-1.png)<!-- -->

```r
# summary(fit$df.prior)

saveRDS(y, file=paste0(dir, "ATAC-seq/results/04_edgeRobject.Rds"))
```

#### Differences across somite trios

We can now compare the chromatin accessibility as somites differentiate (somite I vs II vs III). First, we perform the comparison for each stage. We are under powered in several stages, where we don't have replicates for a particular somite-stage combination.


```r
table(meta$stage, meta$somite)
```

```
##     
##      SI SII SIII
##   8   2   2    2
##   18  3   4    3
##   21  2   4    1
##   25  4   4    3
##   27  3   1    1
##   35  5   4    2
```

We definitely exclude stage27 from the analysis. For all other stages, we define all three pairwise comparisons, and test them all together with an anova-like test.


```r
## contrasts
my.contrasts <- makeContrasts(stage8.IIvsI    = stage8.SII   - stage8.SI, 
                              stage8.IIIvsI   = stage8.SIII  - stage8.SI, 
                              stage8.IIIvsII  = stage8.SIII  - stage8.SII, 
                              stage18.IIvsI   = stage18.SII  - stage18.SI, 
                              stage18.IIIvsI  = stage18.SIII - stage18.SI, 
                              stage18.IIIvsII = stage18.SIII - stage18.SII, 
                              stage21.IIvsI   = stage21.SII  - stage21.SI, 
                              stage21.IIIvsI  = stage21.SIII - stage21.SI, 
                              stage21.IIIvsII = stage21.SIII - stage21.SII, 
                              stage25.IIvsI   = stage25.SII  - stage25.SI, 
                              stage25.IIIvsI  = stage25.SIII - stage25.SI, 
                              stage25.IIIvsII = stage25.SIII - stage25.SII, 
                              stage35.IIvsI   = stage35.SII  - stage35.SI, 
                              stage35.IIIvsI  = stage35.SIII - stage35.SI, 
                              stage35.IIIvsII = stage35.SIII - stage35.SII, 
                              levels=design)

## test
trios.perStage <- list()
trios.perStage[["stage8"]]  <- glmQLFTest(fit, contrast = my.contrasts[,1:3])
trios.perStage[["stage18"]] <- glmQLFTest(fit, contrast = my.contrasts[,4:6])
trios.perStage[["stage21"]] <- glmQLFTest(fit, contrast = my.contrasts[,7:9])
trios.perStage[["stage25"]] <- glmQLFTest(fit, contrast = my.contrasts[,10:12])
trios.perStage[["stage35"]] <- glmQLFTest(fit, contrast = my.contrasts[,13:15])
```

The tests are performed on every window that passed the minimum abundance filter and is within the peak set. However, we are interested in regions of open chromatin as a whole -i.e. *peaks*-, more than the individual windows *per se*. Thus, we merge all adjacent windows that are no more than 150bp apart into regions. We also set an upper limit of 1.5kb to avoid chaining reactions of broadly open areas where a single region could be very long. In these cases, the locus simply gets split into adjacent regions of 1.5kb. For each defined region, we compute a combined p-value against the null hypothesis that none of the windows are differentially accessible.


```r
## merge adjacent windows into regions
merged <- mergeWindows(rowRanges(filtered.data), tol=150, max.width = 1500)

somiteTrios.perStage <- list()
for(contrast in paste0("stage",c(8,18,21,25,35))){
  somiteTrios.perStage[[contrast]] <- merged$region
  mcols(somiteTrios.perStage[[contrast]]) <- combineTests(merged$id, 
                                                          trios.perStage[[contrast]]$table)
  best <- getBestTest(merged$id, trios.perStage[[contrast]]$table)
  mcols(somiteTrios.perStage[[contrast]]) <- cbind(mcols(somiteTrios.perStage[[contrast]]), 
                                                   best[,-c((ncol(best)-2):ncol(best))])
}
saveRDS(somiteTrios.perStage, paste0(dir, 
                                    "ATAC-seq/results/04_diffAccessibility_somiteTrios_perStage.Rds"))
```

Regions are deemed differentially accessible if they have an adjusted p-value < 0.05 and an absolute fold change larger than 1.5.

There is great variability in the number of regions that reach significance for the different stages. This likely reflects the number of replicates and the quality of the samples. Not surprisingly, stage21 is under powered; stage35 should do better but the samples from SIII are not the best quality.

Samples from stage25 were of particularly good quality and, together with many replicates, return the largest number of significant regions.


```r
## retain significant regions
# first, order all test in the same way, and retain the largest fold-change
for(i in 1:length(somiteTrios.perStage)){
  somiteTrios.perStage[[i]] <- somiteTrios.perStage[[i]][order(somiteTrios.perStage[[i]])]
  somiteTrios.perStage[[i]]$maxFC <- apply(mcols(somiteTrios.perStage[[i]])[,11:13], 1, 
                                           function(x) x[which.max(abs(x))])
}
stopifnot(identical(ranges(somiteTrios.perStage[["stage8"]]), 
                    ranges(somiteTrios.perStage[["stage18"]])))
stopifnot(identical(ranges(somiteTrios.perStage[["stage21"]]), 
                    ranges(somiteTrios.perStage[["stage35"]])))

# then, combine into a results table
somiteTrios.perStage.all <- as.data.frame(somiteTrios.perStage[["stage8"]])[,c(1:4,6,14,20)]
colnames(somiteTrios.perStage.all)[(ncol(somiteTrios.perStage.all)-1):ncol(somiteTrios.perStage.all)] <- paste0("stage8", c(".FDR", ".FC"))
somiteTrios.perStage.all$stage18.FDR <- somiteTrios.perStage[["stage18"]]$FDR
somiteTrios.perStage.all$stage18.FC  <- somiteTrios.perStage[["stage18"]]$maxFC
somiteTrios.perStage.all$stage21.FDR <- somiteTrios.perStage[["stage21"]]$FDR
somiteTrios.perStage.all$stage21.FC  <- somiteTrios.perStage[["stage21"]]$maxFC
somiteTrios.perStage.all$stage25.FDR <- somiteTrios.perStage[["stage25"]]$FDR
somiteTrios.perStage.all$stage25.FC  <- somiteTrios.perStage[["stage25"]]$maxFC
somiteTrios.perStage.all$stage35.FDR <- somiteTrios.perStage[["stage35"]]$FDR
somiteTrios.perStage.all$stage35.FC  <- somiteTrios.perStage[["stage35"]]$maxFC

## keep regions with an FDR < 0.05 and an absolute fold-change larger than 1.5 for at least one comparison
keep.8 <- somiteTrios.perStage.all$stage8.FDR < 0.05 & 
  abs(somiteTrios.perStage.all$stage8.FC) > log2(1.5)
keep.18 <- somiteTrios.perStage.all$stage18.FDR < 0.05 & 
  abs(somiteTrios.perStage.all$stage18.FC) > log2(1.5)
keep.21 <- somiteTrios.perStage.all$stage21.FDR < 0.05 & 
  abs(somiteTrios.perStage.all$stage21.FC) > log2(1.5)
keep.25 <- somiteTrios.perStage.all$stage25.FDR < 0.05 & 
  abs(somiteTrios.perStage.all$stage25.FC) > log2(1.5)
keep.35 <- somiteTrios.perStage.all$stage35.FDR < 0.05 & 
  abs(somiteTrios.perStage.all$stage35.FC) > log2(1.5)

somiteTrios.perStage.all <- somiteTrios.perStage.all[
  keep.8 | keep.18 | keep.21 | keep.25 | keep.35, ]
row.names(somiteTrios.perStage.all) <- paste0(somiteTrios.perStage.all$seqnames, ":", 
                                              somiteTrios.perStage.all$start, "-", 
                                              somiteTrios.perStage.all$end)

nDA <- t(data.frame(stage8 = sum(keep.8),
                  stage18 = sum(keep.18),
                  stage21 = sum(keep.21),
                  stage25 = sum(keep.25),
                  stage35 = sum(keep.35)))
colnames(nDA) <- "nDA"
nDA
```

```
##          nDA
## stage8   563
## stage18  161
## stage21    6
## stage25 1028
## stage35   38
```

To overcome the issues of insufficient power, we can also test the average across stages of each somite. This will recover changes that are consistent across stages and miss any stage-specific differences.


```r
## contrasts
my.contrasts <- makeContrasts(
  somiteIIvsI = (stage8.SII + stage18.SII + stage21.SII + 
                   stage25.SII + stage27.SII + stage35.SII)/6 - 
    (stage8.SI + stage18.SI + stage21.SI + 
       stage25.SI + stage27.SI + stage35.SI)/6,
  somiteIIIvsI = (stage8.SIII + stage18.SIII + stage21.SIII + 
                    stage25.SIII + stage27.SIII + stage35.SIII)/6 - 
    (stage8.SI + stage18.SI + stage21.SI + 
       stage25.SI + stage27.SI + stage35.SI)/6,
  somiteIIIvsII = (stage8.SIII + stage18.SIII + stage21.SIII + 
                     stage25.SIII + stage27.SIII + stage35.SIII)/6 - 
    (stage8.SII + stage18.SII + stage21.SII + 
       stage25.SII + stage27.SII + stage35.SII)/6,
  levels=design)

trios.average <- list()
trios.average[["somiteIIvsI"]]   <- glmQLFTest(fit, contrast = my.contrasts[,"somiteIIvsI"])
trios.average[["somiteIIIvsI"]]  <- glmQLFTest(fit, contrast = my.contrasts[,"somiteIIIvsI"])
trios.average[["somiteIIIvsII"]] <- glmQLFTest(fit, contrast = my.contrasts[,"somiteIIIvsII"])

## data.frames
somiteTrios.average <- list()
for(contrast in paste0("somite",c("IIvsI","IIIvsI","IIIvsII"))){
  somiteTrios.average[[contrast]] <- merged$region
  mcols(somiteTrios.average[[contrast]]) <- combineTests(merged$id, 
                                                         trios.average[[contrast]]$table)
  best <- getBestTest(merged$id, trios.average[[contrast]]$table)
  mcols(somiteTrios.average[[contrast]]) <- cbind(mcols(somiteTrios.average[[contrast]]),
                                                  best[,-c((ncol(best)-2):ncol(best))])
}
saveRDS(somiteTrios.average, paste0(dir,
                                    "ATAC-seq/results/04_diffAccessibility_somiteTrios_average.Rds"))
```

We only detect a large number of regions in the comparison between somites I and III, which are the most different.


```r
## same order for all tests
for(i in 1:length(somiteTrios.average)){
  somiteTrios.average[[i]] <- somiteTrios.average[[i]][order(somiteTrios.average[[i]])]
}

## combine into single table
somiteTrios.average.all <- as.data.frame(somiteTrios.average[["somiteIIvsI"]])[,c(1:4,6,10,13)]
colnames(somiteTrios.average.all)[(ncol(somiteTrios.average.all)-1):ncol(somiteTrios.average.all)] <- paste0("somiteIIvsI", c(".FDR", ".FC"))
somiteTrios.average.all$somiteIIIvsII.FDR <- somiteTrios.average[["somiteIIIvsII"]]$FDR
somiteTrios.average.all$somiteIIIvsII.FC <- somiteTrios.average[["somiteIIIvsII"]]$logFC
somiteTrios.average.all$somiteIIIvsI.FDR <- somiteTrios.average[["somiteIIIvsI"]]$FDR
somiteTrios.average.all$somiteIIIvsI.FC <- somiteTrios.average[["somiteIIIvsI"]]$logFC

## keep regions with an FDR < 0.05 and an absolute fold-change larger than 1.5 for at least one comparison
keep.IIvsI <- somiteTrios.average.all$somiteIIvsI.FDR < 0.05 & 
  abs(somiteTrios.average.all$somiteIIvsI.FC) > log2(1.5)
keep.IIIvsII <- somiteTrios.average.all$somiteIIIvsII.FDR < 0.05 & 
  abs(somiteTrios.average.all$somiteIIIvsII.FC) > log2(1.5)
keep.IIIvsI <- somiteTrios.average.all$somiteIIIvsI.FDR < 0.05 & 
  abs(somiteTrios.average.all$somiteIIIvsI.FC) > log2(1.5)

somiteTrios.average.all <- somiteTrios.average.all[keep.IIvsI | keep.IIIvsII | keep.IIIvsI, ]
row.names(somiteTrios.average.all) <- paste0(somiteTrios.average.all$seqnames, ":", 
                                             somiteTrios.average.all$start, "-", 
                                             somiteTrios.average.all$end)

nDA <- t(data.frame(somiteIvsII = sum(keep.IIvsI),
                  somiteIIvsIII = sum(keep.IIIvsII),
                  somiteIvsIII = sum(keep.IIIvsI)))
colnames(nDA) <- "nDA"
nDA
```

```
##                nDA
## somiteIvsII     69
## somiteIIvsIII    1
## somiteIvsIII  1329
```

Finally, we only consider stage-specific changes, those that were not significant in the average test. 


```r
## create a summary table with all DA regions
somiteTrios.all <- rbind(somiteTrios.average.all[,1:5], somiteTrios.perStage.all[,1:5])
somiteTrios.all <- somiteTrios.all[!duplicated(somiteTrios.all),]

somiteTrios.all$average <- ifelse(row.names(somiteTrios.all) %in% 
                                    row.names(somiteTrios.average.all), 1, 0)
somiteTrios.all$stage8 <- ifelse(row.names(somiteTrios.all) %in% 
                          row.names(somiteTrios.perStage.all[
                            somiteTrios.perStage.all$stage8.FDR<0.05,]) & 
                            somiteTrios.all$average == 0, 1, 0)
somiteTrios.all$stage18 <- ifelse(row.names(somiteTrios.all) %in% 
                          row.names(somiteTrios.perStage.all[
                            somiteTrios.perStage.all$stage18.FDR<0.05,]) & 
                            somiteTrios.all$average == 0, 1, 0)
somiteTrios.all$stage21 <- ifelse(row.names(somiteTrios.all) %in% 
                          row.names(somiteTrios.perStage.all[
                            somiteTrios.perStage.all$stage21.FDR<0.05,]) & 
                            somiteTrios.all$average == 0, 1, 0)
somiteTrios.all$stage25 <- ifelse(row.names(somiteTrios.all) %in% 
                          row.names(somiteTrios.perStage.all[
                            somiteTrios.perStage.all$stage25.FDR<0.05,]) & 
                            somiteTrios.all$average == 0, 1, 0)
somiteTrios.all$stage35 <- ifelse(row.names(somiteTrios.all) %in% 
                          row.names(somiteTrios.perStage.all[
                            somiteTrios.perStage.all$stage35.FDR<0.05,]) & 
                            somiteTrios.all$average == 0, 1, 0)
colSums(somiteTrios.all[,6:11])
```

```
## average  stage8 stage18 stage21 stage25 stage35 
##    1366     481      91       2     758      31
```



#### Differences across developmental stages

Next, we test for changes across development. First, test all pairwise comparisons between stages on a per-somite basis. 


```r
## somite I
my.contrasts <- makeContrasts(stage8vs18  = stage18.SI - stage8.SI, 
                              stage8vs21  = stage21.SI - stage8.SI,
                              stage8vs25  = stage25.SI - stage8.SI, 
                              stage8vs27  = stage27.SI - stage8.SI, 
                              stage8vs35  = stage35.SI - stage8.SI, 
                              stage18vs21 = stage21.SI - stage18.SI, 
                              stage18vs25 = stage25.SI - stage18.SI, 
                              stage18vs27 = stage27.SI - stage18.SI, 
                              stage18vs35 = stage35.SI - stage18.SI, 
                              stage21vs25 = stage25.SI - stage21.SI, 
                              stage21vs27 = stage27.SI - stage21.SI, 
                              stage21vs35 = stage35.SI - stage21.SI, 
                              stage25vs27 = stage27.SI - stage25.SI, 
                              stage25vs35 = stage35.SI - stage25.SI, 
                              stage27vs35 = stage35.SI - stage27.SI, 
                              levels=design)

development.perSomite <- list()
development.perSomite[["somiteI"]] <- glmQLFTest(fit, contrast = my.contrasts)

## somite II
my.contrasts <- makeContrasts(stage8vs18  = stage18.SII - stage8.SII, 
                              stage8vs21  = stage21.SII - stage8.SII,
                              stage8vs25  = stage25.SII - stage8.SII, 
                              stage8vs27  = stage27.SII - stage8.SII, 
                              stage8vs35  = stage35.SII - stage8.SII, 
                              stage18vs21 = stage21.SII - stage18.SII, 
                              stage18vs25 = stage25.SII - stage18.SII, 
                              stage18vs27 = stage27.SII - stage18.SII, 
                              stage18vs35 = stage35.SII - stage18.SII, 
                              stage21vs25 = stage25.SII - stage21.SII, 
                              stage21vs27 = stage27.SII - stage21.SII, 
                              stage21vs35 = stage35.SII - stage21.SII, 
                              stage25vs27 = stage27.SII - stage25.SII, 
                              stage25vs35 = stage35.SII - stage25.SII, 
                              stage27vs35 = stage35.SII - stage27.SII, 
                              levels=design)
development.perSomite[["somiteII"]] <- glmQLFTest(fit, contrast = my.contrasts)

## somite III
my.contrasts <- makeContrasts(stage8vs18  = stage18.SIII - stage8.SIII, 
                              stage8vs21  = stage21.SIII - stage8.SIII,
                              stage8vs25  = stage25.SIII - stage8.SIII, 
                              stage8vs27  = stage27.SIII - stage8.SIII, 
                              stage8vs35  = stage35.SIII - stage8.SIII, 
                              stage18vs21 = stage21.SIII - stage18.SIII, 
                              stage18vs25 = stage25.SIII - stage18.SIII, 
                              stage18vs27 = stage27.SIII - stage18.SIII, 
                              stage18vs35 = stage35.SIII - stage18.SIII, 
                              stage21vs25 = stage25.SIII - stage21.SIII, 
                              stage21vs27 = stage27.SIII - stage21.SIII, 
                              stage21vs35 = stage35.SIII - stage21.SIII, 
                              stage25vs27 = stage27.SIII - stage25.SIII, 
                              stage25vs35 = stage35.SIII - stage25.SIII, 
                              stage27vs35 = stage35.SIII - stage27.SIII, 
                              levels=design)
development.perSomite[["somiteIII"]] <- glmQLFTest(fit, contrast = my.contrasts)
```


```r
stage.perSomite <- list()
for(contrast in paste0("somite",c("I","II","III"))){
  stage.perSomite[[contrast]] <- merged$region
  mcols(stage.perSomite[[contrast]]) <- combineTests(merged$id, 
                                                         development.perSomite[[contrast]]$table)
  best <- getBestTest(merged$id, development.perSomite[[contrast]]$table)
  mcols(stage.perSomite[[contrast]]) <- cbind(mcols(stage.perSomite[[contrast]]), 
                                                  best[,-c((ncol(best)-2):ncol(best))])
}
saveRDS(stage.perSomite, paste0(dir, 
                                    "ATAC-seq/results/04_diffAccessibility_stages_perSomite.Rds"))
```

This gives us a few thousand significant regions for somites II and III, but many more are recovered for somite I.


```r
for(i in 1:length(stage.perSomite)){
  stage.perSomite[[i]] <- stage.perSomite[[i]][order(stage.perSomite[[i]])]
  stage.perSomite[[i]]$maxFC <- apply(mcols(stage.perSomite[[i]])[,35:49], 1, 
                                      function(x) x[which.max(abs(x))])
}
stopifnot(identical(ranges(stage.perSomite[["somiteI"]]), ranges(stage.perSomite[["somiteII"]])))
stopifnot(identical(ranges(stage.perSomite[["somiteI"]]), ranges(stage.perSomite[["somiteIII"]])))

stage.perSomite.all <- as.data.frame(stage.perSomite[["somiteI"]])[,c(1:4,6,38,56)]
colnames(stage.perSomite.all)[(ncol(stage.perSomite.all)-1):ncol(stage.perSomite.all)] <- paste0("somiteI", c(".FDR", ".FC"))
stage.perSomite.all$somiteII.FDR <- stage.perSomite[["somiteII"]]$FDR
stage.perSomite.all$somiteII.FC <- stage.perSomite[["somiteII"]]$maxFC
stage.perSomite.all$somiteIII.FDR <- stage.perSomite[["somiteIII"]]$FDR
stage.perSomite.all$somiteIII.FC <- stage.perSomite[["somiteIII"]]$maxFC

## keep only regions that are significant for at least one of the tests
keep.somiteI <- stage.perSomite.all$somiteI.FDR < 0.05 & 
  abs(stage.perSomite.all$somiteI.FC) > log2(1.5)
keep.somiteII <- stage.perSomite.all$somiteII.FDR < 0.05 & 
  abs(stage.perSomite.all$somiteII.FC) > log2(1.5)
keep.somiteIII <- stage.perSomite.all$somiteIII.FDR < 0.05 & 
  abs(stage.perSomite.all$somiteIII.FC) > log2(1.5)

stage.perSomite.all <- stage.perSomite.all[keep.somiteI | keep.somiteII | keep.somiteIII,]
row.names(stage.perSomite.all) <- paste0(stage.perSomite.all$seqnames, ":", 
                                         stage.perSomite.all$start, "-", 
                                         stage.perSomite.all$end)

nDA <- t(data.frame(somiteI   = sum(keep.somiteI),
                    somiteII  = sum(keep.somiteII),
                    somiteIII = sum(keep.somiteIII)))
colnames(nDA) <- "nDA"
nDA
```

```
##             nDA
## somiteI   13891
## somiteII   3822
## somiteIII  4287
```

Now, we instead take the average of all three somites and test for stage differences.


```r
my.contrasts <- makeContrasts(stage8vs18 = (stage18.SIII+stage18.SII+stage18.SI)/3 -
                                (stage8.SIII+stage8.SII+stage8.SI)/3,
                              stage8vs21 = (stage21.SIII+stage21.SII+stage21.SI)/3 -
                                (stage8.SIII+stage8.SII+stage8.SI)/3,
                              stage8vs25 = (stage25.SIII+stage25.SII+stage25.SI)/3 -
                                (stage8.SIII+stage8.SII+stage8.SI)/3,
                              stage8vs27 = (stage27.SIII+stage27.SII+stage27.SI)/3 -
                                (stage8.SIII+stage8.SII+stage8.SI)/3,
                              stage8vs35 = (stage35.SIII+stage35.SII+stage35.SI)/3 -
                                (stage8.SIII+stage8.SII+stage8.SI)/3,
                              stage18vs21 = (stage21.SIII+stage21.SII+stage21.SI)/3 -
                                (stage18.SIII+stage18.SII+stage18.SI)/3,
                              stage18vs25 = (stage25.SIII+stage25.SII+stage25.SI)/3 -
                                (stage18.SIII+stage18.SII+stage18.SI)/3,
                              stage18vs27 = (stage27.SIII+stage27.SII+stage27.SI)/3 -
                                (stage18.SIII+stage18.SII+stage18.SI)/3,
                              stage18vs35 = (stage35.SIII+stage35.SII+stage35.SI)/3 -
                                (stage18.SIII+stage18.SII+stage18.SI)/3,
                              stage21vs25 = (stage25.SIII+stage25.SII+stage25.SI)/3 -
                                (stage21.SIII+stage21.SII+stage21.SI)/3,
                              stage21vs27 = (stage27.SIII+stage27.SII+stage27.SI)/3 -
                                (stage21.SIII+stage21.SII+stage21.SI)/3,
                              stage21vs35 = (stage35.SIII+stage35.SII+stage35.SI)/3 -
                                (stage21.SIII+stage21.SII+stage21.SI)/3,
                              stage25vs27 = (stage27.SIII+stage27.SII+stage27.SI)/3 -
                                (stage25.SIII+stage25.SII+stage25.SI)/3,
                              stage25vs35 = (stage35.SIII+stage35.SII+stage35.SI)/3 -
                                (stage25.SIII+stage25.SII+stage25.SI)/3,
                              stage27vs35 = (stage35.SIII+stage35.SII+stage35.SI)/3 -
                                (stage27.SIII+stage27.SII+stage27.SI)/3, 
                              levels=design)

development.average <- glmQLFTest(fit, contrast = my.contrasts)
```

And we test against the null that none of the stages are different. Therefore, anything with an FDR of 5% or lower indicates a region of open chromatin that is differentially accessible between at least a pair of stages. We also filter out regions with an absolute fold-change lower than 1.5.


```r
stage.average <- merged$region
mcols(stage.average) <- combineTests(merged$id, development.average$table)
best <- getBestTest(merged$id, development.average$table)
mcols(stage.average) <- cbind(mcols(stage.average), best[,-c((ncol(best)-2):ncol(best))])

saveRDS(stage.average, paste0(dir, "ATAC-seq/results/04_diffAccessibility_stages_average.Rds"))
```

Nearly 30 thousand regions are differentially accessible across development.


```r
# somite average
stage.average.all <- as.data.frame(stage.average)
stage.average.all$maxFC <- apply(stage.average.all[,40:54], 1, function(x) x[which.max(abs(x))])
stage.average.all <- stage.average.all[,c(1:4,6,38,56)]
stage.average.all <- stage.average.all[stage.average.all$FDR < 0.05 & 
                                         abs(stage.average.all$maxFC) > log2(1.5),]
row.names(stage.average.all) <- paste0(stage.average.all$seqnames, ":", 
                                       stage.average.all$start, "-", 
                                       stage.average.all$end)

nDA <- t(data.frame(any = nrow(stage.average.all)))
colnames(nDA) <- "nDA"
nDA
```

```
##       nDA
## any 28404
```


Most of the regions identified on a per-somite basis are also significant when taking their average.


```r
ave <- GRanges(stage.average.all$seqnames, IRanges(stage.average.all$start, stage.average.all$end))
perSomite <- GRanges(stage.perSomite.all$seqnames, IRanges(stage.perSomite.all$start, stage.perSomite.all$end))

nDA <- matrix(c(length(perSomite), length(ave), 
                length(unique(subjectHits(findOverlaps(ave, perSomite)))), 
                length(unique(queryHits(findOverlaps(ave, perSomite))))), ncol=2)
colnames(nDA) <- c("total", "overlap")
row.names(nDA) <- c("perSomite", "average")
nDA
```

```
##           total overlap
## perSomite 17073   12870
## average   28404   13576
```

As before, we only consider regions to be somite-specific if they were not identified in the average test.


```r
## create a summary table with all DA regions
stage.all <- rbind(stage.average.all[,1:4], stage.perSomite.all[,1:4])
stage.all <- stage.all[!duplicated(stage.all),]

stage.all$average <- ifelse(row.names(stage.all) %in% row.names(stage.average.all), 1, 0)
stage.all$somiteI <- ifelse(row.names(stage.all) %in% 
                              row.names(stage.perSomite.all[stage.perSomite.all$somiteI.FDR < 0.05,]) & 
                              stage.all$average == 0, 1, 0)
stage.all$somiteII <- ifelse(row.names(stage.all) %in%
                               row.names(stage.perSomite.all[stage.perSomite.all$somiteII.FDR < 0.05,]) & 
                               stage.all$average == 0, 1, 0)
stage.all$somiteIII <- ifelse(row.names(stage.all) %in%
                                row.names(stage.perSomite.all[stage.perSomite.all$somiteIII.FDR < 0.05,]) &
                                stage.all$average == 0, 1, 0)

colSums(stage.all[,-c(1:4)])
```

```
##   average   somiteI  somiteII somiteIII 
##     28404      3997       339       330
```

### Evaluation of the results

In order to evaluate the DA regions and check the somite- and stage-specific changes, we compute the normalised and corrected counts for each peak. We can then use these counts to check the DA regions.


```r
peakCounts <- readRDS(paste0(dir, "ATAC-seq/results/03_peakCounts_csawMerged.Rds"))

# raw counts
data <- log2(assay(peakCounts)+1)
colnames(data) <- substr(peakCounts$bam.files, 83, nchar(peakCounts$bam.files)-14)
tmp <- as.data.frame(rowRanges(peakCounts))
row.names(data) <- paste0(tmp$seqnames, ":", tmp$start, "-", tmp$end)

# loess normalisation
sf.trended <- assay(peakCounts, "offset")
dataNorm <- data - sf.trended/log(2)

## technical effects
design <- model.matrix(~0+group, meta)
colnames(design) <- paste0("stage",levels(meta$group))
dataNorm.corrected <- removeBatchEffect(dataNorm, design = design, covariates = pcs[,1:18])
saveRDS(dataNorm.corrected, file=paste0(dir, "ATAC-seq/results/04_peakCounts_csawMerged.NORM.batchCorrected_18PCs.Rds"))
```


```r
retrieveGeneExpr <- function(dataNorm.corrected, meta, columns=c("stage", "somite", "date"), region){
  d <- as.data.frame(dataNorm.corrected[region,])
  stopifnot(identical(row.names(d), meta$sample))
  d <- cbind(d, meta[,columns])
  d[,2] <- as.factor(d[,2])
  if(ncol(d) <= length(columns)){
    d$cpm <- 0
    d <- d[,c(4,1:3)]
  }else{ colnames(d)[1] <- "cpm" }
  return(d)
}

th <- theme_bw() + theme(axis.text.x = element_text(size=10), axis.title.x = element_text(size=12), axis.text.y = element_text(size=10), axis.title.y = element_text(size=12), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), panel.border = element_blank(), plot.title = element_text(face="bold", hjust = 0.5, size=12))

col.stage <- brewer.pal(n=6, "YlOrRd")

boxplotExpr <- function(group_by=3, colour_by=3, region=NA){
  d <- retrieveGeneExpr(dataNorm.corrected, meta, region = region, columns=c("stage", "somite"))
  if(colour_by == 3){
    ggplot(d, aes(x=d[,group_by], y=cpm)) + 
      geom_boxplot(aes(fill=d[,colour_by])) + 
      scale_fill_manual(values=alpha(rep("orchid4",3), c(0.25, 0.5, 0.75))) +
      xlab(colnames(d)[group_by]) + ylab("log2 CPM") + 
      ggtitle(region) + 
      labs(fill=colnames(d)[colour_by]) + 
      th
  }else{
    ggplot(d, aes(x=d[,group_by], y=cpm)) + 
      geom_boxplot(aes(fill=d[,colour_by])) + 
      scale_fill_manual(values=col.stage) + 
      xlab(colnames(d)[group_by]) + ylab("log2 CPM") + 
      ggtitle(region) + 
      labs(fill=colnames(d)[colour_by]) + 
      th
  }
}
```


#### Somite differentiation {.tabset}

First, we check a few regions that were significant in the average test, selected at random. This shows that:

- Most regions tend to be monotonically up or downregulated.
- The differences in accessibility between adjacent somites are small.

This is consistent with the majority of differences detected between somites I and III only, where the change is big enough to be detected, but not between adjacent somite pairs.


```r
set.seed(297)
sel <- sample(1:sum(somiteTrios.all$average), size = 16)
plots <- list()
for(i in 1:length(sel)){
  plots[[i]] <- boxplotExpr(region = row.names(somiteTrios.all)[sel[i]])
}
ggarrange(plotlist = plots, ncol = 4, nrow = 4, common.legend = TRUE, legend = "bottom")
```

![](04_differentialAccessibility_files/figure-html/somiteTrios_plot-1.png)<!-- -->

Next, we do the same but for the stage-specific differences.

##### Stage 8

The differences between somites at stage 8 are indeed much more pronounced than in other stages; for some, the other stages show a similar behaviour, but with much smaller differences and thus do not reach statistical significance.


```r
sel <- sample(which(somiteTrios.all$stage8==1), size = 16)
plots <- list()
for(i in 1:length(sel)){
  plots[[i]] <- boxplotExpr(region = row.names(somiteTrios.all)[sel[i]], group_by = 2)
}
ggarrange(plotlist = plots, ncol = 4, nrow = 4, common.legend = TRUE, legend = "bottom")
```

![](04_differentialAccessibility_files/figure-html/somiteTrios_plot.s8-1.png)<!-- -->

##### Stage 18

The situation is similar for the stage 18 changes.


```r
sel <- sample(which(somiteTrios.all$stage18==1), size = 16)
plots <- list()
for(i in 1:length(sel)){
  plots[[i]] <- boxplotExpr(region = row.names(somiteTrios.all)[sel[i]], group_by = 2)
}
ggarrange(plotlist = plots, ncol = 4, nrow = 4, common.legend = TRUE, legend = "bottom")
```

![](04_differentialAccessibility_files/figure-html/somiteTrios_plot.s18-1.png)<!-- -->

##### Stage 21

There are only two changes detected with the stage 21 samples that were not picked up by the average test. Both of these are driven by an outlier in somite III, which has only one replicate


```r
sel <- which(somiteTrios.all$stage21==1)
plots <- list()
for(i in 1:length(sel)){
  plots[[i]] <- boxplotExpr(region = row.names(somiteTrios.all)[sel[i]], group_by = 2)
}
ggarrange(plotlist = plots, ncol = 2, nrow = 1, common.legend = TRUE, legend = "bottom")
```

![](04_differentialAccessibility_files/figure-html/somiteTrios_plot.s21-1.png)<!-- -->

##### Stage 25

Again, the differences are much more pronounced in the stage 25 samples, but the other stages tend to show a similar behaviour.


```r
sel <- sample(which(somiteTrios.all$stage25==1), size = 16)
plots <- list()
for(i in 1:length(sel)){
  plots[[i]] <- boxplotExpr(region = row.names(somiteTrios.all)[sel[i]], group_by = 2)
}
ggarrange(plotlist = plots, ncol = 4, nrow = 4, common.legend = TRUE, legend = "bottom")
```

![](04_differentialAccessibility_files/figure-html/somiteTrios_plot.s25-1.png)<!-- -->

##### Stage 35

The regions specifically detected from the samples from stage 35 embryos show quite distinct patterns to the rest of the dataset. This might suggest that the somite differentiation process at later stages of development is different.


```r
sel <- sample(which(somiteTrios.all$stage35==1), size = 16)
plots <- list()
for(i in 1:length(sel)){
  plots[[i]] <- boxplotExpr(region = row.names(somiteTrios.all)[sel[i]], group_by = 2)
}
ggarrange(plotlist = plots, ncol = 4, nrow = 4, common.legend = TRUE, legend = "bottom")
```

![](04_differentialAccessibility_files/figure-html/somiteTrios_plot.s35-1.png)<!-- -->

####

Overall, it seems that the changes detected on a per-somite basis are regions with pronounced differences in one particular stage (which makes it statistically significant), and with a consistent pattern of expression across stages. The exception might be stage 35, which shows lesser consistent behaviour with other stages.

This is evidenced when we compute the correlation of expression values for the somite trios, across stages. Most regions significant in only the stage-specific tests have a high correlation value with at least one other stage, as reflected by the maximum correlation between all pairwise comparisons. However, the median correlation values vary across the whole [-1, 1] interval, indicating that, commonly, at least one stage shows a different pattern; this is consistent with their lack of detection in the average test. 

Also, the median value distributions for regions specific to stages 8 and 35 are biased towards lower and negative correlation values compared to the ones for regions from stages 18 and 25. This might reflect the more dense sampling of the middle stages, which pushes the correlations towards positive values. The regions from stage 35 have consistently lower correlation values, supporting our observation that these regions are not as consistent across development as the rest.


```r
## for each stage-specific region (except for stage21 that only has 2), compute the correlation of expression across stages
stats <- list()
for(s in paste0("stage", c(8,18,25,35))){
  test <- somiteTrios.all[somiteTrios.all[,s]==1,]
  stat <- matrix(ncol=2, nrow=nrow(test))
  for(idx in 1:nrow(test)){
    dat <- dataNorm.corrected[row.names(test)[idx],]
    mat <- matrix(ncol=3, nrow=6)
    i=1
    for(stage in c(8,18,21,25,27,35)){
      j=1
      for(somite in c("SI", "SII", "SIII")){
        mat[i,j] <- mean(dat[which(meta$stage == stage & meta$somite == somite)])
        j <- j+1
      }
      i <- i+1
    }
    stat[idx, 1] <- max(unique(c(cor(t(mat))))[-1])
    stat[idx, 2] <- median(unique(c(cor(t(mat))))[-1])
  }
  stats[[s]] <- stat
}

stats.flat <- as.data.frame(do.call(rbind, stats))
colnames(stats.flat) <- c("max", "median")
stats.flat$stage <- rep(paste0("stage", c(8,18,25,35)), unlist(lapply(stats, nrow)))
stats.flat$stage <- factor(stats.flat$stage, levels=paste0("stage", c(8,18,25,35)))

plots <- list()
plots[[1]] <- ggplot(stats.flat, aes(stage, max, colour=stage)) + 
  geom_boxplot() +
  th + theme(legend.position = "none")
plots[[2]] <- ggplot(stats.flat, aes(stage, median, colour=stage)) + 
  geom_boxplot() + 
  th + theme(legend.position = "none")
ggarrange(plotlist = plots, ncol=2)
```

![](04_differentialAccessibility_files/figure-html/correlations_somite-1.png)<!-- -->

#### Developmental progression {.tabset}

The regions that DA between stages in the average test show a wide range of dynamics.


```r
sel <- sample(1:sum(stage.all$average), size = 16)
plots <- list()
for(i in 1:length(sel)){
  plots[[i]] <- boxplotExpr(region = row.names(stage.all)[sel[i]], group_by = 2, colour_by = 2)
}
ggarrange(plotlist = plots, ncol = 4, nrow = 4, legend = "none")
```

![](04_differentialAccessibility_files/figure-html/stages_plot-1.png)<!-- -->

We saw that the majority of the regions detected in the somite-specific comparisons were also identified in the average test. But there were still some potentially somite-specific changes, especially for somite I.

##### Somite I

In general, we indeed observe much larger differences between somite I samples. Sometimes the pattern is recapitulated by the other somites but with much smaller changes; but some regions seem to be truly somite-specific.


```r
sel <- sample(which(stage.all$somiteI==1), size = 16)
plots <- list()
for(i in 1:length(sel)){
  plots[[i]] <- boxplotExpr(region = row.names(stage.all)[sel[i]], group_by = 3, colour_by = 2)
}
ggarrange(plotlist = plots, ncol = 4, nrow = 4, common.legend = TRUE, legend = "bottom")
```

![](04_differentialAccessibility_files/figure-html/stages_plot.sI-1.png)<!-- -->

##### Somite II

A similar situation is observed for the somite II samples, although the other somites seem to have a similar pattern much more often.


```r
sel <- sample(which(stage.all$somiteII==1), size = 16)
plots <- list()
for(i in 1:length(sel)){
  plots[[i]] <- boxplotExpr(region = row.names(stage.all)[sel[i]], group_by = 3, colour_by = 2)
}
ggarrange(plotlist = plots, ncol = 4, nrow = 4, common.legend = TRUE, legend = "bottom")
```

![](04_differentialAccessibility_files/figure-html/stages_plot.sII-1.png)<!-- -->

##### Somite III

And the same is true for somite III samples.


```r
sel <- sample(which(stage.all$somiteIII==1), size = 16)
plots <- list()
for(i in 1:length(sel)){
  plots[[i]] <- boxplotExpr(region = row.names(stage.all)[sel[i]], group_by = 3, colour_by = 2)
}
ggarrange(plotlist = plots, ncol = 4, nrow = 4, common.legend = TRUE, legend = "bottom")
```

![](04_differentialAccessibility_files/figure-html/stages_plot.sIII-1.png)<!-- -->

####


```r
write.table(somiteTrios.all, paste0(dir, "ATAC-seq/results/04_DAregions_somiteTrios.tsv"), 
            quote = FALSE, sep="\t")
write.table(stage.all, paste0(dir, "ATAC-seq/results/04_DAregions_stages.tsv"), 
            quote = FALSE, sep="\t")
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
##  [1] RColorBrewer_1.1-2          ggpubr_0.4.0               
##  [3] ggplot2_3.3.2               UpSetR_1.4.0               
##  [5] edgeR_3.28.1                limma_3.42.2               
##  [7] csaw_1.20.0                 SummarizedExperiment_1.16.1
##  [9] DelayedArray_0.12.3         BiocParallel_1.20.1        
## [11] matrixStats_0.56.0          Biobase_2.46.0             
## [13] GenomicRanges_1.38.0        GenomeInfoDb_1.22.1        
## [15] IRanges_2.20.2              S4Vectors_0.24.4           
## [17] BiocGenerics_0.32.0        
## 
## loaded via a namespace (and not attached):
##  [1] bitops_1.0-6             bit64_4.0.2              progress_1.2.2          
##  [4] httr_1.4.2               tools_3.6.1              backports_1.1.8         
##  [7] R6_2.4.1                 DBI_1.1.0                colorspace_1.4-1        
## [10] withr_2.2.0              tidyselect_1.1.0         gridExtra_2.3           
## [13] prettyunits_1.1.1        bit_4.0.3                curl_4.3                
## [16] compiler_3.6.1           labeling_0.3             rtracklayer_1.46.0      
## [19] scales_1.1.1             askpass_1.1              rappdirs_0.3.1          
## [22] stringr_1.4.0            digest_0.6.25            Rsamtools_2.2.3         
## [25] foreign_0.8-72           rmarkdown_2.3            rio_0.5.16              
## [28] XVector_0.26.0           pkgconfig_2.0.3          htmltools_0.5.0         
## [31] dbplyr_1.4.4             readxl_1.3.1             rlang_0.4.7             
## [34] rstudioapi_0.11          RSQLite_2.2.0            farver_2.0.3            
## [37] generics_0.0.2           zip_2.0.4                dplyr_1.0.1             
## [40] car_3.0-8                RCurl_1.98-1.2           magrittr_1.5            
## [43] GenomeInfoDbData_1.2.2   Matrix_1.2-18            Rcpp_1.0.5              
## [46] munsell_0.5.0            abind_1.4-5              lifecycle_0.2.0         
## [49] stringi_1.4.6            yaml_2.2.1               carData_3.0-4           
## [52] zlibbioc_1.32.0          plyr_1.8.6               BiocFileCache_1.10.2    
## [55] grid_3.6.1               blob_1.2.1               forcats_0.5.0           
## [58] crayon_1.3.4             lattice_0.20-41          cowplot_1.0.0           
## [61] splines_3.6.1            Biostrings_2.54.0        haven_2.3.1             
## [64] GenomicFeatures_1.38.2   hms_0.5.3                locfit_1.5-9.4          
## [67] knitr_1.29               pillar_1.4.6             ggsignif_0.6.0          
## [70] biomaRt_2.42.1           XML_3.99-0.3             glue_1.4.1              
## [73] evaluate_0.14            data.table_1.12.8        vctrs_0.3.2             
## [76] cellranger_1.1.0         gtable_0.3.0             openssl_1.4.2           
## [79] purrr_0.3.4              tidyr_1.1.1              assertthat_0.2.1        
## [82] openxlsx_4.1.5           xfun_0.16                broom_0.7.0             
## [85] rstatix_0.6.0            tibble_3.0.3             GenomicAlignments_1.22.1
## [88] AnnotationDbi_1.48.0     memoise_1.1.0            statmod_1.4.34          
## [91] ellipsis_0.3.1
```


