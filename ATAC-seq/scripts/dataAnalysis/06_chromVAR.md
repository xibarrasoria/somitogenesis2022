---
title: "chromVAR analysis"
date: '10 April, 2021'
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



`chromVAR` is a method that exploits the aggregated signal from peaks containing a particular motif, k-mer or annotation, to uncover those that are correlated with variability in accessibility across samples. Since the method relies on identifying matches to motif models, and these occur by chance often, the size of the regions that are analysed is important. Thus, we will start from the window counts, instead of the peak counts, since peaks are of different sizes.

The windows we used for normalisation were 150 bp in length, sliding 50bp to cover the whole genome. We then filtered to retain only those with enough counts and overlapping peaks. Counts are then normalised and batch corrected as done previously.


```r
## counts from 150 bp windows, that slide by 50 bp, along the whole genome
## filtered to retain only those overlapping peaks (as called by MACS2) and of certain abundance
windowCounts <- readRDS(paste0(dir, "ATAC-seq/results/03_windowCounts_filteredWindows.Rds"))
# add sample names
colnames(windowCounts) <- substr(windowCounts$bam.files, 83, nchar(windowCounts$bam.files)-14)

## add metadata
meta <- read.table(paste0(dir, "ATAC-seq/data/metadata_ATACseq_GQ.tsv"), header = TRUE)
# make stage and somite factors with correct levels
meta$somite <- factor(meta$somite, levels = paste0("S", c("III","II","I")))
meta$stage <- paste0("stage", meta$stage)
meta$stage <- factor(meta$stage, levels = paste0("stage", c(8,18,21,25,27,35)))

stopifnot(identical(meta$sample, colnames(windowCounts)))
colData(windowCounts) <- cbind(colData(windowCounts), meta[,c(3:6)])
## remove unwanted columns
colData(windowCounts) <- colData(windowCounts)[,c(2,5:8)]

## chromVar expects the total library size in a column named 'depth'
colnames(colData(windowCounts))[1] <- "depth"
## add GC content information
windowCounts <- addGCBias(windowCounts, genome = BSgenome.Mmusculus.UCSC.mm10)


## add normalised counts, by using the `offsets` to do trended normalisation
# offsets are in natural log scale; convert everything to log2
offsets <- assay(windowCounts, "offset")
assay(windowCounts, 'normalised') <- log2(assay(windowCounts)+0.5) - offsets/log(2)


## regress out batch effects
# PCs for correction 
pcs <- read.table(paste0(dir, "ATAC-seq/results/03_pcs_residuals.tab"))

# design matrix
meta$group <- factor(paste(meta$stage, meta$somite, sep="."))
design <- model.matrix(~0+group, meta)
colnames(design) <- paste0("stage",levels(meta$group))

## regress out PCs
normCounts <- assay(windowCounts, 'normalised')
normCounts.corr <- removeBatchEffect(normCounts, 
                                     design = design,
                                     covariates = pcs[,1:18])
## assign to windowCounts and save for later use
assay(windowCounts, 'corrected') <- normCounts.corr
saveRDS(windowCounts, paste0(dir, "ATAC-seq/results/06_windowCounts_filteredWindows_normCorr.Rds"))

## compute PCA for visualisation later on
vars <- rowVars(normCounts.corr)
tmp <- normCounts.corr[order(vars, decreasing=TRUE)[1:5000],]
pca <- prcomp(t(tmp))
pca <- as.data.frame(pca$x)
pca <- cbind(pca, meta[match(row.names(pca), meta$sample),])

ggplot(pca, aes(PC1, PC2, colour=stage)) +
  geom_point() +
  scale_color_manual(values = cols.stage) +
  th
```

![](06_chromVAR_files/figure-html/prepare_data-1.png)<!-- -->


```
##             used   (Mb) gc trigger   (Mb) limit (Mb)   max used   (Mb)
## Ncells  10535639  562.7   19005835 1015.1         NA   13354272  713.2
## Vcells 502330624 3832.5 1168127324 8912.2     102400 1138989755 8689.9
```

We start by computing deviation scores for mouse transcription factor (TF) motifs provided in the `chromVARmotifs` package that includes a non-redundant collection curated from [cisBP](http://cisbp.ccbr.utoronto.ca/). 


```r
# filter overlapping windows; the most abundant is retained
counts_filtered <- filterPeaks(counts, non_overlapping = TRUE)

# match which peaks contain each motif
data("mouse_pwms_v2") # mouse collection from cisBP; non-redundant
motif_ix <- matchMotifs(mouse_pwms_v2, counts_filtered, 
                        genome = BSgenome.Mmusculus.UCSC.mm10)

# compute deviation scores
dev <- computeDeviations(object = counts_filtered, annotations = motif_ix)
# deviation scores are very highly correlated if using only normalised data
# scores between raw and normalised data are almost identical
# proceed with batch-corrected for consistency with the rest of the analyses
```


Based on the deviation scores, we can compute the pairwise correlations between samples. Samples tend to cluster by stage, although segregation isn't perfect.


```r
## sample correlation based on deviation scores
sample_cor <- getSampleCorrelation(dev)

pheatmap(as.dist(sample_cor),
         annotation_col = as.data.frame(colData(dev)[,c(1,3,4)]),
         clustering_distance_rows = as.dist(1-sample_cor), 
         clustering_distance_cols = as.dist(1-sample_cor))
```

![](06_chromVAR_files/figure-html/sample_correlation-1.png)<!-- -->

Next, we compute the variability for each TF across samples. There are many TFs with large variability, and these might have important functions regulating somite differentiation and development.


```r
## variability for each motif across samples
variability <- computeVariability(dev)
# plot(density(variability$variability))

variability[order(variability$variability, decreasing = TRUE),][1:10,]
```

```
##                                         name variability bootstrap_lower_bound
## ENSMUSG00000005698_LINE295_Ctcf_D       Ctcf    14.13876             11.911461
## ENSMUSG00000063972_LINE1636_Nr6a1_I    Nr6a1    12.92274             10.896669
## ENSMUSG00000021255_LINE1599_Esrrb_D_N3 Esrrb    12.65683             10.761224
## ENSMUSG00000024955_LINE1608_Esrra_D_N1 Esrra    12.49177             10.606439
## ENSMUSG00000026610_LINE1613_Esrrg_D    Esrrg    12.04931             10.263966
## ENSMUSG00000026398_LINE1611_Nr5a2_D_N1 Nr5a2    11.87091             10.124178
## ENSMUSG00000026751_LINE6719_Nr5a1_I_N5 Nr5a1    11.87091             10.124178
## ENSMUSG00000015053_LINE1112_Gata2_I    Gata2    11.48068              9.330171
## ENSMUSG00000021944_LINE1116_Gata4_D_N1 Gata4    11.47710              9.228610
## ENSMUSG00000031162_LINE1118_Gata1_D    Gata1    11.18169              8.863822
##                                        bootstrap_upper_bound p_value
## ENSMUSG00000005698_LINE295_Ctcf_D                   15.84936       0
## ENSMUSG00000063972_LINE1636_Nr6a1_I                 14.49801       0
## ENSMUSG00000021255_LINE1599_Esrrb_D_N3              14.14793       0
## ENSMUSG00000024955_LINE1608_Esrra_D_N1              13.97100       0
## ENSMUSG00000026610_LINE1613_Esrrg_D                 13.51851       0
## ENSMUSG00000026398_LINE1611_Nr5a2_D_N1              13.22330       0
## ENSMUSG00000026751_LINE6719_Nr5a1_I_N5              13.22330       0
## ENSMUSG00000015053_LINE1112_Gata2_I                 13.28908       0
## ENSMUSG00000021944_LINE1116_Gata4_D_N1              13.49940       0
## ENSMUSG00000031162_LINE1118_Gata1_D                 13.17649       0
##                                        p_value_adj
## ENSMUSG00000005698_LINE295_Ctcf_D                0
## ENSMUSG00000063972_LINE1636_Nr6a1_I              0
## ENSMUSG00000021255_LINE1599_Esrrb_D_N3           0
## ENSMUSG00000024955_LINE1608_Esrra_D_N1           0
## ENSMUSG00000026610_LINE1613_Esrrg_D              0
## ENSMUSG00000026398_LINE1611_Nr5a2_D_N1           0
## ENSMUSG00000026751_LINE6719_Nr5a1_I_N5           0
## ENSMUSG00000015053_LINE1112_Gata2_I              0
## ENSMUSG00000021944_LINE1116_Gata4_D_N1           0
## ENSMUSG00000031162_LINE1118_Gata1_D              0
```
We can use PCA to visualise the samples based on their accessibility data; and then overlay the deviation scores for variable TFs. For example, Esrra, one of the top scoring TFs, is very well correlated with somite stage, suggesting that sites bound by this family of TFs increase in accessibility as development progresses.


```r
palette <- colorRamp2(breaks = c(-1,0,1), colors = c("steelblue3","white", "indianred3"))

## plot TF z-scores on PCA
tf <- "Esrra"
i <- grep(tf, row.names(dev))
df <- cbind(pca, z=assay(dev, 'z')[i,])

plots <- list()
plots[[1]] <- ggplot(df, aes(PC1, PC2, colour=stage)) +
  geom_point() +
  scale_color_manual(values = cols.stage) +
  th
plots[[2]] <- ggplot(df, aes(PC1, PC2)) +
  geom_point(aes(fill=z), colour="grey50", pch=21) +
  scale_fill_gradientn(colours = palette(seq(-1,1,0.1))) +
  ggtitle(tf) +
  th
ggarrange(plotlist = plots, ncol=2, nrow=1, legend = "bottom", align = "h")
```

![](06_chromVAR_files/figure-html/plot_stage-1.png)<!-- -->

Other TFs show variability between the somite trios, although there are many fewer instances of these.


```r
## plot TF z-scores on PCA
tf <- "Cdx2"
i <- grep(tf, row.names(dev))
df <- cbind(pca, z=assay(dev, 'z')[i,])

plots <- list()
plots[[1]] <- ggplot(df, aes(PC1, PC2, colour=stage)) +
  geom_point() +
  scale_color_manual(values = cols.stage) +
  th
plots[[2]] <- ggplot(df, aes(PC1, PC2)) +
  geom_point(aes(fill=z), colour="grey50", pch=21) +
  scale_fill_gradientn(colours = palette(seq(-1,1,0.1))) +
  ggtitle(tf) +
  th
plots[[3]] <- ggplot(df, aes(stage, z, colour=stage)) +
  geom_boxplot() +
  scale_color_manual(values = cols.stage) +
  ggtitle(tf) +
  th
plots[[4]] <- ggplot(df, aes(stage, z, colour=somite)) +
  geom_boxplot() +
  scale_color_manual(values = cols.somite) +
  ggtitle(tf) +
  th
ggarrange(plotlist = plots, ncol=2, nrow=2, legend = "none", align = "h")
```

![](06_chromVAR_files/figure-html/plot_somite-1.png)<!-- -->
To better identify TFs that vary with respect to our groups of interest, we directly test for differences in deviation scores between somites and stages.


```r
## determine if there is a significant difference between deviation scores among groups
diff_acc.somite <- differentialDeviations(dev, "somite")
diff_acc.somite <- diff_acc.somite[order(diff_acc.somite$p_value_adjusted),]
# diff_acc.somite[diff_acc.somite$p_value_adjusted<0.05,]

diff_acc.stage <- differentialDeviations(dev, "stage")
diff_acc.stage <- diff_acc.stage[order(diff_acc.stage$p_value_adjusted),]
# diff_acc.stage[diff_acc.stage$p_value_adjusted < 0.05,]
```

In total, 173 TFs are significantly different between the somite trios, and 754 between stages.



#### Cooperative binding

By looking at the difference in variability in accessibility between sites with pairs of TFs, versus when these occur alone, can suggest cooperative/competitive binding. We look at this *additional variability* present when both proteins bind together for the Hox proteins and their reported cofactors: Pbx1-4, Meis1-3 and Prep1-2 (Pknox1-2). Some patterns are evident but nothing particularly clear between Hox proteins and their cofactors, perhaps because these often occur in tripartite complexes, not just as heterodimers.


```r
select <- c(grep("Hox", colnames(motif_ix)),
            grep("Pbx", colnames(motif_ix)),
            grep("Meis", colnames(motif_ix)),
            grep("Pknox", colnames(motif_ix)))

synergy_hox <- getAnnotationSynergy(counts_filtered, motif_ix[,select])
row.names(synergy_hox) <- unlist(lapply(strsplit(row.names(synergy_hox),"_"), '[[', 3))
colnames(synergy_hox) <- unlist(lapply(strsplit(colnames(synergy_hox),"_"), '[[', 3))

pheatmap(synergy_hox, breaks = seq(-5, 5, length.out=101))
```

![](06_chromVAR_files/figure-html/synergy-1.png)<!-- -->

```r
correlation_hox <- getAnnotationCorrelation(counts_filtered, motif_ix[,select])
row.names(correlation_hox) <- unlist(lapply(strsplit(row.names(correlation_hox),"_"), '[[', 3))
colnames(correlation_hox) <- unlist(lapply(strsplit(colnames(correlation_hox),"_"), '[[', 3))
# pheatmap(correlation_hox)
```

#### Hox genes

Below are the deviation scores for all Hox genes across stages.


```r
## plot TF z-scores on PCA
plots <- list()
for(hox in grep("Hox", row.names(variability), value = TRUE)){
  i <- grep(hox, row.names(dev))
  df <- cbind(pca, z=assay(dev, 'z')[i,])
  
  plots[[hox]] <- ggplot(df, aes(stage, z, colour=stage)) +
    geom_boxplot() +
    scale_color_manual(values = cols.stage) +
    ggtitle(unlist(strsplit(hox, "_"))[3]) +
    th
}

idx <- grep("Hoxa",names(plots))
idx <- idx[c(4,3,11,2,9,10,8,7,1,6,5)]
ggarrange(plotlist = plots[idx], ncol=4, nrow=3, legend = "none", align = "hv")
```

![](06_chromVAR_files/figure-html/plot_hox-1.png)<!-- -->

```r
idx <- grep("Hoxb",names(plots))
idx <- idx[c(2,10,7,4,5,1,6,9,3,8)]
ggarrange(plotlist = plots[idx], ncol=4, nrow=3, legend = "none", align = "hv")
```

![](06_chromVAR_files/figure-html/plot_hox-2.png)<!-- -->

```r
idx <- grep("Hoxc",names(plots))
idx <- idx[c(9,6,4,3,7,5,2,8,1)]
ggarrange(plotlist = plots[idx], ncol=4, nrow=3, legend = "none", align = "hv")
```

![](06_chromVAR_files/figure-html/plot_hox-3.png)<!-- -->

```r
idx <- grep("Hoxd",names(plots))
idx <- idx[c(4,9,5,3,7,8,6,2,1)]
ggarrange(plotlist = plots[idx], ncol=4, nrow=3, legend = "none", align = "hv")
```

![](06_chromVAR_files/figure-html/plot_hox-4.png)<!-- -->


### PBM-HOMEO motif collection

A study by [Berger et al. Cell, 2008](https://doi.org/10.1016/j.cell.2008.05.024) profiled the binding affinity of mouse homeodomain TFs by using microarrays containing all possible 10nt sequences. Although binding doesn't occur under physiological conditions, and no cofactors are available, these data represent a well-controlled assay to define sequence affinity for each TF, including all Hox proteins.

We can use these motifs instead of the cisBP collection used above. In this case, results will be limited only to homeodomain TFs, but these play many important roles during development and somitogenesis. Based on the deviation scores, samples cluster by both stage and somite.


```r
homeo <- getJasparMotifs(species = "Mus musculus", collection = "PBM_HOMEO")
homeo_ix <- matchMotifs(homeo, counts_filtered, 
                        genome = BSgenome.Mmusculus.UCSC.mm10)

# compute deviation scores
dev_homeo <- computeDeviations(object = counts_filtered, annotations = homeo_ix)

# sample correaltion
sample_cor_homeo <- getSampleCorrelation(dev_homeo)
pheatmap(as.dist(sample_cor_homeo),
         annotation_col = as.data.frame(colData(dev_homeo)[,c(1,3,4)]),
         clustering_distance_rows = as.dist(1-sample_cor_homeo), 
         clustering_distance_cols = as.dist(1-sample_cor_homeo))
```

![](06_chromVAR_files/figure-html/homeo-1.png)<!-- -->

And in terms of variability in accessibility deviations, late Hox proteins show much higher variation compared to mid or early Hox members, or other homeodomain TFs. This is driven by the unique accessibility of these Hox genes in the stage35 samples. 

Cdx1 and Cdx2 are also very highly ranked, perhaps due to their consistent variation across somite trios.


```r
## variability for each motif across samples
variability_homeo <- computeVariability(dev_homeo)
variability_homeo <- variability_homeo[order(variability_homeo$variability, decreasing = TRUE),]

## plot
variability_homeo$Hox <- ifelse(grepl("Hox", variability_homeo$name), "Hox", "non-Hox")
variability_homeo$label <- ifelse(variability_homeo$variability > 3.5, variability_homeo$name, "")

ggplot(variability_homeo, aes(1:nrow(variability_homeo), variability, colour=Hox, label=label)) +
  geom_point() +
  scale_color_manual(values = c('Hox'="indianred3", 'non-Hox'="grey60")) +
  geom_text_repel(show.legend=FALSE) +
  xlab("homeodomain TFs") +
  labs(colour="") +
  th + theme(axis.text.x = element_blank(),
             axis.ticks.x = element_blank())
```

![](06_chromVAR_files/figure-html/variability_homoe-1.png)<!-- -->

The picture for these TFs when using the cisBP motif database is pretty similar, with the main differences being:

- Hoxc9, Hoxc10 and Hoxd11 are more variable when using the homeodomain-specific motifs.
- Hoxc8 and Hoxd8 are more variable with the cisBP motifs.


```r
## compare variability estimates with the two sets of motifs
# same plot with the cisBP data
tmp <- variability[variability$name %in% variability_homeo$name,]
tmp <- tmp[order(tmp$variability, decreasing=TRUE),]
tmp$label <- ifelse(tmp$variability>4, tmp$name, "")
tmp$Hox <- ifelse(grepl("Hox", tmp$name), "Hox", "non-Hox")

plots <- list()
plots[[1]] <- ggplot(tmp, aes(1:nrow(tmp), variability, colour=Hox, label=label)) +
  geom_point() +
  scale_color_manual(values = c('Hox'="indianred3", 'non-Hox'="grey60")) +
  geom_text_repel(size=4, show.legend=FALSE, max.overlaps = 50) +
  ggtitle("using cisBP motifs") +
  xlab("homeodomain TFs") +
  labs(colour="") +
  th + theme(axis.text.x = element_blank(),
             axis.ticks.x = element_blank(),
             legend.position = "none")

# direct comparison
tfs <- intersect(tmp$name, variability_homeo$name)
tmp <- data.frame(tf = tfs, 
                  cisBP = variability[match(tfs, variability$name),]$variability,
                  homeo = variability_homeo[match(tfs, variability_homeo$name),]$variability)
plots[[2]] <- ggplot(tmp, aes(cisBP, homeo, label=tf)) +
  geom_point() +
  geom_text_repel() +
  ggtitle("variability across motif databases") +
  # geom_abline(slope = 1, intercept = 0, col="grey60", lty=2) +
  th
ggarrange(plotlist = plots, ncol=2, align = "hv")
```

![](06_chromVAR_files/figure-html/compare-1.png)<!-- -->

### K-mer analysis

Repeating the analysis using k-mers instead of motif PWMs returns similar results. The most variable k-mers correspond to the motifs identified before.


```r
kmer_ix <- matchKmers(8, counts_filtered, genome = BSgenome.Mmusculus.UCSC.mm10)
dev_kmer <- computeDeviations(counts_filtered, kmer_ix)

# kmer_cov <- deviationsCovariability(kmer_dev)
# plotKmerMismatch("CTAATTCA", kmer_cov)

variability_kmer <- computeVariability(dev_kmer)
variability_kmer <- variability_kmer[order(variability_kmer$variability, decreasing = TRUE),]
variability_kmer$label <- ifelse(variability_kmer$variability > 7.5, variability_kmer$name, "")

ggplot(variability_kmer, aes(1:nrow(variability_kmer), variability, label=label)) +
  geom_point() +
  geom_text_repel(size=3, max.overlaps = 100, show.legend=FALSE) +
  xlab("k-mers (k=8)") +
  labs(colour="") +
  geom_hline(yintercept = 1.5, lty=2) +
  th + theme(axis.text.x = element_blank(),
             axis.ticks.x = element_blank())
```

![](06_chromVAR_files/figure-html/kmer_deviations-1.png)<!-- -->


```r
## plot TF z-scores on PCA
tf <- "CATCAATC"
i <- grep(tf, row.names(dev_kmer))
df <- cbind(pca, z=assay(dev_kmer, 'z')[i,])

plots <- list()
plots[[1]] <- ggplot(df, aes(PC1, PC2, colour=stage)) +
  geom_point() +
  scale_color_manual(values = cols.stage) +
  th
plots[[2]] <- ggplot(df, aes(PC1, PC2)) +
  geom_point(aes(fill=z), colour="grey50", pch=21) +
  scale_fill_gradientn(colours = palette(seq(-1,1,0.1))) +
  ggtitle(tf) +
  th
plots[[3]] <- ggplot(df, aes(stage, z, colour=stage)) +
  geom_boxplot() +
  scale_color_manual(values = cols.stage) +
  ggtitle(tf) +
  th
plots[[4]] <- ggplot(df, aes(stage, z, colour=somite)) +
  geom_boxplot() +
  scale_color_manual(values = cols.somite) +
  ggtitle(tf) +
  th
ggarrange(plotlist = plots, ncol=2, nrow=2, legend = "none", align = "h")
```

![](06_chromVAR_files/figure-html/plot_kmer-1.png)<!-- -->

However, the scores for specific k-mers could be useful later on.




### Signal by chromosomal location instead of TF motifs

Another way of finding loci that are variable across samples is to look at the accessibility signal for windows across the genome, based on chromosomal location. With this approach, sets of 25 windows are grouped as the *regions of interest* and deviation and variability scores are computed. The range slides down by 10 windows, to cover the whole genome.

A very large number of regions capture variability across samples; some of the top hits are, not surprisingly, loci containing Hox genes from late paralogue groups, which are only accessible in stage35.


```r
# match which peaks are in each chromosomal window
cis_ix <- getCisGroups(counts_filtered, grpsize = 25, stepsize = 10) 

# compute deviation scores
dev_chr <- computeDeviations(object = counts_filtered, annotations = cis_ix)

# variability per TF
variability_chr <- computeVariability(dev_chr)
variability_chr <- variability_chr[order(variability_chr$variability, decreasing = TRUE),]
variability_chr$label <- ifelse(variability_chr$variability>7, variability_chr$name, "")

ggplot(variability_chr, aes(1:nrow(variability_chr), variability, label=label)) +
  geom_point() +
  geom_text_repel(size=3, max.overlaps = 100) +
  xlab("transcription factors") +
  labs(colour="") +
  th + theme(axis.text.x = element_blank(),
             axis.ticks.x = element_blank())
```

![](06_chromVAR_files/figure-html/chr_loc-1.png)<!-- -->

And the correlation between deviation profiles separates samples by stage well, although the thoracic samples mix somewhat.


```r
# correlation between samples
sample_cor_chr <- getSampleCorrelation(dev_chr)
pheatmap(as.dist(sample_cor_chr),
         annotation_col = as.data.frame(colData(dev_chr)[,c(3,4)]),
         clustering_distance_rows = as.dist(1-sample_cor_chr), 
         clustering_distance_cols = as.dist(1-sample_cor_chr))
```

![](06_chromVAR_files/figure-html/chr_correlation-1.png)<!-- -->

As before, we can test for differential deviation scores between somite trios and between stages. This results in over 2K regions between somite trios and nearly 30K between stages.


```r
## determine if there is a significant difference between deviation scores among groups
diff_acc.chrLoc.somite <- differentialDeviations(dev_chr, "somite")
diff_acc.chrLoc.somite <- diff_acc.chrLoc.somite[order(diff_acc.chrLoc.somite$p_value_adjusted),]
# diff_acc.chrLoc.somite[diff_acc.chrLoc.somite$p_value_adjusted<0.05,]

diff_acc.chrLoc.stage <- differentialDeviations(dev_chr, "stage")
diff_acc.chrLoc.stage <- diff_acc.chrLoc.stage[order(diff_acc.chrLoc.stage$p_value_adjusted),]
# diff_acc.chrLoc.stage[diff_acc.chrLoc.stage$p_value_adjusted < 0.05,]

## annotate loci
# restrict to most significant hits
# somite trios
diff_acc.chrLoc.somite.diff <- diff_acc.chrLoc.somite[diff_acc.chrLoc.somite$p_value_adjusted < 0.05,]
diff_acc.chrLoc.somite.diff$region <- row.names(diff_acc.chrLoc.somite.diff)
# convert to genomic ranges
diff_acc.chrLoc.somite.diff.gr <- t(apply(diff_acc.chrLoc.somite.diff, 1, function(x){
  windows <- rowRanges(cis_ix)[which(assay(cis_ix[,x[3]]))]
  c(chr = as.character(unique(seqnames(windows))), start=min(start(windows)), end=max(end(windows)),
    fdr = x[2])
}))
diff_acc.chrLoc.somite.diff.gr <- GRanges(diff_acc.chrLoc.somite.diff.gr[,1],
                                          IRanges(as.numeric(diff_acc.chrLoc.somite.diff.gr[,2]), 
                                                  as.numeric(diff_acc.chrLoc.somite.diff.gr[,3])),
                                          region = row.names(diff_acc.chrLoc.somite.diff.gr),
                                          FDR = diff_acc.chrLoc.somite.diff.gr[,4])
# annotate
mcols(diff_acc.chrLoc.somite.diff.gr) <- cbind(mcols(diff_acc.chrLoc.somite.diff.gr),
                                               detailRanges(diff_acc.chrLoc.somite.diff.gr, 
                                                            txdb=TxDb.Mmusculus.UCSC.mm10.knownGene,
                                                            orgdb=org.Mm.eg.db, 
                                                            promoter=c(0, 0), dist=0))

# stages
diff_acc.chrLoc.stage.diff <- diff_acc.chrLoc.stage[diff_acc.chrLoc.stage$p_value_adjusted < 0.05,]
diff_acc.chrLoc.stage.diff$region <- row.names(diff_acc.chrLoc.stage.diff)
# convert to genomic ranges
diff_acc.chrLoc.stage.diff.gr <- t(apply(diff_acc.chrLoc.stage.diff, 1, function(x){
  windows <- rowRanges(cis_ix)[which(assay(cis_ix[,x[3]]))]
  c(chr = as.character(unique(seqnames(windows))), start=min(start(windows)), end=max(end(windows)),
    fdr = x[2])
}))
diff_acc.chrLoc.stage.diff.gr <- GRanges(diff_acc.chrLoc.stage.diff.gr[,1],
                                          IRanges(as.numeric(diff_acc.chrLoc.stage.diff.gr[,2]), 
                                                  as.numeric(diff_acc.chrLoc.stage.diff.gr[,3])),
                                          region = row.names(diff_acc.chrLoc.stage.diff.gr),
                                          FDR = diff_acc.chrLoc.stage.diff.gr[,4])
# annotate
mcols(diff_acc.chrLoc.stage.diff.gr) <- cbind(mcols(diff_acc.chrLoc.stage.diff.gr),
                                              detailRanges(diff_acc.chrLoc.stage.diff.gr, 
                                                           txdb=TxDb.Mmusculus.UCSC.mm10.knownGene,
                                                           orgdb=org.Mm.eg.db, 
                                                           promoter=c(0, 0), dist=0))
```

These deviation scores based on genomic location will be useful to infer the accessibility landscape across samples for particular genes.



```r
saveRDS(counts_filtered, paste0(dir, "ATAC-seq/results/06_windowsFilteredCHROMVAR_counts_norm_14PCs.Rds"))

## cisBP motifs
saveRDS(motif_ix, paste0(dir, "ATAC-seq/results/06_motif_matches.Rds"))
saveRDS(dev, paste0(dir, "ATAC-seq/results/06_motif_deviationScores.Rds"))
write.table(variability, paste0(dir, "ATAC-seq/results/06_motif_variability.tsv"),
            quote = FALSE, sep="\t")
write.table(diff_acc.somite, paste0(dir, "ATAC-seq/results/06_motif_diffAcc_somiteTrios.tsv"),
            quote = FALSE, sep="\t")
write.table(diff_acc.stage, paste0(dir, "ATAC-seq/results/06_motif_diffAcc_stages.tsv"),
            quote = FALSE, sep="\t")
saveRDS(synergy_hox, paste0(dir, "ATAC-seq/results/06_motif_synergyScores.Rds"))
saveRDS(correlation_hox, paste0(dir, "ATAC-seq/results/06_motif_correlationScores.Rds"))

## k-mers
saveRDS(kmer_ix, paste0(dir, "ATAC-seq/results/06_8kmer_matches.Rds"))
saveRDS(dev_kmer, paste0(dir, "ATAC-seq/results/06_8kmer_deviationScores.Rds"))
write.table(variability_kmer, paste0(dir, "ATAC-seq/results/06_8kmer_variability.tsv"),
            quote = FALSE, sep="\t")

## chromosomal location
saveRDS(cis_ix, paste0(dir, "ATAC-seq/results/06_chrLocation_matches.Rds"))
saveRDS(dev_chr, paste0(dir, "ATAC-seq/results/06_chrLocation_deviationScores.Rds"))
write.table(variability_chr, paste0(dir, "ATAC-seq/results/06_chrLocation_variability.tsv"),
            quote = FALSE, sep="\t")
saveRDS(diff_acc.chrLoc.somite.diff.gr, 
        paste0(dir, "ATAC-seq/results/06_chrLocation_diffAcc_somiteTrios.Rds"))
saveRDS(diff_acc.chrLoc.stage.diff.gr, 
        paste0(dir, "ATAC-seq/results/06_chrLocation_diffAcc_stages.Rds"))
```



```r
sessionInfo()
```

```
## R version 4.0.3 (2020-10-10)
## Platform: x86_64-apple-darwin17.0 (64-bit)
## Running under: macOS Catalina 10.15.7
## 
## Matrix products: default
## BLAS:   /Library/Frameworks/R.framework/Versions/4.0/Resources/lib/libRblas.dylib
## LAPACK: /Library/Frameworks/R.framework/Versions/4.0/Resources/lib/libRlapack.dylib
## 
## locale:
## [1] en_GB.UTF-8/en_GB.UTF-8/en_GB.UTF-8/C/en_GB.UTF-8/en_GB.UTF-8
## 
## attached base packages:
## [1] parallel  stats4    stats     graphics  grDevices utils     datasets 
## [8] methods   base     
## 
## other attached packages:
##  [1] TxDb.Mmusculus.UCSC.mm10.knownGene_3.10.0
##  [2] GenomicFeatures_1.42.3                   
##  [3] org.Mm.eg.db_3.12.0                      
##  [4] AnnotationDbi_1.52.0                     
##  [5] csaw_1.24.3                              
##  [6] circlize_0.4.12                          
##  [7] pheatmap_1.0.12                          
##  [8] RColorBrewer_1.1-2                       
##  [9] ggrepel_0.9.1                            
## [10] ggpubr_0.4.0                             
## [11] ggplot2_3.3.3                            
## [12] BiocParallel_1.24.1                      
## [13] BSgenome.Mmusculus.UCSC.mm10_1.4.0       
## [14] BSgenome_1.58.0                          
## [15] rtracklayer_1.50.0                       
## [16] Biostrings_2.58.0                        
## [17] XVector_0.30.0                           
## [18] motifmatchr_1.12.0                       
## [19] chromVARmotifs_0.2.0                     
## [20] chromVAR_1.12.0                          
## [21] limma_3.46.0                             
## [22] SummarizedExperiment_1.20.0              
## [23] Biobase_2.50.0                           
## [24] GenomicRanges_1.42.0                     
## [25] GenomeInfoDb_1.26.7                      
## [26] IRanges_2.24.1                           
## [27] S4Vectors_0.28.1                         
## [28] BiocGenerics_0.36.0                      
## [29] MatrixGenerics_1.2.1                     
## [30] matrixStats_0.58.0                       
## 
## loaded via a namespace (and not attached):
##   [1] readxl_1.3.1                backports_1.2.1            
##   [3] BiocFileCache_1.14.0        plyr_1.8.6                 
##   [5] lazyeval_0.2.2              TFBSTools_1.28.0           
##   [7] digest_0.6.27               htmltools_0.5.1.1          
##   [9] GO.db_3.12.1                fansi_0.4.2                
##  [11] JASPAR2016_1.18.0           magrittr_2.0.1             
##  [13] memoise_2.0.0               openxlsx_4.2.3             
##  [15] readr_1.4.0                 annotate_1.68.0            
##  [17] R.utils_2.10.1              askpass_1.1                
##  [19] prettyunits_1.1.1           colorspace_2.0-0           
##  [21] rappdirs_0.3.3              blob_1.2.1                 
##  [23] haven_2.3.1                 xfun_0.22                  
##  [25] dplyr_1.0.5                 crayon_1.4.1               
##  [27] RCurl_1.98-1.3              jsonlite_1.7.2             
##  [29] TFMPvalue_0.0.8             glue_1.4.2                 
##  [31] gtable_0.3.0                zlibbioc_1.36.0            
##  [33] DelayedArray_0.16.3         car_3.0-10                 
##  [35] shape_1.4.5                 abind_1.4-5                
##  [37] scales_1.1.1                edgeR_3.32.1               
##  [39] DBI_1.1.1                   rstatix_0.7.0              
##  [41] miniUI_0.1.1.1              Rcpp_1.0.6                 
##  [43] progress_1.2.2              viridisLite_0.3.0          
##  [45] xtable_1.8-4                foreign_0.8-81             
##  [47] bit_4.0.4                   DT_0.17                    
##  [49] htmlwidgets_1.5.3           httr_1.4.2                 
##  [51] nabor_0.5.0                 ellipsis_0.3.1             
##  [53] farver_2.1.0                pkgconfig_2.0.3            
##  [55] XML_3.99-0.6                R.methodsS3_1.8.1          
##  [57] dbplyr_2.1.1                sass_0.3.1                 
##  [59] locfit_1.5-9.4              utf8_1.2.1                 
##  [61] labeling_0.4.2              tidyselect_1.1.0           
##  [63] rlang_0.4.10                reshape2_1.4.4             
##  [65] later_1.1.0.1               munsell_0.5.0              
##  [67] cellranger_1.1.0            tools_4.0.3                
##  [69] cachem_1.0.4                DirichletMultinomial_1.32.0
##  [71] generics_0.1.0              RSQLite_2.2.5              
##  [73] broom_0.7.6                 evaluate_0.14              
##  [75] stringr_1.4.0               fastmap_1.1.0              
##  [77] yaml_2.2.1                  knitr_1.31                 
##  [79] bit64_4.0.5                 zip_2.1.1                  
##  [81] caTools_1.18.2              purrr_0.3.4                
##  [83] KEGGREST_1.30.1             mime_0.10                  
##  [85] R.oo_1.24.0                 poweRlaw_0.70.6            
##  [87] xml2_1.3.2                  pracma_2.3.3               
##  [89] biomaRt_2.46.3              compiler_4.0.3             
##  [91] plotly_4.9.3                curl_4.3                   
##  [93] png_0.1-7                   ggsignif_0.6.1             
##  [95] tibble_3.1.0                bslib_0.2.4                
##  [97] stringi_1.5.3               highr_0.8                  
##  [99] forcats_0.5.1               lattice_0.20-41            
## [101] CNEr_1.26.0                 Matrix_1.3-2               
## [103] vctrs_0.3.7                 pillar_1.5.1               
## [105] lifecycle_1.0.0             jquerylib_0.1.3            
## [107] GlobalOptions_0.1.2         cowplot_1.1.1              
## [109] data.table_1.14.0           bitops_1.0-6               
## [111] httpuv_1.5.5                R6_2.5.0                   
## [113] promises_1.2.0.1            rio_0.5.26                 
## [115] codetools_0.2-18            gtools_3.8.2               
## [117] assertthat_0.2.1            seqLogo_1.56.0             
## [119] openssl_1.4.3               withr_2.4.1                
## [121] GenomicAlignments_1.26.0    Rsamtools_2.6.0            
## [123] GenomeInfoDbData_1.2.4      hms_1.0.0                  
## [125] grid_4.0.3                  tidyr_1.1.3                
## [127] rmarkdown_2.7               carData_3.0-4              
## [129] shiny_1.6.0
```

