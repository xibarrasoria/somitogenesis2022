---
title: "Somite maturation"
date: '10 August, 2022'
output:
  html_document:
    keep_md: true
    fig_width: 5
    fig_height: 5
    fig_caption: yes
    code_folding: hide
    df_print: paged
    toc: true
    toc_float: 
      collapsed: false
---



We have processed, QCed, normalised and batch-corrected the RNA- and ATAC-seq data from mouse somites. We have also identified differentially expressed genes and differentially accessible chromatin regions across somite trios and embryonic development.


```r
### RNA-seq
## metadata
meta.rna <- read.table(paste0(dir, "RNA-seq/data/metadata_RNAseq.tsv"), 
                       stringsAsFactors = FALSE, header = TRUE)
meta.rna <- meta.rna[meta.rna$QC == 1 & meta.rna$wrongStage == 0,]

## expression data
geneCounts <- read.table(paste0(dir, "RNA-seq/data/geneCounts.NORM_batchCorrected_14PCs.tsv"), 
                       check.names = FALSE, stringsAsFactors = FALSE)
geneCounts <- geneCounts[,which(colnames(geneCounts) %in% c("gene", meta.rna$sample)),]


### ATAC-seq
## metadata
meta.atac <- read.table(paste0(dir, "ATAC-seq/data/metadata_ATACseq_GQ.tsv"), 
                   stringsAsFactors = FALSE, header = TRUE)
meta.atac <- meta.atac[meta.atac$QCpass==1,]

## accessibility data
peakCounts <- readRDS(paste0(dir,
                         "ATAC-seq/results/04_peakCounts_csawMerged.NORM.batchCorrected_18PCs.Rds"))
```

We now turn our attention to the genomic loci that change as somites mature and differentiate into the different somite derivatives. By comparing somites I, II and III, we can characterise the changes that occur after segmentation, up to ~6 hours of maturation.

### Transcriptional changes of somite maturation

We have identified just under three thousand genes that are significantly differentially expressed between the somite trios. 


```r
degs <- read.table(paste0(dir, "RNA-seq/results/04_DEresults_summary_somiteTrios.tsv"),
                   stringsAsFactors = FALSE)
degs$gene <- geneCounts[match(row.names(degs), row.names(geneCounts)),1]
degs <- degs[,c(ncol(degs), 1:(ncol(degs)-1))]

c(average = sum(degs$ave), stage_wise = sum(degs$stageSpecific))
```

```
##    average stage_wise 
##       1235       1742
```

```r
## keep only significant genes but save the universe
universe <- data.frame(gene_id=row.names(degs), gene_name=degs$gene)
degs <- degs[rowSums(degs[,-1])>0,]
```

We have observed that most genes show small changes in expression that increase as somites differentiate, such that many changes only become detectable when comparing somites I and III. To get a better indication of the expression dynamics of this set of DE genes, we cluster their expression profiles, using the average across stages as a proxy.

We obtain nine different clusters, that represent all the expected patterns of expression. Up- and down-regulated genes are the most prevalent, with a much smaller number showing distinct expression in somite II. However, there are differences on the kinetics of the changes.


```r
## remove log transformation
data <- 2^geneCounts[,-1]

## compute mean expression for each somite across stages
means.rna <- matrix(ncol=3, nrow=nrow(degs))
colnames(means.rna) <- c("SI", "SII", "SIII")
row.names(means.rna) <- row.names(degs)
for(somite in c("SI", "SII", "SIII")){
  samples <- meta.rna[meta.rna$somite == somite,]$sample
  means.rna[,somite] <- rowMeans(data[row.names(degs), samples])
}

## standardise
means.rna <- t(apply(means.rna, 1, function(x) x/max(x)))

## use correlation as measure of dissimilarity
corrs <- cor(t(means.rna), method = "pearson")
dissim <- sqrt(0.5*(1-corrs))
dist <- as.dist(dissim)
hc <- hclust(dist, method = "average")

clusters <- cutreeDynamic(hc, distM = as.matrix(dist), 
                     minClusterSize=100, method="hybrid", 
                     deepSplit = 1, verbose = 0)
names(clusters) <- row.names(means.rna)

plots <- list()
for(c in 1:max(clusters)){
  genes <- names(clusters[clusters==c])
  dat <- means.rna[genes,]
  means <- colMeans(dat)
  std <- apply(dat, 2, sd)
  
  df <- data.frame(x = 1:3, y = means)
  df.sd <- data.frame(x = 1:3, ymin = means-std, ymax = means+std)
  
  plots[[c]] <- ggplot() + 
    geom_line(data = df, aes(x, y), size=1) +
    geom_ribbon(data = df.sd, aes(x=x, ymin=ymin, ymax=ymax), alpha=0.2, fill="black") +
    xlab(expression("somite")) + ylab(expression("expression")) + 
    ggtitle(paste0("Cluster ", c, " (n = ", nrow(dat), ")")) + 
    ylim(0,1.1) + 
    th + 
    scale_x_continuous(breaks=1:3, labels=c("SI", "SII", "SIII"))
}
ggarrange(plotlist = plots[c(1,2,5, 4,6,3, 7,9,8)], ncol=3, nrow=3)
```

![](01_somiteMaturation_files/figure-html/dynamics_degs-1.png)<!-- -->

To better define the patterns of change, we classify each DE gene according to rules that capture each of the patterns.


```r
## group genes based on defined behaviours
means.rna <- as.data.frame(means.rna)

dynamics.rna <- list()
dynamics.rna[['down.gradual']] <- means.rna[means.rna$SI==1 & 
                                          means.rna$SII<0.85 & 
                                          means.rna$SIII<0.75 & 
                                          means.rna$SII - means.rna$SIII > 0.15,]
dynamics.rna[['down.fast']] <- means.rna[means.rna$SI==1 & 
                                       means.rna$SII<0.85 & 
                                       abs(means.rna$SII - means.rna$SIII) < 0.15,]
dynamics.rna[['down.slow']] <- means.rna[abs(means.rna$SI - means.rna$SII) < 0.15 & 
                                       means.rna$SIII<0.85,]

dynamics.rna[['up.gradual']] <- means.rna[means.rna$SIII==1 & 
                                        means.rna$SII<0.85 & 
                                        means.rna$SI<0.75 & 
                                        means.rna$SII - means.rna$SI > 0.15,]
dynamics.rna[['up.slow']] <- means.rna[means.rna$SIII==1 & 
                                     means.rna$SII<0.85 & 
                                     abs(means.rna$SI - means.rna$SII) < 0.15,]
dynamics.rna[['up.fast']] <- means.rna[abs(means.rna$SII - means.rna$SIII) < 0.15 & 
                                     means.rna$SI<0.85,]

dynamics.rna[['peak']] <- means.rna[means.rna$SII==1 & 
                                  means.rna$SII - means.rna$SI > 0.15 & 
                                  means.rna$SII - means.rna$SIII > 0.15,]
dynamics.rna[['dip']] <- means.rna[means.rna$SII < means.rna$SI & 
                                 means.rna$SII < means.rna$SIII &
                                 means.rna$SII < 0.85 & 
                                 means.rna$SI - means.rna$SII > 0.15 & 
                                 means.rna$SIII - means.rna$SII > 0.15,]

## anything not classified so far has very small changes, between 0.85 and 1
covered <- do.call('c', lapply(dynamics.rna, row.names))
dynamics.rna[['flat']] <- means.rna[setdiff(row.names(means.rna), covered),]
covered <- do.call('c', lapply(dynamics.rna, row.names))

profiles.trios <- names(covered)
names(profiles.trios) <- covered
profiles.trios[grep("dip", profiles.trios)] <- "dip"
profiles.trios[grep("peak", profiles.trios)] <- "peak"
profiles.trios[grep("flat", profiles.trios)] <- "flat"
profiles.trios[grep("down.g", profiles.trios)] <- "down.gradual"
profiles.trios[grep("down.f", profiles.trios)] <- "down.fast"
profiles.trios[grep("down.s", profiles.trios)] <- "down.slow"
profiles.trios[grep("up.g", profiles.trios)] <- "up.gradual"
profiles.trios[grep("up.f", profiles.trios)] <- "up.fast"
profiles.trios[grep("up.s", profiles.trios)] <- "up.slow"

# par(mfrow=c(3,3))
# for(dyn in names(dynamics.rna)){
#   tmp <- dynamics.rna[[dyn]]
#   plot(as.numeric(tmp[1,]), type="l", ylim=c(0,1), main=paste0(dyn, " (n=", nrow(tmp), ")"), xlab="", ylab="relative expression", axes=FALSE)
#   box(bty="l"); axis(2, at=c(0,0.5,1), las=2); axis(1, at=1:3, labels = c("SI", "SII", "SIII"))
#   for(i in 2:nrow(tmp)){
#     lines(as.numeric(tmp[i,]))  
#   }
# }
```


```r
## Use the defined profiles to classify the 'flat' genes, by correlation
# mean cluster behaviour
profile.means <- list()
for(c in setdiff(unique(profiles.trios), 'flat')){
  g <- names(profiles.trios[profiles.trios==c])
  test <- means.rna[g,]
  profile.means[[c]] <- colMeans(test)
}
profile.means <- do.call('rbind', profile.means)

# classify each flat gene
genes <- names(profiles.trios[profiles.trios=='flat'])
dat <- means.rna[genes,]
profiles.trios.flat <- c()
for(i in 1:nrow(dat)){
  profiles.trios.flat <- c(profiles.trios.flat, 
                           names(which.max(
                             apply(profile.means, 1, function(x) cor(x, as.numeric(dat[i,])) )
                             )) )
}
names(profiles.trios.flat) <- genes

# check that the classification was successful
# par(mfrow=c(3,3))
# for(c in names(dynamics.rna)[-9]){
#   genes <- names(profiles.trios.flat[profiles.trios.flat==c])
#   test <- means.rna[genes,]
#   plot(as.numeric(test[1,]), type="l", ylim=c(0.8,1), main=paste0(c, " (n=", nrow(test), ")"), xlab="", ylab="relative expression", axes=FALSE)
#   box(bty="l"); axis(2, at=c(0.85,1), las=2); axis(1, at=1:3, labels = c("SI", "SII", "SIII"))
#   for(i in 2:nrow(test)){
#     lines(as.numeric(test[i,]))
#   }
# }
```

This results in cleaner clusters.


```r
## add to rest of the data
profiles.trios[names(profiles.trios.flat)] <- profiles.trios.flat

## add to summary df
degs$profile <- profiles.trios[match(row.names(degs), names(profiles.trios))]
degs$profile.gral <- ifelse(grepl("up", degs$profile), "up",
                            ifelse(degs$profile == "dip", "dip",
                                   ifelse(degs$profile == "peak", "peak", "down")))

## plot
classes <- c("gradual", "fast", "slow")
classes <- c(paste("up", classes, sep="."), 
             paste("down", classes, sep="."),
             "peak", "dip")
plots <- list()
for(c in classes){
  genes <- names(profiles.trios[profiles.trios==c])
  dat <- means.rna[genes,]
  means <- colMeans(dat)
  std <- apply(dat, 2, sd)
  
  df <- data.frame(x = 1:3, y = means)
  df.sd <- data.frame(x = 1:3, ymin = means-std, ymax = means+std)
  
  plots[[c]] <- ggplot() + 
    geom_line(data = df, aes(x, y), size=1) +
    geom_ribbon(data = df.sd, aes(x=x, ymin=ymin, ymax=ymax), alpha=0.2, fill="black") +
    xlab(expression("somite")) + ylab(expression("expression")) + 
    ggtitle(paste0("Cluster ", c, " (n = ", nrow(dat), ")")) + 
    ylim(0,1.1) + 
    th + 
    scale_x_continuous(breaks=1:3, labels=c("SI", "SII", "SIII"))
}
ggarrange(plotlist = plots, ncol=3, nrow=3)
```

![](01_somiteMaturation_files/figure-html/clustered_degs-1.png)<!-- -->

Finally, we classify gene expression according to these patterns, but considering each stage separately, to assess whether the dynamics are consistent across development. For all those genes showing up or downregulation in at least two stages, the vast majority show the same direction of change across development.


```r
## get data fro DE genes
dat <- data[rownames(degs),]

## compute average across replicates, for each stage; annotate profile
means.rna.stage <- list()
dynamics.rna.stage <- list()
for(stage in c(8,18,21,25,27,35)){
  samples <- meta.rna[meta.rna$stage==stage,]$sample
  means.rna.stage[[paste0("stage",stage)]] <- sapply(paste0("S", c("I-","II-","III")), function(s) 
    rowMeans(dat[,grep(s, colnames(dat[,samples]), value = TRUE)]))
  colnames(means.rna.stage[[paste0("stage",stage)]]) <- paste0("S", c("I","II","III"))
  ## standardise
  means.rna.stage[[paste0("stage",stage)]] <- t(apply(means.rna.stage[[paste0("stage",stage)]], 1, function(x) x/max(x)))

  dyn <- list()
  tmp <- as.data.frame(means.rna.stage[[paste0("stage",stage)]])
  dyn[['down_gradual']] <- tmp[tmp$SI==1 & tmp$SII<0.85 & tmp$SIII<0.75 & tmp$SII - tmp$SIII > 0.15,]
  dyn[['down_fast']] <- tmp[tmp$SI==1 & tmp$SII<0.85 & abs(tmp$SII - tmp$SIII) < 0.15,]
  dyn[['down_slow']] <- tmp[abs(tmp$SI - tmp$SII) < 0.15 & tmp$SIII<0.85,]
  dyn[['up_gradual']] <- tmp[tmp$SIII==1 & tmp$SII<0.85 & tmp$SI<0.75 & tmp$SII - tmp$SI > 0.15,]
  dyn[['up_slow']] <- tmp[tmp$SIII==1 & tmp$SII<0.85 & abs(tmp$SI - tmp$SII) < 0.15,]
  dyn[['up_fast']] <- tmp[abs(tmp$SII - tmp$SIII) < 0.15 & tmp$SI<0.85,]
  dyn[['peak']] <- tmp[tmp$SII==1 & tmp$SII - tmp$SI > 0.15 & tmp$SII - tmp$SIII > 0.15,]
  dyn[['dip']] <- tmp[tmp$SII < tmp$SI & tmp$SII < tmp$SIII & tmp$SII < 0.85 & tmp$SI - tmp$SII > 0.15 & tmp$SIII - tmp$SII > 0.15,]
  
  ## anything not classified so far has very small changes, between 0.85 and 1
  covered <- do.call('c', lapply(dyn, row.names))
  dyn[['flat']] <- tmp[setdiff(row.names(tmp), covered),]
  covered <- do.call('c', lapply( dyn, row.names))
  
  ## join
  dyn <- do.call(rbind, dyn)
  dyn$gene <- unlist(lapply(strsplit(row.names(dyn), ".", fixed=TRUE),'[[', 2))
  dyn$profile <- unlist(lapply(strsplit(row.names(dyn), ".", fixed=TRUE),'[[', 1))
  
  dynamics.rna.stage[[paste0("stage",stage)]] <- dyn
}

## integrate across stages
compare <- do.call(rbind, dynamics.rna.stage)
compare$stage <- unlist(lapply(strsplit(row.names(compare), ".", fixed=TRUE),'[[', 1))

## simplify profile
compare$trend <- compare$profile
compare[grep("up", compare$profile),'trend'] <- "up"
compare[grep("down", compare$profile),'trend'] <- "down"

## count stages with up/down profile
n_profiles_per_gene <- compare %>%
  group_by(gene) %>%
  summarise(n_up = sum(trend=="up"),
            n_down = sum(trend=="down"))

## classify consistent vs discordant genes
# since classification is not perfect, when the majority of stages show consistent behaviour, keep as consistent
n_profiles_per_gene$direction <- paste(n_profiles_per_gene$n_up, n_profiles_per_gene$n_down, sep=".")
n_profiles_per_gene$consistent <- ifelse(n_profiles_per_gene$direction %in% 
                                           c(paste(0,2:6,sep="."), paste(1,3:5,sep="."), paste(2,c(0,4),sep="."),
                                             paste(3,0:1,sep="."), paste(4,0:2,sep="."), paste(5,0:1,sep="."), "6.0"),
                                         TRUE, FALSE)
# ignore genes where less than two stages show up/down behaviour
n_profiles_per_gene[n_profiles_per_gene$direction %in% c("0.0","0.1","1.0"),'consistent'] <- NA
n_profiles_per_gene$gene <- geneCounts[n_profiles_per_gene$gene,1]
# table(n_profiles_per_gene$direction, n_profiles_per_gene$consistent)
n_profiles_per_gene$direction_str <- ifelse(n_profiles_per_gene$consistent & 
                                              n_profiles_per_gene$n_up > n_profiles_per_gene$n_down, "up",
                                            ifelse(n_profiles_per_gene$consistent, "down", "NA"))

round(prop.table(table(consistent_across_stages=n_profiles_per_gene$consistent))*100, 2)
```

```
## consistent_across_stages
## FALSE  TRUE 
##  24.2  75.8
```

And from these genes with consistent behaviour, there are equal number of genes up and downregulated.


```r
table(n_profiles_per_gene[n_profiles_per_gene$consistent,]$direction_str)
```

```
## 
## down   up 
## 1120 1029
```



#### Overrepresented biological processes

To gain insight into the function of the DE genes along somite maturation we use Gene Ontology enrichment analysis. First, we look for significantly enriched terms when considering **all** DE genes.

The list of DE genes is highly enriched for many terms; many relate to somitogenesis, patterning and regionalisation; others refer to cell motility and epithelial to mesenchymal transition (EMT) processes; some more are involved in the development of the musculo-skeletal system.


```r
## GO enrichment
GO.somite <- topGOtable(DEgenes = degs$gene,
                        BGgenes = universe$gene_name,
                        topGO_method2 = "elim",
                        ontology = "BP",
                        geneID = "symbol",
                        addGeneToTerms = TRUE,
                        mapping = "org.Mm.eg.db",
                        topTablerows = 500)
## significant terms
GO.somite[c(1,2,4:7,16,17,19,20,23,28,31:34,37,40,43:46,48,
            61:62,71,73,77,83,98,99,102,103,107,119,121,213), c(2:5,7)]
```

<div data-pagedtable="false">
  <script data-pagedtable-source type="application/json">
{"columns":[{"label":[""],"name":["_rn_"],"type":[""],"align":["left"]},{"label":["Term"],"name":[1],"type":["chr"],"align":["left"]},{"label":["Annotated"],"name":[2],"type":["int"],"align":["right"]},{"label":["Significant"],"name":[3],"type":["int"],"align":["right"]},{"label":["Expected"],"name":[4],"type":["dbl"],"align":["right"]},{"label":["p.value_elim"],"name":[5],"type":["dbl"],"align":["right"]}],"data":[{"1":"somitogenesis","2":"74","3":"35","4":"11.51","5":"1.80e-08","_rn_":"1"},{"1":"anterior/posterior pattern specification","2":"222","3":"85","4":"34.52","5":"1.90e-08","_rn_":"2"},{"1":"negative regulation of epithelial cell differentiation","2":"41","3":"20","4":"6.38","5":"6.10e-07","_rn_":"4"},{"1":"extracellular matrix organization","2":"256","3":"79","4":"39.81","5":"6.70e-07","_rn_":"5"},{"1":"embryonic skeletal system development","2":"136","3":"51","4":"21.15","5":"1.40e-06","_rn_":"6"},{"1":"cell adhesion","2":"1073","3":"284","4":"166.85","5":"1.50e-06","_rn_":"7"},{"1":"negative regulation of canonical Wnt signaling pathway","2":"121","3":"37","4":"18.81","5":"2.30e-05","_rn_":"16"},{"1":"skeletal system development","2":"476","3":"148","4":"74.02","5":"2.40e-05","_rn_":"17"},{"1":"cell fate specification","2":"86","3":"34","4":"13.37","5":"2.60e-05","_rn_":"19"},{"1":"cartilage development","2":"183","3":"50","4":"28.46","5":"2.90e-05","_rn_":"20"},{"1":"positive regulation of cell-substrate adhesion","2":"119","3":"42","4":"18.50","5":"5.00e-05","_rn_":"23"},{"1":"embryonic pattern specification","2":"63","3":"26","4":"9.80","5":"9.20e-05","_rn_":"28"},{"1":"osteoblast differentiation","2":"181","3":"50","4":"28.14","5":"1.20e-04","_rn_":"31"},{"1":"regulation of Rho protein signal transduction","2":"121","3":"35","4":"18.81","5":"1.30e-04","_rn_":"32"},{"1":"regulation of non-canonical Wnt signaling pathway","2":"25","3":"12","4":"3.89","5":"1.40e-04","_rn_":"33"},{"1":"regulation of cellular response to growth factor stimulus","2":"254","3":"62","4":"39.50","5":"1.40e-04","_rn_":"34"},{"1":"morphogenesis of an epithelium","2":"522","3":"145","4":"81.17","5":"1.80e-04","_rn_":"37"},{"1":"smooth muscle tissue development","2":"26","3":"12","4":"4.04","5":"2.20e-04","_rn_":"40"},{"1":"cell-matrix adhesion","2":"194","3":"58","4":"30.17","5":"2.40e-04","_rn_":"43"},{"1":"negative regulation of ERK1 and ERK2 cascade","2":"66","3":"22","4":"10.26","5":"2.60e-04","_rn_":"44"},{"1":"cardiac epithelial to mesenchymal transition","2":"30","3":"13","4":"4.66","5":"2.60e-04","_rn_":"45"},{"1":"muscle tissue morphogenesis","2":"75","3":"24","4":"11.66","5":"2.80e-04","_rn_":"46"},{"1":"positive regulation of cell differentiation","2":"793","3":"181","4":"123.31","5":"2.90e-04","_rn_":"48"},{"1":"cell migration","2":"1273","3":"292","4":"197.95","5":"4.60e-04","_rn_":"61"},{"1":"negative regulation of cell migration","2":"239","3":"57","4":"37.16","5":"4.80e-04","_rn_":"62"},{"1":"mesenchyme development","2":"271","3":"85","4":"42.14","5":"7.50e-04","_rn_":"71"},{"1":"negative regulation of cell differentiation","2":"585","3":"144","4":"90.96","5":"7.80e-04","_rn_":"73"},{"1":"skeletal muscle cell differentiation","2":"67","3":"21","4":"10.42","5":"8.90e-04","_rn_":"77"},{"1":"proximal/distal pattern formation","2":"34","3":"13","4":"5.29","5":"1.11e-03","_rn_":"83"},{"1":"mesenchymal cell development","2":"88","3":"32","4":"13.68","5":"1.44e-03","_rn_":"98"},{"1":"positive regulation of cell migration","2":"486","3":"117","4":"75.57","5":"1.47e-03","_rn_":"99"},{"1":"Notch signaling pathway","2":"165","3":"45","4":"25.66","5":"1.53e-03","_rn_":"102"},{"1":"somite rostral/caudal axis specification","2":"13","3":"7","4":"2.02","5":"1.56e-03","_rn_":"103"},{"1":"striated muscle tissue development","2":"390","3":"103","4":"60.64","5":"1.62e-03","_rn_":"107"},{"1":"fibroblast growth factor receptor signaling pathway","2":"66","3":"20","4":"10.26","5":"1.85e-03","_rn_":"119"},{"1":"positive regulation of cell population proliferation","2":"788","3":"166","4":"122.53","5":"2.04e-03","_rn_":"121"},{"1":"regionalization","2":"351","3":"119","4":"54.58","5":"5.74e-03","_rn_":"213"}],"options":{"columns":{"min":{},"max":[10]},"rows":{"min":[10],"max":[10]},"pages":{}}}
  </script>
</div>

These are all terms we would expect from somite maturation, and capture the progression from segmentation to EMT. Additionally, when we look at the proportion of genes associated with each enriched term that are down or upregulated, we can see that: 

- Terms with a majority of downregulated genes relate to segmentation and patterning, which we know is prevalent in the PSM.
- Terms with a majority of upregulated genes are instead concerned with EMT, which is being established in the more mature somites (III), to initiate migration of the derivative lineages.


```r
## create object to use with `GeneTonic`
## expression data as dds object - from edgeR
y <- readRDS(paste0(dir, "RNA-seq/results/04_edgeRobject_pca.Rds"))
dds <- as.DESeqDataSet(y)
dds <- estimateSizeFactors(dds)
colnames(rowData(dds))[2] <- "SYMBOL"
# add normalised and batch-corrected data
stopifnot(identical(row.names(dds), row.names(geneCounts)))
assay(dds, 'cpm') <- geneCounts[,-1]

## DE results
# integrate all tests, keeping the most significant
# fold-changes are consistent, such that negative values indicate expression is going down from one somite to the next.
somite_perStage <- list()
for(stage in paste0("stage",c(8,18,21,25,27,35))){
  somite_perStage[[stage]] <- read.table(
    paste0(dir, "RNA-seq/results/04_DEresults_somiteTrios_", stage, ".tsv"))
}
somite_average <- list()
for(contrast in c("somiteIvsII", "somiteIIvsIII", "somiteIvsIII")){
  somite_average[[contrast]] <- read.table(
    paste0(dir, "RNA-seq/results/04_DEresults_somiteTrios_", contrast, ".tsv"))
}
# for each DE gene, collect stats for the most significant result
# genes detected in the average test first
degs.somite_ave <- t(sapply(row.names(degs[degs$ave==1,]), function(x){
  tmp <- rbind(somite_average[[1]][x,],
               somite_average[[2]][x,],
               somite_average[[3]][x,])
  tmp$average <- 1
  return(tmp[which.min(tmp$FDR),-7])
}, simplify = FALSE))
degs.somite_ave <- do.call(rbind, degs.somite_ave)
# genes only detected in one/some stage(s)
degs.somite_stage <- t(sapply(row.names(degs[degs$stageSpecific==1,]), function(x){
  tmp <- rbind(somite_perStage[[1]][x,c(1,9,5:8)],
               somite_perStage[[2]][x,c(1,9,5:8)],
               somite_perStage[[3]][x,c(1,9,5:8)],
               somite_perStage[[4]][x,c(1,9,5:8)],
               somite_perStage[[5]][x,c(1,9,5:8)],
               somite_perStage[[6]][x,c(1,9,5:8)])
  tmp$average <- 0
  return(tmp[which.min(tmp$FDR),])
}, simplify = FALSE))
degs.somite_stage <- do.call(rbind, degs.somite_stage)
colnames(degs.somite_stage)[2] <- "logFC"
# join - keep info of which are stage specific
degs.somite_all <- rbind(degs.somite_ave, degs.somite_stage)
# rename columns to conform with geneTonic
colnames(degs.somite_all)[colnames(degs.somite_all) == "PValue"] <- "pvalue"
colnames(degs.somite_all)[colnames(degs.somite_all) == "logFC"] <- "log2FoldChange"
colnames(degs.somite_all)[colnames(degs.somite_all) == "logCPM"] <- "baseMean"
colnames(degs.somite_all)[colnames(degs.somite_all) == "FDR"] <- "padj"
# get from the logCPM to something on the scale of the baseMean
degs.somite_all$baseMean <- (2^degs.somite_all$baseMean) * mean(y$samples$lib.size) / 1e6
# use the constructor for DESeqResults
degs_res <- DESeqResults(DataFrame(degs.somite_all))
colnames(degs_res)[1] <- "SYMBOL"
# recover ensembl IDs for rownames
row.names(degs_res) <- row.names(geneCounts[match(degs_res$SYMBOL, geneCounts$gene),])

## reformat GO enrichment results
GO.somite_res <- shake_topGOtableResult(GO.somite)
# add aggregated scores (z-score for each gene set based on the logFC from DE genes)
GO.somite_res <- get_aggrscores(GO.somite_res, degs_res, universe)

##
geneTon <- GeneTonic_list(dds = dds,
                          res_de = degs_res,
                          res_enrich = GO.somite_res,
                          annotation_obj = universe)

terms <- GO.somite[c(1,2,4:7,16,17,19,20,23,28,31:34,37,40,43:46,48,
                     61:62,71,73,77,83,98,99,102,103,107,119,121,213),1]

## define outliers based on z scores
thr <- c(median(GO.somite_res$z_score)+mad(GO.somite_res$z_score),
         median(GO.somite_res$z_score)-mad(GO.somite_res$z_score))

ggplot(GO.somite_res, aes(z_score, -log10(gs_pvalue),
                          size = DE_count, label=gs_description)) +
  geom_point(colour="grey50", shape=21,
             aes(fill = aggr_score)) + 
  scale_fill_gradientn(colours = rev(brewer.pal(n=11, "RdYlBu")), 
                       limits=c(-1,1), labels = c("< -1", seq(-0.5,0.5,0.5), "> 1")) +
  xlab("z score") +
  ylab("-log10 p-value") +
  geom_vline(xintercept = 0, lty=2, colour="grey70") +
  # geom_text_repel(data = GO.somite_res[(GO.somite_res$gs_id %in% terms & 
  #                                        abs(GO.somite_res$z_score) > 1) |
  #                                        -log10(GO.somite_res$gs_pvalue) > 5.5,],
  geom_text_repel(data = GO.somite_res[GO.somite_res$gs_pvalue < 0.01 &
                                         GO.somite_res$DE_count > 10 &
                                         (GO.somite_res$z_score > thr[1] |
                                            GO.somite_res$z_score < thr[2]),],
                  colour = "grey50",
                  min.segment.length = unit(0, 'lines'),
                  nudge_y = 0.3,
                  size=3, max.overlaps = 200) +
  scale_size(range=c(0,10), breaks=c(30,50,100,500,1000,1500)) +
  th 
```

![](01_somiteMaturation_files/figure-html/gene_tonic-1.png)<!-- -->

```r
# GO.somite_res_simplified <- gs_simplify(GO.somite_res, gs_overlap = 0.7)
# gs_volcano(GO.somite_res_simplified,
#            color_by = "aggr_score",
#            gs_ids = terms,
#            scale_circles = 0.1)

# gs_summary_overview(GO.somite_res,
#                     n_gs = 30,
#                     p_value_column = "gs_pvalue",
#                     color_by = "z_score")
```


```r
## compute number fo DE genes up/downregulated
tmp <- GO.somite_res[GO.somite_res$gs_pvalue < 0.01 &
                       GO.somite_res$DE_count > 10 &
                       (GO.somite_res$z_score > thr[1] |
                          GO.somite_res$z_score < thr[2]),]
tmp <- tmp[order(tmp$z_score),]
tmp$DE_down <- sapply(tmp$gs_genes, function(x)
  sum(degs_res[degs_res$SYMBOL %in% unlist(strsplit(x, ",")),]$log2FoldChange < 0))
tmp$DE_up <- sapply(tmp$gs_genes, function(x)
  sum(degs_res[degs_res$SYMBOL %in% unlist(strsplit(x, ",")),]$log2FoldChange > 0))
stopifnot(all(tmp$DE_down+tmp$DE_up == tmp$DE_count))

## plot
tmp <- reshape2::melt(tmp[,c('gs_description', 'DE_down', 'DE_up')])
tmp$gs_description <- factor(tmp$gs_description, levels = tmp$gs_description[1:57])

ggplot(tmp, aes(x=gs_description, y=value, fill=variable, alpha=0.25)) + 
  geom_bar(position="fill", stat="identity") +
  scale_fill_manual(values = c(DE_down = "steelblue", DE_up = "indianred")) +
  xlab("") + 
  ylab("proportion of DE genes") +
  geom_hline(yintercept = 0.5) +
  th + theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5))
```

![](01_somiteMaturation_files/figure-html/Goterms_props-1.png)<!-- -->

Thus, these all validate we are capturing meaningful gene expression changes.


### Chromatin remodelling during somite maturation

We identify a similar number of chromatin regions that are differentially accessible as somites mature.


```r
dars <- read.table(paste0(dir, "ATAC-seq/results/04_DAregions_somiteTrios.tsv"), 
                   stringsAsFactors = FALSE)

dars$stageSpecific <- ifelse(rowSums(dars[,7:11])>0, 1, 0)
c(average = sum(dars$ave), stage_wise = sum(dars$stageSpecific))
```

```
##    average stage_wise 
##       1366       1335
```

We expect these to show the same dynamics observed for gene expression, with loci either opening or closing, albeit with different kinetics. Indeed, by taking the mean accessibility across stages, we observe the same clusters as with gene expression, but many fewer loci show *peak* or *dip* behaviour, which is expected since this would require very quick chromatin changes. 


```r
## remove log transformation
data <- 2^peakCounts

## compute mean expression for each somite across stages
means.atac <- matrix(ncol=3, nrow=nrow(dars))
colnames(means.atac) <- c("SI", "SII", "SIII")
row.names(means.atac) <- row.names(dars)
for(somite in c("SI", "SII", "SIII")){
  samples <- meta.atac[meta.atac$somite == somite,]$sample
  means.atac[,somite] <- rowMeans(data[row.names(dars), samples])
}

## standardise
means.atac <- t(apply(means.atac, 1, function(x) x/max(x)))

## use correlation as measure of dissimilarity
corrs <- cor(t(means.atac), method = "pearson")
dissim <- sqrt(0.5*(1-corrs))
dist <- as.dist(dissim)
hc <- hclust(dist, method = "average")

clusters <- cutreeDynamic(hc, distM = as.matrix(dist), 
                     minClusterSize=30, method="hybrid", 
                     deepSplit = 1, verbose = 0)
names(clusters) <- row.names(means.atac)

plots <- list()
for(c in 1:max(clusters)){
  genes <- names(clusters[clusters==c])
  dat <- means.atac[genes,]
  means <- colMeans(dat)
  std <- apply(dat, 2, sd)
  
  df <- data.frame(x = 1:3, y = means)
  df.sd <- data.frame(x = 1:3, ymin = means-std, ymax = means+std)
  
  plots[[c]] <- ggplot() + 
    geom_line(data = df, aes(x, y), size=1) +
    geom_ribbon(data = df.sd, aes(x=x, ymin=ymin, ymax=ymax), alpha=0.2, fill="black") +
    xlab(expression("somite")) + ylab(expression("accessibility")) + 
    ggtitle(paste0("Cluster ", c, " (n = ", nrow(dat), ")")) + 
    ylim(0,1.1) + 
    th + 
    scale_x_continuous(breaks=1:3, labels=c("SI", "SII", "SIII"))
}
ggarrange(plotlist = plots[c(2,3,4, 1,5,6, 7, 8)], ncol=3, nrow=3)
```

![](01_somiteMaturation_files/figure-html/dynamics_dars-1.png)<!-- -->

As with the RNA expression, we refine the clusters by using predefined rules that characterise each pattern.


```r
## group genes based on defined behaviours
means.atac <- as.data.frame(means.atac)

dynamics.atac <- list()
dynamics.atac[['down.gradual']] <- means.atac[means.atac$SI==1 & 
                                          means.atac$SII<0.85 & 
                                          means.atac$SIII<0.75 & 
                                          means.atac$SII - means.atac$SIII > 0.15,]
dynamics.atac[['down.fast']] <- means.atac[means.atac$SI==1 & 
                                       means.atac$SII<0.85 & 
                                       abs(means.atac$SII - means.atac$SIII) < 0.15,]
dynamics.atac[['down.slow']] <- means.atac[abs(means.atac$SI - means.atac$SII) < 0.15 & 
                                       means.atac$SIII<0.85,]

dynamics.atac[['up.gradual']] <- means.atac[means.atac$SIII==1 & 
                                        means.atac$SII<0.85 & 
                                        means.atac$SI<0.75 & 
                                        means.atac$SII - means.atac$SI > 0.15,]
dynamics.atac[['up.slow']] <- means.atac[means.atac$SIII==1 & 
                                     means.atac$SII<0.85 & 
                                     abs(means.atac$SI - means.atac$SII) < 0.15,]
dynamics.atac[['up.fast']] <- means.atac[abs(means.atac$SII - means.atac$SIII) < 0.15 & 
                                     means.atac$SI<0.85,]

dynamics.atac[['peak']] <- means.atac[means.atac$SII==1 & 
                                  means.atac$SII - means.atac$SI > 0.15 & 
                                  means.atac$SII - means.atac$SIII > 0.15,]
dynamics.atac[['dip']] <- means.atac[means.atac$SII < means.atac$SI & 
                                 means.atac$SII < means.atac$SIII &
                                 means.atac$SII < 0.85 & 
                                 means.atac$SI - means.atac$SII > 0.15 & 
                                 means.atac$SIII - means.atac$SII > 0.15,]

## anything not classified so far has very small changes, between 0.85 and 1
covered <- do.call('c', lapply(dynamics.atac, row.names))
dynamics.atac[['flat']] <- means.atac[setdiff(row.names(means.atac), covered),]
covered <- do.call('c', lapply(dynamics.atac, row.names))

profiles.trios.atac <- names(covered)
names(profiles.trios.atac) <- covered
profiles.trios.atac[grep("dip", profiles.trios.atac)] <- "dip"
profiles.trios.atac[grep("peak", profiles.trios.atac)] <- "peak"
profiles.trios.atac[grep("flat", profiles.trios.atac)] <- "flat"
profiles.trios.atac[grep("down.g", profiles.trios.atac)] <- "down.gradual"
profiles.trios.atac[grep("down.f", profiles.trios.atac)] <- "down.fast"
profiles.trios.atac[grep("down.s", profiles.trios.atac)] <- "down.slow"
profiles.trios.atac[grep("up.g", profiles.trios.atac)] <- "up.gradual"
profiles.trios.atac[grep("up.f", profiles.trios.atac)] <- "up.fast"
profiles.trios.atac[grep("up.s", profiles.trios.atac)] <- "up.slow"

# par(mfrow=c(3,3))
# for(dyn in names(dynamics.atac)){
#   tmp <- dynamics.atac[[dyn]]
#   plot(as.numeric(tmp[1,]), type="l", ylim=c(0,1),
#        main=paste0(dyn, " (n=", nrow(tmp), ")"),
#        xlab="", ylab="relative expression", axes=FALSE)
#   box(bty="l")
#   axis(2, at=c(0,0.5,1), las=2)
#   axis(1, at=1:3, labels = c("SI", "SII", "SIII"))
#   for(i in 2:nrow(tmp)){
#     lines(as.numeric(tmp[i,]))
#   }
# }
```


```r
## Use the defined profiles to classify the 'flat' genes, by correlation
# mean cluster behaviour
profile.means <- list()
for(c in setdiff(unique(profiles.trios.atac), 'flat')){
  g <- names(profiles.trios.atac[profiles.trios.atac==c])
  test <- means.atac[g,]
  profile.means[[c]] <- colMeans(test)
}
profile.means <- do.call('rbind', profile.means)

# classify each flat gene
regions <- names(profiles.trios.atac[profiles.trios.atac=='flat'])
dat <- means.atac[regions,]
profiles.trios.atac.flat <- c()
for(i in 1:nrow(dat)){
  profiles.trios.atac.flat <- c(profiles.trios.atac.flat, 
                                names(which.max(
                                  apply(profile.means, 1, function(x) cor(x, as.numeric(dat[i,])) )
                             )) )
}
names(profiles.trios.atac.flat) <- regions

# check that the classification was successful
# par(mfrow=c(3,3))
# for(c in names(dynamics.atac)[-9]){
#   regions <- names(profiles.trios.atac.flat[profiles.trios.atac.flat==c])
#   test <- means.atac[regions,]
#   plot(as.numeric(test[1,]), type="l", ylim=c(0.8,1),
#        main=paste0(c, " (n=", nrow(test), ")"),
#        xlab="", ylab="relative expression", axes=FALSE)
#   box(bty="l")
#   axis(2, at=c(0.85,1), las=2)
#   axis(1, at=1:3, labels = c("SI", "SII", "SIII"))
#   for(i in 2:nrow(test)){
#     lines(as.numeric(test[i,]))
#   }
# }
```

This shifts many regions from the *gradual* pattern to *slow* or *fast* kinetics instead; but, overall, the clusters defined by both methods agree well.


```r
## add to rest of the data
profiles.trios.atac[names(profiles.trios.atac.flat)] <- profiles.trios.atac.flat

## add to summary df
dars$profile <- profiles.trios.atac[match(row.names(dars), names(profiles.trios.atac))]
dars$profile.gral <- ifelse(grepl("up", dars$profile), "up",
                            ifelse(dars$profile == "dip", "dip",
                                   ifelse(dars$profile == "peak", "peak", "down")))

## plot
classes <- c("gradual", "fast", "slow")
classes <- c(paste("up", classes, sep="."), 
             paste("down", classes, sep="."),
             "peak", "dip")
plots <- list()
for(c in classes){
  regions <- names(profiles.trios.atac[profiles.trios.atac==c])
  dat <- means.atac[regions,]
  means <- colMeans(dat)
  std <- apply(dat, 2, sd)
  
  df <- data.frame(x = 1:3, y = means)
  df.sd <- data.frame(x = 1:3, ymin = means-std, ymax = means+std)
  
  plots[[c]] <- ggplot() + 
    geom_line(data = df, aes(x, y), size=1) +
    geom_ribbon(data = df.sd, aes(x=x, ymin=ymin, ymax=ymax), alpha=0.2, fill="black") +
    xlab(expression("somite")) + ylab(expression("accessibility")) + 
    ggtitle(paste0("Cluster ", c, " (n = ", nrow(dat), ")")) + 
    ylim(0,1.1) + 
    th + 
    scale_x_continuous(breaks=1:3, labels=c("SI", "SII", "SIII"))
}
ggarrange(plotlist = plots, ncol=3, nrow=3)
```

![](01_somiteMaturation_files/figure-html/clustered_dars-1.png)<!-- -->

To better understand how these DA regions are regulating somite maturation we first examine their genomic location, using the definitions from previous analyses (see `ATAC-seq/scripts/dataAnalysis/05_annotation_openChromatinRegions.Rmd` for details).

The majority of DA regions overlap or are close to genes, but few overlap the promoter region.


```r
## annotation of regions based on genomic context
peakAnn <- read.table(paste0(dir, "ATAC-seq/results/05_peaks_classAnnotation.tsv"), sep="\t")

## subset to DA regions
dars.ann <- peakAnn[substr(row.names(dars), 4, 30),]
rbind(n = colSums(dars.ann[,7:11]),
      pct = round(colSums(dars.ann[,7:11])/nrow(dars.ann)*100, 2))
```

```
##     promoter   genic proximal distal intergenic
## n     238.00 1358.00   720.00 334.00      51.00
## pct     8.81   50.28    26.66  12.37       1.89
```

Compared to the distribution of **all** open chromatin regions, we observe a depletion of promoters, with a corresponding increase of distal and intergenic peaks.


```r
t(data.frame(allPeaks.pct=round(colSums(peakAnn[,7:11])/nrow(peakAnn)*100, 2)))
```

```
##              promoter genic proximal distal intergenic
## allPeaks.pct    18.72 49.27    24.13   6.89       0.99
```

```r
df <- data.frame(all = round(colSums(peakAnn[,7:11])/nrow(peakAnn)*100, 2),
                 DA = round(colSums(dars.ann[,7:11])/nrow(dars.ann)*100, 2))
df <- data.frame(class=rep(row.names(df),2),
                 DA = c(rep(FALSE, 5), rep(TRUE, 5)),
                 prop = c(df$all, df$DA))

ggplot(df, aes(DA, prop, fill=class)) + 
  geom_bar(stat="identity") +
  scale_fill_manual(values = brewer.pal(n=5, "Spectral")) +
  ylab("proportion of all peaks") +
  th
```

![](01_somiteMaturation_files/figure-html/genomic_context_allPeaks-1.png)<!-- -->

This suggests that enhancer elements play an important role in mediating somite maturation. 

#### Overrepresented biological processes

We expect the DA regions that are changing along the somite trios to be regulating the biological processes captured by the DE genes. To perform a similar enrichment analysis but using the DA peaks instead, we use the Genomic Regions Enrichment of Annotations Tool (GREAT).

All open chromatin regions identified in the somites are our background, and we test for enrichments in the DA regions specifically. A number of GO terms are significantly enriched, and these agree well with those observed when using DE genes directly. The terms that are significant based on the DA regions tend to have much lower p-values in the corresponding enrichment test using DE genes. 


```r
## background is all open chromatin regions
background <- GRanges(seqnames = paste0("chr",peakAnn$seqnames),
                      ranges = IRanges(peakAnn$start, peakAnn$end))
## DA regions in somite trios
da <- GRanges(seqnames = dars$seqnames, 
              ranges = IRanges(dars$start, dars$end),
              profile = dars$profile.gral,
              stageSpecific = dars$stageSpecific)

## run GREAT
job = submitGreatJob(gr = da,
                     bg = background,
                     species = "mm10")

# retrieve which genes are associated with each peak
great.peak_geneAnn <- plotRegionGeneAssociationGraphs(job, plot=FALSE)
# length(intersect(degs$gene, unique(great.peak_geneAnn$gene)))
# 678 (22.78%)  # DE genes that are the 'closest' gene to a DA peak

## enrichment results
go.dars <- getEnrichmentTables(job, ontology = "GO Biological Process")
pheno.dars <- getEnrichmentTables(job, category = "Phenotype")

## GO biological process
# go.dars[[1]]
# epithelium dev/morph; skeletal system dev/morph.
# regulation of (canonical) WNT signalling
# chemotaxis; cell adhesion; extracellular matrix organisation
# somite dev; segmentation; AP patterning

GO.somite$p.value_great <- go.dars[[1]][match(GO.somite$GO.ID, go.dars[[1]]$ID),]$Hyper_Raw_PValue
GO.somite$sig_great <- GO.somite$p.value_great<0.05
ggplot(GO.somite[!is.na(GO.somite$p.value_great),], aes(sig_great, p.value_classic, colour=sig_great)) +
  geom_boxplot() +
  scale_color_manual(values=c('TRUE'="red", 'FALSE'="black")) +
  xlab("term is significant with GREAT") +
  ylab("p-value of enrichment with DE genes") +
  th + theme(legend.position = "none")
```

![](01_somiteMaturation_files/figure-html/great-1.png)<!-- -->

Interestingly, GREAT also performs enrichment analysis of mouse gene KO phenotypes, and there is a very strong enrichment of skeleton abnormalities:


```r
## mouse phenotype
pheno.dars[['Mouse Phenotype Single KO']][c(1:3,12:14,28,30,63,64,68,81,103,104,121,133), 
                                          c(1,2,6,11:9,13)]
```

<div data-pagedtable="false">
  <script data-pagedtable-source type="application/json">
{"columns":[{"label":[""],"name":["_rn_"],"type":[""],"align":["left"]},{"label":["ID"],"name":[1],"type":["chr"],"align":["left"]},{"label":["name"],"name":[2],"type":["chr"],"align":["left"]},{"label":["Hyper_Fold_Enrichment"],"name":[3],"type":["dbl"],"align":["right"]},{"label":["Total_Genes_Annotated"],"name":[4],"type":["dbl"],"align":["right"]},{"label":["Hyper_Background_Gene_Hits"],"name":[5],"type":["dbl"],"align":["right"]},{"label":["Hyper_Foreground_Gene_Hits"],"name":[6],"type":["dbl"],"align":["right"]},{"label":["Hyper_Adjp_BH"],"name":[7],"type":["dbl"],"align":["right"]}],"data":[{"1":"MP:0005508","2":"abnormal skeleton morphology","3":"1.238860","4":"1832","5":"1815","6":"518","7":"9.914109e-09","_rn_":"1"},{"1":"MP:0005390","2":"skeleton phenotype","3":"1.226806","4":"1928","5":"1908","6":"535","7":"1.757182e-08","_rn_":"2"},{"1":"MP:0002114","2":"abnormal axial skeleton morphology","3":"1.291961","4":"1157","5":"1150","6":"356","7":"1.757182e-08","_rn_":"3"},{"1":"MP:0004624","2":"abnormal thoracic cage morphology","3":"1.430837","4":"432","5":"432","6":"138","7":"3.492609e-06","_rn_":"12"},{"1":"MP:0000137","2":"abnormal vertebrae morphology","3":"1.402894","4":"509","5":"506","6":"155","7":"3.647828e-06","_rn_":"13"},{"1":"MP:0008272","2":"abnormal endochondral bone ossification","3":"1.879190","4":"93","5":"92","6":"45","7":"7.553657e-06","_rn_":"14"},{"1":"MP:0008271","2":"abnormal bone ossification","3":"1.409746","4":"341","5":"337","6":"115","7":"3.108301e-04","_rn_":"28"},{"1":"MP:0004625","2":"abnormal rib joint","3":"1.834788","4":"73","5":"73","6":"32","7":"3.541209e-04","_rn_":"30"},{"1":"MP:0002113","2":"abnormal skeleton development","3":"1.309800","4":"382","5":"379","6":"137","7":"2.622729e-03","_rn_":"63"},{"1":"MP:0000155","2":"asymmetric rib joints","3":"2.012126","4":"39","5":"39","6":"19","7":"2.776540e-03","_rn_":"64"},{"1":"MP:0000150","2":"abnormal rib morphology","3":"1.358321","4":"328","5":"328","6":"99","7":"3.613192e-03","_rn_":"68"},{"1":"MP:0002759","2":"abnormal caudal vertebrae morphology","3":"1.581086","4":"113","5":"111","6":"39","7":"8.390097e-03","_rn_":"81"},{"1":"MP:0004986","2":"abnormal osteoblast morphology","3":"1.561818","4":"112","5":"111","6":"43","7":"1.281535e-02","_rn_":"103"},{"1":"MP:0003098","2":"decreased tendon stiffness","3":"4.239318","4":"5","5":"5","6":"3","7":"1.281535e-02","_rn_":"104"},{"1":"MP:0012181","2":"increased somite number","3":"3.976006","4":"3","5":"3","6":"2","7":"1.857033e-02","_rn_":"121"},{"1":"MP:0000459","2":"abnormal presacral vertebrae morphology","3":"1.324361","4":"280","5":"279","6":"83","7":"2.398504e-02","_rn_":"133"}],"options":{"columns":{"min":{},"max":[10]},"rows":{"min":[10],"max":[10]},"pages":{}}}
  </script>
</div>

These indicate the the DA regions are in close proximity to genes that, when inactivated, lead to abnormalities in the establishment and development of the body plan.

In summary, the DA regions of open chromatin are nearby genes that are involved in physiologically relevant processes, consistent with those overrepresented in the set of DE genes. Furthermore, some of the genes nearby DA regions, when knocked out, result in phenotypes involving anomalies in the body plan and skeleton.

#### Motif enrichment analysis

Given that many of the DA regions are distal to genes, it is likely that they are regulatory elements that mediate their actions through transcription factor (TF) binding. To investigate this, we use motif enrichment analysis to search for TFs that are overrepresented in the set of DA peaks, compared to non-DA peaks.

To identify motifs that occur more often in the DA vs nonDA peaks, we use `Analysis of Motif Enrichment` (AME) from the `MEME` suite (http://meme-suite.org/tools/ame). The sequences of the DA regions are compared to an equivalent set of non-DA regions that have the same length distribution.


```r
## retrieve the FASTA sequences of the DA regions
# all DA
write.table(paste0("samtools faidx ", dir, "REFERENCE/Mus_musculus.GRCm38.dna.primary_assembly.fa ",
                   substr(row.names(dars), 4, 50), " >> ", dir,
                   "RNA+ATAC/results/01_DAregions_somiteTrios.fasta"), 
            paste0(dir, "RNA+ATAC/scripts/dataProcessing/01_getDNAseq_DAregions_somiteTrios.sh"), 
            quote = FALSE, row.names = FALSE, col.names = FALSE)

## similar set of nonDA to use as a control instead of shuffled seqs (should have the same length distribution)
## since 2701 is a very small subset of all nonDA regions, we select 25 times as many, for a more representative background and enhanced statistical power
l <- dars$end-dars$start+1 ## DA regions lengths
freqs <- table(cut(l,6)) ## define 10 intervals and count number of seqs in each
intervals <- gsub("[(](.+),(.+)[]]","\\1-\\2",names(freqs))

peaks_nonDA <- peaks[peaks$name %in% setdiff(peaks$name, row.names(dars)),]
sel <- c() ## randomly sample nonDA regions from all peaks that have the same length distribution
for(i in 1:length(intervals)){
  limits <- as.numeric(unlist(strsplit(intervals[i], "-")))
  sel <- c(sel, sample(peaks_nonDA[width(peaks_nonDA) > limits[1] & width(peaks_nonDA) <= limits[2],],
                       freqs[i]*25))
}
control <- do.call('c', sel) # 40% the total peaks
# plot(density(l))
# lines(density(width(control)), col="red")

write.table(paste0("samtools faidx ", dir, "REFERENCE/Mus_musculus.GRCm38.dna.primary_assembly.fa ",
                   substr(control$name, 4, 50), " >> ", dir,
                   "RNA+ATAC/results/01_nonDA_controlRegions_x25_somiteTrios.fasta"), 
            paste0(dir, "RNA+ATAC/scripts/dataProcessing/01_getDNAseq_nonDA_controlRegions_somiteTrios.sh"), 
            quote = FALSE, row.names = FALSE, col.names = FALSE)

## submit to AME (http://meme-suite.org/tools/ame) to identify enriched TF motifs
## ran against the human and mouse HOCOMOCOv11_full motif databases
## all other parameters default
## downloaded results and saved in dir/RNA+ATAC/results/01_AME_DAregions_somiteTrios.tsv
```

Significantly enriched motifs are filtered to remove non-expressed TFs. Overall, a couple hundred different TFs are identified, although many of these are paralogs likely to identify very similar motifs.


```r
## motif annotation downloaded from http://hocomoco11.autosome.ru/downloads_v11
motif_ann <- read.table(paste0(dir, "REFERENCE/HOCOMOCOv11_full_annotation_MOUSE_mono.tsv"), 
                        sep="\t", header = TRUE, stringsAsFactors = FALSE)
tmp <- read.table(paste0(dir, "REFERENCE/HOCOMOCOv11_full_annotation_HUMAN_mono.tsv"), 
                  sep="\t", header = TRUE, stringsAsFactors = FALSE)
colnames(tmp) <- colnames(motif_ann)
## for human TFs, transform gene name to lower case to match mouse gene names (won't work for all, but good enough approximation)
tmp$Transcription.factor <- paste0(substr(tmp$Transcription.factor, 1, 1), 
                                   tolower(substr(tmp$Transcription.factor, 2, 10)))
motif_ann <- rbind(motif_ann, tmp)

## significantly enriched motifs in DA vs nonDA sequences
motifs_da <- read.table(paste0(dir, "RNA+ATAC/results/01_AME_DAregions_somiteTrios.tsv"),
                        sep = "\t", header = TRUE, stringsAsFactors = FALSE)[,-c(2,4)]

## remove TFs that are not expressed
not_expr <- which(motif_ann[motif_ann$Model %in% motifs_da$motif_ID,]$Transcription.factor %notin% geneCounts$gene)
not_expr <- motif_ann[motif_ann$Model %in% motifs_da$motif_ID,][not_expr,]$Model

motifs_da <- motifs_da[motifs_da$motif_ID %notin% not_expr,]

## annotate TF gene
motifs_da$TF <- motif_ann[match(motifs_da$motif_ID, motif_ann$Model),]$Transcription.factor
motifs_da$family <- motif_ann[match(motifs_da$motif_ID, motif_ann$Model),]$TF.family
motifs_da$subfamily <- motif_ann[match(motifs_da$motif_ID, motif_ann$Model),]$TF.subfamily
motifs_da$species <- ifelse(grepl("MOUSE", motifs_da$motif_ID), "mouse", "human")

# a few motifs are missing from the annotation
tmp <- substr(motifs_da[is.na(motifs_da$TF),]$motif_ID, 1,4)
tmp <- paste0(substr(tmp,1,1), tolower(substr(tmp,2,3)), substr(tmp,4,4))
motifs_da[is.na(motifs_da$TF),]$TF <- tmp

# when both human and mouse exist, remove human
duplicates <- motifs_da[duplicated(motifs_da$TF),]$TF
duplicates <- motifs_da[motifs_da$TF %in% duplicates & motifs_da$species == "human", ]
motifs_da <- motifs_da[!(motifs_da$rank %in% duplicates$rank),]

length(unique(motifs_da$TF)) # 201
```

```
## [1] 201
```

```r
# table(motifs_da$subfamily)
```

The top hits include many **Hox** and **Cdx** genes, as well as TFs from the **Forkhead box (Fox)**, **Sox-related**, **NK-related** and other **homeodomain** families.


```r
motifs_da$label <- ifelse(motifs_da$rank < 25 | motifs_da$X.TP/motifs_da$X.FP > 1.6,
                          motifs_da$TF, "")
ggplot(motifs_da, aes(log2(X.TP/X.FP), -log10(adj_p.value), label=label)) + 
  geom_point(aes(colour=species)) + 
  xlab("log2 fold-change (DA / nonDA)") +
  ylab("-log10(adj p-value)") +
  geom_text_repel(size = 3, max.overlaps = 200) +
  th
```

![](01_somiteMaturation_files/figure-html/motifs_enriched-1.png)<!-- -->

- **Hox** and **Cdx** factors are fundamental in the antero-posterior patterning of the embryo, and confer distinct identities to somites from different axial levels.
- [**Twist1**](https://doi.org/10.1038/cr.2011.144): an essential factor for mesoderm specification and differentiation, with well documented expression in the PSM and somites. It is embryonic lethal, and KO embryos show somite defects (among many others). It has well-known functions in the regulation of **epithelial-mesenchymal transition**, and is involved in maintaining mesenchymal cells in an undifferentiated state, thus blocking premature myogenesis and osteogenesis. Thus, it is not surprising that it is involved in the regulation of the somite maturation process.
- [**Msgn1**](https://doi.org/10.1242/dev.110908) referred to as a master regulator of paraxial presomitic mesoderm differentiation. It is involved in the differentiation of the PSM, but also in regulating the expression of [cyclic genes from the Notch pathway](https://doi.org/10.1038/ncomms1381). Thus, we would expect the DA regions with *Msgn1* binding sites to be closing during somite maturation.

Some of these are specifically found in DA regions that open or close as somites mature. As expected, the majority of *Twist1* matches occur in opening loci, whereas *Msgn1* sites are much more prevalent in closing regions. Other factors observed preferentially in opening regions are related to retinoic-acid signalling (COT1-2 correspond to TFs *Nr2f1-2*, which are activated by RA).


```r
motifs_da_hits <- read.table(paste0(dir, "RNA+ATAC/results/01_AME_DAregions_somiteTrios_hits.tsv"),
                             sep = "\t", header = TRUE, stringsAsFactors = FALSE)[,-c(1)]
motifs_da_hits <- motifs_da_hits[motifs_da_hits$motif_ID %in% motifs_da$motif_ID,]

motif_dyn_pct <- matrix(ncol=3, nrow=nrow(motifs_da))
row.names(motif_dyn_pct) <- motifs_da$motif_ID
for(i in motifs_da$motif_ID){
  tmp <- motifs_da_hits[motifs_da_hits$motif_ID ==  i,]
  n <- table(profiles.trios.atac[paste0("chr", tmp[tmp$class == "tp",]$seq_ID)])
  motif_dyn_pct[i,] <- c(sum(n[grep("down",names(n))]),
                         sum(n[grep("up",names(n))]), sum(tmp$class=="tp"))
}
motif_dyn_pct <- as.data.frame(motif_dyn_pct)
colnames(motif_dyn_pct) <- c("down", "up", "total")
motif_dyn_pct$down_pct <- motif_dyn_pct$down/motif_dyn_pct$total*100
motif_dyn_pct$up_pct <- motif_dyn_pct$up/motif_dyn_pct$total*100

## TFs that are biased to closing/opening regions
tmp <- motif_dyn_pct[motif_dyn_pct$down_pct > 60 | motif_dyn_pct$up_pct > 60,]
tmp[order(tmp$up_pct),4:5]
```

<div data-pagedtable="false">
  <script data-pagedtable-source type="application/json">
{"columns":[{"label":[""],"name":["_rn_"],"type":[""],"align":["left"]},{"label":["down_pct"],"name":[1],"type":["dbl"],"align":["right"]},{"label":["up_pct"],"name":[2],"type":["dbl"],"align":["right"]}],"data":[{"1":"68.18182","2":"28.95623","_rn_":"MSGN1_MOUSE.H11MO.0.C"},{"1":"33.99716","2":"60.45519","_rn_":"COT1_MOUSE.H11MO.1.C"},{"1":"34.33396","2":"60.60038","_rn_":"COT1_MOUSE.H11MO.0.B"},{"1":"31.97492","2":"64.57680","_rn_":"RXRG_HUMAN.H11MO.0.B"},{"1":"27.96935","2":"66.66667","_rn_":"COT2_MOUSE.H11MO.2.B"},{"1":"28.75536","2":"67.38197","_rn_":"COT2_MOUSE.H11MO.1.B"},{"1":"27.41433","2":"67.91277","_rn_":"NR2C2_MOUSE.H11MO.0.A"},{"1":"9.65251","2":"88.03089","_rn_":"TWST1_MOUSE.H11MO.0.B"}],"options":{"columns":{"min":{},"max":[10]},"rows":{"min":[10],"max":[10]},"pages":{}}}
  </script>
</div>


```r
motif_dyn_pct$z <- (motif_dyn_pct$up-motif_dyn_pct$down)/sqrt(sum(motif_dyn_pct$up, motif_dyn_pct$down))
motif_dyn_pct$pval <- motifs_da[match(row.names(motif_dyn_pct), motifs_da$motif_ID),]$adj_p.value
motif_dyn_pct$TF <- motifs_da[match(row.names(motif_dyn_pct), motifs_da$motif_ID),]$TF
motif_dyn_pct$da <- motifs_da[match(row.names(motif_dyn_pct), motifs_da$motif_ID),]$X.TP
motif_dyn_pct$nonDa <- motifs_da[match(row.names(motif_dyn_pct), motifs_da$motif_ID),]$X.FP
tmp$TF <- motifs_da[match(row.names(tmp), motifs_da$motif_ID),]$TF

ggplot(motif_dyn_pct, aes(z, da/nonDa)) +
  geom_point(pch=21, colour="grey70") +
  geom_point(data=motif_dyn_pct[abs(motif_dyn_pct$z)>0.3,], aes(fill=z),
             pch=21, colour="grey70") +
  scale_fill_gradientn(colours = rev(brewer.pal(n=9, "RdYlBu"))) +
  geom_vline(xintercept = 0, lty=2, colour="grey80") +
  geom_text_repel(data=motif_dyn_pct[motif_dyn_pct$da/motif_dyn_pct$nonDa>1.6 | 
                                       abs(motif_dyn_pct$z)>0.3,], 
                  aes(label=TF),
                  max.overlaps = 50) +
  xlab("z score") +
  ylab("peaks with TFBS\nDA / nonDA") +
  ylim(c(1,2)) +
  th 
```

![](01_somiteMaturation_files/figure-html/motifs_direction_plot-1.png)<!-- -->

Overall, these data suggest some of the important regulators driving somite maturation.

#### Chromatin accessibility variation

Using `chromVar` we have observed a few thousand genomic regions that show different accessibility between the somite trios, and the results from that analysis are consistent with the motif enrichments observed above.

For example, we identified *Cdx* factors as some of the most variable across samples, and this variability was strongly correlated with somite differentiation. Many of the strongly enriched motifs indeed show different accessibility levels between somite trios at sites predicted to be bound by the corresponding TFs. The motifs observed preferentially within opening regions indeed show this behaviour, and most of the others are instead in closing chromatin loci.


```r
dev <- readRDS(paste0(dir, "ATAC-seq/results/06_motif_deviationScores.Rds"))

## plot deviation scores for enriched TFs
df <- as.data.frame(colData(dev))
df$stage <- substr(df$stage, 6,7)
df$somite <- factor(df$somite, levels = c("SI", "SII", "SIII"))
plots <- list()
for(tf in c("Hoxb13", "Evx2", "Cdx2", "Nkx61", "Lhx9",
            "Meis1", "Twist1", "Rxra", "Nr2c1", "Nr2f2")){
  i <- grep(tf, row.names(dev))
  df$z <- assay(dev, 'z')[i,]
  
  plots[[tf]] <- ggplot(df, aes(somite, z, colour=somite)) +
    geom_boxplot() +
    scale_color_manual(values = cols.somite) +
    ggtitle(tf) +
    xlab("") +
    th
}
ggarrange(plotlist = plots, ncol=5, nrow=2, legend = "none", align = "h")
```

![](01_somiteMaturation_files/figure-html/chromVar-1.png)<!-- -->

Overall, we observe consistent results between the enrichment analysis and the accessibility pattern aggregated over motif instances in open chromatin.



```r
write.table(degs, paste0(dir, "RNA+ATAC/results/01_DEresults_summary_somiteTrios_withClusters.tsv"),
            quote = FALSE, sep="\t")
write.table(dars, paste0(dir, "RNA+ATAC/results/01_DAregions_summary_somiteTrios_withClusters.tsv"),
            quote = FALSE, sep="\t")

write.table(GO.somite, paste0(dir, "RNA+ATAC/results/01_GO_enrichment_DEGs_somiteTrios.tsv"),
            quote = FALSE, sep = "\t", row.names = FALSE)
save(geneTon, dds, degs_res, GO.somite_res, universe, 
     file = paste0(dir, "RNA+ATAC/results/01_geneTonic_objects.RData"))

write.table(go.dars, paste0(dir, "RNA+ATAC/results/01_GO_enrichment_DARs_somiteTrios.tsv"),
            quote = FALSE, sep = "\t", row.names = FALSE)
saveRDS(pheno.dars, paste0(dir, "RNA+ATAC/results/01_mousePhenotype_enrichment_DARs_somiteTrios.Rds"))

write.table(motif_dyn_pct, paste0(dir, "RNA+ATAC/results/01_motifEnrichment_somiteTrios.tsv"),
            quote = FALSE, sep="\t")
```



```r
sessionInfo()
```

```
## R version 4.1.0 (2021-05-18)
## Platform: x86_64-apple-darwin17.0 (64-bit)
## Running under: macOS Big Sur 10.16
## 
## Matrix products: default
## BLAS:   /Library/Frameworks/R.framework/Versions/4.1/Resources/lib/libRblas.dylib
## LAPACK: /Library/Frameworks/R.framework/Versions/4.1/Resources/lib/libRlapack.dylib
## 
## locale:
## [1] en_GB.UTF-8/en_GB.UTF-8/en_GB.UTF-8/C/en_GB.UTF-8/en_GB.UTF-8
## 
## attached base packages:
## [1] parallel  stats4    stats     graphics  grDevices utils     datasets 
## [8] methods   base     
## 
## other attached packages:
##  [1] dplyr_1.0.8                        DESeq2_1.32.0                     
##  [3] DEFormats_1.20.0                   GeneTonic_1.4.1                   
##  [5] pcaExplorer_2.18.0                 qvalue_2.24.0                     
##  [7] BSgenome.Mmusculus.UCSC.mm10_1.4.0 BSgenome_1.60.0                   
##  [9] rtracklayer_1.52.1                 Biostrings_2.60.2                 
## [11] XVector_0.32.0                     org.Mm.eg.db_3.13.0               
## [13] SummarizedExperiment_1.22.0        MatrixGenerics_1.4.3              
## [15] matrixStats_0.62.0                 chromVAR_1.14.0                   
## [17] dynamicTreeCut_1.63-1              RColorBrewer_1.1-3                
## [19] ggrepel_0.9.1                      ggpubr_0.4.0                      
## [21] ggplot2_3.3.5                      rGREAT_1.24.0                     
## [23] topGO_2.44.0                       SparseM_1.81                      
## [25] GO.db_3.13.0                       AnnotationDbi_1.54.1              
## [27] Biobase_2.52.0                     graph_1.70.0                      
## [29] GenomicRanges_1.44.0               GenomeInfoDb_1.28.4               
## [31] IRanges_2.26.0                     S4Vectors_0.30.2                  
## [33] BiocGenerics_0.38.0               
## 
## loaded via a namespace (and not attached):
##   [1] rappdirs_0.3.3              AnnotationForge_1.34.1     
##   [3] R.methodsS3_1.8.1           pkgmaker_0.32.2            
##   [5] tidyr_1.2.0                 bit64_4.0.5                
##   [7] knitr_1.38                  DelayedArray_0.18.0        
##   [9] R.utils_2.11.0              data.table_1.14.2          
##  [11] KEGGREST_1.32.0             TFBSTools_1.30.0           
##  [13] RCurl_1.98-1.6              doParallel_1.0.17          
##  [15] generics_0.1.2              cowplot_1.1.1              
##  [17] RSQLite_2.2.12              bit_4.0.4                  
##  [19] tzdb_0.3.0                  webshot_0.5.3              
##  [21] xml2_1.3.3                  httpuv_1.6.5               
##  [23] assertthat_0.2.1            DirichletMultinomial_1.34.0
##  [25] viridis_0.6.2               xfun_0.30                  
##  [27] hms_1.1.1                   jquerylib_0.1.4            
##  [29] evaluate_0.15               promises_1.2.0.1           
##  [31] TSP_1.2-0                   fansi_1.0.3                
##  [33] restfulr_0.0.13             progress_1.2.2             
##  [35] caTools_1.18.2              dendextend_1.15.2          
##  [37] dbplyr_2.1.1                Rgraphviz_2.36.0           
##  [39] igraph_1.3.0                DBI_1.1.2                  
##  [41] geneplotter_1.70.0          htmlwidgets_1.5.4          
##  [43] purrr_0.3.4                 ellipsis_0.3.2             
##  [45] crosstalk_1.2.0             backports_1.4.1            
##  [47] annotate_1.70.0             gridBase_0.4-7             
##  [49] biomaRt_2.48.3              vctrs_0.4.1                
##  [51] Cairo_1.5-15                abind_1.4-5                
##  [53] cachem_1.0.6                withr_2.5.0                
##  [55] ggforce_0.3.3               checkmate_2.0.0            
##  [57] GenomicAlignments_1.28.0    prettyunits_1.1.1          
##  [59] cluster_2.1.3               backbone_2.0.3             
##  [61] lazyeval_0.2.2              seqLogo_1.58.0             
##  [63] crayon_1.5.1                genefilter_1.74.1          
##  [65] labeling_0.4.2              edgeR_3.34.1               
##  [67] pkgconfig_2.0.3             tweenr_1.0.2               
##  [69] seriation_1.3.5             rlang_1.0.6                
##  [71] lifecycle_1.0.1             miniUI_0.1.1.1             
##  [73] colourpicker_1.1.1          registry_0.5-1             
##  [75] filelock_1.0.2              BiocFileCache_2.0.0        
##  [77] GOstats_2.58.0              polyclip_1.10-0            
##  [79] rngtools_1.5.2              Matrix_1.4-1               
##  [81] carData_3.0-5               base64enc_0.1-3            
##  [83] GlobalOptions_0.1.2         pheatmap_1.0.12            
##  [85] png_0.1-7                   viridisLite_0.4.0          
##  [87] rjson_0.2.21                bitops_1.0-7               
##  [89] shinydashboard_0.7.2        R.oo_1.24.0                
##  [91] visNetwork_2.1.0            blob_1.2.3                 
##  [93] shape_1.4.6                 rintrojs_0.3.0             
##  [95] stringr_1.4.0               readr_2.1.2                
##  [97] rstatix_0.7.0               shinyAce_0.4.1             
##  [99] ggsignif_0.6.3              CNEr_1.28.0                
## [101] scales_1.2.0                memoise_2.0.1              
## [103] GSEABase_1.54.0             magrittr_2.0.3             
## [105] plyr_1.8.7                  zlibbioc_1.38.0            
## [107] threejs_0.3.3               compiler_4.1.0             
## [109] BiocIO_1.2.0                clue_0.3-60                
## [111] Rsamtools_2.8.0             cli_3.2.0                  
## [113] Category_2.58.0             MASS_7.3-56                
## [115] tidyselect_1.1.2            stringi_1.7.6              
## [117] shinyBS_0.61.1              highr_0.9                  
## [119] yaml_2.3.5                  locfit_1.5-9.5             
## [121] grid_4.1.0                  sass_0.4.1                 
## [123] tools_4.1.0                 circlize_0.4.14            
## [125] rstudioapi_0.13             TFMPvalue_0.0.8            
## [127] foreach_1.5.2               gridExtra_2.3              
## [129] farver_2.1.0                digest_0.6.29              
## [131] shiny_1.7.1                 pracma_2.3.8               
## [133] Rcpp_1.0.8.3                car_3.0-12                 
## [135] broom_0.8.0                 later_1.3.0                
## [137] shinyWidgets_0.6.4          httr_1.4.2                 
## [139] ComplexHeatmap_2.8.0        colorspace_2.0-3           
## [141] XML_3.99-0.9                splines_4.1.0              
## [143] RBGL_1.68.0                 expm_0.999-6               
## [145] plotly_4.10.0               xtable_1.8-4               
## [147] jsonlite_1.8.0              poweRlaw_0.70.6            
## [149] heatmaply_1.3.0             R6_2.5.1                   
## [151] pillar_1.7.0                htmltools_0.5.2            
## [153] mime_0.12                   NMF_0.24.0                 
## [155] glue_1.6.2                  fastmap_1.1.0              
## [157] DT_0.22                     BiocParallel_1.26.2        
## [159] bs4Dash_2.0.3               codetools_0.2-18           
## [161] utf8_1.2.2                  lattice_0.20-45            
## [163] bslib_0.3.1                 tibble_3.1.6               
## [165] curl_4.3.2                  gtools_3.9.2               
## [167] survival_3.3-1              limma_3.48.3               
## [169] rmarkdown_2.13              munsell_0.5.0              
## [171] GetoptLong_1.0.5            GenomeInfoDbData_1.2.6     
## [173] iterators_1.0.14            shinycssloaders_1.0.0      
## [175] reshape2_1.4.4              gtable_0.3.0
```

