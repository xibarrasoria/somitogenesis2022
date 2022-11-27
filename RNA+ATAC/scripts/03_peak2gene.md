---
title: "Regulatory control of somitogenesis"
date: '18 August, 2022'
output:
  html_document:
    keep_md: true
    fig_width: 5
    fig_height: 5
    fig_caption: yes
    df_print: paged
    code_folding: hide
    toc: true
    toc_float: 
      collapsed: false
---



Having matched RNA and ATAC-seq profiles provides an opportunity to characterise the regulatory mechanisms that control somitogenesis.


```r
## expression data
meta.rna <- read.table(paste0(dir, "RNA-seq/data/metadata_RNAseq.tsv"), 
                       stringsAsFactors = FALSE, header = TRUE)
meta.rna <- meta.rna[meta.rna$QC == 1 & meta.rna$wrongStage == 0,]
meta.rna$stage <- paste0("stage", meta.rna$stage)
meta.rna$somite <- factor(meta.rna$somite, levels=c("SIII","SII","SI"))

data.rna <- read.table(paste0(dir, "RNA-seq/data/geneCounts.NORM_batchCorrected_14PCs.tsv"),
                       check.names = FALSE)
data.rna <- data.rna[,c('gene', meta.rna$sample)]
universe <- data.rna$gene
names(universe) <- row.names(data.rna)

## accessibility data
meta.atac <- read.table(paste0(dir, "ATAC-seq/data/metadata_ATACseq_GQ.tsv"), 
                   stringsAsFactors = FALSE, header = TRUE)
meta.atac <- meta.atac[meta.atac$QCpass==1,]
meta.atac$stage <- paste0("stage", meta.atac$stage)
meta.atac$somite <- factor(meta.atac$somite, levels=c("SIII","SII","SI"))

data.atac <- readRDS(paste0(dir, "ATAC-seq/results/04_peakCounts_csawMerged.NORM.batchCorrected_18PCs.Rds"))
```

### Peak to gene linking

One way to do this is to identify chromatin loci in the vicinity of DE genes that show coordinated activity patterns; these loci can be postulated as putative regulatory elements driving the changes in gene expression. We can exploit the variation across development, which usually involves large changes in activity, to find gene-peak pairs using correlation analysis. To determine whether a given correlation score is significantly larger than expected by chance, we compare correlation values of the same gene against peaks that are in other chromosomes, which do not contribute to *cis* regulatory mechanisms. 

For this analysis, we use the method implemented in `FigR` (with minor modifications), which uses `chromVAR` to determine a set of background peaks matched for accessibility and GC content to model the null distribution of correlation scores for each gene-peak pair. Naturally, the analysis is restricted to the 43 samples where we generated good-quality libraries for both the RNA and ATAC datasets. We look for correlated peaks within 100kb of a DE gene.



There are 15734 gene-peak pairs that are significantly correlated. These involve nearly 6K DE genes, linked to ~13.5K peaks.


```r
c(n_genes = length(unique(links$Gene)),
  n_peaks = length(unique(links$PeakRanges)))
```

```
## n_genes n_peaks 
##    5909   13657
```

Most gene-peak pairs show moderate correlation values (median=0.382547). However, a few pairs with very low correlation coefficients are deemed significant by the method. These correspond to cases where gene expression and accessibility are not associated (correlation values close to 0) but the set of background peaks used to compute significance are biased towards negative correlations. This leads to a significant p-value even if the correlation of the pair in question is not significant itself.


```r
ggplot(links, aes(1, rObs)) +
  geom_violin() +
  geom_boxplot(width=0.1) +
  xlab("") +
  ylab("Spearman correlation") +
  th + theme(axis.text.x = element_blank(),
             axis.ticks.x = element_blank())
```

![](03_peak2gene_files/figure-html/links_corr-1.png)<!-- -->

Thus, to consider a gene-peak pair as significantly correlated we require a minimum correlation coefficient of 0.3 on top of a significant p-value, taking the number of links to 12.8K.


```r
links <- links[links$rObs>0.3,]
nrow(links)
```

```
## [1] 12803
```

Although most genes are linked to a single or a few peaks, a small proportion are connected to more than 5 peaks (349 genes). In contrast, the vast majority of peaks (88.91%) regulate a single gene.


```r
par(mfrow=c(1,2))
plot(table(table(links$Gene)), bty="l",
     xlab="number of peaks per gene",
     ylab="number of genes")
plot(table(table(links$PeakRanges)), bty="l",
     xlab="number of genes per peak",
     ylab="number of peaks")
```

![](03_peak2gene_files/figure-html/links_per_feature-1.png)<!-- -->

```r
# summary(as.numeric(table(links$Gene)))
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 1.00    1.00    2.00    2.452   3.000  30.000  
# summary(as.numeric(table(links$PeakRanges)))
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 1.000   1.000   1.000   1.147   1.000   8.000 
```

The set of defined pairs involve around half of all the DE genes, for both those that are differential across somite trios, and across stages. 


```r
## check how many of the correlated regions are differential
dars.somite <- read.table(paste0(dir, "RNA+ATAC/results/01_DAregions_summary_somiteTrios_withClusters.tsv"))
dars.stage <- read.table(paste0(dir, "RNA+ATAC/results/02_DAregions_summary_stage_fate.tsv"))
dars <- union(row.names(dars.somite), row.names(dars.stage))

links$DE <- ifelse(links$Gene %in% degs.somite$gene & links$Gene %in% degs.stage$gene, "both",
                   ifelse(links$Gene %in% degs.somite$gene, "somite",
                          ifelse(links$Gene %in% degs.stage$gene, "stage", "nonDE")))
links$DE_stage <- degs.stage[match(links$Gene, degs.stage$gene),]$fate 
# the same proportion of DE genes in each fate are linked

links$DA <- ifelse(links$PeakRanges %in% row.names(dars.somite) & links$PeakRanges %in% row.names(dars.stage), "both",
                   ifelse(links$PeakRanges %in% row.names(dars.somite), "somite",
                          ifelse(links$PeakRanges %in% row.names(dars.stage), "stage", "nonDA")))
links$DA_stage <- dars.stage[match(links$PeakRanges, row.names(dars.stage)),]$fate 

df <- data.frame(n_DE_linked = c(length(unique(links[links$DE %in% c("both", "somite"),]$Gene)),
                                 length(unique(links[links$DE %in% c("both", "stage"),]$Gene))),
                 DE_pct = c(round(length(unique(links[links$DE %in% c("both", "somite"),]$Gene))/nrow(degs.somite)*100,2),
                            round(length(unique(links[links$DE %in% c("both", "stage"),]$Gene))/nrow(degs.stage)*100,2)),
                 n_DA_linked = c(length(unique(links[links$DA %in% c("both", "somite"),]$PeakRanges)),
                                 length(unique(links[links$DA %in% c("both", "stage"),]$PeakRanges))),
                 DA_pct = c(round(length(unique(links[links$DA %in% c("both", "somite"),]$PeakRanges))/nrow(dars.somite)*100,2),
                            round(length(unique(links[links$DA %in% c("both", "stage"),]$PeakRanges))/nrow(dars.stage)*100,2)))
row.names(df) <- c("somite", "stage")
df
```

<div data-pagedtable="false">
  <script data-pagedtable-source type="application/json">
{"columns":[{"label":[""],"name":["_rn_"],"type":[""],"align":["left"]},{"label":["n_DE_linked"],"name":[1],"type":["int"],"align":["right"]},{"label":["DE_pct"],"name":[2],"type":["dbl"],"align":["right"]},{"label":["n_DA_linked"],"name":[3],"type":["int"],"align":["right"]},{"label":["DA_pct"],"name":[4],"type":["dbl"],"align":["right"]}],"data":[{"1":"1530","2":"51.39","3":"294","4":"10.88","_rn_":"somite"},{"1":"5059","2":"47.32","3":"3317","4":"10.05","_rn_":"stage"}],"options":{"columns":{"min":{},"max":[10]},"rows":{"min":[10],"max":[10]},"pages":{}}}
  </script>
</div>

Although some peaks are nearby their putative target, the majority are dozens of kilobases away.


```r
## TSS coordinates
TSSg <- FigR::mm10TSSRanges

# distance between gene and peak
tmp <- unlist(lapply(strsplit(links$PeakRanges, ":"), '[[', 2))
links.gr <- GenomicInteractions(GRanges(seqnames = unlist(lapply(strsplit(links$PeakRanges, ":"), '[[', 1)),
                                        IRanges(start = as.numeric(unlist(lapply(strsplit(tmp, "-"), '[[', 1))),
                                                end = as.numeric(unlist(lapply(strsplit(tmp, "-"), '[[', 2))))),
                                TSSg[match(links$Gene, TSSg$gene_name),])
distance <- distance(anchorOne(links.gr), anchorTwo(links.gr))
links$distance <- distance

bins <- data.frame(distance=c('0-200','200-2e3','2e3-5e3','5e3-10e3','10e3-50e3','50e3-100e3'),
                   n=c(sum(distance <= 200),
                       sum(distance > 200 & distance <= 2e3),
                       sum(distance > 2e3 & distance <= 5e3),
                       sum(distance > 5e3 & distance <= 10e3),
                       sum(distance > 10e3 & distance <= 50e3),
                       sum(distance > 50e3 & distance <= 100e3)))
bins$n_pct <- round(bins$n/length(distance)*100,2)
bins
```

<div data-pagedtable="false">
  <script data-pagedtable-source type="application/json">
{"columns":[{"label":["distance"],"name":[1],"type":["chr"],"align":["left"]},{"label":["n"],"name":[2],"type":["int"],"align":["right"]},{"label":["n_pct"],"name":[3],"type":["dbl"],"align":["right"]}],"data":[{"1":"0-200","2":"736","3":"5.75"},{"1":"200-2e3","2":"573","3":"4.48"},{"1":"2e3-5e3","2":"571","3":"4.46"},{"1":"5e3-10e3","2":"850","3":"6.64"},{"1":"10e3-50e3","2":"5027","3":"39.26"},{"1":"50e3-100e3","2":"5046","3":"39.41"}],"options":{"columns":{"min":{},"max":[10]},"rows":{"min":[10],"max":[10]},"pages":{}}}
  </script>
</div>

The median distance between linked peaks and genes is of 37.845kb.

Supporting the regulatory potential of the linked peaks, a larger proportion overlap ENCODE cCREs characterised as enhancers (H3K27ac-high, H3K4me3-low) compared to all peaks detected in somites, and all peaks within 100kb of a DE gene. Furthermore, the overlap greatly increases if we restrict to links with the strongest correlations.


```r
encode <- read.table(paste0(dir, "ATAC-seq/results/05_peaks_overlapExternalData.tsv"), sep="\t")
row.names(encode) <- paste0("chr", row.names(encode))

encode$encode <- ifelse(encode$pELS | encode$dELS, "enhancer",
                        ifelse(encode$PLS, "promoter",
                               ifelse(encode$CTCF_bound, "CTCFonly",
                                      ifelse(encode$H3K4me3, "H3K4me3only", "not_in_encode"))))
links$encode <- encode[match(links$PeakRanges, row.names(encode)),]$encode
links$fantom <- encode[match(links$PeakRanges, row.names(encode)),]$FANTOM
links$vista <- encode[match(links$PeakRanges, row.names(encode)),]$VISTA

links$encode <- as.factor(links$encode)

# round(rbind(all_peaks = prop.table(table(encode$encode)[1:4])*100,
#             peaks_proximal_to_DEGs = prop.table(table(encode[unique(cisCor$PeakRanges),]$encode)[1:4])*100,
#             peaks_linked_to_DEGs = prop.table(table(links$encode)[1:4])*100,
#             linked_0.5 = prop.table(table(links[links$rObs>0.5,]$encode)[1:4])*100,
#             linked_0.65 = prop.table(table(links[links$rObs>0.65,]$encode)[1:4])*100), 2)[,c(4,2,1,3)]

df <- data.frame(type = rep(c("enhancer", "not_in_encode"), each=7),
                 prop = c(sum(encode$encode=="enhancer")/nrow(encode)*100,
                          sum(encode[unique(cisCor$PeakRanges),]$encode=="enhancer")/nrow(encode[unique(cisCor$PeakRanges),])*100,
                          sum(links[links$rObs>0.3,]$encode=="enhancer")/nrow(links[links$rObs>0.3,])*100,
                          sum(links[links$rObs>0.4,]$encode=="enhancer")/nrow(links[links$rObs>0.4,])*100,
                          sum(links[links$rObs>0.5,]$encode=="enhancer")/nrow(links[links$rObs>0.5,])*100,
                          sum(links[links$rObs>0.6,]$encode=="enhancer")/nrow(links[links$rObs>0.6,])*100,
                          sum(links[links$rObs>0.7,]$encode=="enhancer")/nrow(links[links$rObs>0.7,])*100,
                          sum(encode$encode=="not_in_encode")/nrow(encode)*100,
                          sum(encode[unique(cisCor$PeakRanges),]$encode=="not_in_encode")/nrow(encode[unique(cisCor$PeakRanges),])*100,
                          sum(links[links$rObs>0.3,]$encode=="not_in_encode")/nrow(links[links$rObs>0.3,])*100,
                          sum(links[links$rObs>0.4,]$encode=="not_in_encode")/nrow(links[links$rObs>0.4,])*100,
                          sum(links[links$rObs>0.5,]$encode=="not_in_encode")/nrow(links[links$rObs>0.5,])*100,
                          sum(links[links$rObs>0.6,]$encode=="not_in_encode")/nrow(links[links$rObs>0.6,])*100,
                          sum(links[links$rObs>0.7,]$encode=="not_in_encode")/nrow(links[links$rObs>0.7,])*100))

ggplot(df, aes(rep(1:7,2), prop, colour=type)) +
  geom_point(size=2) +
  geom_line() +
  xlab("") +
  ylab("% of peaks") +
  ylim(c(0,100)) +
  scale_x_continuous(breaks=1:7,
                     labels=c("all peaks", "peaks near DEGs", paste0("linked peaks > 0.", 3:7))) +
  th + theme(axis.text.x = element_text(angle=45, hjust=1),
             panel.grid.major.y = element_line(size=0.5, colour="grey90"),
             panel.grid.minor.y = element_line(size=0.25, colour="grey95"))
```

![](03_peak2gene_files/figure-html/encode-1.png)<!-- -->

Similarly, while only ~10% of all peaks overlap with FANTOM5 enhancers, the number increases to 14% for linked peaks and this goes up as the correlation threshold becomes more stringent.


```r
df <- data.frame(type = "in_FANTOM",
                 prop = c(sum(encode$FANTOM)/nrow(encode)*100,
                          sum(encode[unique(cisCor$PeakRanges),]$FANTOM)/nrow(encode[unique(cisCor$PeakRanges),])*100,
                          sum(links[links$rObs>0.3,]$fantom)/nrow(links[links$rObs>0.3,])*100,
                          sum(links[links$rObs>0.4,]$fantom)/nrow(links[links$rObs>0.4,])*100,
                          sum(links[links$rObs>0.5,]$fantom)/nrow(links[links$rObs>0.5,])*100,
                          sum(links[links$rObs>0.6,]$fantom)/nrow(links[links$rObs>0.6,])*100,
                          sum(links[links$rObs>0.7,]$fantom)/nrow(links[links$rObs>0.7,])*100))

ggplot(df, aes(1:7, prop, colour=type)) +
  geom_point(size=2) +
  geom_line() +
  xlab("") +
  ylab("% of peaks") +
  ylim(c(0,100)) +
  scale_x_continuous(breaks=1:7,
                     labels=c("all peaks", "peaks near DEGs", paste0("linked peaks > 0.", 3:7))) +
  th + theme(axis.text.x = element_text(angle=45, hjust=1),
             panel.grid.major.y = element_line(size=0.5, colour="grey90"),
             panel.grid.minor.y = element_line(size=0.25, colour="grey95"))
```

![](03_peak2gene_files/figure-html/fantom-1.png)<!-- -->

### Highly-regulated genes

As mentioned above, although most genes are linked only to a single or few peaks, a few hundred have many more connections. The authors of `FigR` coined the term [*domains of regulatory chromatin (DORCs)*](https://doi.org/10.1016/j.cell.2020.09.056) to refer to these loci of high interactivity, and showed that they "are enriched for lineage-determining genes and overlap with known super-enhancers". In our case, many of the genes linked to dozens of peaks correspond to late Hox genes, which is expected since it has been shown that Hox expression regulation relies on coordinated chromatin remodeling. Many other genes are identified as well. 


```r
# Determine DORC genes
dorcGenes <- dorcJPlot(dorcTab = links,
                       cutoff = 6,
                       returnGeneList = TRUE)
```

![](03_peak2gene_files/figure-html/dorcs-1.png)<!-- -->

This restricted set of genes with a large number of linked peaks are still enriched for key processes relevant for somitogenesis, such as patterning and the development and morphogenesis of the different somite derivative lineages. These data suggest that the expression of genes that are crucial in the development and differentiation of somites are under intricate regulatory control.


```r
go <- topGOtable(DEgenes = dorcGenes,
                 BGgenes = universe,
                 topGO_method2 = "elim",
                 ontology = "BP",
                 geneID = "symbol",
                 fullNamesInRows = FALSE,
                 addGeneToTerms = TRUE,
                 mapping = "org.Mm.eg.db",
                 topTablerows = 100)
go[,c(2,4,7,9)]
```

<div data-pagedtable="false">
  <script data-pagedtable-source type="application/json">
{"columns":[{"label":[""],"name":["_rn_"],"type":[""],"align":["left"]},{"label":["Term"],"name":[1],"type":["chr"],"align":["left"]},{"label":["Significant"],"name":[2],"type":["int"],"align":["right"]},{"label":["p.value_elim"],"name":[3],"type":["dbl"],"align":["right"]},{"label":["genes"],"name":[4],"type":["chr"],"align":["left"]}],"data":[{"1":"anterior/posterior pattern specification","2":"33","3":"4.10e-14","4":"Alx4,Btg2,Cdx1,Cdx2,Celsr1,Dll1,Foxb1,Gbx2,Hey1,Hoxa10,Hoxa11,Hoxa2,Hoxa4,Hoxa9,Hoxb1,Hoxb3,Hoxc10,Hoxc11,Hoxc13,Hoxc4,Hoxc5,Hoxc6,Hoxc8,Hoxd10,Hoxd11,Hoxd13,Hoxd8,Lef1,Mesp1,Tbx3,Tcf15,Wnt5a,Wt1","_rn_":"1"},{"1":"embryonic digit morphogenesis","2":"11","3":"4.80e-07","4":"Alx4,Bcl2l11,Hoxa11,Hoxc11,Hoxd11,Hoxd12,Hoxd13,Sall1,Tbx2,Tbx3,Wnt5a","_rn_":"2"},{"1":"proximal/distal pattern formation","2":"8","3":"7.30e-07","4":"Dll1,Hoxa10,Hoxa11,Hoxa9,Hoxc10,Hoxc11,Hoxd10,Hoxd11","_rn_":"3"},{"1":"positive regulation of transcription by ...","2":"47","3":"6.30e-06","4":"Ablim1,Acvrl1,Alx4,Atf3,Bcl11b,Bmp7,Cdx1,Crtc1,Dab2ip,Dll1,E2f2,Ebf2,Eng,Etv4,Fgfr1,Gbx2,Hey1,Hoxa10,Hoxa13,Hoxa2,Hoxa9,Hoxb1,Hoxb3,Hoxc10,Hoxc11,Hoxc13,Hoxc4,Hoxd10,Hoxd13,Hoxd8,Jund,Klf1,Lef1,Meis2,Mesp1,Nfic,Plagl1,Sall1,Sall4,Slc38a3,Smad5,Smyd3,Ssbp4,Tcf15,Tox2,Wnt5a,Wt1","_rn_":"4"},{"1":"embryonic forelimb morphogenesis","2":"7","3":"1.20e-05","4":"Alx4,Crabp2,Hoxa11,Hoxa13,Hoxa9,Hoxd11,Tbx3","_rn_":"5"},{"1":"epithelial cell development","2":"19","3":"1.50e-05","4":"Bcl11b,Cdh5,Col18a1,Dact2,Dll1,Eng,Epha2,Fgfr1,Flna,Hoxa13,Hoxb13,Id1,Rab25,S1pr3,Src,Tjp3,Whrn,Wnt5a,Wt1","_rn_":"6"},{"1":"organ induction","2":"6","3":"4.10e-05","4":"Fgfr1,Hoxa11,Hoxc11,Hoxd11,Mesp1,Wnt5a","_rn_":"7"},{"1":"vasculogenesis","2":"10","3":"4.20e-05","4":"Ackr3,Aplnr,Egfl7,Eng,Epha2,Fgfr1,Hdac7,Hey1,Hoxa13,Wt1","_rn_":"8"},{"1":"axon guidance","2":"18","3":"6.00e-05","4":"Ablim1,Bcl11b,Bmp7,Crmp1,Cxcl12,Draxin,Efna2,Epha2,Ephb2,Etv4,Gbx1,Gbx2,Hoxa2,Prtg,Sema3b,Sema6c,Vstm2l,Wnt5a","_rn_":"9"},{"1":"embryonic limb morphogenesis","2":"21","3":"7.90e-05","4":"Alx4,Bcl2l11,Bmp7,Crabp2,Fgfr1,Hoxa10,Hoxa11,Hoxa13,Hoxa9,Hoxc10,Hoxc11,Hoxd10,Hoxd11,Hoxd12,Hoxd13,Lef1,Sall1,Sall4,Tbx2,Tbx3,Wnt5a","_rn_":"10"},{"1":"negative regulation of Notch signaling p...","2":"6","3":"8.90e-05","4":"Bmp7,Dlk1,Dll1,Egfl7,Gdpd5,Hey1","_rn_":"11"},{"1":"uterus development","2":"5","3":"9.10e-05","4":"Hoxa10,Hoxa11,Hoxa9,Src,Wnt5a","_rn_":"12"},{"1":"olfactory bulb interneuron differentiati...","2":"4","3":"1.10e-04","4":"Fgfr1,Sall1,Sall3,Wnt5a","_rn_":"13"},{"1":"embryonic skeletal joint morphogenesis","2":"4","3":"1.60e-04","4":"Bmp7,Hoxa11,Hoxc11,Hoxd11","_rn_":"14"},{"1":"embryonic skeletal system morphogenesis","2":"13","3":"1.90e-04","4":"Alx4,Bmp7,Hoxa11,Hoxa2,Hoxa4,Hoxa9,Hoxb1,Hoxb3,Hoxc11,Hoxc4,Hoxd10,Hoxd11,Hspg2","_rn_":"15"},{"1":"branching involved in mammary gland duct...","2":"5","3":"2.20e-04","4":"Epha2,Etv4,Src,Tbx3,Wnt5a","_rn_":"16"},{"1":"branching involved in prostate gland mor...","2":"4","3":"2.20e-04","4":"Bmp7,Hoxa13,Hoxb13,Hoxd13","_rn_":"17"},{"1":"skeletal system development","2":"39","3":"2.70e-04","4":"Adamts12,Alx4,Bmp7,Cdx1,Dlk1,Epha2,Fgf18,Fgfr1,Hotair,Hoxa10,Hoxa11,Hoxa2,Hoxa4,Hoxa9,Hoxb1,Hoxb3,Hoxc10,Hoxc11,Hoxc4,Hoxc5,Hoxc6,Hoxc8,Hoxd1,Hoxd10,Hoxd11,Hoxd12,Hoxd13,Hoxd8,Hspg2,Ltbp3,Mmp9,Nab2,Pkdcc,Smad5,Snai1,Src,Tbx3,Tcf15,Wnt5a","_rn_":"18"},{"1":"somitogenesis","2":"8","3":"2.80e-04","4":"Cdx1,Cdx2,Dll1,Foxb1,Lef1,Mesp1,Tcf15,Wnt5a","_rn_":"19"},{"1":"heart looping","2":"7","3":"3.00e-04","4":"Dll1,Eng,Lbx1,Mesp1,Tbx2,Tbx3,Wnt5a","_rn_":"20"},{"1":"cell surface receptor signaling pathway ...","2":"5","3":"3.20e-04","4":"Dll1,Hey1,Mesp1,Snai1,Wnt5a","_rn_":"21"},{"1":"negative regulation of transcription by ...","2":"34","3":"3.60e-04","4":"Atf3,Bach2,Bcl11a,Btg2,Cbfa2t3,Cdx2,Chd3,Dab2ip,Ehmt2,Eng,Fgfr1,Foxp4,Hdac7,Hey1,Hoxa2,Hoxb13,Hoxb3,Hoxc8,Hoxd8,Id1,Id3,Lef1,Mad2l2,Mxi1,Nfic,Prdm6,Sall1,Sall4,Smad5,Snai1,Tbx2,Tbx3,Wt1,Zfp217","_rn_":"22"},{"1":"angiogenesis","2":"31","3":"4.50e-04","4":"Ackr3,Acvrl1,Adgra2,Angpt2,Aplnr,Card10,Cdh5,Col18a1,Cxcl12,Dab2ip,Dll1,E2f2,Egfl7,Emilin1,Eng,Epha2,Ephb2,Fgf18,Fgfr1,Flna,Gbx2,Hdac7,Hspg2,Id1,Lef1,Mmp9,Shb,Smad5,Tbxa2r,Tnfaip2,Wnt5a","_rn_":"23"},{"1":"regulation of morphogenesis of an epithe...","2":"7","3":"4.60e-04","4":"Adamts12,Bmp7,Fgfr1,Hoxd13,Sall1,Tbx2,Wnt5a","_rn_":"24"},{"1":"hindlimb morphogenesis","2":"6","3":"5.20e-04","4":"Alx4,Hotair,Hoxd10,Sall1,Sall3,Tbx3","_rn_":"25"},{"1":"muscle organ development","2":"17","3":"5.30e-04","4":"Alx4,Atf3,Btg2,Dll1,Eng,Hdac7,Hey1,Hoxd10,Id3,Lbx1,Lef1,Lmna,Meg3,Plagl1,Tcf15,Wnt5a,Wt1","_rn_":"26"},{"1":"forelimb morphogenesis","2":"10","3":"5.70e-04","4":"Alx4,Crabp2,Hoxa11,Hoxa13,Hoxa9,Hoxd10,Hoxd11,Sall1,Sall3,Tbx3","_rn_":"27"},{"1":"negative regulation of epithelial cell d...","2":"8","3":"6.00e-04","4":"Acvrl1,Dll1,Foxp4,Id1,Mmp9,S1pr3,Sall1,Tbx3","_rn_":"28"},{"1":"metanephros development","2":"8","3":"6.50e-04","4":"Bmp7,Hoxa11,Hoxc11,Hoxd11,Id3,Kif26b,Sall1,Wt1","_rn_":"29"},{"1":"BMP signaling pathway","2":"13","3":"6.80e-04","4":"Acvrl1,Bmp7,Cdh5,Eng,Hoxa13,Htra1,Id1,Itga3,Lef1,Rgmb,Slc39a5,Smad5,Wnt5a","_rn_":"30"},{"1":"embryonic skeletal system development","2":"18","3":"7.30e-04","4":"Alx4,Bmp7,Dlk1,Hoxa11,Hoxa2,Hoxa4,Hoxa9,Hoxb1,Hoxb3,Hoxc11,Hoxc4,Hoxc5,Hoxc6,Hoxd1,Hoxd10,Hoxd11,Hspg2,Wnt5a","_rn_":"31"},{"1":"spinal cord motor neuron differentiation","2":"5","3":"7.30e-04","4":"Gbx1,Gdpd5,Hoxc10,Hoxd10,Lbx1","_rn_":"32"},{"1":"negative regulation of osteoblast differ...","2":"6","3":"7.40e-04","4":"Dlk1,Fgfr1,Hdac7,Hoxa2,Id1,Id3","_rn_":"33"},{"1":"negative regulation of ERK1 and ERK2 cas...","2":"7","3":"7.40e-04","4":"Atf3,Csk,Dab2ip,Dusp3,Emilin1,Ephb2,Spry4","_rn_":"34"},{"1":"positive regulation of chondrocyte diffe...","2":"4","3":"7.90e-04","4":"Fgf18,Hoxa11,Hoxd11,Pkdcc","_rn_":"35"},{"1":"odontogenesis of dentin-containing tooth","2":"7","3":"9.70e-04","4":"Bcl11b,Bcl2l11,Bmp7,Dll1,Dlx3,Lef1,Nfic","_rn_":"36"},{"1":"digestive tract development","2":"9","3":"9.90e-04","4":"Alx4,Cdx2,Foxp4,Hoxd13,Pkdcc,Sall1,Src,Tbx2,Wnt5a","_rn_":"37"},{"1":"sprouting angiogenesis","2":"9","3":"9.90e-04","4":"Adgra2,Aplnr,Card10,Dll1,E2f2,Epha2,Hdac7,Lef1,Tbxa2r","_rn_":"38"},{"1":"negative regulation of cell population p...","2":"26","3":"1.01e-03","4":"Ackr3,Acvrl1,Bcl11b,Bmp7,Bnipl,Btg2,Cbfa2t3,Cdh5,Csk,Cxcl12,Dab2ip,Dlk1,Dll1,Eng,Etv4,H19,Ifi30,Lbx1,Lmna,Mmp9,Nfe2,Rian,Spint1,Timp2,Wnt5a,Wt1","_rn_":"39"},{"1":"transforming growth factor beta receptor...","2":"11","3":"1.01e-03","4":"Acvrl1,Cdh5,Emilin1,Eng,Htra1,Itga3,Ltbp3,Pmepa1,Smad5,Src,Zfp703","_rn_":"40"},{"1":"male genitalia development","2":"4","3":"1.17e-03","4":"Hoxa13,Hoxd13,Tbx3,Wt1","_rn_":"41"},{"1":"paraxial mesoderm development","2":"4","3":"1.17e-03","4":"Fgfr1,Lef1,Tcf15,Wnt5a","_rn_":"42"},{"1":"dorsal aorta development","2":"3","3":"1.26e-03","4":"Acvrl1,Eng,Hey1","_rn_":"43"},{"1":"negative regulation of endothelial cell ...","2":"3","3":"1.26e-03","4":"Acvrl1,Id1,S1pr3","_rn_":"44"},{"1":"cell fate determination","2":"5","3":"1.27e-03","4":"Ebf2,Hoxa2,Lbx1,Mesp1,Tbx2","_rn_":"45"},{"1":"hair follicle development","2":"8","3":"1.28e-03","4":"Alx4,Celsr1,Hoxc13,Mfsd12,Prss8,Snai1,Tradd,Wnt5a","_rn_":"46"},{"1":"anatomical structure formation involved ...","2":"63","3":"1.46e-03","4":"Ackr3,Acvrl1,Adgra2,Angpt2,Aplnr,Bcl2l11,Bmp7,Card10,Casq2,Cdh5,Cdx1,Cdx2,Celsr1,Col18a1,Cxcl12,Dab2ip,Dll1,Dusp5,E2f2,Egfl7,Emilin1,Eng,Epha2,Ephb2,Fgf18,Fgfr1,Flna,Foxb1,Gbx2,Hdac7,Hey1,Hoxa11,Hoxb1,Hoxc11,Hoxd11,Hspg2,Id1,Kif26b,Krt8,Lef1,Mesp1,Mfng,Mmp9,Mxi1,Ncmap,Nfe2,Pnldc1,Prickle1,Sall1,Sall4,Shb,Smad5,Snai1,Spint1,Tanc1,Tbx2,Tbx3,Tbxa2r,Tcf15,Tnfaip2,Whrn,Wnt5a,Wt1","_rn_":"47"},{"1":"roof of mouth development","2":"8","3":"1.47e-03","4":"Alx4,Ephb2,Lef1,Pkdcc,Snai1,Tbx2,Tbx3,Wnt5a","_rn_":"48"},{"1":"labyrinthine layer development","2":"6","3":"1.53e-03","4":"Bmp7,Cdx2,Hey1,Lef1,Nfe2,Spint1","_rn_":"49"},{"1":"positive regulation of epithelial cell m...","2":"10","3":"1.56e-03","4":"Adgra2,Fgf18,Fgfr1,Hdac7,Itga3,Mmp9,Plpp3,Rab25,Src,Wnt5a","_rn_":"50"},{"1":"negative regulation of myoblast differen...","2":"4","3":"1.67e-03","4":"Dll1,Id3,Prickle1,Tbx3","_rn_":"51"},{"1":"blood vessel development","2":"44","3":"1.69e-03","4":"Ackr3,Acvrl1,Adgra2,Angpt2,Aplnr,Bmp7,Card10,Cdh5,Cdx2,Col18a1,Cxcl12,Dab2ip,Dll1,Dlx3,E2f2,Egfl7,Emilin1,Eng,Epha2,Ephb2,Fgf18,Fgfr1,Flna,Gbx2,Hdac7,Hey1,Hoxa13,Hspg2,Id1,Lef1,Mesp1,Mmp9,Nfe2,Plpp3,Prickle1,Shb,Smad5,Spint1,Tbx2,Tbx3,Tbxa2r,Tnfaip2,Wnt5a,Wt1","_rn_":"52"},{"1":"mesodermal cell fate specification","2":"3","3":"1.70e-03","4":"Fgfr1,Hoxa11,Mesp1","_rn_":"53"},{"1":"mammary gland development","2":"13","3":"1.70e-03","4":"Aprt,Atp2c2,Bcl2l11,Epha2,Etv4,Foxb1,Hoxa9,Lef1,Src,Tbx2,Tbx3,Wnt5a,Zfp703","_rn_":"54"},{"1":"respiratory system development","2":"13","3":"1.80e-03","4":"Celsr1,Fgf18,Fgfr1,Foxp4,Gng8,Id1,Itga3,Lef1,Ltbp3,Meg3,Pkdcc,Wnt5a,Wt1","_rn_":"55"},{"1":"negative regulation of canonical Wnt sig...","2":"9","3":"1.82e-03","4":"Dab2ip,Draxin,Mad2l2,Mesp1,Nkd2,Prickle1,Shisa3,Tmem88,Wnt5a","_rn_":"56"},{"1":"negative regulation of developmental pro...","2":"46","3":"1.88e-03","4":"Acvrl1,Adamts12,Angpt2,Bcl11a,Bmp7,Dab2ip,Dlk1,Dll1,Draxin,E2f2,Emilin1,Epha2,Ephb2,Fgfr1,Foxp4,H19,Hdac7,Hey1,Hoxa2,Hoxa9,Id1,Id3,Lbx1,Lef1,Lmna,Ltbp3,Mad2l2,Meis2,Mesp1,Mmp9,Nfe2,Prdm6,Prickle1,Prtg,S1pr3,Sall1,Sema3b,Sema6c,Shb,Snai1,Spry4,Tbx2,Tbx3,Tnpo2,Wnt5a,Wt1","_rn_":"57"},{"1":"positive regulation of angiogenesis","2":"10","3":"1.92e-03","4":"Acvrl1,Angpt2,Aplnr,Cdh5,Dll1,Eng,Fgf18,Mmp9,Tbxa2r,Wnt5a","_rn_":"58"},{"1":"cell-cell adhesion mediated by cadherin","2":"4","3":"1.97e-03","4":"Cdh5,Mad2l2,Plekha7,Wnt5a","_rn_":"59"},{"1":"regulation of multicellular organismal d...","2":"59","3":"2.07e-03","4":"Acvrl1,Adamts12,Adgra2,Angpt2,Aplnr,Bcl11a,Bmp7,Ccr7,Cdh5,Crabp2,Cxcl12,Dab2ip,Dlk1,Dll1,Draxin,E2f2,Eef2k,Emilin1,Eng,Epha2,Ephb2,Etv4,Fgf18,Fgfr1,Foxp4,Gbx1,Gbx2,Golga4,Hey1,Hoxa11,Hoxa9,Hoxb3,Hoxd11,Id1,Lef1,Lpar3,Ltbp3,Meis2,Mesp1,Mmp9,Ncmap,Nfe2,Nlgn2,Pkdcc,Prtg,Rbm19,S1pr3,Sall1,Sema3b,Sema6c,Shb,Spint1,Tbx2,Tbx3,Tbxa2r,Tradd,Twf2,Wnt5a,Wt1","_rn_":"60"},{"1":"branching morphogenesis of an epithelial...","2":"15","3":"2.12e-03","4":"Bmp7,Celsr1,Cxcl12,Eng,Epha2,Etv4,Gbx2,Hoxa11,Hoxd11,Lef1,Sall1,Src,Tbx3,Wnt5a,Wt1","_rn_":"61"},{"1":"regulation of morphogenesis of a branchi...","2":"6","3":"2.20e-03","4":"Bcl11a,Bmp7,Fgfr1,Hoxd13,Sall1,Wnt5a","_rn_":"62"},{"1":"anterior/posterior axis specification, e...","2":"3","3":"2.22e-03","4":"Tbx3,Wnt5a,Wt1","_rn_":"63"},{"1":"positive regulation of Notch signaling p...","2":"5","3":"2.30e-03","4":"Dll1,Enho,Mesp1,Mfng,Src","_rn_":"64"},{"1":"post-anal tail morphogenesis","2":"4","3":"2.31e-03","4":"Epha2,Hotair,Tcf15,Wnt5a","_rn_":"65"},{"1":"negative regulation of multicellular org...","2":"40","3":"2.64e-03","4":"Acvrl1,Adamts12,Angpt2,Atp1a2,Bcl11a,Bmp7,Card10,Ccr7,Csk,Dab2ip,Dlk1,Dll1,Draxin,Dusp5,E2f2,Emilin1,Epha2,Ephb2,Etv4,Fgfr1,Foxp4,Gnai2,H19,Hdac7,Hey1,Hoxc10,Id1,Lef1,Lmna,Ltbp3,Mad2l2,Mesp1,Nfe2,Prickle1,Prtg,Sema3b,Sema6c,Tbxa2r,Wnt5a,Wt1","_rn_":"66"},{"1":"columnar/cuboidal epithelial cell differ...","2":"7","3":"2.66e-03","4":"Cdx2,Dll1,Dlx3,Foxp4,Lef1,Src,Wnt5a","_rn_":"67"},{"1":"lateral sprouting from an epithelium","2":"3","3":"2.84e-03","4":"Bmp7,Celsr1,Wnt5a","_rn_":"68"},{"1":"cardiac cell fate commitment","2":"3","3":"2.84e-03","4":"Mesp1,Tbx3,Wt1","_rn_":"69"},{"1":"anatomical structure regression","2":"3","3":"2.84e-03","4":"Cd248,Lef1,Smad5","_rn_":"70"},{"1":"muscle cell fate commitment","2":"3","3":"2.84e-03","4":"Tbx2,Tbx3,Wt1","_rn_":"71"},{"1":"negative regulation of pathway-restricte...","2":"3","3":"2.84e-03","4":"Emilin1,Eng,Pmepa1","_rn_":"72"},{"1":"bone remodeling","2":"7","3":"2.85e-03","4":"Csk,Dlk1,Efna2,Epha2,Ltbp3,Src,Syt7","_rn_":"73"},{"1":"negative regulation of cellular response...","2":"7","3":"3.05e-03","4":"Adamts12,Adgra2,Dab2ip,Emilin1,Htra1,Spry4,Wnt5a","_rn_":"74"},{"1":"regulation of BMP signaling pathway","2":"7","3":"3.26e-03","4":"Acvrl1,Cdh5,Eng,Hoxa13,Htra1,Itga3,Wnt5a","_rn_":"75"},{"1":"positive regulation of protein localizat...","2":"6","3":"3.33e-03","4":"Epb41l2,Epha2,Ephb2,Itga3,Nkd2,Zdhhc5","_rn_":"76"},{"1":"negative regulation of protein kinase ac...","2":"11","3":"3.38e-03","4":"Bmp7,Dab2ip,Dusp3,Dusp5,Ephb2,Gprc5a,Pcp4,Sh3bp5,Shb,Smyd3,Spry4","_rn_":"77"},{"1":"regulation of RNA metabolic process","2":"110","3":"3.45e-03","4":"Ablim1,Acvrl1,Ahnak,Alx4,Atf3,Bach2,Bcl11a,Bcl11b,Bmp7,Btaf1,Btg2,Cbfa2t3,Cdx1,Cdx2,Chd3,Crtc1,Dab2ip,Dhx34,Dll1,Dlx3,Dusp5,E2f2,Ebf2,Ehmt2,Eng,Etv4,Fgfr1,Flna,Foxb1,Foxp4,Gbx1,Gbx2,Hdac7,Hey1,Hotair,Hoxa10,Hoxa11,Hoxa11os,Hoxa13,Hoxa2,Hoxa4,Hoxa9,Hoxb1,Hoxb13,Hoxb3,Hoxc10,Hoxc11,Hoxc12,Hoxc13,Hoxc4,Hoxc5,Hoxc6,Hoxc8,Hoxd1,Hoxd10,Hoxd11,Hoxd12,Hoxd13,Hoxd8,Id1,Id3,Igf2bp1,Jund,Klf1,Lbhd1,Lbx1,Lef1,Mad2l2,Meis2,Mesp1,Mex3d,Mxi1,Mybbp1a,Nab2,Nacc1,Nanos3,Nfe2,Nfic,Pknox2,Plagl1,Plpp3,Prdm6,Prickle1,Rgmb,Sall1,Sall3,Sall4,Siva1,Slc38a3,Slc39a5,Smad5,Smyd3,Snai1,Src,Ssbp4,Tbx2,Tbx3,Tcf15,Tox2,Tradd,Trim8,Upf1,Wnt5a,Wt1,Ybx2,Zfhx2,Zfp217,Zfp362,Zfp703,Zfp931","_rn_":"78"},{"1":"positive regulation of epithelial to mes...","2":"5","3":"3.47e-03","4":"Bmp7,Eng,Lef1,Snai1,Zfp703","_rn_":"79"},{"1":"male gonad development","2":"8","3":"3.49e-03","4":"Bcl2l11,Dhx37,Flna,Hoxa10,Hoxa11,Hoxa9,Wnt5a,Wt1","_rn_":"80"},{"1":"negative regulation of cellular macromol...","2":"54","3":"3.53e-03","4":"Abca7,Acvrl1,Atf3,Bach2,Bcl11a,Bmp7,Btaf1,Btg2,Cbfa2t3,Cdx2,Chd3,Dab2ip,Dusp5,Ehmt2,Eng,Fgfr1,Flna,Foxp4,Hdac7,Hey1,Hotair,Hoxa2,Hoxb13,Hoxb3,Hoxc8,Hoxd8,Id1,Id3,Igf2bp1,Lef1,Mad2l2,Mesp1,Mex3d,Mxi1,Mybbp1a,Nab2,Nacc1,Nanos3,Nfic,Prdm6,Prickle1,Sall1,Sall4,Smad5,Snai1,Src,Tbx2,Tbx3,Upf1,Wnt5a,Wt1,Ybx2,Zfp217,Zfp703","_rn_":"81"},{"1":"neuron fate specification","2":"4","3":"3.54e-03","4":"Dll1,Ehmt2,Hoxc10,Hoxd10","_rn_":"82"},{"1":"positive regulation of collateral sprout...","2":"3","3":"3.56e-03","4":"Bcl11a,Crabp2,Lpar3","_rn_":"83"},{"1":"regulation of mesenchymal cell apoptotic...","2":"3","3":"3.56e-03","4":"Bmp7,Hoxa13,Wt1","_rn_":"84"},{"1":"fibroblast activation","2":"3","3":"3.56e-03","4":"Cygb,Fgfr1,Ulk3","_rn_":"85"},{"1":"cell adhesion","2":"48","3":"3.72e-03","4":"Ackr3,Acvrl1,Adamts12,Angpt2,BC024139,Bcl2l11,Bmp7,Ccr7,Cd59b,Cdh15,Cdh5,Celsr1,Cercam,Col18a1,Cxcl12,Dact2,Dll1,Dusp3,Egfl7,Emilin1,Eng,Epha2,Flna,Itga3,Itga9,Kif26b,Lef1,Lgals1,Lrrc4,Mad2l2,Nlgn2,Pcdh1,Plekha7,Plpp3,Podxl,Podxl2,Prtg,Rgmb,Rnd1,Shb,Slc4a1,Spry4,Src,Tgfbi,Tjp3,Vstm2l,Wnt5a,Zfp703","_rn_":"86"},{"1":"negative regulation of reproductive proc...","2":"5","3":"3.83e-03","4":"Bmp7,Shb,Snai1,Wnt5a,Wt1","_rn_":"87"},{"1":"negative regulation of axonogenesis","2":"6","3":"3.89e-03","4":"Bcl11a,Draxin,Ephb2,Sema3b,Sema6c,Wnt5a","_rn_":"88"},{"1":"mammary gland epithelial cell proliferat...","2":"4","3":"4.03e-03","4":"Epha2,Etv4,Wnt5a,Zfp703","_rn_":"89"},{"1":"purine ribonucleotide catabolic process","2":"4","3":"4.03e-03","4":"Nudt18,Nudt3,Nudt4,Pde9a","_rn_":"90"},{"1":"cartilage development","2":"14","3":"4.11e-03","4":"Adamts12,Bmp7,Fgf18,Fgfr1,Hoxa11,Hoxb3,Hoxc4,Hoxd11,Hspg2,Ltbp3,Pkdcc,Smad5,Snai1,Wnt5a","_rn_":"91"},{"1":"regulation of cell development","2":"27","3":"4.15e-03","4":"Adgra2,Bcl11a,Bmp7,Cdh5,Crabp2,Cxcl12,Dll1,Draxin,Eef2k,Ephb2,Flna,Golga4,Hey1,Hoxa11,Hoxb3,Hoxd11,Id1,Lpar3,Prtg,S1pr3,Sema3b,Sema6c,Shb,Spint1,Spry4,Twf2,Wnt5a","_rn_":"92"},{"1":"ureteric bud morphogenesis","2":"6","3":"4.19e-03","4":"Hey1,Hoxa11,Hoxd11,Kif26b,Sall1,Wt1","_rn_":"93"},{"1":"T cell chemotaxis","2":"3","3":"4.37e-03","4":"Ccr7,Cxcl12,Wnt5a","_rn_":"94"},{"1":"regulation of cell fate specification","2":"3","3":"4.37e-03","4":"Fgfr1,Mesp1,Wnt5a","_rn_":"95"},{"1":"locomotory exploration behavior","2":"3","3":"4.37e-03","4":"Atp1a2,Gad1,Nlgn2","_rn_":"96"},{"1":"mesoderm formation","2":"8","3":"4.45e-03","4":"Bmp7,Epha2,Fgfr1,Hoxa11,Lef1,Mesp1,Snai1,Wnt5a","_rn_":"97"},{"1":"negative regulation of chemotaxis","2":"5","3":"4.61e-03","4":"Angpt2,Dusp3,Sema3b,Sema6c,Wnt5a","_rn_":"98"},{"1":"cell-cell junction organization","2":"10","3":"4.63e-03","4":"Cdh5,Csk,Epha2,Flna,Hdac7,Plekha7,Snai1,Tjp3,Whrn,Zfp703","_rn_":"99"},{"1":"positive regulation of cell migration","2":"26","3":"4.70e-03","4":"Ackr3,Adgra2,Bmp7,Ccr7,Cdh5,Col18a1,Cxcl12,Dab2ip,Fgf18,Fgfr1,Flna,Gnai2,Hdac7,Itga3,Lef1,Mmp9,Plpp3,Plvap,Rab25,Sema3b,Sema6c,Snai1,Src,Tradd,Wnt5a,Zfp703","_rn_":"100"}],"options":{"columns":{"min":{},"max":[10]},"rows":{"min":[10],"max":[10]},"pages":{}}}
  </script>
</div>

Interestingly, while genes most active at different fates are all as likely to be associated to at least one peak, there is a strong depletion of genes expressed in thoracic somites among DORC genes.


```r
## check proportion per stage
tmp <- unique(links[,c('Gene','DE_stage')])
df <- rbind(all_DEGs = prop.table(table(degs.stage$fate))*100,
            all_linked_DEGs = prop.table(table(tmp$DE_stage))*100,
            DORC_DEGs  = prop.table(table(tmp[tmp$Gene %in% dorcGenes,]$DE_stage))*100)
df[,c(1,4,2,3)]
```

```
##                 cervical  thoracic   lumbar   sacral
## all_DEGs        34.49304 17.902584 17.58449 30.01988
## all_linked_DEGs 32.70559 15.861059 17.99540 33.43796
## DORC_DEGs       38.05310  8.259587 15.63422 38.05310
```


### Regulatory control of DORC genes

Next, we can ask whether the peaks linked to the highly-regulated genes are enriched for particular TF binding sites, compared to a set of (matched) background peaks. These TFs have the potential of regulating the expression of the gene in question. This hypothesis becomes stronger if the expression of the TF itself is correlated with the accessibility values of the peaks. `FigR` performs these two calculations, and combines the results into a *regulation score*. Positive/negative scores correspond to activating/repressing interactions.



We observe a set of DORCs with very strong TFBS enrichment scores, all corresponding to enrichment for *Nr6a1* binding sites; DORCs regulated by *Nr6a1* show negative scores, indicating TF expression is negatively correlated with chromatin accessibility. Other transcription factors result in more modest motif enrichment scores, with both repressing and enhancing activities.


```r
ggplot(fig.d, aes(Corr.log10P, Enrichment.log10P, color=Score, shape=as.factor(Motif=="Nr6a1"))) +
  geom_point_rast(size=0.75) + 
  scale_color_gradientn(colours = jdb_palette("solar_extra"), 
                        limits=c(-2,2),
                        oob = scales::squish,
                        breaks=scales::breaks_pretty(n=3)) +
  geom_vline(xintercept = 0, lty=2, colour="grey") +
  xlab(expression('(signed) -log'[10]*' p-value TF expression - accessibility correlation')) +
  ylab(expression('-log'[10]*' p-value TF motif enrichment')) +
  labs(shape = "Nr6a1 enrichment") +
  th
```

![](03_peak2gene_files/figure-html/figr_plot-1.png)<!-- -->

Across all DORC genes, some TFs consistently show large regulation scores, as shown below, with *Nr6a1* again standing out. Several of the TFs highlighted have been already picked up in previous analyses, given the enrichment of their binding sites in DA peaks.


```r
rankDrivers(fig.d, rankBy = "meanScore", interactive = FALSE)
```

![](03_peak2gene_files/figure-html/figr_important_tfs-1.png)<!-- -->

Generally, TFs are associated with a handful to a dozen different targets, but some factors show more widespread activity (again, *Nr6a1*, *Ebf1* and *Sox6*, all acting as repressors, and several Hox genes).


```r
rankDrivers(fig.d, score.cut = 1, rankBy = "nTargets", interactive = TRUE)
```

```
## Ranking TFs by total number of associated DORCs ..
```

```
## Using absolute score cut-off of: 1 ..
```

```
## Warning: `gather_()` was deprecated in tidyr 1.2.0.
## Please use `gather()` instead.
## This warning is displayed once every 8 hours.
## Call `lifecycle::last_lifecycle_warnings()` to see where this warning was generated.
```

```{=html}
<div id="htmlwidget-3f6c7e599dbc8e36b9f2" style="width:768px;height:384px;" class="plotly html-widget"></div>
<script type="application/json" data-for="htmlwidget-3f6c7e599dbc8e36b9f2">{"x":{"data":[{"orientation":"v","width":[0.9,0.9,0.9,0.9,0.9,0.9,0.9,0.899999999999999,0.899999999999999,0.899999999999999,0.899999999999999,0.899999999999999,0.899999999999999,0.899999999999999,0.899999999999999,0.899999999999999,0.899999999999999,0.899999999999999,0.899999999999999,0.899999999999999,0.899999999999999,0.899999999999999,0.899999999999999,0.899999999999999,0.899999999999999,0.899999999999999,0.899999999999999,0.899999999999999,0.899999999999999,0.899999999999999,0.899999999999999,0.900000000000002,0.900000000000006,0.900000000000006,0.900000000000006,0.900000000000006,0.900000000000006,0.900000000000006,0.900000000000006,0.900000000000006,0.900000000000006,0.900000000000006,0.900000000000006,0.900000000000006,0.900000000000006,0.900000000000006,0.900000000000006,0.900000000000006,0.900000000000006,0.900000000000006,0.900000000000006,0.900000000000006,0.900000000000006,0.900000000000006,0.900000000000006,0.900000000000006,0.900000000000006,0.900000000000006,0.900000000000006,0.900000000000006,0.900000000000006,0.900000000000006,0.900000000000006,0.900000000000006,0.900000000000006,0.900000000000006,0.900000000000006,0.900000000000006,0.900000000000006,0.900000000000006,0.900000000000006,0.900000000000006,0.900000000000006,0.900000000000006,0.900000000000006,0.900000000000006,0.900000000000006,0.900000000000006,0.900000000000006,0.900000000000006,0.900000000000006,0.900000000000006,0.900000000000006,0.900000000000006,0.900000000000006,0.900000000000006,0.900000000000006,0.900000000000006,0.900000000000006,0.900000000000006,0.900000000000006,0.900000000000006,0.900000000000006,0.900000000000006,0.900000000000006,0.900000000000006,0.900000000000006,0.900000000000006,0.900000000000006,0.900000000000006,0.900000000000006,0.900000000000006,0.900000000000006,0.900000000000006,0.900000000000006,0.900000000000006,0.900000000000006,0.900000000000006,0.900000000000006,0.900000000000006,0.900000000000006,0.900000000000006,0.900000000000006,0.900000000000006,0.900000000000006,0.900000000000006,0.900000000000006,0.900000000000006,0.900000000000006,0.900000000000006,0.900000000000006,0.900000000000006,0.900000000000006,0.900000000000006,0.900000000000006,0.900000000000006,0.900000000000006,0.899999999999991,0.899999999999977,0.899999999999977,0.899999999999977,0.899999999999977,0.899999999999977,0.899999999999977,0.899999999999977,0.899999999999977,0.899999999999977,0.899999999999977,0.899999999999977,0.899999999999977,0.899999999999977,0.899999999999977,0.899999999999977,0.899999999999977,0.899999999999977,0.899999999999977,0.899999999999977,0.899999999999977,0.899999999999977,0.899999999999977,0.899999999999977,0.899999999999977,0.899999999999977,0.899999999999977,0.899999999999977,0.899999999999977,0.899999999999977,0.899999999999977,0.899999999999977,0.899999999999977,0.899999999999977,0.899999999999977,0.899999999999977,0.899999999999977,0.899999999999977,0.899999999999977,0.899999999999977,0.899999999999977,0.899999999999977,0.899999999999977,0.899999999999977,0.899999999999977,0.899999999999977,0.899999999999977,0.899999999999977,0.899999999999977,0.899999999999977,0.899999999999977,0.899999999999977,0.899999999999977,0.899999999999977,0.899999999999977,0.899999999999977],"base":[0,0,2,0,0,2,0,1,2,1,5,0,0,0,0,0,0,0,0,1,0,1,2,0,1,10,0,1,0,0,0,0,0,0,0,0,0,1,2,0,0,0,8,1,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,1,0,-2,-2,0,0,-3,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-6,0,0,0,0,0,0,0,-5,-15,-5,0,0,0,0,0,-12,-8,0,0,0,0,0,0,0,-9,-11,0,0,-9,-10,0,0,-27,-23,0,-13,-16,0,0,-13,0,0,0,0,0,0,0,0,-20,-28,0,0,-33,-52,0,-82,-94,0],"x":[1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63,64,65,66,67,68,69,70,71,72,73,74,75,76,77,78,79,80,81,82,83,84,85,86,87,88,89,90,91,92,93,94,95,96,97,98,99,100,101,102,103,104,105,106,107,108,109,110,111,112,113,114,115,116,117,118,119,120,121,122,123,124,125,126,127,128,129,130,131,132,133,134,135,136,137,138,139,140,141,142,143,144,145,146,147,148,149,150,151,152,153,154,155,156,157,158,159,160,161,162,163,164,165,166,167,168,169,170,171,172,173,174,175,176,177,178,179,180,181,182,183],"y":[53,45,42,26,20,22,18,18,18,14,17,11,11,10,9,9,8,8,7,8,7,8,9,7,7,15,5,6,5,5,5,4,4,4,4,3,3,4,5,3,3,3,11,4,3,2,2,2,2,2,2,2,2,2,3,1,1,1,1,1,1,1,1,2,1,1,1,1,1,1,1,1,1,1,1,1,2,1,1,1,2,1,1,1,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,2,0,0,0,0,0,0,0,1,6,1,0,0,0,0,0,4,2,0,0,0,0,0,0,0,1,2,0,0,1,1,0,0,9,7,0,2,3,0,0,1,0,0,0,0,0,0,0,0,1,5,0,0,2,11,0,2,1,0],"text":["Motif: Plagl1<br />numActivated: 53<br />numRepressed:   0","Motif: Hoxc10<br />numActivated: 45<br />numRepressed:   0","Motif: Hoxa10<br />numActivated: 42<br />numRepressed:   2","Motif: Smad5<br />numActivated: 26<br />numRepressed:   0","Motif: Hoxa11<br />numActivated: 20<br />numRepressed:   0","Motif: Prdm11<br />numActivated: 22<br />numRepressed:   2","Motif: Hoxd11<br />numActivated: 18<br />numRepressed:   0","Motif: Mxi1<br />numActivated: 18<br />numRepressed:   1","Motif: Sall4<br />numActivated: 18<br />numRepressed:   2","Motif: Sox9<br />numActivated: 14<br />numRepressed:   1","Motif: Hoxd10<br />numActivated: 17<br />numRepressed:   5","Motif: Foxp1<br />numActivated: 11<br />numRepressed:   0","Motif: Tcf3<br />numActivated: 11<br />numRepressed:   0","Motif: Rxra<br />numActivated: 10<br />numRepressed:   0","Motif: Gbx2<br />numActivated:  9<br />numRepressed:   0","Motif: Pbx4<br />numActivated:  9<br />numRepressed:   0","Motif: Foxc2<br />numActivated:  8<br />numRepressed:   0","Motif: Snai2<br />numActivated:  8<br />numRepressed:   0","Motif: Gli1<br />numActivated:  7<br />numRepressed:   0","Motif: Hoxa7<br />numActivated:  8<br />numRepressed:   1","Motif: Mef2c<br />numActivated:  7<br />numRepressed:   0","Motif: Nfib<br />numActivated:  8<br />numRepressed:   1","Motif: Pax9<br />numActivated:  9<br />numRepressed:   2","Motif: Zfp287<br />numActivated:  7<br />numRepressed:   0","Motif: Ets2<br />numActivated:  7<br />numRepressed:   1","Motif: Cdx1<br />numActivated: 15<br />numRepressed:  10","Motif: E2f7<br />numActivated:  5<br />numRepressed:   0","Motif: Hoxb8<br />numActivated:  6<br />numRepressed:   1","Motif: Ikzf1<br />numActivated:  5<br />numRepressed:   0","Motif: Mybl2<br />numActivated:  5<br />numRepressed:   0","Motif: Tcf12<br />numActivated:  5<br />numRepressed:   0","Motif: Foxl1<br />numActivated:  4<br />numRepressed:   0","Motif: Hoxd13<br />numActivated:  4<br />numRepressed:   0","Motif: Phf21a<br />numActivated:  4<br />numRepressed:   0","Motif: Twist1<br />numActivated:  4<br />numRepressed:   0","Motif: Creb3l2<br />numActivated:  3<br />numRepressed:   0","Motif: Esrrg<br />numActivated:  3<br />numRepressed:   0","Motif: Foxn3<br />numActivated:  4<br />numRepressed:   1","Motif: Hey1<br />numActivated:  5<br />numRepressed:   2","Motif: Hmga2<br />numActivated:  3<br />numRepressed:   0","Motif: Hoxb7<br />numActivated:  3<br />numRepressed:   0","Motif: Hoxc11<br />numActivated:  3<br />numRepressed:   0","Motif: Sp8<br />numActivated: 11<br />numRepressed:   8","Motif: Srebf1<br />numActivated:  4<br />numRepressed:   1","Motif: Tet1<br />numActivated:  3<br />numRepressed:   0","Motif: Foxb1<br />numActivated:  2<br />numRepressed:   0","Motif: Hes6<br />numActivated:  2<br />numRepressed:   0","Motif: Hoxa2<br />numActivated:  2<br />numRepressed:   0","Motif: Klf1<br />numActivated:  2<br />numRepressed:   0","Motif: Myc<br />numActivated:  2<br />numRepressed:   0","Motif: Rest<br />numActivated:  2<br />numRepressed:   0","Motif: Sp4<br />numActivated:  2<br />numRepressed:   0","Motif: Uncx<br />numActivated:  2<br />numRepressed:   0","Motif: Zic3<br />numActivated:  2<br />numRepressed:   0","Motif: Zic5<br />numActivated:  3<br />numRepressed:   1","Motif: Cdx4<br />numActivated:  1<br />numRepressed:   0","Motif: Cxxc4<br />numActivated:  1<br />numRepressed:   0","Motif: Foxc1<br />numActivated:  1<br />numRepressed:   0","Motif: Foxk1<br />numActivated:  1<br />numRepressed:   0","Motif: Gfi1b<br />numActivated:  1<br />numRepressed:   0","Motif: Hesx1<br />numActivated:  1<br />numRepressed:   0","Motif: Hoxa6<br />numActivated:  1<br />numRepressed:   0","Motif: Hoxb2<br />numActivated:  1<br />numRepressed:   0","Motif: Hoxb6<br />numActivated:  2<br />numRepressed:   1","Motif: Hoxd12<br />numActivated:  1<br />numRepressed:   0","Motif: Hoxd4<br />numActivated:  1<br />numRepressed:   0","Motif: Hoxd8<br />numActivated:  1<br />numRepressed:   0","Motif: Irf3<br />numActivated:  1<br />numRepressed:   0","Motif: Jun<br />numActivated:  1<br />numRepressed:   0","Motif: Klf6<br />numActivated:  1<br />numRepressed:   0","Motif: Maf<br />numActivated:  1<br />numRepressed:   0","Motif: Nr2f1<br />numActivated:  1<br />numRepressed:   0","Motif: Osr1<br />numActivated:  1<br />numRepressed:   0","Motif: Pax1<br />numActivated:  1<br />numRepressed:   0","Motif: Pbx3<br />numActivated:  1<br />numRepressed:   0","Motif: Sox12<br />numActivated:  1<br />numRepressed:   0","Motif: Sp5<br />numActivated:  2<br />numRepressed:   1","Motif: Stat3<br />numActivated:  1<br />numRepressed:   0","Motif: Thra<br />numActivated:  1<br />numRepressed:   0","Motif: Zfp13<br />numActivated:  1<br />numRepressed:   0","Motif: Zfp692<br />numActivated:  2<br />numRepressed:   1","Motif: Zfp711<br />numActivated:  1<br />numRepressed:   0","Motif: Cux1<br />numActivated:  1<br />numRepressed:   1","Motif: Ybx3<br />numActivated:  1<br />numRepressed:   1","Motif: Arid5b<br />numActivated:  0<br />numRepressed:   1","Motif: Arnt2<br />numActivated:  0<br />numRepressed:   1","Motif: Dnajc21<br />numActivated:  1<br />numRepressed:   2","Motif: Dpf3<br />numActivated:  0<br />numRepressed:   1","Motif: Esrra<br />numActivated:  0<br />numRepressed:   1","Motif: Foxf1<br />numActivated:  0<br />numRepressed:   1","Motif: Foxo1<br />numActivated:  0<br />numRepressed:   1","Motif: Foxp2<br />numActivated:  0<br />numRepressed:   1","Motif: Hbp1<br />numActivated:  0<br />numRepressed:   1","Motif: Heyl<br />numActivated:  0<br />numRepressed:   1","Motif: Hic2<br />numActivated:  0<br />numRepressed:   1","Motif: Hoxa1<br />numActivated:  0<br />numRepressed:   1","Motif: Hoxa4<br />numActivated:  0<br />numRepressed:   1","Motif: Hoxd9<br />numActivated:  0<br />numRepressed:   1","Motif: Hsf2<br />numActivated:  0<br />numRepressed:   1","Motif: Jund<br />numActivated:  0<br />numRepressed:   1","Motif: Lhx1<br />numActivated:  0<br />numRepressed:   1","Motif: Meis1<br />numActivated:  0<br />numRepressed:   1","Motif: Pax8<br />numActivated:  0<br />numRepressed:   1","Motif: Rfx7<br />numActivated:  0<br />numRepressed:   1","Motif: Setbp1<br />numActivated:  0<br />numRepressed:   1","Motif: Six4<br />numActivated:  0<br />numRepressed:   1","Motif: Snai1<br />numActivated:  0<br />numRepressed:   1","Motif: Srf<br />numActivated:  0<br />numRepressed:   1","Motif: Tbx18<br />numActivated:  0<br />numRepressed:   1","Motif: Tead2<br />numActivated:  0<br />numRepressed:   1","Motif: Tigd2<br />numActivated:  0<br />numRepressed:   1","Motif: Twist2<br />numActivated:  0<br />numRepressed:   1","Motif: Zbtb1<br />numActivated:  0<br />numRepressed:   1","Motif: Zbtb18<br />numActivated:  0<br />numRepressed:   1","Motif: Zeb1<br />numActivated:  0<br />numRepressed:   1","Motif: Zfp46<br />numActivated:  0<br />numRepressed:   1","Motif: Zfp558<br />numActivated:  0<br />numRepressed:   1","Motif: Arid3b<br />numActivated:  0<br />numRepressed:   2","Motif: Glis1<br />numActivated:  0<br />numRepressed:   2","Motif: Rarb<br />numActivated:  0<br />numRepressed:   2","Motif: Rora<br />numActivated:  0<br />numRepressed:   2","Motif: Six2<br />numActivated:  0<br />numRepressed:   2","Motif: Sox13<br />numActivated:  0<br />numRepressed:   2","Motif: T<br />numActivated:  0<br />numRepressed:   2","Motif: Tbx6<br />numActivated:  2<br />numRepressed:   4","Motif: Tcf4<br />numActivated:  0<br />numRepressed:   2","Motif: Tead1<br />numActivated:  0<br />numRepressed:   2","Motif: Zfp317<br />numActivated:  0<br />numRepressed:   2","Motif: Zfp454<br />numActivated:  0<br />numRepressed:   2","Motif: Zfp770<br />numActivated:  0<br />numRepressed:   2","Motif: Zscan20<br />numActivated:  0<br />numRepressed:   2","Motif: Arid3a<br />numActivated:  0<br />numRepressed:   3","Motif: Ddit3<br />numActivated:  1<br />numRepressed:   4","Motif: Hoxb9<br />numActivated:  6<br />numRepressed:   9","Motif: Hoxc9<br />numActivated:  1<br />numRepressed:   4","Motif: Nr2c2<br />numActivated:  0<br />numRepressed:   3","Motif: Tcf7<br />numActivated:  0<br />numRepressed:   3","Motif: Zfp189<br />numActivated:  0<br />numRepressed:   3","Motif: Zfp422<br />numActivated:  0<br />numRepressed:   3","Motif: Zfp580<br />numActivated:  0<br />numRepressed:   3","Motif: Ebf2<br />numActivated:  4<br />numRepressed:   8","Motif: Ebf3<br />numActivated:  2<br />numRepressed:   6","Motif: Foxp4<br />numActivated:  0<br />numRepressed:   4","Motif: Hoxb4<br />numActivated:  0<br />numRepressed:   4","Motif: Pknox2<br />numActivated:  0<br />numRepressed:   4","Motif: Zfp248<br />numActivated:  0<br />numRepressed:   4","Motif: Zfp263<br />numActivated:  0<br />numRepressed:   4","Motif: Srebf2<br />numActivated:  0<br />numRepressed:   5","Motif: Meis3<br />numActivated:  0<br />numRepressed:   6","Motif: Hand1<br />numActivated:  1<br />numRepressed:   8","Motif: Hoxa5<br />numActivated:  2<br />numRepressed:   9","Motif: Prdm6<br />numActivated:  0<br />numRepressed:   7","Motif: Tbx22<br />numActivated:  0<br />numRepressed:   7","Motif: Trps1<br />numActivated:  1<br />numRepressed:   8","Motif: Gli2<br />numActivated:  1<br />numRepressed:   9","Motif: Rara<br />numActivated:  0<br />numRepressed:   8","Motif: Thap11<br />numActivated:  0<br />numRepressed:   8","Motif: Hmga1<br />numActivated:  9<br />numRepressed:  18","Motif: Hoxc8<br />numActivated:  7<br />numRepressed:  16","Motif: Zfp143<br />numActivated:  0<br />numRepressed:   9","Motif: Zfp932<br />numActivated:  2<br />numRepressed:  11","Motif: Hoxc4<br />numActivated:  3<br />numRepressed:  13","Motif: Prdm5<br />numActivated:  0<br />numRepressed:  10","Motif: Zfp637<br />numActivated:  0<br />numRepressed:  10","Motif: Hoxb5<br />numActivated:  1<br />numRepressed:  12","Motif: Mycn<br />numActivated:  0<br />numRepressed:  11","Motif: Tbx2<br />numActivated:  0<br />numRepressed:  11","Motif: Zeb2<br />numActivated:  0<br />numRepressed:  11","Motif: Wt1<br />numActivated:  0<br />numRepressed:  13","Motif: Hoxc5<br />numActivated:  0<br />numRepressed:  14","Motif: Klf3<br />numActivated:  0<br />numRepressed:  14","Motif: Tcf7l2<br />numActivated:  0<br />numRepressed:  14","Motif: Klf10<br />numActivated:  0<br />numRepressed:  17","Motif: Glis2<br />numActivated:  1<br />numRepressed:  19","Motif: Sox11<br />numActivated:  5<br />numRepressed:  23","Motif: Lef1<br />numActivated:  0<br />numRepressed:  21","Motif: Nr1d2<br />numActivated:  0<br />numRepressed:  24","Motif: Rreb1<br />numActivated:  2<br />numRepressed:  31","Motif: Hoxa9<br />numActivated: 11<br />numRepressed:  41","Motif: Zfp647<br />numActivated:  0<br />numRepressed:  43","Motif: Ebf1<br />numActivated:  2<br />numRepressed:  80","Motif: Sox6<br />numActivated:  1<br />numRepressed:  93","Motif: Nr6a1<br />numActivated:  0<br />numRepressed: 100"],"type":"bar","textposition":"none","marker":{"autocolorscale":false,"color":"rgba(205,38,38,1)","line":{"width":0.377952755905512,"color":"rgba(211,211,211,1)"}},"name":"numActivatedY","legendgroup":"numActivatedY","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"orientation":"v","width":[0.9,0.9,0.9,0.9,0.9,0.9,0.9,0.899999999999999,0.899999999999999,0.899999999999999,0.899999999999999,0.899999999999999,0.899999999999999,0.899999999999999,0.899999999999999,0.899999999999999,0.899999999999999,0.899999999999999,0.899999999999999,0.899999999999999,0.899999999999999,0.899999999999999,0.899999999999999,0.899999999999999,0.899999999999999,0.899999999999999,0.899999999999999,0.899999999999999,0.899999999999999,0.899999999999999,0.899999999999999,0.900000000000002,0.900000000000006,0.900000000000006,0.900000000000006,0.900000000000006,0.900000000000006,0.900000000000006,0.900000000000006,0.900000000000006,0.900000000000006,0.900000000000006,0.900000000000006,0.900000000000006,0.900000000000006,0.900000000000006,0.900000000000006,0.900000000000006,0.900000000000006,0.900000000000006,0.900000000000006,0.900000000000006,0.900000000000006,0.900000000000006,0.900000000000006,0.900000000000006,0.900000000000006,0.900000000000006,0.900000000000006,0.900000000000006,0.900000000000006,0.900000000000006,0.900000000000006,0.900000000000006,0.900000000000006,0.900000000000006,0.900000000000006,0.900000000000006,0.900000000000006,0.900000000000006,0.900000000000006,0.900000000000006,0.900000000000006,0.900000000000006,0.900000000000006,0.900000000000006,0.900000000000006,0.900000000000006,0.900000000000006,0.900000000000006,0.900000000000006,0.900000000000006,0.900000000000006,0.900000000000006,0.900000000000006,0.900000000000006,0.900000000000006,0.900000000000006,0.900000000000006,0.900000000000006,0.900000000000006,0.900000000000006,0.900000000000006,0.900000000000006,0.900000000000006,0.900000000000006,0.900000000000006,0.900000000000006,0.900000000000006,0.900000000000006,0.900000000000006,0.900000000000006,0.900000000000006,0.900000000000006,0.900000000000006,0.900000000000006,0.900000000000006,0.900000000000006,0.900000000000006,0.900000000000006,0.900000000000006,0.900000000000006,0.900000000000006,0.900000000000006,0.900000000000006,0.900000000000006,0.900000000000006,0.900000000000006,0.900000000000006,0.900000000000006,0.900000000000006,0.900000000000006,0.900000000000006,0.900000000000006,0.900000000000006,0.900000000000006,0.900000000000006,0.899999999999991,0.899999999999977,0.899999999999977,0.899999999999977,0.899999999999977,0.899999999999977,0.899999999999977,0.899999999999977,0.899999999999977,0.899999999999977,0.899999999999977,0.899999999999977,0.899999999999977,0.899999999999977,0.899999999999977,0.899999999999977,0.899999999999977,0.899999999999977,0.899999999999977,0.899999999999977,0.899999999999977,0.899999999999977,0.899999999999977,0.899999999999977,0.899999999999977,0.899999999999977,0.899999999999977,0.899999999999977,0.899999999999977,0.899999999999977,0.899999999999977,0.899999999999977,0.899999999999977,0.899999999999977,0.899999999999977,0.899999999999977,0.899999999999977,0.899999999999977,0.899999999999977,0.899999999999977,0.899999999999977,0.899999999999977,0.899999999999977,0.899999999999977,0.899999999999977,0.899999999999977,0.899999999999977,0.899999999999977,0.899999999999977,0.899999999999977,0.899999999999977,0.899999999999977,0.899999999999977,0.899999999999977,0.899999999999977,0.899999999999977],"base":[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-1,-1,-1,-1,-2,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-2,-2,-2,-2,-2,-2,-2,-4,-2,-2,-2,-2,-2,-2,-3,-4,-9,-4,-3,-3,-3,-3,-3,-8,-6,-4,-4,-4,-4,-4,-5,-6,-8,-9,-7,-7,-8,-9,-8,-8,-18,-16,-9,-11,-13,-10,-10,-12,-11,-11,-11,-13,-14,-14,-14,-17,-19,-23,-21,-24,-31,-41,-43,-80,-93,-100],"x":[1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63,64,65,66,67,68,69,70,71,72,73,74,75,76,77,78,79,80,81,82,83,84,85,86,87,88,89,90,91,92,93,94,95,96,97,98,99,100,101,102,103,104,105,106,107,108,109,110,111,112,113,114,115,116,117,118,119,120,121,122,123,124,125,126,127,128,129,130,131,132,133,134,135,136,137,138,139,140,141,142,143,144,145,146,147,148,149,150,151,152,153,154,155,156,157,158,159,160,161,162,163,164,165,166,167,168,169,170,171,172,173,174,175,176,177,178,179,180,181,182,183],"y":[0,0,2,0,0,2,0,1,2,1,5,0,0,0,0,0,0,0,0,1,0,1,2,0,1,10,0,1,0,0,0,0,0,0,0,0,0,1,2,0,0,0,8,1,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,1,0,1,1,1,1,2,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,4,2,2,2,2,2,2,3,4,9,4,3,3,3,3,3,8,6,4,4,4,4,4,5,6,8,9,7,7,8,9,8,8,18,16,9,11,13,10,10,12,11,11,11,13,14,14,14,17,19,23,21,24,31,41,43,80,93,100],"text":["Motif: Plagl1<br />numActivated: 53<br />numRepressed:   0","Motif: Hoxc10<br />numActivated: 45<br />numRepressed:   0","Motif: Hoxa10<br />numActivated: 42<br />numRepressed:   2","Motif: Smad5<br />numActivated: 26<br />numRepressed:   0","Motif: Hoxa11<br />numActivated: 20<br />numRepressed:   0","Motif: Prdm11<br />numActivated: 22<br />numRepressed:   2","Motif: Hoxd11<br />numActivated: 18<br />numRepressed:   0","Motif: Mxi1<br />numActivated: 18<br />numRepressed:   1","Motif: Sall4<br />numActivated: 18<br />numRepressed:   2","Motif: Sox9<br />numActivated: 14<br />numRepressed:   1","Motif: Hoxd10<br />numActivated: 17<br />numRepressed:   5","Motif: Foxp1<br />numActivated: 11<br />numRepressed:   0","Motif: Tcf3<br />numActivated: 11<br />numRepressed:   0","Motif: Rxra<br />numActivated: 10<br />numRepressed:   0","Motif: Gbx2<br />numActivated:  9<br />numRepressed:   0","Motif: Pbx4<br />numActivated:  9<br />numRepressed:   0","Motif: Foxc2<br />numActivated:  8<br />numRepressed:   0","Motif: Snai2<br />numActivated:  8<br />numRepressed:   0","Motif: Gli1<br />numActivated:  7<br />numRepressed:   0","Motif: Hoxa7<br />numActivated:  8<br />numRepressed:   1","Motif: Mef2c<br />numActivated:  7<br />numRepressed:   0","Motif: Nfib<br />numActivated:  8<br />numRepressed:   1","Motif: Pax9<br />numActivated:  9<br />numRepressed:   2","Motif: Zfp287<br />numActivated:  7<br />numRepressed:   0","Motif: Ets2<br />numActivated:  7<br />numRepressed:   1","Motif: Cdx1<br />numActivated: 15<br />numRepressed:  10","Motif: E2f7<br />numActivated:  5<br />numRepressed:   0","Motif: Hoxb8<br />numActivated:  6<br />numRepressed:   1","Motif: Ikzf1<br />numActivated:  5<br />numRepressed:   0","Motif: Mybl2<br />numActivated:  5<br />numRepressed:   0","Motif: Tcf12<br />numActivated:  5<br />numRepressed:   0","Motif: Foxl1<br />numActivated:  4<br />numRepressed:   0","Motif: Hoxd13<br />numActivated:  4<br />numRepressed:   0","Motif: Phf21a<br />numActivated:  4<br />numRepressed:   0","Motif: Twist1<br />numActivated:  4<br />numRepressed:   0","Motif: Creb3l2<br />numActivated:  3<br />numRepressed:   0","Motif: Esrrg<br />numActivated:  3<br />numRepressed:   0","Motif: Foxn3<br />numActivated:  4<br />numRepressed:   1","Motif: Hey1<br />numActivated:  5<br />numRepressed:   2","Motif: Hmga2<br />numActivated:  3<br />numRepressed:   0","Motif: Hoxb7<br />numActivated:  3<br />numRepressed:   0","Motif: Hoxc11<br />numActivated:  3<br />numRepressed:   0","Motif: Sp8<br />numActivated: 11<br />numRepressed:   8","Motif: Srebf1<br />numActivated:  4<br />numRepressed:   1","Motif: Tet1<br />numActivated:  3<br />numRepressed:   0","Motif: Foxb1<br />numActivated:  2<br />numRepressed:   0","Motif: Hes6<br />numActivated:  2<br />numRepressed:   0","Motif: Hoxa2<br />numActivated:  2<br />numRepressed:   0","Motif: Klf1<br />numActivated:  2<br />numRepressed:   0","Motif: Myc<br />numActivated:  2<br />numRepressed:   0","Motif: Rest<br />numActivated:  2<br />numRepressed:   0","Motif: Sp4<br />numActivated:  2<br />numRepressed:   0","Motif: Uncx<br />numActivated:  2<br />numRepressed:   0","Motif: Zic3<br />numActivated:  2<br />numRepressed:   0","Motif: Zic5<br />numActivated:  3<br />numRepressed:   1","Motif: Cdx4<br />numActivated:  1<br />numRepressed:   0","Motif: Cxxc4<br />numActivated:  1<br />numRepressed:   0","Motif: Foxc1<br />numActivated:  1<br />numRepressed:   0","Motif: Foxk1<br />numActivated:  1<br />numRepressed:   0","Motif: Gfi1b<br />numActivated:  1<br />numRepressed:   0","Motif: Hesx1<br />numActivated:  1<br />numRepressed:   0","Motif: Hoxa6<br />numActivated:  1<br />numRepressed:   0","Motif: Hoxb2<br />numActivated:  1<br />numRepressed:   0","Motif: Hoxb6<br />numActivated:  2<br />numRepressed:   1","Motif: Hoxd12<br />numActivated:  1<br />numRepressed:   0","Motif: Hoxd4<br />numActivated:  1<br />numRepressed:   0","Motif: Hoxd8<br />numActivated:  1<br />numRepressed:   0","Motif: Irf3<br />numActivated:  1<br />numRepressed:   0","Motif: Jun<br />numActivated:  1<br />numRepressed:   0","Motif: Klf6<br />numActivated:  1<br />numRepressed:   0","Motif: Maf<br />numActivated:  1<br />numRepressed:   0","Motif: Nr2f1<br />numActivated:  1<br />numRepressed:   0","Motif: Osr1<br />numActivated:  1<br />numRepressed:   0","Motif: Pax1<br />numActivated:  1<br />numRepressed:   0","Motif: Pbx3<br />numActivated:  1<br />numRepressed:   0","Motif: Sox12<br />numActivated:  1<br />numRepressed:   0","Motif: Sp5<br />numActivated:  2<br />numRepressed:   1","Motif: Stat3<br />numActivated:  1<br />numRepressed:   0","Motif: Thra<br />numActivated:  1<br />numRepressed:   0","Motif: Zfp13<br />numActivated:  1<br />numRepressed:   0","Motif: Zfp692<br />numActivated:  2<br />numRepressed:   1","Motif: Zfp711<br />numActivated:  1<br />numRepressed:   0","Motif: Cux1<br />numActivated:  1<br />numRepressed:   1","Motif: Ybx3<br />numActivated:  1<br />numRepressed:   1","Motif: Arid5b<br />numActivated:  0<br />numRepressed:   1","Motif: Arnt2<br />numActivated:  0<br />numRepressed:   1","Motif: Dnajc21<br />numActivated:  1<br />numRepressed:   2","Motif: Dpf3<br />numActivated:  0<br />numRepressed:   1","Motif: Esrra<br />numActivated:  0<br />numRepressed:   1","Motif: Foxf1<br />numActivated:  0<br />numRepressed:   1","Motif: Foxo1<br />numActivated:  0<br />numRepressed:   1","Motif: Foxp2<br />numActivated:  0<br />numRepressed:   1","Motif: Hbp1<br />numActivated:  0<br />numRepressed:   1","Motif: Heyl<br />numActivated:  0<br />numRepressed:   1","Motif: Hic2<br />numActivated:  0<br />numRepressed:   1","Motif: Hoxa1<br />numActivated:  0<br />numRepressed:   1","Motif: Hoxa4<br />numActivated:  0<br />numRepressed:   1","Motif: Hoxd9<br />numActivated:  0<br />numRepressed:   1","Motif: Hsf2<br />numActivated:  0<br />numRepressed:   1","Motif: Jund<br />numActivated:  0<br />numRepressed:   1","Motif: Lhx1<br />numActivated:  0<br />numRepressed:   1","Motif: Meis1<br />numActivated:  0<br />numRepressed:   1","Motif: Pax8<br />numActivated:  0<br />numRepressed:   1","Motif: Rfx7<br />numActivated:  0<br />numRepressed:   1","Motif: Setbp1<br />numActivated:  0<br />numRepressed:   1","Motif: Six4<br />numActivated:  0<br />numRepressed:   1","Motif: Snai1<br />numActivated:  0<br />numRepressed:   1","Motif: Srf<br />numActivated:  0<br />numRepressed:   1","Motif: Tbx18<br />numActivated:  0<br />numRepressed:   1","Motif: Tead2<br />numActivated:  0<br />numRepressed:   1","Motif: Tigd2<br />numActivated:  0<br />numRepressed:   1","Motif: Twist2<br />numActivated:  0<br />numRepressed:   1","Motif: Zbtb1<br />numActivated:  0<br />numRepressed:   1","Motif: Zbtb18<br />numActivated:  0<br />numRepressed:   1","Motif: Zeb1<br />numActivated:  0<br />numRepressed:   1","Motif: Zfp46<br />numActivated:  0<br />numRepressed:   1","Motif: Zfp558<br />numActivated:  0<br />numRepressed:   1","Motif: Arid3b<br />numActivated:  0<br />numRepressed:   2","Motif: Glis1<br />numActivated:  0<br />numRepressed:   2","Motif: Rarb<br />numActivated:  0<br />numRepressed:   2","Motif: Rora<br />numActivated:  0<br />numRepressed:   2","Motif: Six2<br />numActivated:  0<br />numRepressed:   2","Motif: Sox13<br />numActivated:  0<br />numRepressed:   2","Motif: T<br />numActivated:  0<br />numRepressed:   2","Motif: Tbx6<br />numActivated:  2<br />numRepressed:   4","Motif: Tcf4<br />numActivated:  0<br />numRepressed:   2","Motif: Tead1<br />numActivated:  0<br />numRepressed:   2","Motif: Zfp317<br />numActivated:  0<br />numRepressed:   2","Motif: Zfp454<br />numActivated:  0<br />numRepressed:   2","Motif: Zfp770<br />numActivated:  0<br />numRepressed:   2","Motif: Zscan20<br />numActivated:  0<br />numRepressed:   2","Motif: Arid3a<br />numActivated:  0<br />numRepressed:   3","Motif: Ddit3<br />numActivated:  1<br />numRepressed:   4","Motif: Hoxb9<br />numActivated:  6<br />numRepressed:   9","Motif: Hoxc9<br />numActivated:  1<br />numRepressed:   4","Motif: Nr2c2<br />numActivated:  0<br />numRepressed:   3","Motif: Tcf7<br />numActivated:  0<br />numRepressed:   3","Motif: Zfp189<br />numActivated:  0<br />numRepressed:   3","Motif: Zfp422<br />numActivated:  0<br />numRepressed:   3","Motif: Zfp580<br />numActivated:  0<br />numRepressed:   3","Motif: Ebf2<br />numActivated:  4<br />numRepressed:   8","Motif: Ebf3<br />numActivated:  2<br />numRepressed:   6","Motif: Foxp4<br />numActivated:  0<br />numRepressed:   4","Motif: Hoxb4<br />numActivated:  0<br />numRepressed:   4","Motif: Pknox2<br />numActivated:  0<br />numRepressed:   4","Motif: Zfp248<br />numActivated:  0<br />numRepressed:   4","Motif: Zfp263<br />numActivated:  0<br />numRepressed:   4","Motif: Srebf2<br />numActivated:  0<br />numRepressed:   5","Motif: Meis3<br />numActivated:  0<br />numRepressed:   6","Motif: Hand1<br />numActivated:  1<br />numRepressed:   8","Motif: Hoxa5<br />numActivated:  2<br />numRepressed:   9","Motif: Prdm6<br />numActivated:  0<br />numRepressed:   7","Motif: Tbx22<br />numActivated:  0<br />numRepressed:   7","Motif: Trps1<br />numActivated:  1<br />numRepressed:   8","Motif: Gli2<br />numActivated:  1<br />numRepressed:   9","Motif: Rara<br />numActivated:  0<br />numRepressed:   8","Motif: Thap11<br />numActivated:  0<br />numRepressed:   8","Motif: Hmga1<br />numActivated:  9<br />numRepressed:  18","Motif: Hoxc8<br />numActivated:  7<br />numRepressed:  16","Motif: Zfp143<br />numActivated:  0<br />numRepressed:   9","Motif: Zfp932<br />numActivated:  2<br />numRepressed:  11","Motif: Hoxc4<br />numActivated:  3<br />numRepressed:  13","Motif: Prdm5<br />numActivated:  0<br />numRepressed:  10","Motif: Zfp637<br />numActivated:  0<br />numRepressed:  10","Motif: Hoxb5<br />numActivated:  1<br />numRepressed:  12","Motif: Mycn<br />numActivated:  0<br />numRepressed:  11","Motif: Tbx2<br />numActivated:  0<br />numRepressed:  11","Motif: Zeb2<br />numActivated:  0<br />numRepressed:  11","Motif: Wt1<br />numActivated:  0<br />numRepressed:  13","Motif: Hoxc5<br />numActivated:  0<br />numRepressed:  14","Motif: Klf3<br />numActivated:  0<br />numRepressed:  14","Motif: Tcf7l2<br />numActivated:  0<br />numRepressed:  14","Motif: Klf10<br />numActivated:  0<br />numRepressed:  17","Motif: Glis2<br />numActivated:  1<br />numRepressed:  19","Motif: Sox11<br />numActivated:  5<br />numRepressed:  23","Motif: Lef1<br />numActivated:  0<br />numRepressed:  21","Motif: Nr1d2<br />numActivated:  0<br />numRepressed:  24","Motif: Rreb1<br />numActivated:  2<br />numRepressed:  31","Motif: Hoxa9<br />numActivated: 11<br />numRepressed:  41","Motif: Zfp647<br />numActivated:  0<br />numRepressed:  43","Motif: Ebf1<br />numActivated:  2<br />numRepressed:  80","Motif: Sox6<br />numActivated:  1<br />numRepressed:  93","Motif: Nr6a1<br />numActivated:  0<br />numRepressed: 100"],"type":"bar","textposition":"none","marker":{"autocolorscale":false,"color":"rgba(54,100,139,1)","line":{"width":0.377952755905512,"color":"rgba(211,211,211,1)"}},"name":"numRepressedY","legendgroup":"numRepressedY","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"x":[0.4,183.6],"y":[0,0],"text":"","type":"scatter","mode":"lines","line":{"width":1.88976377952756,"color":"rgba(0,0,0,1)","dash":"solid"},"hoveron":"points","showlegend":false,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null}],"layout":{"margin":{"t":30.6118721461187,"r":7.30593607305936,"b":32.8767123287671,"l":43.1050228310502},"plot_bgcolor":"rgba(255,255,255,1)","paper_bgcolor":"rgba(255,255,255,1)","font":{"color":"rgba(0,0,0,1)","family":"","size":14.6118721461187},"xaxis":{"domain":[0,1],"automargin":true,"type":"linear","autorange":false,"range":[0.4,183.6],"tickmode":"array","ticktext":["Plagl1","Hoxc10","Hoxa10","Smad5","Hoxa11","Prdm11","Hoxd11","Mxi1","Sall4","Sox9","Hoxd10","Foxp1","Tcf3","Rxra","Gbx2","Pbx4","Foxc2","Snai2","Gli1","Hoxa7","Mef2c","Nfib","Pax9","Zfp287","Ets2","Cdx1","E2f7","Hoxb8","Ikzf1","Mybl2","Tcf12","Foxl1","Hoxd13","Phf21a","Twist1","Creb3l2","Esrrg","Foxn3","Hey1","Hmga2","Hoxb7","Hoxc11","Sp8","Srebf1","Tet1","Foxb1","Hes6","Hoxa2","Klf1","Myc","Rest","Sp4","Uncx","Zic3","Zic5","Cdx4","Cxxc4","Foxc1","Foxk1","Gfi1b","Hesx1","Hoxa6","Hoxb2","Hoxb6","Hoxd12","Hoxd4","Hoxd8","Irf3","Jun","Klf6","Maf","Nr2f1","Osr1","Pax1","Pbx3","Sox12","Sp5","Stat3","Thra","Zfp13","Zfp692","Zfp711","Cux1","Ybx3","Arid5b","Arnt2","Dnajc21","Dpf3","Esrra","Foxf1","Foxo1","Foxp2","Hbp1","Heyl","Hic2","Hoxa1","Hoxa4","Hoxd9","Hsf2","Jund","Lhx1","Meis1","Pax8","Rfx7","Setbp1","Six4","Snai1","Srf","Tbx18","Tead2","Tigd2","Twist2","Zbtb1","Zbtb18","Zeb1","Zfp46","Zfp558","Arid3b","Glis1","Rarb","Rora","Six2","Sox13","T","Tbx6","Tcf4","Tead1","Zfp317","Zfp454","Zfp770","Zscan20","Arid3a","Ddit3","Hoxb9","Hoxc9","Nr2c2","Tcf7","Zfp189","Zfp422","Zfp580","Ebf2","Ebf3","Foxp4","Hoxb4","Pknox2","Zfp248","Zfp263","Srebf2","Meis3","Hand1","Hoxa5","Prdm6","Tbx22","Trps1","Gli2","Rara","Thap11","Hmga1","Hoxc8","Zfp143","Zfp932","Hoxc4","Prdm5","Zfp637","Hoxb5","Mycn","Tbx2","Zeb2","Wt1","Hoxc5","Klf3","Tcf7l2","Klf10","Glis2","Sox11","Lef1","Nr1d2","Rreb1","Hoxa9","Zfp647","Ebf1","Sox6","Nr6a1"],"tickvals":[1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63,64,65,66,67,68,69,70,71,72,73,74,75,76,77,78,79,80,81,82,83,84,85,86,87,88,89,90,91,92,93,94,95,96,97,98,99,100,101,102,103,104,105,106,107,108,109,110,111,112,113,114,115,116,117,118,119,120,121,122,123,124,125,126,127,128,129,130,131,132,133,134,135,136,137,138,139,140,141,142,143,144,145,146,147,148,149,150,151,152,153,154,155,156,157,158,159,160,161,162,163,164,165,166,167,168,169,170,171,172,173,174,175,176,177,178,179,180,181,182,183],"categoryorder":"array","categoryarray":["Plagl1","Hoxc10","Hoxa10","Smad5","Hoxa11","Prdm11","Hoxd11","Mxi1","Sall4","Sox9","Hoxd10","Foxp1","Tcf3","Rxra","Gbx2","Pbx4","Foxc2","Snai2","Gli1","Hoxa7","Mef2c","Nfib","Pax9","Zfp287","Ets2","Cdx1","E2f7","Hoxb8","Ikzf1","Mybl2","Tcf12","Foxl1","Hoxd13","Phf21a","Twist1","Creb3l2","Esrrg","Foxn3","Hey1","Hmga2","Hoxb7","Hoxc11","Sp8","Srebf1","Tet1","Foxb1","Hes6","Hoxa2","Klf1","Myc","Rest","Sp4","Uncx","Zic3","Zic5","Cdx4","Cxxc4","Foxc1","Foxk1","Gfi1b","Hesx1","Hoxa6","Hoxb2","Hoxb6","Hoxd12","Hoxd4","Hoxd8","Irf3","Jun","Klf6","Maf","Nr2f1","Osr1","Pax1","Pbx3","Sox12","Sp5","Stat3","Thra","Zfp13","Zfp692","Zfp711","Cux1","Ybx3","Arid5b","Arnt2","Dnajc21","Dpf3","Esrra","Foxf1","Foxo1","Foxp2","Hbp1","Heyl","Hic2","Hoxa1","Hoxa4","Hoxd9","Hsf2","Jund","Lhx1","Meis1","Pax8","Rfx7","Setbp1","Six4","Snai1","Srf","Tbx18","Tead2","Tigd2","Twist2","Zbtb1","Zbtb18","Zeb1","Zfp46","Zfp558","Arid3b","Glis1","Rarb","Rora","Six2","Sox13","T","Tbx6","Tcf4","Tead1","Zfp317","Zfp454","Zfp770","Zscan20","Arid3a","Ddit3","Hoxb9","Hoxc9","Nr2c2","Tcf7","Zfp189","Zfp422","Zfp580","Ebf2","Ebf3","Foxp4","Hoxb4","Pknox2","Zfp248","Zfp263","Srebf2","Meis3","Hand1","Hoxa5","Prdm6","Tbx22","Trps1","Gli2","Rara","Thap11","Hmga1","Hoxc8","Zfp143","Zfp932","Hoxc4","Prdm5","Zfp637","Hoxb5","Mycn","Tbx2","Zeb2","Wt1","Hoxc5","Klf3","Tcf7l2","Klf10","Glis2","Sox11","Lef1","Nr1d2","Rreb1","Hoxa9","Zfp647","Ebf1","Sox6","Nr6a1"],"nticks":null,"ticks":"outside","tickcolor":"rgba(51,51,51,1)","ticklen":3.65296803652968,"tickwidth":0.66417600664176,"showticklabels":false,"tickfont":{"color":null,"family":null,"size":0},"tickangle":-0,"showline":true,"linecolor":"rgba(0,0,0,1)","linewidth":0.66417600664176,"showgrid":false,"gridcolor":null,"gridwidth":0,"zeroline":false,"anchor":"y","title":{"text":"Ranked TF Motifs","font":{"color":"rgba(0,0,0,1)","family":"","size":14.6118721461187}},"hoverformat":".2f"},"yaxis":{"domain":[0,1],"automargin":true,"type":"linear","autorange":false,"range":[-107.65,60.65],"tickmode":"array","ticktext":[100,50,0,50],"tickvals":[-100,-50,0,50],"categoryorder":"array","categoryarray":[100,50,0,50],"nticks":null,"ticks":"outside","tickcolor":"rgba(51,51,51,1)","ticklen":3.65296803652968,"tickwidth":0.66417600664176,"showticklabels":true,"tickfont":{"color":"rgba(0,0,0,1)","family":"","size":11.689497716895},"tickangle":-0,"showline":true,"linecolor":"rgba(0,0,0,1)","linewidth":0.66417600664176,"showgrid":false,"gridcolor":null,"gridwidth":0,"zeroline":false,"anchor":"x","title":{"text":"# Associated genes <br />abs(Score) >= 1","font":{"color":"rgba(0,0,0,1)","family":"","size":14.6118721461187}},"hoverformat":".2f"},"shapes":[{"type":"rect","fillcolor":null,"line":{"color":null,"width":0,"linetype":[]},"yref":"paper","xref":"paper","x0":0,"x1":1,"y0":0,"y1":1}],"showlegend":false,"legend":{"bgcolor":"rgba(255,255,255,1)","bordercolor":"transparent","borderwidth":1.88976377952756,"font":{"color":"rgba(0,0,0,1)","family":"","size":11.689497716895}},"hovermode":"closest","barmode":"relative"},"config":{"doubleClick":"reset","modeBarButtonsToAdd":["hoverclosest","hovercompare"],"showSendToCloud":false},"source":"A","attrs":{"af1579b199fb":{"x":{},"y":{},"fill":{},"numActivated":{},"numRepressed":{},"type":"bar"},"af157af46c27":{"yintercept":{}}},"cur_data":"af1579b199fb","visdat":{"af1579b199fb":["function (y) ","x"],"af157af46c27":["function (y) ","x"]},"highlight":{"on":"plotly_click","persistent":false,"dynamic":false,"selectize":false,"opacityDim":0.2,"selected":{"opacity":1},"debounce":0},"shinyEvents":["plotly_hover","plotly_click","plotly_selected","plotly_relayout","plotly_brushed","plotly_brushing","plotly_clickannotation","plotly_doubleclick","plotly_deselect","plotly_afterplot","plotly_sunburstclick"],"base_url":"https://plot.ly"},"evals":[],"jsHooks":[]}</script>
```

We can then group DORCs (*rows*) that are regulated in a concerted manner by groups of TFs (*columns*). 


```r
score.cut <- 1.25
set.seed(2948)
regulons <- plotfigRHeatmap(figR.d = fig.d,
                            score.cut = score.cut,
                            row_names_gp = gpar(fontsize=8, fontface="italic"),
                            column_names_gp = gpar(fontsize=10),
                            show_row_dend = TRUE, k=4)
regulons <- draw(regulons)
```

![](03_peak2gene_files/figure-html/figr_heatmap-1.png)<!-- -->

```r
## get cluster membership
DORCsToKeep <- fig.d %>% dplyr::filter(abs(Score) >= score.cut) %>% pull(DORC) %>% unique()
TFsToKeep <- fig.d %>% dplyr::filter(abs(Score) >= score.cut) %>% pull(Motif) %>% unique()
net.d <- fig.d %>% dplyr::filter(DORC %in% DORCsToKeep & Motif %in% TFsToKeep) %>%
    reshape2::dcast(DORC ~ Motif) %>%
    tibble::column_to_rownames("DORC") %>% as.matrix()

clusters <- as.data.frame(unlist(row_order(regulons)))
clusters$cluster <- paste0("cluster", substr(row.names(clusters), 1, 1))
clusters$gene <- row.names(net.d[clusters[,1],])
clusters$cluster <- factor(clusters$cluster, levels=paste0("cluster", names(row_order(regulons))))
```

And visualise the expression of the DORC genes, to link each module to its activity pattern across stages.


```r
data <- data.rna[match(clusters$gene, row.names(data.rna)),]
data <- t(apply(data, 1, function(x) (x-mean(x))/sd(x) ))

m <- meta.rna[match(colnames(data), meta.rna$sample),]
m$stage <- factor(m$stage, levels=paste0("stage", c(8,18,21,25,27,35)))
m$somite <- factor(m$somite, levels=c("SIII","SII","SI"))
order <- order(m$stage, m$somite)

## heatmap annotation
ha  <- HeatmapAnnotation(df = data.frame(stage = m[order,]$stage, 
                                         somite = m[order,]$somite), 
                         col = list(stage = cols.stage, somite = cols.somite))
Heatmap(data[,order], 
        cluster_columns = FALSE, 
        cluster_rows = FALSE,
        col=colorRamp2(breaks = c(-4,-2,0,2,3), 
                       colors = c("steelblue4","steelblue","white", 
                                  "indianred3","indianred4")),
        row_names_gp = gpar(fontsize=8, fontface="italic"),
        column_names_gp = gpar(fontsize=10),
        name = "z-score", 
        show_row_names = TRUE, 
        show_column_names = FALSE, 
        top_annotation = ha,
        row_split = clusters$cluster,
        cluster_row_slices = FALSE)
```

![](03_peak2gene_files/figure-html/figr_expr_heatmap-1.png)<!-- -->
Two main modules are evident, with genes expressed either early or late in development. The role of *Nr6a1* is prominent, and Hox genes also show regulatory effects consistent with developmental progression. 


```r
saveRDS(cisCor, paste0(dir, "RNA+ATAC/results/03_gene_peak_correlations.Rds"))
write.table(links, paste0(dir, "RNA+ATAC/results/03_gene_peak_links.tsv"), quote = FALSE, sep = "\t", row.names = FALSE)
saveRDS(fig.d, paste0(dir, "RNA+ATAC/results/03_FigR_network.Rds"))
saveRDS(net.d, paste0(dir, "RNA+ATAC/results/03_regulationScores_thr1.25.Rds"))
write.table(clusters, paste0(dir, "RNA+ATAC/results/03_regulons_clusters.tsv"), quote = FALSE, sep="\t", row.names = FALSE)
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
##  [1] grid      parallel  stats4    stats     graphics  grDevices utils    
##  [8] datasets  methods   base     
## 
## other attached packages:
##  [1] BSgenome.Mmusculus.UCSC.mm10_1.4.0 BSgenome_1.60.0                   
##  [3] rtracklayer_1.52.1                 Biostrings_2.60.2                 
##  [5] XVector_0.32.0                     org.Mm.eg.db_3.13.0               
##  [7] FigR_0.1.0                         motifmatchr_1.14.0                
##  [9] chromVAR_1.14.0                    cowplot_1.1.1                     
## [11] Matrix_1.4-1                       pbmcapply_1.5.0                   
## [13] doParallel_1.0.17                  iterators_1.0.14                  
## [15] foreach_1.5.2                      dynamicTreeCut_1.63-1             
## [17] pcaExplorer_2.18.0                 topGO_2.44.0                      
## [19] SparseM_1.81                       GO.db_3.13.0                      
## [21] AnnotationDbi_1.54.1               graph_1.70.0                      
## [23] networkD3_0.4                      circlize_0.4.14                   
## [25] ComplexHeatmap_2.8.0               RColorBrewer_1.1-3                
## [27] BuenColors_0.5.6                   MASS_7.3-56                       
## [29] GGally_2.1.2                       ggrastr_1.0.1                     
## [31] ggpubr_0.4.0                       ggrepel_0.9.1                     
## [33] ggplot2_3.3.5                      dplyr_1.0.8                       
## [35] GenomicInteractions_1.26.0         InteractionSet_1.20.0             
## [37] SummarizedExperiment_1.22.0        Biobase_2.52.0                    
## [39] GenomicRanges_1.44.0               GenomeInfoDb_1.28.4               
## [41] IRanges_2.26.0                     S4Vectors_0.30.2                  
## [43] BiocGenerics_0.38.0                MatrixGenerics_1.4.3              
## [45] matrixStats_0.62.0                
## 
## loaded via a namespace (and not attached):
##   [1] rappdirs_0.3.3              AnnotationForge_1.34.1     
##   [3] R.methodsS3_1.8.1           nabor_0.5.0                
##   [5] pkgmaker_0.32.2             tidyr_1.2.0                
##   [7] bit64_4.0.5                 knitr_1.38                 
##   [9] R.utils_2.11.0              DelayedArray_0.18.0        
##  [11] data.table_1.14.2           rpart_4.1.16               
##  [13] KEGGREST_1.32.0             TFBSTools_1.30.0           
##  [15] RCurl_1.98-1.6              AnnotationFilter_1.16.0    
##  [17] generics_0.1.2              snow_0.4-4                 
##  [19] GenomicFeatures_1.44.2      RSQLite_2.2.12             
##  [21] tzdb_0.3.0                  bit_4.0.4                  
##  [23] webshot_0.5.3               xml2_1.3.3                 
##  [25] httpuv_1.6.5                assertthat_0.2.1           
##  [27] DirichletMultinomial_1.34.0 viridis_0.6.2              
##  [29] xfun_0.30                   hms_1.1.1                  
##  [31] jquerylib_0.1.4             evaluate_0.15              
##  [33] promises_1.2.0.1            TSP_1.2-0                  
##  [35] fansi_1.0.3                 restfulr_0.0.13            
##  [37] progress_1.2.2              caTools_1.18.2             
##  [39] dendextend_1.15.2           dbplyr_2.1.1               
##  [41] Rgraphviz_2.36.0            igraph_1.3.0               
##  [43] DBI_1.1.2                   geneplotter_1.70.0         
##  [45] htmlwidgets_1.5.4           reshape_0.8.9              
##  [47] Rmpfr_0.8-7                 purrr_0.3.4                
##  [49] ellipsis_0.3.2              crosstalk_1.2.0            
##  [51] backports_1.4.1             annotate_1.70.0            
##  [53] gridBase_0.4-7              biomaRt_2.48.3             
##  [55] vctrs_0.4.1                 ensembldb_2.16.4           
##  [57] Cairo_1.5-15                abind_1.4-5                
##  [59] cachem_1.0.6                withr_2.5.0                
##  [61] Gviz_1.36.2                 checkmate_2.0.0            
##  [63] GenomicAlignments_1.28.0    prettyunits_1.1.1          
##  [65] cluster_2.1.3               seqLogo_1.58.0             
##  [67] lazyeval_0.2.2              crayon_1.5.1               
##  [69] genefilter_1.74.1           labeling_0.4.2             
##  [71] pkgconfig_2.0.3             vipor_0.4.5                
##  [73] ProtGenerics_1.24.0         seriation_1.3.5            
##  [75] nnet_7.3-17                 rlang_1.0.2                
##  [77] lifecycle_1.0.1             miniUI_0.1.1.1             
##  [79] registry_0.5-1              filelock_1.0.2             
##  [81] BiocFileCache_2.0.0         doSNOW_1.0.20              
##  [83] GOstats_2.58.0              dichromat_2.0-0            
##  [85] rngtools_1.5.2              carData_3.0-5              
##  [87] base64enc_0.1-3             beeswarm_0.4.0             
##  [89] GlobalOptions_0.1.2         pheatmap_1.0.12            
##  [91] png_0.1-7                   viridisLite_0.4.0          
##  [93] rjson_0.2.21                bitops_1.0-7               
##  [95] shinydashboard_0.7.2        R.oo_1.24.0                
##  [97] blob_1.2.3                  shape_1.4.6                
##  [99] stringr_1.4.0               readr_2.1.2                
## [101] jpeg_0.1-9                  rstatix_0.7.0              
## [103] shinyAce_0.4.1              ggsignif_0.6.3             
## [105] CNEr_1.28.0                 scales_1.2.0               
## [107] memoise_2.0.1               GSEABase_1.54.0            
## [109] magrittr_2.0.3              plyr_1.8.7                 
## [111] zlibbioc_1.38.0             threejs_0.3.3              
## [113] compiler_4.1.0              BiocIO_1.2.0               
## [115] clue_0.3-60                 DESeq2_1.32.0              
## [117] Rsamtools_2.8.0             cli_3.2.0                  
## [119] Category_2.58.0             htmlTable_2.4.0            
## [121] Formula_1.2-4               tidyselect_1.1.2           
## [123] stringi_1.7.6               shinyBS_0.61.1             
## [125] highr_0.9                   yaml_2.3.5                 
## [127] locfit_1.5-9.5              latticeExtra_0.6-29        
## [129] sass_0.4.1                  VariantAnnotation_1.38.0   
## [131] tools_4.1.0                 rstudioapi_0.13            
## [133] TFMPvalue_0.0.8             foreign_0.8-82             
## [135] gridExtra_2.3               farver_2.1.0               
## [137] digest_0.6.29               FNN_1.1.3                  
## [139] pracma_2.3.8                shiny_1.7.1                
## [141] Rcpp_1.0.8.3                car_3.0-12                 
## [143] broom_0.8.0                 later_1.3.0                
## [145] httr_1.4.2                  biovizBase_1.40.0          
## [147] colorspace_2.0-3            XML_3.99-0.9               
## [149] splines_4.1.0               RBGL_1.68.0                
## [151] plotly_4.10.0               gmp_0.6-5                  
## [153] xtable_1.8-4                poweRlaw_0.70.6            
## [155] jsonlite_1.8.0              heatmaply_1.3.0            
## [157] R6_2.5.1                    Hmisc_4.7-0                
## [159] pillar_1.7.0                htmltools_0.5.2            
## [161] mime_0.12                   NMF_0.24.0                 
## [163] glue_1.6.2                  fastmap_1.1.0              
## [165] DT_0.22                     BiocParallel_1.26.2        
## [167] codetools_0.2-18            utf8_1.2.2                 
## [169] lattice_0.20-45             bslib_0.3.1                
## [171] tibble_3.1.6                curl_4.3.2                 
## [173] ggbeeswarm_0.6.0            gtools_3.9.2               
## [175] magick_2.7.3                survival_3.3-1             
## [177] limma_3.48.3                rmarkdown_2.13             
## [179] munsell_0.5.0               GetoptLong_1.0.5           
## [181] GenomeInfoDbData_1.2.6      reshape2_1.4.4             
## [183] gtable_0.3.0
```

