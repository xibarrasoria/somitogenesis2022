---
title: "Open chromatin regions"
date: '06 September, 2020'
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



We have defined a set of open chromatin loci, some of which show changes in accessibility across somites and development. Now we characterise these open chromatin regions.


```r
## metadata
meta <- read.table(paste0(dir, "ATAC-seq/data/metadata_ATACseq_GQ.tsv"), 
                   stringsAsFactors = FALSE, header = TRUE)
meta <- meta[meta$QCpass==1,]

## peak counts
peaks <- readRDS(paste0(dir, "ATAC-seq/results/03_peakCounts_csawMerged.Rds"))
peaks <- rowRanges(peaks)

## ENSEMBL annotation doesn't prefix 'chr' to chromosome names; remove from peaks
seqlevelsStyle(peaks) <- "NCBI"
names(peaks) <- paste0(seqnames(peaks), ":", start(peaks), "-", end(peaks))
```


### Annotation of open chromatin regions

The open chromatin regions are considered putative regulatory elements, since they allow the interaction of regulatory proteins with the underlying DNA. These proteins can in turn modify gene expression. The two major, and best understood, regulatory element types are promoters and enhancers. Promoters are straight-forward to identify because they lie within genes' transcription start sites (TSS), and the gene catalog is very well annotated in the mouse genome. Enhancers, on the other hand, are more challenging since they can be located anywhere in the genome, although they tend to be *close* to the genes they regulate. Furthermore, other elements exist, such as repressors and insulators, that are not well-characterised and cannot be differentiated from enhancers without additional data.

Thus, we can first separate promoters from everything else, by identifying peaks that lie in TSS regions, using Ensembl's gene annotation. We define a *promoter* as the interval ±200bp around the TSS. This strategy results in almost a fifth of all peaks annotated as promoters:


```r
## annotate open chromatin regions
## we use csaw function detailRanges to annotate the peaks
## but we want to use specifically Ensembl's annotation, version 96, to keep compatibility with the RNA-seq results
## thus, create a transcriptDB object from Ensembl v96
mmusculusEnsembl.v96 <- makeTxDbFromBiomart(biomart = "ENSEMBL_MART_ENSEMBL", 
                                            dataset = "mmusculus_gene_ensembl", 
                                            host="apr2019.archive.ensembl.org")
mcols(peaks) <- detailRanges(peaks, 
                             txdb = mmusculusEnsembl.v96, 
                             orgdb = org.Mm.eg.db, 
                             dist = 5000, # flanking distance
                             promoter = c(200, 200),
                             key.field = "ENSEMBL") 

c(nPeaks = length(peaks[grep(":P", peaks$overlap),]), 
  pct = round(length(peaks[grep(":P", peaks$overlap),])/length(peaks)*100, 2))
```

```
##   nPeaks      pct 
## 32367.00    18.72
```

From the remaining peaks, the majority are within gene bodies, and only a fifth are distal to annotated genes. The remaining ~10% are within 5kb of annotated genes.


```r
peaks$promoter <- ifelse(grepl(":P", peaks$overlap), 1, 0)
distal <- peaks[peaks$overlap == "" & peaks$left=="" & peaks$right == "",]

dist <- distanceToNearest(peaks[peaks$promoter==0,], genes(mmusculusEnsembl.v96))

plot(density(log10(mcols(dist)$distance+1)), 
     main = "peaks outside promoters",
     xlab="log10 distance to nearest gene",
     bty="l")
text(4.5, 0.35, cex = 0.75,
     labels = paste0(length(distal), " peaks (", round(length(distal)/length(peaks)*100, 2), "%)"))
text(1.5, 1.35, cex=0.75,
     labels = paste0(sum(mcols(dist)$distance==0), " peaks (",
                     round(sum(mcols(dist)$distance==0)/length(peaks)*100, 2), "%)"))
```

![](05_annotation_openChromatinRegions_files/figure-html/dist_nearest-1.png)<!-- -->

To better annotate the peaks that do not overlap TSSs we need additional data. We exploit available data from ENCODE3 [https://screen.encodeproject.org/], FANTOM5 [https://fantom.gsc.riken.jp/5/] and the VISTA enhancer browser [https://enhancer.lbl.gov/].

### Annotation with external data

The ENCODE3 project has combined data on chromatin accessibility from DNAse-seq with H3K4me3, H3K27ac and CTCF ChIP-seq tracks to define promoter and enhancer elements. Data comes from a variety of cell lines and tissues and is combined into a consensus set of candidate cis-Regulatory Elements (cCREs).


```r
encode <- read.table(paste0(dir, "ATAC-seq/data/externalData/mm10-ccREs.bed"))
encode <- GRanges(seqnames = substr(encode$V1,4,15),
                  IRanges(start = encode$V2, end = encode$V3),
                  name = encode$V5,
                  ann = encode$V6)

## find overlaps
overlap.encode <- as.data.frame(peaks[queryHits(findOverlaps(peaks, encode)),], row.names = NULL)
colnames(overlap.encode) <- paste0("peaks_", colnames(overlap.encode))
overlap.encode <- cbind(overlap.encode, 
                        as.data.frame(encode[subjectHits(findOverlaps(peaks, encode)),]))
overlap.encode$locus <- paste0(overlap.encode$peaks_seqnames, ":", 
                               overlap.encode$peaks_start, "-", 
                               overlap.encode$peaks_end)

## concatenate all ENCODE annotations for each peak
peaks.externalAnn <- as.data.frame(peaks)
peaks.externalAnn$PLS <- ifelse(row.names(peaks.externalAnn) %in% 
                                  overlap.encode[grep("PLS", overlap.encode$ann),]$locus, 1, 0)
peaks.externalAnn$pELS <- ifelse(row.names(peaks.externalAnn) %in% 
                                   overlap.encode[grep("pELS", overlap.encode$ann),]$locus, 1, 0)
peaks.externalAnn$dELS <- ifelse(row.names(peaks.externalAnn) %in% 
                                   overlap.encode[grep("dELS", overlap.encode$ann),]$locus, 1, 0)
peaks.externalAnn$CTCF_bound <- ifelse(row.names(peaks.externalAnn) %in%
                                         overlap.encode[grep("CTCF-bound", overlap.encode$ann),]$locus, 1, 0)
peaks.externalAnn$H3K4me3 <- ifelse(row.names(peaks.externalAnn) %in%
                                      overlap.encode[grep("DNase-H3K4me3", overlap.encode$ann),]$locus, 1, 0)
```

From the 339815 annotated cCREs in the mouse genome, 174090 (51.23%) overlap ATAC-seq peaks in the somites. Conversely, out of the 172901 peaks identified in the somites, 101377 (58.63%) overlap at least one cCRE. This indicates overall congruence between datasets, but also highlights many regulatory elements that haven't been captured by the cell populations used in ENCODE3 and, perhaps, somite-specific elements. 

FANTOM5 provides enhancer annotations based on CAGE data; bidirectional transcription initiation events have been show to mark enhancer elements. 


```r
## blocks correspond to merged regions of transcription initiation evens in the plus and minus strand
fantom <- read.table(paste0(dir, "ATAC-seq/data/externalData/F5.mm10.enhancers.bed"))
colnames(fantom) <- c("chr", "start", "end", "name", "score", "strand", 
                      "thickStart", "thickEnd", "col", 
                      "blockCount", "blockSizes", "blockStarts")
fantom <- GRanges(seqnames = substr(fantom$chr, 4, 15),
                  IRanges(start = fantom$start, end = fantom$end))

## find overlaps
overlap.fantom <- as.data.frame(peaks[queryHits(findOverlaps(peaks, fantom)),], row.names = NULL)
colnames(overlap.fantom) <- paste0("peaks_", colnames(overlap.fantom))
overlap.fantom <- cbind(overlap.fantom, 
                        as.data.frame(fantom[subjectHits(findOverlaps(peaks, fantom)),]))
overlap.fantom$locus <- paste0(overlap.fantom$peaks_seqnames, ":", 
                               overlap.fantom$peaks_start, "-", 
                               overlap.fantom$peaks_end)

peaks.externalAnn$FANTOM <- ifelse(row.names(peaks.externalAnn) %in% overlap.fantom$locus, 1, 0)
```

Out of the 49797 annotated enhancers in the mouse genome, 17968 (36.08%) overlap somite ATAC-seq peaks. This represents 10.39% of the ATAC-seq peaks.

Finally, the VISTA enhancer browser contains a collection of human and mouse putative regulatory elements tested in mouse transgenic models to assess whether these sequences are able to drive expression of a reporter gene, and determine tissue-specific expression patterns. For the human sequences, the orthologous mouse sequence has been defined. 


```r
vista <- read.table(paste0(dir, "ATAC-seq/data/externalData/VISTAenhancers_mm10.tsv"), 
                    sep="\t", header = TRUE)
vista <- GRanges(seqnames = vista$chr_mm10,
                 IRanges(start = vista$start_mm10, end = vista$end_mm10),
                 ID = vista$ID,
                 genes = vista$bracketing_genes,
                 expression = vista$expression)

## find overlaps
overlap.vista <- as.data.frame(peaks[queryHits(findOverlaps(peaks, vista)),], row.names = NULL)
colnames(overlap.vista) <- paste0("peaks_", colnames(overlap.vista))
overlap.vista <- cbind(overlap.vista, 
                       as.data.frame(vista[subjectHits(findOverlaps(peaks, vista)),]))
overlap.vista$locus <- paste0(overlap.vista$peaks_seqnames, ":", 
                              overlap.vista$peaks_start, "-",
                              overlap.vista$peaks_end)

peaks.externalAnn$VISTA <- ifelse(row.names(peaks.externalAnn) %in% overlap.vista$locus, 1, 0)

## how many annotation categories are present for each peak
peaks.externalAnn$sum <- rowSums(peaks.externalAnn[,9:16])

## capture whether the enhancer was positive or negative in embryos
peaks.externalAnn[peaks.externalAnn$VISTA==1,]$VISTA <- as.character(overlap.vista[
  match(row.names(peaks.externalAnn[peaks.externalAnn$VISTA==1,]), overlap.vista$locus),]$expression)
```

There are many fewer of these VISTA enhancers, with only 3124 tested elements, of which 1723 (55.15%) overlap ATAC-seq peaks. Many of the sequences tested are quite large, so 2616 peaks overlap VISTA enhancers. From these, 1050 VISTA enhancers have positive expression in the reporter assay.


```r
## retrieve VISTA enhancer with expression annotated in somites
vista.somite <- read.table(paste0(dir, "ATAC-seq/data/externalData/VISTAenhancers_expr_somite.tsv"),
                           stringsAsFactors = FALSE, sep="\t", fill=TRUE)
vista.somite <- subset(vista.somite, vista.somite$V1 == "Human" | vista.somite$V1 == "Mouse")

ids <- unlist(lapply(strsplit(vista.somite$V3, " "), '[[', 3))
ids <- paste0(ifelse(vista.somite$V1=="Human", "hs", "mm"), ids)

## save which peaks are enhancers with *proved* expression in somites
peaks_VISTAsomite <- overlap.vista[overlap.vista$ID %in% ids,]
## and add this information to the general table
peaks.externalAnn$VISTAsomite <- ifelse(row.names(peaks.externalAnn) %in% peaks_VISTAsomite$locus, 
                                        1, 0)
```

Furthermore, 67 of the VISTA enhancers with positive expression have positive signal in the somites; 52 of these overlap 110 ATAC-seq peaks.


```r
tmp <- peaks.externalAnn[peaks.externalAnn$sum>0,9:16]
tmp$VISTA <- ifelse(tmp$VISTA == "Positive", 1,
                    ifelse(tmp$VISTA == "Negative", -1, 0))
tmp <- as.matrix(tmp)

tmp <- tmp[order(tmp[,1], tmp[,2], tmp[,3], tmp[,4], tmp[,5], tmp[,6], tmp[,7], tmp[,8], decreasing = TRUE),]

Heatmap(tmp,
        cluster_rows = FALSE, 
        cluster_columns = FALSE, 
        show_row_names = FALSE,
        show_heatmap_legend = FALSE)
```

### Classification of ATAC-seq peaks {.tabset}

We classify the peaks based on their genomic context:

- Promoter: if the peak is within 200bp of an annotated TSS.
- Genic: if the peak overlaps an annotated gene, but not the TSS.
- Proximal: if the peak is within 25kb of an annotated gene.
- Distal: if the peak is farther than 25kb from an annotated gene, but within 100kb.
- Intergenic: any other peaks.

Recent studies linking putative enhancers to their target genes have shown that most enhancers regulate genes that are located *close-by*. Gasperini et al. (Cell, 2019) found that the median distance between enhancers and target genes was 24.1kb, and the vast majority of pairs were within 100kb. Thus, peaks that are very far from any genes are less likely to be regulators of gene expression.

Apart from promoters, most peaks fall within genes, mostly in introns. Very few peaks are in gene desserts (*intergenic* class, which doesn't have a gene within 100kb), with the vast majority within 25kb from an annotated gene.


```r
## compute distance to nearest gene for all peaks
dist <- distanceToNearest(peaks, genes(mmusculusEnsembl.v96))

## annotate each peak
peaks.classes <- peaks.externalAnn[,c(1:4,6)]
peaks.classes$nearest_gene <- mcols(dist)$distance

peaks.classes$promoter <- ifelse(grepl(":P", peaks.classes$overlap), 1, 0)
peaks.classes$genic <- ifelse(peaks.classes$overlap != "" & 
                                peaks.classes$promoter == 0, 
                              1, 0)
peaks.classes$proximal <- ifelse(peaks.classes$overlap == "" & 
                                   peaks.classes$nearest_gene <= 25e3, 
                                 1, 0)
peaks.classes$distal <- ifelse(peaks.classes$nearest_gene > 25e3 & 
                                 peaks.classes$nearest_gene <= 100e3, 
                               1, 0)
peaks.classes$intergenic <- ifelse(peaks.classes$nearest_gene > 100e3, 1, 0)

## for peaks in promoters, recover the gene
peaks.classes$promoter_gene <- NA
promoters <- promoters(mmusculusEnsembl.v96, 
                       upstream = 200, downstream = 200,
                       columns = "GENEID")
# assign gene name based on ensembl ID
gene_ann <- read.table(paste0(dir, "RNA-seq/data/Mus_musculus.GRCm38.96.ann"), 
                       stringsAsFactors = FALSE, row.names = 1)
promoters$geneName <- gene_ann[unlist(promoters$GENEID),1]
tmp <- peaks.classes[peaks.classes$promoter==1,]
for(i in 1:nrow(tmp)){
  # if(i %% 1000 == 0){ print(i) }
  peaks.classes[row.names(tmp)[i], 'promoter_gene'] <- paste(
    unique(unlist(
      promoters[subjectHits(findOverlaps(peaks[row.names(tmp)[i]], promoters))]$geneName)),
    collapse = ";")
}
```

```r
## summarise
df <- data.frame(colSums(peaks.classes[,7:11]))
colnames(df) <- "count"
df$class <- row.names(df)
df$class <- factor(df$class, levels=c("promoter", "genic", "proximal", "distal", "intergenic"))

ggplot(df, aes(class, count, fill=class)) + 
  geom_bar(stat="identity") +
  geom_text(aes(label=count)) +
  xlab("") + ylab("number of peaks") +
  th + theme(legend.position = "none")
```

![](05_annotation_openChromatinRegions_files/figure-html/classify-1.png)<!-- -->

```r
df$pct <- round(df$count/length(peaks)*100, 2)
df[,c(2,1,3)]
```

```
##                 class count   pct
## promoter     promoter 32367 18.72
## genic           genic 85194 49.27
## proximal     proximal 41724 24.13
## distal         distal 11912  6.89
## intergenic intergenic  1704  0.99
```

Therefore, the vast majority of the regions of open chromatin we identify in the somites have regulatory potential. 

Some of these have supporting evidence from the external databases, which makes them stronger candidates.

#### Promoter peaks

Not surprisingly, there is very good agreement between our classified promoters and the ENCODE2 cCREs with promoter-like signatures (PLS). Of note is that a large number of promoter peaks are bound by CTCF. As expected, only a small number of promoters overlap with FANTOM enhancers.


```r
cols <- gg_color_hue(n=5)

upset(peaks.externalAnn[row.names(peaks.classes[peaks.classes$promoter==1,]),10:16],
                    main.bar.color = cols[1], order.by = "freq", nintersects = 10)
```

![](05_annotation_openChromatinRegions_files/figure-html/promoters-1.png)<!-- -->

#### Genic peaks

Peaks that overlap genes but not the TSS are instead overlapped by distal enhancer-like elements and a larger proportion of FANTOM enhancers. A much smaller proportion are bound by CTCF, compared to promoter elements.


```r
upset(peaks.externalAnn[row.names(peaks.classes[peaks.classes$genic==1,]),10:16],
                    main.bar.color = cols[2], order.by = "freq", nintersects = 10)
```

![](05_annotation_openChromatinRegions_files/figure-html/genic-1.png)<!-- -->

#### Proximal peaks

For peaks that do not overlap genes but are within 25kb, again we see a large overlap with ENCODE's dELS elements, and many are also identified in the FANTOM5 data.


```r
upset(peaks.externalAnn[row.names(peaks.classes[peaks.classes$proximal==1,]),10:16],
                    main.bar.color = cols[3], order.by = "freq", nintersects = 10)
```

![](05_annotation_openChromatinRegions_files/figure-html/proximal-1.png)<!-- -->

####  Distal peaks

A similar picture is observed for distal peaks, which are those that are farther from genes, but within 100kb.


```r
upset(peaks.externalAnn[row.names(peaks.classes[peaks.classes$distal==1,]),10:16],
                    main.bar.color = cols[4], order.by = "freq", nintersects = 10)
```

![](05_annotation_openChromatinRegions_files/figure-html/distal-1.png)<!-- -->

####  Intergenic peaks

And only a few of the peaks that are away from genes have overlaps with the external data.


```r
upset(peaks.externalAnn[row.names(peaks.classes[peaks.classes$intergenic==1,]),10:16],
                    main.bar.color = cols[5], order.by = "freq", nintersects = 10)
```

![](05_annotation_openChromatinRegions_files/figure-html/intergenic-1.png)<!-- -->

###

Overall, promoter peaks are very well represented in the external databases, with over 85% annotated in ENCODE3. In contrast, over half of all peaks falling close to genes are not captured in the external catalogues. This is a reflection of the incomplete nature of these compendia, since only a limited set of cell types and tissues has been considered. 


```r
encode_ccres <- row.names(peaks.externalAnn[peaks.externalAnn$PLS == 1 | 
                                              peaks.externalAnn$pELS == 1 |
                                            peaks.externalAnn$dELS == 1, ])
ctcf_only <- setdiff(row.names(peaks.externalAnn[peaks.externalAnn$CTCF_bound == 1,]), 
                     encode_ccres)
fantom <- row.names(peaks.externalAnn[peaks.externalAnn$FANTOM == 1,])
vista <- row.names(peaks.externalAnn[peaks.externalAnn$VISTA != 0,])

peaks.classes$in_encode <- ifelse(row.names(peaks.classes) %in% encode_ccres, 1, 0)
peaks.classes$in_fantom <- ifelse(row.names(peaks.classes) %in% fantom, 1, 0)
peaks.classes$in_vista <- ifelse(row.names(peaks.classes) %in% vista, 1, 0)
peaks.classes$CTCF_only <- ifelse(row.names(peaks.classes) %in% ctcf_only, 1, 0)


df <- data.frame(promoters = colSums(peaks.classes[peaks.classes$promoter == 1, 13:16]),
           genic = colSums(peaks.classes[peaks.classes$genic == 1, 13:16]),
           proximal = colSums(peaks.classes[peaks.classes$proximal == 1, 13:16]),
           distal = colSums(peaks.classes[peaks.classes$distal == 1, 13:16]),
           intergenic = colSums(peaks.classes[peaks.classes$intergenic == 1, 13:16]))

df
```

```
##           promoters genic proximal distal intergenic
## in_encode     28050 37432    17336   4057        425
## in_fantom      3341  8271     5370    933         55
## in_vista        184  1345      778    261         48
## CTCF_only       658  6538     4538   1074        102
```

```r
df2 <- round(t(t(df)/colSums(peaks.classes[,7:11]))*100,2)
colnames(df2) <- paste0(colnames(df2), "%")
df2
```

```
##           promoters% genic% proximal% distal% intergenic%
## in_encode      86.66  43.94     41.55   34.06       24.94
## in_fantom      10.32   9.71     12.87    7.83        3.23
## in_vista        0.57   1.58      1.86    2.19        2.82
## CTCF_only       2.03   7.67     10.88    9.02        5.99
```

Nonetheless, putative enhancer elements with available experimental data supporting their role as enhancers will be very useful. 


```r
write.table(peaks.classes, paste0(dir, "ATAC-seq/results/05_peaks_classAnnotation.tsv"), 
            quote = FALSE, sep="\t")
write.table(peaks.externalAnn, paste0(dir, "ATAC-seq/results/05_peaks_overlapExternalData.tsv"), 
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
##  [1] UpSetR_1.4.0                RColorBrewer_1.1-2         
##  [3] ggpubr_0.4.0                ggplot2_3.3.2              
##  [5] org.Mm.eg.db_3.10.0         GenomicFeatures_1.38.2     
##  [7] AnnotationDbi_1.48.0        csaw_1.20.0                
##  [9] SummarizedExperiment_1.16.1 DelayedArray_0.12.3        
## [11] BiocParallel_1.20.1         matrixStats_0.56.0         
## [13] Biobase_2.46.0              GenomicRanges_1.38.0       
## [15] GenomeInfoDb_1.22.1         IRanges_2.20.2             
## [17] S4Vectors_0.24.4            BiocGenerics_0.32.0        
## 
## loaded via a namespace (and not attached):
##  [1] bitops_1.0-6             bit64_4.0.2              progress_1.2.2          
##  [4] httr_1.4.2               tools_3.6.1              backports_1.1.8         
##  [7] R6_2.4.1                 DBI_1.1.0                colorspace_1.4-1        
## [10] withr_2.2.0              gridExtra_2.3            tidyselect_1.1.0        
## [13] prettyunits_1.1.1        bit_4.0.3                curl_4.3                
## [16] compiler_3.6.1           labeling_0.3             rtracklayer_1.46.0      
## [19] scales_1.1.1             askpass_1.1              rappdirs_0.3.1          
## [22] stringr_1.4.0            digest_0.6.25            Rsamtools_2.2.3         
## [25] foreign_0.8-72           rmarkdown_2.3            rio_0.5.16              
## [28] XVector_0.26.0           pkgconfig_2.0.3          htmltools_0.5.0         
## [31] dbplyr_1.4.4             limma_3.42.2             rlang_0.4.7             
## [34] readxl_1.3.1             RSQLite_2.2.0            farver_2.0.3            
## [37] generics_0.0.2           dplyr_1.0.1              zip_2.0.4               
## [40] car_3.0-8                RCurl_1.98-1.2           magrittr_1.5            
## [43] GenomeInfoDbData_1.2.2   Matrix_1.2-18            Rcpp_1.0.5              
## [46] munsell_0.5.0            abind_1.4-5              lifecycle_0.2.0         
## [49] stringi_1.4.6            yaml_2.2.1               edgeR_3.28.1            
## [52] carData_3.0-4            zlibbioc_1.32.0          plyr_1.8.6              
## [55] BiocFileCache_1.10.2     grid_3.6.1               blob_1.2.1              
## [58] forcats_0.5.0            crayon_1.3.4             lattice_0.20-41         
## [61] Biostrings_2.54.0        haven_2.3.1              hms_0.5.3               
## [64] locfit_1.5-9.4           knitr_1.29               pillar_1.4.6            
## [67] ggsignif_0.6.0           biomaRt_2.42.1           XML_3.99-0.3            
## [70] glue_1.4.1               evaluate_0.14            data.table_1.12.8       
## [73] vctrs_0.3.2              cellranger_1.1.0         gtable_0.3.0            
## [76] openssl_1.4.2            purrr_0.3.4              tidyr_1.1.1             
## [79] assertthat_0.2.1         xfun_0.16                openxlsx_4.1.5          
## [82] broom_0.7.0              rstatix_0.6.0            tibble_3.0.3            
## [85] GenomicAlignments_1.22.1 memoise_1.1.0            ellipsis_0.3.1
```

