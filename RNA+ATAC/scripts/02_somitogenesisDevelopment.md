---
title: "Somitogenesis across development"
date: '11 August, 2022'
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



We have explored the genes and chromatin loci that change as somites mature and prepare to undergo EMT. Now we can turn our attention to the changes that happen at different stages of development.

We have data from six different developmental stages. Although the segmentation process is conserved, the somites generated at different stages will give rise to different structures, depending on their position in the antero-posterior axis.

Thus, we can use these data to explore what genes and regulatory regions are different in somites that will develop into different structures.


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


### Transcriptional changes across development

Differential expression analysis identified over ten thousand genes that change in expression between at least a pair of stages. Most show the same behaviour irrespective of the maturation level of each somite, which is not surprising given that the differences between somite trios were very small and restricted to ~3K genes. There are still a few thousand genes that are identified only when performing the analysis in each somite level independently.


```r
degs <- read.table(paste0(dir, "RNA-seq/results/04_DEresults_summary_stage.tsv"),
                   stringsAsFactors = FALSE)
degs$gene <- geneCounts[match(row.names(degs), row.names(geneCounts)),1]
degs <- degs[,c(ncol(degs), 1:(ncol(degs)-1))]

c(average = sum(degs$ave), somite_wise = sum(degs$somiteSpecific))
```

```
##     average somite_wise 
##        7714        2977
```

```r
## keep only significant genes but save the universe
universe <- data.frame(gene_id=row.names(degs), gene_name=degs$gene)
degs <- degs[rowSums(degs[,-1])>0,]
```

But, in the original analysis we observed that most of these *somite-specific* changes show consistent patterns across all three somite levels, but reach significance in only one of the tests. Thus, very few truly somite-specific changes are observed.

Now, let's evaluate the expression dynamics of these DE genes. We can roughly split the DE genes into 12 clusters that show different expression dynamics.


```r
## count data for DEGs
data <- geneCounts[row.names(degs),-1]
## z-score
data <- 2^data-1
data <- t(apply(data, 1, function(x) (x-mean(x))/sd(x) ))

## split genes into clusters
dist <- dist(data)
hc <- hclust(dist)
cut <- cutreeDynamic(hc, distM = as.matrix(dist), 
                     minClusterSize = 500, method = "hybrid", 
                     deepSplit = 1, verbose = 0)

## order clusters to follow developmental progression
clusters.degs <- factor(paste0("cluster",cut), levels = paste0("cluster",c(8,2,6,7, 4,5, 3, 11,12,9,1, 10)))
names(clusters.degs) <- row.names(data)
table(clusters.degs)
```

```
## clusters.degs
##  cluster8  cluster2  cluster6  cluster7  cluster4  cluster5  cluster3 cluster11 
##       725      1254       752       739       988       813      1184       585 
## cluster12  cluster9  cluster1 cluster10 
##       555       675      1790       631
```

Most of these clusters show highest expression in one of the four types of somites (i.e., cervical, thoracic, lumbar, or sacral).


```r
## use loess to model the expression pattern of each gene
tmp <- 2^geneCounts[row.names(degs),-1]
stage <- factor(paste0("stage", meta.rna$stage), levels=paste0("stage",c(8,18,21,25,27,35)))
fits <- apply(tmp, 1, function(x) stats::lowess(x~stage))

curves <- unique(as.data.frame(lapply(fits, '[', 2)))
colnames(curves) <- row.names(tmp)
row.names(curves) <- paste0("stage",c(8,18,21,25,27,35))
curvesStd <- t(apply(t(curves), 1, function(x) x/max(x)))

plots <- list()
for(c in levels(clusters.degs)){
  g <- names(clusters.degs[clusters.degs==c])
  test <- curvesStd[g,]
  mean <- colMeans(test)
  std <- apply(test, 2, sd)
  df.sub <- data.frame(x=1:6, ymin=mean-std, ymax=mean+std)
  df.sub.avg <- data.frame(x=1:6, y=mean)
  plots[[c]] <- ggplot() + 
    geom_line(data = df.sub.avg, aes(x,y), colour="black" , size=1) +
    geom_ribbon(data = df.sub, aes(x=x, ymin=ymin, ymax=ymax), alpha=0.2, fill="black") +
    xlab("stage") + ylab("expression") + 
    ggtitle(paste(c, "-", nrow(test), "genes")) + 
    ylim(0,1.1) +
    scale_x_discrete(limits=c("8","18","21","25","27","35")) + 
    th
}
ggarrange(plotlist = plots, ncol = 4, nrow = 3)
```

![](02_somitogenesisDevelopment_files/figure-html/patterns-1.png)<!-- -->

Better visualised with a heatmap:


```r
## order columns based on somite number
m <- meta.rna
m$somite <- factor(m$somite, levels=c("SIII","SII","SI"))
stopifnot(identical(colnames(data), m$sample)) # make sure the metadata corresponds with count data matrix
order <- order(m$stage, m$somite)

## heatmap annotation
ha  <- HeatmapAnnotation(df = data.frame(stage = paste0("stage", m[order,]$stage), 
                                         somite = m[order,]$somite), 
                         col = list(stage = cols.stage, somite = cols.somite))

h.stage <- Heatmap(data[,order], 
                   cluster_columns = FALSE, 
                   col=colorRamp2(breaks = c(-4,-2,0,2,3), 
                                  colors = c("steelblue4","steelblue","white", 
                                             "indianred3","indianred4")),
                   name = "z-score", 
                   show_row_names = FALSE, 
                   show_column_names = FALSE, 
                   top_annotation = ha, 
                   row_split = clusters.degs,
                   cluster_row_slices = FALSE)
h.stage <- draw(h.stage)
```

![](02_somitogenesisDevelopment_files/figure-html/heatmap-1.png)<!-- -->

```r
## save to use later
data.stage.degs <- data
saveRDS(data.stage.degs, paste0(dir, "RNA+ATAC/results/02_stage_DEGs_data_for_heatmap.Rds"))
```


```r
## add cluster identity to DEGs
degs$cluster <- clusters.degs
degs$fate <- ifelse(degs$cluster %in% paste0("cluster", c(8,2,6,7)), "cervical",
                    ifelse(degs$cluster %in% paste0("cluster", c(4,5)), "thoracic",
                           ifelse(degs$cluster %in% paste0("cluster", c(3,11)), "lumbar", "sacral")))
degs[degs$cluster=="cluster10",]$fate <- NA
```


#### Overrepresented biological processes {.tabset}

As expected, the set of DE genes is highly enriched for terms associated with developmental processes, including the musculo-skeletal system.


```r
## GO enrichment
GO.all <- topGOtable(DEgenes = degs$gene,
                     BGgenes = universe$gene_name,
                     topGO_method2 = "elim",
                     ontology = "BP",
                     geneID = "symbol",
                     addGeneToTerms = TRUE,
                     mapping = "org.Mm.eg.db",
                     topTablerows = 1e3)
## significant terms
GO.all[c(2,3,6,8,11,12,17,18,20,27,28,31,34,38,39,42,43,49,71,81,84,87), c(2:5,7)]
```

<div data-pagedtable="false">
  <script data-pagedtable-source type="application/json">
{"columns":[{"label":[""],"name":["_rn_"],"type":[""],"align":["left"]},{"label":["Term"],"name":[1],"type":["chr"],"align":["left"]},{"label":["Annotated"],"name":[2],"type":["int"],"align":["right"]},{"label":["Significant"],"name":[3],"type":["int"],"align":["right"]},{"label":["Expected"],"name":[4],"type":["dbl"],"align":["right"]},{"label":["p.value_elim"],"name":[5],"type":["dbl"],"align":["right"]}],"data":[{"1":"positive regulation of ERK1 and ERK2 cascade","2":"160","3":"123","4":"87.07","5":"2.6e-09","_rn_":"2"},{"1":"angiogenesis","2":"467","3":"342","4":"254.13","5":"8.1e-09","_rn_":"3"},{"1":"embryonic skeletal system morphogenesis","2":"103","3":"85","4":"56.05","5":"8.1e-08","_rn_":"6"},{"1":"regulation of cell migration","2":"802","3":"565","4":"436.42","5":"3.3e-07","_rn_":"8"},{"1":"epithelial cell differentiation","2":"501","3":"359","4":"272.63","5":"8.8e-07","_rn_":"11"},{"1":"anterior/posterior pattern specification","2":"222","3":"163","4":"120.80","5":"8.9e-07","_rn_":"12"},{"1":"proximal/distal pattern formation","2":"34","3":"31","4":"18.50","5":"3.9e-06","_rn_":"17"},{"1":"cell adhesion","2":"1073","3":"755","4":"583.89","5":"4.0e-06","_rn_":"18"},{"1":"negative regulation of canonical Wnt signaling pathway","2":"121","3":"90","4":"65.84","5":"4.3e-06","_rn_":"20"},{"1":"positive regulation of epithelial cell proliferation","2":"180","3":"125","4":"97.95","5":"2.4e-05","_rn_":"27"},{"1":"cartilage development","2":"183","3":"145","4":"99.58","5":"2.8e-05","_rn_":"28"},{"1":"smooth muscle cell differentiation","2":"72","3":"56","4":"39.18","5":"3.2e-05","_rn_":"31"},{"1":"smooth muscle tissue development","2":"26","3":"24","4":"14.15","5":"3.3e-05","_rn_":"34"},{"1":"cell fate determination","2":"36","3":"31","4":"19.59","5":"5.6e-05","_rn_":"38"},{"1":"endothelial cell migration","2":"196","3":"144","4":"106.66","5":"6.5e-05","_rn_":"39"},{"1":"Notch signaling pathway","2":"165","3":"114","4":"89.79","5":"7.6e-05","_rn_":"42"},{"1":"positive regulation of epithelial to mesenchymal transition","2":"45","3":"37","4":"24.49","5":"8.5e-05","_rn_":"43"},{"1":"positive regulation of osteoblast differentiation","2":"66","3":"51","4":"35.91","5":"9.7e-05","_rn_":"49"},{"1":"response to retinoic acid","2":"75","3":"56","4":"40.81","5":"2.3e-04","_rn_":"71"},{"1":"positive regulation of epithelial cell differentiation","2":"45","3":"36","4":"24.49","5":"3.1e-04","_rn_":"81"},{"1":"myotube differentiation","2":"101","3":"72","4":"54.96","5":"3.6e-04","_rn_":"84"},{"1":"skin development","2":"217","3":"151","4":"118.08","5":"4.2e-04","_rn_":"87"}],"options":{"columns":{"min":{},"max":[10]},"rows":{"min":[10],"max":[10]},"pages":{}}}
  </script>
</div>

We can also perform the analysis for each set of clusters that show highest expression in a given fate.

##### Cervical

There are 3470 genes with highest expression in cervical somites (stage 8): the first four clusters in the heatmap above.

These genes are enriched for broad developmental and metabolic terms, with epithelium development being highly significant.


```r
## cervical genes
GO.cervical <- topGOtable(DEgenes = degs[degs$fate=="cervical",]$gene,
                     BGgenes = universe$gene_name,
                     topGO_method2 = "elim",
                     ontology = "BP",
                     geneID = "symbol",
                     addGeneToTerms = TRUE,
                     mapping = "org.Mm.eg.db",
                     topTablerows = 1e3)
## significant terms
GO.cervical[c(1,11,25,28,30,32,36,40,53,60,79,86,87), c(2:5,7,9)]
```

<div data-pagedtable="false">
  <script data-pagedtable-source type="application/json">
{"columns":[{"label":[""],"name":["_rn_"],"type":[""],"align":["left"]},{"label":["Term"],"name":[1],"type":["chr"],"align":["left"]},{"label":["Annotated"],"name":[2],"type":["int"],"align":["right"]},{"label":["Significant"],"name":[3],"type":["int"],"align":["right"]},{"label":["Expected"],"name":[4],"type":["dbl"],"align":["right"]},{"label":["p.value_elim"],"name":[5],"type":["dbl"],"align":["right"]},{"label":["genes"],"name":[6],"type":["chr"],"align":["left"]}],"data":[{"1":"epithelial cell differentiation","2":"501","3":"134","4":"77.87","5":"1.70e-08","6":"Abcb1b,Acvr2b,Ajap1,Akap9,Alg10b,Alox8,Apc,Arhgef26,Atoh8,Barx1,Basp1,Bcr,Bmp2,Cav1,Ccnd1,Cdh1,Cdhr2,Cdx2,Ceacam1,Cgn,Cldn1,Clock,Crb3,Crhr2,Ctnnd1,Cul7,Cyp26b1,Dab2,Dact2,Deup1,Dlx3,Dlx5,E2f4,E2f7,Elf5,Epb41l5,Esrp1,Etv2,Exph5,Ezr,F2rl1,Fam20c,Fat4,Fgf8,Fgfr2,Fgfr3,Flna,Fndc3a,Foxa3,Foxb1,Foxl2,Foxp4,Fzd2,Gata4,Gli1,Gli2,Grem1,Grhl2,Hes5,Hey2,Hnf1b,Hnf4a,Hoxa5,Hoxb5,Ihh,Jag1,Jmjd1c,Klf15,Klf5,Klf7,Lhx1,Maf,Mafg,Map2k1,Map7,Marveld2,Myadm,Myo3b,Myo6,Ncor2,Nkx2-2,Osr1,Pck2,Pdzd7,Pelo,Pgr,Plaur,Pls1,Pof1b,Pou3f1,Pphln1,Ppl,Prdm1,Prlr,Ptch1,Ptprs,Rab25,Rara,Rarg,Rest,Rfx6,Rxra,Safb2,Scrib,Sh3bp1,Sipa1l3,Slc4a5,Slc4a7,Sox2,Sox3,Sox4,Spint2,Spred2,St14,Stat5b,Sult2b1,Tfap2c,Tfcp2l1,Tgm1,Thra,Tjp1,Tjp2,Tjp3,Tmem132e,Tprn,Triobp,Tst,Upk1a,Upk1b,Whrn,Yipf6,Zeb2,Zfas1,Zfp36l1","_rn_":"1"},{"1":"anterior/posterior pattern specification","2":"222","3":"54","4":"34.50","5":"4.00e-04","6":"Acvr2b,Aldh1a2,Apc,Axin2,Barx1,Bmi1,Bmp2,Cdx1,Cdx2,Dll3,Epb41l5,Fgf8,Foxb1,Foxh1,Gata4,Gbx2,Gli2,Gpc3,Helt,Hes3,Hes5,Hes7,Hey2,Hnf1b,Hoxa1,Hoxa2,Hoxa4,Hoxa5,Hoxb1,Hoxb2,Hoxb3,Hoxb4,Hoxb5,Hoxb6,Hoxb8,Hoxc4,Hoxc5,Ifitm1,Lhx1,Msx1,Nog,Osr1,Otx1,Pcgf2,Pofut1,Rarg,Ror2,Sfrp1,Tbx6,Tdrd5,Tshz1,Wls,Xrcc2,Zeb2","_rn_":"11"},{"1":"positive regulation of ERK1 and ERK2 cascade","2":"160","3":"40","4":"24.87","5":"1.22e-03","6":"Acta2,Adcyap1,Akap12,Apoe,App,Arrb2,Bmp2,Braf,Camk2d,Ccl5,Cd44,Cysltr2,Ddt,F2rl1,Fga,Fgb,Fgf15,Fgf8,Fgfbp3,Fgfr2,Fgfr3,Fgfr4,Fgg,Gata4,Gcnt2,Lrp1,Map2k1,Map2k7,Mapk3,Necab2,Pdgfra,Prkca,Prkcz,Pten,Ptk2b,Pycard,Rara,Serpinf2,Shc1,Trpv4","_rn_":"25"},{"1":"morphogenesis of embryonic epithelium","2":"181","3":"44","4":"28.13","5":"1.33e-03","6":"Aldh1a2,Cdk20,Cited2,Dlc1,Epb41l5,Fgf8,Fgfr2,Folr1,Fzd2,Gli2,Grem1,Grhl2,Hand1,Hes5,Hnf1b,Lama5,Llgl2,Lrp2,Mthfd1,Nog,Osr1,Prickle1,Ptch1,Rara,Rarg,Rdh10,Ret,Rgma,Scrib,Sema4c,Setd2,Sfrp1,Sfrp5,Shank3,Sox4,Spint1,Spint2,St14,Sufu,Tead2,Trim71,Vegfc,Wnt6,Zeb2","_rn_":"28"},{"1":"retinoic acid metabolic process","2":"16","3":"8","4":"2.49","5":"1.34e-03","6":"Aldh1a2,Aldh8a1,Bco2,Crabp1,Cyp26b1,Cyp2s1,Rbp1,Rdh10","_rn_":"30"},{"1":"endochondral bone morphogenesis","2":"56","3":"18","4":"8.70","5":"1.47e-03","6":"Alpl,Axin2,Col27a1,Dlx5,Ext2,Ihh,Inppl1,Mmp14,Mmp16,Nab1,Rara,Rarg,Runx2,Serpinh1,Smpd3,Thbs3,Trpv4,Zmpste24","_rn_":"32"},{"1":"embryonic skeletal system morphogenesis","2":"103","3":"28","4":"16.01","5":"1.69e-03","6":"Asxl2,Bmi1,Dlg1,Fgf8,Fgfr2,Grhl2,Hoxa1,Hoxa2,Hoxa4,Hoxa5,Hoxb1,Hoxb2,Hoxb3,Hoxb4,Hoxb5,Hoxb6,Hoxb8,Hoxc4,Lhx1,Mmp14,Mmp16,Nog,Osr1,Pcgf2,Pdgfra,Rdh10,Runx2,Setd2","_rn_":"36"},{"1":"response to retinoic acid","2":"75","3":"22","4":"11.66","5":"1.80e-03","6":"Abca1,Aldh1a2,Cyp26b1,Dnaaf2,Gata4,Hoxa1,Hoxa2,Hsd17b2,Map7,Osr1,Ptk2b,Rara,Rarg,Rbp4,Ret,Rxra,Rxrg,Scamp3,Sox2,Synj1,Tead1,Tead2","_rn_":"40"},{"1":"regulation of cell migration","2":"802","3":"153","4":"124.65","5":"3.17e-03","6":"Acta2,Actn4,Adamts1,Adgra2,Adora2b,Adra2a,Ajuba,Akap12,Amot,Amotl2,Apc,Apex1,Apoe,App,Arf6,Arhgdia,Atoh8,Bag4,Bcar1,Bcr,Bex4,Bmp2,Braf,Bst1,Camk2b,Camk2d,Card10,Carmil1,Cav1,Ccbe1,Ccl5,Cdh1,Cdh13,Ceacam1,Cited2,Cldn1,Cln3,Cmklr1,Ctnna2,Dab2,Dab2ip,Dag1,Dapk3,Ddt,Ddx58,Diaph1,Dlc1,Dnaja4,Dock4,Dpp4,Efemp1,Emc10,Emp2,Epb41l5,Epha1,Evl,F2rl1,F7,Fbln1,Fga,Fgf15,Fgf8,Flna,Fut4,Fut9,Gab1,Gcnt2,Gli1,Gpsm3,Grb7,Grem1,Hdac4,Hnf4a,Hspa5,Igf2,Iqgap1,Iqsec1,Irs2,Itga5,Jag1,Jam3,Jcad,Lama5,Lgmn,Lrp1,Map2k1,Mapk3,Mien1,Mmp14,Mmp2,Myadm,Myc,Myo1c,Nbl1,Nog,Nrg1,Ntng1,Olfm1,Padi2,Pawr,Pbld1,Pdgfra,Pgr,Phlda2,Pik3cb,Pik3r2,Pkp2,Pla2g7,Plcg2,Plxna1,Plxna3,Plxnc1,Prkca,Prkce,Prr5,Pten,Ptk2b,Ptpn23,Ptprz1,Pycard,Rab25,Rcc2,Reln,Ret,Rffl,Ripk3,Ror2,Rras2,Rreb1,Sema3e,Sema4c,Sema4d,Sema5b,Sfrp1,Sh3bp1,Smim22,Smpd3,Sox14,Srcin1,Srpx2,Stard13,Tacstd2,Tbx5,Tjp1,Tmsb10,Tradd,Trib1,Trpv4,Unc5d,Vegfc,Vil1,Wnt5b,Zswim6","_rn_":"53"},{"1":"BMP signaling pathway","2":"133","3":"33","4":"20.67","5":"3.55e-03","6":"Acvr2b,Bmp2,Cav1,Dand5,Dlx5,Etv2,Fam83g,Foxd1,Fst,Gata4,Gdf6,Gpc3,Grem1,Hes5,Htra1,Lrp2,Mapk3,Msx1,Nbl1,Nog,Pelo,Rgma,Rnf165,Ror2,Runx2,Sfrp1,Sfrp5,Slc39a5,Smad1,Smad9,Smpd3,Tcf7l2,Ube2o","_rn_":"60"},{"1":"SMAD protein signal transduction","2":"66","3":"19","4":"10.26","5":"4.49e-03","6":"Afp,Atoh8,Bmp2,Btbd11,Dab2,Gata4,Gdf6,Hnf4a,Jade2,Nceh1,Pbld1,Rbm14,Rbpms,Ror2,Smad1,Smad9,Sptbn1,Trf,Veph1","_rn_":"79"},{"1":"ossification involved in bone maturation","2":"19","3":"8","4":"2.95","5":"5.08e-03","6":"Asxl2,Bmp2,Fat4,Grem1,Rflna,Sema4d,Snx10,Thbs3","_rn_":"86"},{"1":"epithelial cell morphogenesis","2":"44","3":"14","4":"6.84","5":"5.29e-03","6":"Arhgef26,Cdh1,Cgn,Dact2,Epb41l5,Grhl2,Ihh,Jmjd1c,Pof1b,Rab25,Scrib,Sipa1l3,Spint2,St14","_rn_":"87"}],"options":{"columns":{"min":{},"max":[10]},"rows":{"min":[10],"max":[10]},"pages":{}}}
  </script>
</div>

We can also compute a dendrogram that captures the similarity between the top 50 enriched terms, to identify groups that represent different gene sets.


```r
## GO term relationship
geneTon <- shake_topGOtableResult(GO.cervical)
gs_dendro(geneTon,
          n_gs = 50,
          gs_dist_type = "kappa", 
          clust_method = "ward.D2",
          color_leaves_by = "gs_pvalue",
          color_branches_by = "clusters",
          create_plot = TRUE)
```

![](02_somitogenesisDevelopment_files/figure-html/cervical_dend-1.png)<!-- -->


##### Thoracic

There are 1801 genes with highest expression in thoracic somites (stages 18, 21, and 25): the next two clusters in the heatmap.

These genes are highly enriched for relevant terms, despite being a much smaller set. We observe many of the terms identified when using the complete set of DE genes. We find many terms related to the development of the musculo-skeletal system and somitogenesis/patterning.


```r
## thoracic genes
GO.thoracic <- topGOtable(DEgenes = degs[degs$fate=="thoracic",]$gene,
                     BGgenes = universe$gene_name,
                     topGO_method2 = "elim",
                     ontology = "BP",
                     geneID = "symbol",
                     addGeneToTerms = TRUE,
                     mapping = "org.Mm.eg.db",
                     topTablerows = 1e3)
## significant terms
GO.thoracic[c(1,2,5,7,11,12,14,16:18,21,25:27,30,39,45,47,54,57,70,78,79,88,91,99), c(2:5,7,9)]
```

<div data-pagedtable="false">
  <script data-pagedtable-source type="application/json">
{"columns":[{"label":[""],"name":["_rn_"],"type":[""],"align":["left"]},{"label":["Term"],"name":[1],"type":["chr"],"align":["left"]},{"label":["Annotated"],"name":[2],"type":["int"],"align":["right"]},{"label":["Significant"],"name":[3],"type":["int"],"align":["right"]},{"label":["Expected"],"name":[4],"type":["dbl"],"align":["right"]},{"label":["p.value_elim"],"name":[5],"type":["dbl"],"align":["right"]},{"label":["genes"],"name":[6],"type":["chr"],"align":["left"]}],"data":[{"1":"cholesterol biosynthetic process","2":"44","3":"17","4":"3.45","5":"1.30e-08","6":"Cyp51,Dhcr24,Dhcr7,Fdft1,Fdps,Hmgcr,Hmgcs1,Insig1,Lss,Msmo1,Mvd,Mvk,Nsdhl,Pmvk,Sc5d,Srebf1,Srebf2","_rn_":"1"},{"1":"somitogenesis","2":"74","3":"18","4":"5.81","5":"1.20e-05","6":"Abi1,Dll1,Dmrt2,Foxc1,Foxc2,Lef1,Meox1,Meox2,Myf5,Myf6,Nkx3-1,Notch1,Pax3,Ripply1,Sema3c,Smad3,Tbx18,Tcf15","_rn_":"2"},{"1":"negative regulation of canonical Wnt signaling pathway","2":"121","3":"22","4":"9.50","5":"1.60e-04","6":"Amer1,Ccdc88c,Cdh2,Dact1,Ddit3,Frzb,Fzd4,Invs,Lrp4,Lzts2,Mcc,Mllt3,Nkd1,Notch1,Nphp3,Sfrp4,Slc9a3r1,Tbx18,Tle1,Tle4,Tmem131l,Tmem170b","_rn_":"5"},{"1":"cartilage development","2":"183","3":"29","4":"14.37","5":"2.10e-04","6":"Bbs2,Bmp3,Carm1,Ctnnb1,Esrra,Fgf18,Fgf2,Fgf9,Fgfr1,Frzb,Gnas,Hoxa3,Hoxd3,Lnpk,Matn3,Mboat2,Mycn,Myf5,Prrx2,Six2,Smad3,Snx19,Sox6,Srf,Stc1,Sulf1,Sulf2,Thbs1,Uncx","_rn_":"7"},{"1":"cell-matrix adhesion","2":"194","3":"30","4":"15.23","5":"2.50e-04","6":"Cask,Ccl25,Cd36,Col13a1,Ctnnb1,Eda,Efna5,Epdr1,Fn1,Hoxa7,Hoxd3,Itgb5,Macf1,Myf5,Nf1,Nid2,Nrp1,Pcsk5,Ptprk,Rin2,Smad3,Srf,Tesk2,Thbs1,Tiam1,Trpm7,Vcl,Vtn,Vwa2,Zyx","_rn_":"11"},{"1":"regulation of cellular response to growth factor stimulus","2":"254","3":"40","4":"19.94","5":"3.70e-04","6":"Bmper,Cask,Ctnnb1,Dll1,Fbn2,Fgf18,Fgf2,Fgf9,Fgfr1,Fgfrl1,Fstl1,Fzd4,Hes1,Hhip,Htra3,Il17rd,Lox,Ltbp1,Myocd,Notch1,Notch2,Rasl11b,Sema6a,Sfrp4,Shisa2,Skil,Smad3,Snx25,Sorl1,Spry1,Sulf1,Sulf2,Tgfb3,Thbs1,Tspan12,Vasn,Vegfb,Wasf1,Zdhhc17,Zfp423","_rn_":"12"},{"1":"embryonic skeletal system morphogenesis","2":"103","3":"22","4":"8.09","5":"6.40e-04","6":"Ctnnb1,Eya1,Foxc2,Gnas,Hoxa3,Hoxa6,Hoxa7,Hoxb7,Hoxb9,Hoxc9,Hoxd3,Hoxd4,Mycn,Myf5,Prrx2,Six1,Six2,Six4,Smad3,Tgfb3,Twist1,Twist2","_rn_":"14"},{"1":"bone development","2":"206","3":"30","4":"16.17","5":"7.20e-04","6":"Abi1,Akap13,Amer1,Ano6,Cadm1,Carm1,Col13a1,Fgf18,Fli1,Foxc1,Foxn3,Gnas,Has2,Insig1,Lox,Map2k6,Mbtps2,Meis1,Mpig6b,Notch2,Ppib,Ripply1,Sfrp4,Slc9b2,Srf,Stc1,Sulf1,Sulf2,Thbs1,Twist1","_rn_":"16"},{"1":"skin development","2":"217","3":"31","4":"17.04","5":"8.20e-04","6":"Ash1l,Col5a1,Ctnnb1,Dhcr24,Dll1,Eda,Edaradd,Etv4,Evpl,Flnb,Foxc1,Gnas,Hoxa7,Krt36,Lor,Lrp4,Mysm1,Nf1,Nme2,Notch1,Nsdhl,Opn3,Palld,Scd1,Slitrk5,Srf,Tcf15,Tcf7l1,Tnfrsf19,Wdr48,Zdhhc21","_rn_":"17"},{"1":"regulation of Wnt signaling pathway","2":"294","3":"48","4":"23.08","5":"8.40e-04","6":"Amer1,Aspm,Ccdc88c,Cdh2,Cxxc4,Dact1,Ddit3,Depdc1b,Eda,Fgf2,Fgf9,Foxl1,Frzb,Fzd4,Fzd7,Hm629797,Hmga2,Invs,Kank1,Lef1,Lrp4,Lrrk2,Lzts2,Macf1,Mcc,Mdfic,Mllt3,Nfatc4,Nkd1,Notch1,Nphp3,Nxn,Plpp3,Sall1,Sfrp4,Shisa2,Slc9a3r1,Smad3,Sulf1,Sulf2,Tbl1xr1,Tbx18,Tiam1,Tle1,Tle4,Tmem131l,Tmem170b,Tmem237","_rn_":"18"},{"1":"skeletal muscle tissue development","2":"163","3":"25","4":"12.80","5":"9.10e-04","6":"Asb2,Btg2,Ccnt2,Ctnnb1,Dll1,Eln,Flnb,Foxp1,Hdac5,Heyl,Hmgcr,Meox2,Myf5,Myf6,Myl6b,Myocd,Nf1,Notch1,Nupr1,Pax3,Six1,Six4,Skil,Svil,Twist1","_rn_":"21"},{"1":"negative regulation of epithelial cell differentiation","2":"41","3":"10","4":"3.22","5":"1.00e-03","6":"Cited1,Ctnnb1,Dll1,Foxp1,Frzb,Hes1,Hoxa7,Notch1,Sall1,Spry1","_rn_":"25"},{"1":"anterior/posterior pattern specification","2":"222","3":"41","4":"17.43","5":"1.03e-03","6":"Abi1,Bhlhe41,Btg2,Ctnnb1,Ddit3,Dll1,Dmrt2,Foxc1,Foxc2,Hes1,Hey1,Heyl,Hoxa3,Hoxa6,Hoxa7,Hoxb7,Hoxb9,Hoxc6,Hoxc8,Hoxc9,Hoxd3,Hoxd4,Lef1,Meox1,Meox2,Mllt3,Myf5,Myf6,Nkx3-1,Notch1,Pax3,Pcsk5,Ripply1,Sema3c,Six2,Smad3,Srf,Tbx18,Tcf15,Tcf7l1,Wnt2b","_rn_":"26"},{"1":"embryonic pattern specification","2":"63","3":"13","4":"4.95","5":"1.03e-03","6":"Cited1,Ctnnb1,Dll1,Meis1,Meox1,Meox2,Nrp1,Ooep,Ripply1,Sema3f,Smad3,Stil,Tcf7l1","_rn_":"27"},{"1":"mesenchymal cell development","2":"88","3":"16","4":"6.91","5":"1.23e-03","6":"Acvr1,Cdh2,Ednra,Fn1,Foxc1,Foxc2,Hes1,Hey1,Heyl,Notch1,Nrp1,Pax3,Sema3c,Sema3f,Sema6a,Twist1","_rn_":"30"},{"1":"negative regulation of fibroblast growth factor receptor signaling pathway","2":"18","3":"6","4":"1.41","5":"1.87e-03","6":"Fgfrl1,Shisa2,Spry1,Sulf1,Sulf2,Thbs1","_rn_":"39"},{"1":"muscle cell development","2":"177","3":"25","4":"13.90","5":"2.92e-03","6":"Adra1a,Akap13,Akap6,Csrp2,Fdps,Fhl2,Fhod3,Foxp1,Hdac5,Hes1,Lox,Myf5,Myf6,Myocd,Nfatc4,Nkx2-6,Notch1,Pdgfrb,Sgcb,Six1,Six4,Srf,Tbx18,Tmtc3,Wdr1","_rn_":"45"},{"1":"regulation of myoblast differentiation","2":"47","3":"10","4":"3.69","5":"3.02e-03","6":"Boc,Ddit3,Dll1,Eid2b,Mkx,Myf5,Myf6,Myocd,Notch1,Nr2c2","_rn_":"47"},{"1":"axis specification","2":"89","3":"15","4":"6.99","5":"3.70e-03","6":"Cited1,Ctnnb1,Ddit3,Dll1,Hey1,Hspb11,Notch1,Notch2,Ripply1,Six2,Srf,Stc1,Stil,Tbx18,Tcf7l1","_rn_":"54"},{"1":"negative regulation of striated muscle cell differentiation","2":"34","3":"8","4":"2.67","5":"4.03e-03","6":"Bhlhe41,Dll1,Foxp1,Fzd7,Hdac5,Myocd,Notch1,Sox6","_rn_":"57"},{"1":"negative regulation of angiogenesis","2":"92","3":"15","4":"7.22","5":"5.10e-03","6":"Adgrb3,Cd36,Col4a2,Ctnnb1,Efna3,Foxc1,Hhip,Ism1,Krit1,Nf1,Pde3b,Sema6a,Slc12a2,Sulf1,Thbs1","_rn_":"70"},{"1":"skeletal system morphogenesis","2":"234","3":"41","4":"18.37","5":"5.79e-03","6":"Carm1,Col13a1,Ctnnb1,Eya1,Fbn2,Fgf18,Foxc1,Foxc2,Foxn3,Gnas,Has2,Hhip,Hoxa3,Hoxa6,Hoxa7,Hoxb7,Hoxb9,Hoxc8,Hoxc9,Hoxd3,Hoxd4,Insig1,Mycn,Myf5,Pdgfrb,Prrx2,Ripply1,Sfrp4,Six1,Six2,Six4,Smad3,Sox6,Stc1,Tcf15,Tgfb3,Thbs1,Twist1,Twist2,Uncx,Wdr48","_rn_":"78"},{"1":"transforming growth factor beta receptor signaling pathway","2":"158","3":"22","4":"12.40","5":"5.97e-03","6":"Acvr1,Cited1,Fbn2,Fut8,Htra3,Il17rd,Itgb5,Lox,Lrrc32,Ltbp1,Ltbp4,Myocd,Ptprk,Rasl11b,Skil,Smad3,Snx25,Spry1,Tgfb3,Thbs1,Vasn,Zyx","_rn_":"79"},{"1":"positive regulation of mesenchymal cell proliferation","2":"37","3":"8","4":"2.90","5":"6.95e-03","6":"Ctnnb1,Fgf9,Fgfr1,Foxp1,Mycn,Prrx2,Six1,Tbx18","_rn_":"88"},{"1":"positive regulation of epithelial cell proliferation","2":"180","3":"24","4":"14.13","5":"7.25e-03","6":"Akt3,Arg1,Ctnnb1,Eya1,Fgf18,Fgf2,Fgf9,Fgfr1,Foxp1,Fzd7,Has2,Hmga2,Iqgap3,Itpr1,Nf1,Nme2,Notch1,Notch2,Prkd1,Tbx18,Twist1,Twist2,Vegfb,Wdr48","_rn_":"91"},{"1":"osteoblast differentiation","2":"181","3":"24","4":"14.21","5":"7.76e-03","6":"Acvr1,Bmp3,Cited1,Ctnnb1,Esrra,Fbn2,Fgf2,Fgf9,Fgfr1,Fhl2,Gabbr1,Gnas,Hdac5,Hey1,Lox,Map2k6,Nf1,Notch1,Pdlim7,Prkd1,Smad3,Smoc1,Twist1,Twist2","_rn_":"99"}],"options":{"columns":{"min":{},"max":[10]},"rows":{"min":[10],"max":[10]},"pages":{}}}
  </script>
</div>

```r
# cholesterol biosynthesis is specific to stage 25
```



```r
## GO term relationship
geneTon <- shake_topGOtableResult(GO.thoracic)
gs_dendro(geneTon,
          n_gs = 50,
          gs_dist_type = "kappa", 
          clust_method = "ward.D2",
          color_leaves_by = "gs_pvalue",
          color_branches_by = "clusters",
          create_plot = TRUE)
```

![](02_somitogenesisDevelopment_files/figure-html/thoracic_dend-1.png)<!-- -->

##### Lumbar

There are 1769 genes with highest expression in lumbar somites (stage 27): the next two clusters in the heatmap.

In this case, terms related to erythrocyte development, angiogenesis and vasculogenesis are prevalent, together with differentiation and cell migration.


```r
## lumbar genes
GO.lumbar <- topGOtable(DEgenes = degs[degs$fate=="lumbar",]$gene,
                     BGgenes = universe$gene_name,
                     topGO_method2 = "elim",
                     ontology = "BP",
                     geneID = "symbol",
                     addGeneToTerms = TRUE,
                     mapping = "org.Mm.eg.db",
                     topTablerows = 1e3)
## significant terms
GO.lumbar[c(1,4,5,8,21,22,25,29,34,36:38,51,61), c(2:5,7,9)]
```

<div data-pagedtable="false">
  <script data-pagedtable-source type="application/json">
{"columns":[{"label":[""],"name":["_rn_"],"type":[""],"align":["left"]},{"label":["Term"],"name":[1],"type":["chr"],"align":["left"]},{"label":["Annotated"],"name":[2],"type":["int"],"align":["right"]},{"label":["Significant"],"name":[3],"type":["int"],"align":["right"]},{"label":["Expected"],"name":[4],"type":["dbl"],"align":["right"]},{"label":["p.value_elim"],"name":[5],"type":["dbl"],"align":["right"]},{"label":["genes"],"name":[6],"type":["chr"],"align":["left"]}],"data":[{"1":"erythrocyte development","2":"39","3":"20","4":"4.13","5":"7.90e-08","6":"Alas1,Alas2,Ank1,Bpgm,Dmtn,Epb42,Fam210b,Flvcr1,Gata1,Hba-a1,Hba-a2,Hba-x,Klf1,Ptbp3,Rhag,Rhd,Slc11a2,Slc4a1,Tal1,Tmod3","_rn_":"1"},{"1":"negative regulation of mesenchymal cell apoptotic process","2":"11","3":"7","4":"1.17","5":"3.30e-05","6":"Bmp7,Pax2,Pax8,Pou3f4,Shh,Sox9,Wt1","_rn_":"4"},{"1":"cell fate determination","2":"36","3":"13","4":"3.81","5":"4.40e-05","6":"Ascl1,Bmp4,Chrdl1,Gata2,Gata3,Gata6,Isl1,Mesp1,Pax2,Pax6,Sox17,Wnt1,Wnt7a","_rn_":"5"},{"1":"embryonic hindlimb morphogenesis","2":"35","3":"12","4":"3.71","5":"1.60e-04","6":"Alx3,Bmp4,Mecom,Osr2,Pitx1,Pitx2,Rarb,Rspo2,Shh,Tbx4,Wnt7a,Zbtb16","_rn_":"8"},{"1":"angiogenesis","2":"467","3":"79","4":"49.49","5":"1.01e-03","6":"Adam12,Adgrg1,Adra2b,Angpt4,Angptl6,Anxa1,Anxa3,Apela,Apln,Atp2b4,Bmp4,Ccn1,Cd40,Cd59a,Cela1,Chrna7,Cib1,Clec14a,Col8a2,Cxcl12,Cyp1b1,Dcn,E2f2,Emilin1,Ets1,Flt4,Gadd45a,Gata2,Gata6,Gpx1,Hand2,Hbegf,Hgf,Hhex,Hif3a,Htatip2,Id1,Isl1,Itgb8,Jak1,Jmjd8,Jun,Jup,Kctd10,Kdr,Lgals3,Mfge8,Mmp9,Ndnf,Nrxn3,Pdcd10,Pdgfa,Pitx2,Pknox1,Plau,Rasip1,Rhob,Sema4a,Serpinf1,Sfrp2,Shh,Sirt6,Sox17,Sparc,Stat3,Tal1,Tbx20,Tbx4,Tbxa2r,Tcf21,Tek,Tgfa,Thsd7a,Tie1,Tmem100,Tnfaip2,Unc5b,Wnt7a,Wnt7b","_rn_":"21"},{"1":"positive regulation of osteoblast differentiation","2":"66","3":"16","4":"6.99","5":"1.16e-03","6":"Acvr2a,Bmp4,Bmp6,Bmp7,Ccn1,Ccn4,Cebpa,Cebpb,Ddr2,Nell1,Npnt,Nppc,Scube2,Sfrp2,Tmem119,Wnt7b","_rn_":"22"},{"1":"regulation of morphogenesis of an epithelium","2":"61","3":"15","4":"6.46","5":"1.40e-03","6":"Ar,Bmp4,Bmp7,Fgf10,Gata3,Gdnf,Hgf,Ntn4,Pax2,Pax8,Pdgfa,Shh,Sirt6,Sox9,Wnt2","_rn_":"25"},{"1":"positive regulation of epithelial cell migration","2":"142","3":"27","4":"15.05","5":"1.83e-03","6":"Angpt4,Anxa1,Anxa3,Bmp4,Capn7,Cd40,Cib1,Ets1,Fgf10,Flt4,Gata2,Gata3,Glipr2,Hbegf,Itga3,Jun,Kdr,Map2k3,Mmp9,Rhob,Shh,Sox9,Sparc,Stat5a,Tac1,Tek,Wnt7a","_rn_":"29"},{"1":"regulation of Notch signaling pathway","2":"83","3":"18","4":"8.80","5":"2.31e-03","6":"Adam10,Arrdc1,Ascl1,Bmp2k,Bmp7,Dlk1,Dlk2,Fgf10,Gata2,Gdpd5,Kctd10,Mesp1,Mfng,Pdcd10,Postn,Robo2,Stat3,Wnt1","_rn_":"34"},{"1":"cartilage development involved in endochondral bone morphogenesis","2":"35","3":"10","4":"3.71","5":"2.64e-03","6":"Anxa6,Col1a1,Col2a1,Nppc,Poc1a,Rarb,Scube2,Shox2,Sik3,Sox9","_rn_":"36"},{"1":"positive regulation of angiogenesis","2":"146","3":"27","4":"15.47","5":"2.76e-03","6":"Adam12,Angpt4,Anxa3,Cd40,Cela1,Chrna7,Cyp1b1,Ets1,Gata2,Gata6,Hgf,Isl1,Itgb8,Jak1,Jmjd8,Jup,Kdr,Lgals3,Mmp9,Rhob,Sfrp2,Shh,Sirt6,Stat3,Tbxa2r,Tek,Tie1","_rn_":"37"},{"1":"vasculogenesis","2":"91","3":"19","4":"9.64","5":"2.80e-03","6":"Apela,Asb4,Epor,Foxf1,Glmn,Hhex,Itgb8,Kdr,Pitx2,Rasip1,Shh,Sox17,Tbx20,Tek,Tie1,Tmem100,Wnt7a,Wnt7b,Wt1","_rn_":"38"},{"1":"striated muscle cell differentiation","2":"270","3":"43","4":"28.61","5":"4.17e-03","6":"Alpk3,Bmp4,Bvht,Cacnb4,Ccn4,Col14a1,Csrp1,Csrp3,Cxcl12,Ehd2,Flot1,Gata6,Gpx1,Hand2,Hdac3,Igfbp5,Isl1,Kel,Krt19,Mesp1,Myl2,Myof,Nln,Ntn3,Pitx2,Pmp22,Popdc2,Rarb,Rbm38,Rgs2,Rgs4,Shh,Shox2,Sorbs2,Tmem119,Tmod1,Tmod2,Tmod3,Tnnt2,Tnnt3,Wnt1,Wt1,Xk","_rn_":"51"},{"1":"axis specification","2":"89","3":"18","4":"9.43","5":"5.09e-03","6":"Bmp4,Fgf10,Foxa2,Hhex,Lefty1,Mesp1,Mesp2,Pax6,Pcsk6,Peg12,Pitx2,Ripply2,Shh,Ttc8,Vax2,Wnt1,Wnt7a,Wt1","_rn_":"61"}],"options":{"columns":{"min":{},"max":[10]},"rows":{"min":[10],"max":[10]},"pages":{}}}
  </script>
</div>


```r
## GO term relationship
geneTon <- shake_topGOtableResult(GO.lumbar)
gs_dendro(geneTon,
          n_gs = 50,
          gs_dist_type = "kappa", 
          clust_method = "ward.D2",
          color_leaves_by = "gs_pvalue",
          color_branches_by = "clusters",
          create_plot = TRUE)
```

![](02_somitogenesisDevelopment_files/figure-html/lumbar_dend-1.png)<!-- -->


##### Sacral

Finally, there are 3020 genes with highest expression in sacral somites (stage 35): the next three clusters in the heatmap.

We still see an enrichment for angiogenesis and blood vessel development. Immune-related terms are also common.


```r
## sacral genes
GO.sacral <- topGOtable(DEgenes = degs[degs$fate=="sacral",]$gene,
                     BGgenes = universe$gene_name,
                     topGO_method2 = "elim",
                     ontology = "BP",
                     geneID = "symbol",
                     addGeneToTerms = TRUE,
                     mapping = "org.Mm.eg.db",
                     topTablerows = 1e3)
## significant terms
GO.sacral[c(1,4,6,7,14,31,33,45,47,70), c(2:5,7,9)]
```

<div data-pagedtable="false">
  <script data-pagedtable-source type="application/json">
{"columns":[{"label":[""],"name":["_rn_"],"type":[""],"align":["left"]},{"label":["Term"],"name":[1],"type":["chr"],"align":["left"]},{"label":["Annotated"],"name":[2],"type":["int"],"align":["right"]},{"label":["Significant"],"name":[3],"type":["int"],"align":["right"]},{"label":["Expected"],"name":[4],"type":["dbl"],"align":["right"]},{"label":["p.value_elim"],"name":[5],"type":["dbl"],"align":["right"]},{"label":["genes"],"name":[6],"type":["chr"],"align":["left"]}],"data":[{"1":"inflammatory response","2":"475","3":"129","4":"81.98","5":"1.50e-07","6":"Acer3,Adamts12,Adcy7,Adora2a,Ager,Aldh2,Aoah,B4galt1,Bcl6b,Btk,C3ar1,C5ar1,Calcrl,Casp1,Casp4,Casp6,Ccl24,Ccl3,Ccl6,Ccl9,Cd200r1,Cd68,Cdh5,Cdk19,Chst1,Csf1r,Cspg4,Ctla2a,Ctsc,Ctss,Cuedc2,Cx3cl1,Cx3cr1,Cyba,Cybb,Dpep1,Ecm1,Egfr,Epha2,F8,Fcer1g,Fcgr1,Fcgr3,Ffar4,Fndc4,Gja1,Grn,Gsdmd,H2-T23,Hdac7,Hmgb1,Hnrnpa0,Hp,Hyal3,Icam1,Igf1,Il17d,Il17rc,Il18rap,Il1rl2,Il4ra,Irak2,Irf3,Irf5,Itga2,Itgav,Itgb1,Itgb2,Kit,Lat,Lgals9,Loxl3,Lpl,Lxn,Ly86,Ly96,Mapk8,Mgll,Naip2,Ncf1,Ndfip1,Ndst1,Ndufc2,Ndufs4,Nfkbid,Nlrc3,Nmi,Nrros,Parp4,Pf4,Pik3ap1,Pik3cg,Pla2g4a,Pld4,Plgrkt,Ppbp,Prcp,Prdx2,Prkcq,Ptger3,Ptgfr,Ptgs1,Ptn,Rhbdd3,S1pr3,Scnn1b,Serpinb1a,Slit2,Smpdl3b,Socs5,Stab1,Syk,Syt11,Tcirg1,Tgfb1,Themis2,Timp1,Tlr4,Tlr7,Tnfaip6,Tnfaip8l2,Tnfrsf11a,Tnfrsf1a,Trem2,Tril,Tyro3,Usp18,Wnt5a,Zfp580","_rn_":"1"},{"1":"positive regulation of angiogenesis","2":"146","3":"44","4":"25.20","5":"8.40e-05","6":"Acvrl1,Adm,Angpt2,Aqp1,Btg1,C3ar1,C5ar1,Ccl24,Cd34,Cdh5,Ctsh,Cx3cl1,Cx3cr1,Cybb,Ecm1,Eng,Erap1,Flt1,Grn,Hmgb1,Hspb1,Hspb6,Itgb1,Itgb2,Itgb3,Klf4,Nos3,Notch4,Pgf,Pik3r6,Prkd2,Ramp2,Rapgef3,Rras,Runx1,Sash1,Sema5a,Smoc2,Stim1,Tert,Tgfbr2,Tnfrsf1a,Vegfd,Wnt5a","_rn_":"4"},{"1":"angiogenesis","2":"467","3":"125","4":"80.60","5":"1.30e-04","6":"Ackr3,Acvrl1,Adam15,Adm,Ahr,Angpt1,Angpt2,Angptl1,Angptl2,Angptl4,Apold1,Aqp1,Atp5b,B4galt1,Btg1,C3ar1,C5ar1,Calcrl,Casp8,Ccl24,Ccn2,Cd34,Cdh5,Cemip2,Cldn5,Col15a1,Col8a1,Cspg4,Ctsh,Cx3cl1,Cx3cr1,Cybb,Dll4,Dysf,Ecm1,Ecscr,Efna1,Egfl7,Eng,Enpp2,Epas1,Epha2,Ephb2,Ephb4,Erap1,Esm1,Fbln5,Fgfbp1,Flt1,Foxj2,Grn,H2-M3,Hdac7,Hmgb1,Hpse,Hspb1,Hspb6,Hspg2,Itgav,Itgb1,Itgb2,Itgb3,Klf2,Klf4,Lemd3,Mcam,Mmp19,Mmrn2,Naxe,Nos3,Notch3,Notch4,Nrxn1,Parva,Pecam1,Pf4,Pgf,Pik3c2a,Pik3cg,Pik3r3,Pik3r6,Plk2,Plxnd1,Ppp1r16b,Prcp,Prkd2,Prkx,Ptn,Ptprb,Ptprm,Ramp1,Ramp2,Ramp3,Rapgef3,Rgcc,Rhoj,Rnh1,Robo4,Rras,Rspo3,Runx1,S1pr1,Sash1,Sema5a,Slit2,Smad5,Smoc2,Sox18,Stab1,Stat1,Stim1,Syk,Tbx1,Tcf4,Tert,Tgfb2,Tgfbr1,Tgfbr2,Tgfbr3,Thbs2,Tnfrsf1a,Vav3,Vegfd,Vhl,Wnt5a","_rn_":"6"},{"1":"cell adhesion","2":"1073","3":"271","4":"185.18","5":"1.30e-04","6":"Abat,Ackr3,Acvrl1,Adam15,Adam2,Adam23,Adamts12,Adamts18,Adgre5,Adgrv1,Adora2a,Afdn,Ager,Aif1,Akip1,Alox12,Angpt1,Angpt2,Antxr1,Anxa9,Apbb1ip,Aplp1,Arhgap6,Ass1,Atp5b,B4galt1,Bcl2,Bcl2l11,Bmx,Btn1a1,Ccm2l,Ccn2,Cd200r1,Cd33,Cd34,Cd59b,Cd63,Cd83,Cd9,Cd93,Cd99l2,Cdh11,Cdh12,Cdh15,Cdh5,Cdk5,Cdk5r1,Cdk6,Cdon,Chl1,Chrd,Clca2,Cldn3,Cldn34c1,Cldn4,Cldn5,Cldn8,Cntn2,Col15a1,Col6a2,Col8a1,Coro1a,Csf3r,Ctnna3,Ctnnd2,Cx3cl1,Cx3cr1,Cyth3,Cytip,Dchs1,Dlg5,Dmd,Dusp22,Dusp3,Dysf,Edil3,Efemp2,Efna1,Efnb3,Egfl6,Egfl7,Egfr,Eng,Enpp2,Epcam,Epha2,Epha3,Epha4,Epha7,Ephb4,Erf,Esam,Fbln2,Fbln5,Fbln7,Fermt1,Fermt3,Fes,Fgl2,Flrt2,Fndc3b,Gas6,Gp5,Gp9,Gpc4,Gpm6b,Gsk3b,H2-M3,H2-T23,Hfe,Hmcn1,Hpse,Hsd17b12,Hspb1,Hyal3,Icam1,Icam2,Igf1,Igfbp2,Igsf21,Il15,Il1rl2,Il4ra,Irf1,Itga1,Itga2,Itga4,Itgal,Itgav,Itgb1,Itgb2,Itgb3,Itpkb,Jag2,Jam2,Kirrel3,Kit,Kitl,Klf4,L1cam,Lama1,Lama3,Lama4,Lamc2,Laptm5,Lgals1,Lgals8,Lgals9,Limch1,Lmo7,Loxl3,Lrrc4,Mad2l2,Madcam1,Magi1,Magi2,Mcam,Megf11,Megf9,Mertk,Mfap4,Milr1,Mpzl2,Myo10,Myo1g,Mypn,Ncan,Nckap1l,Ndfip1,Nectin2,Nectin4,Nexmif,Nfkbid,Ninj1,Nlgn3,Notch4,Nrxn1,Ntm,Olfm4,Omd,Onecut2,Pard3b,Parva,Parvb,Parvg,Pcdh1,Pcdh10,Pcdh12,Pcdh17,Pdgfb,Pdpn,Pecam1,Perp,Pgm5,Phldb2,Pik3r6,Plek,Plxna2,Plxnd1,Podxl2,Ppm1f,Ppp1r12a,Ppp3ca,Prdx2,Prkar1a,Prkcq,Prkd2,Prkg1,Prkx,Ptger3,Ptk7,Ptn,Ptprc,Ptprm,Ptpru,Rac2,Rac3,Rasal3,Rgcc,Rhod,Rras,Runx1,Runx3,S100a10,S1pr1,Sash3,Scarf1,Selenok,Selplg,Sema5a,Serpine2,Serpini1,Skap1,Slc7a11,Slit2,Smad7,Smoc2,Snai2,Sned1,Socs5,Sorbs3,Sparcl1,Spock2,Src,St3gal4,Stab1,Syk,Tcam1,Tenm2,Tgfb1,Tgfb2,Tgfbi,Tgfbr2,Thbs2,Thsd1,Tmem8b,Tnc,Tnfaip6,Tnfaip8l2,Tnr,Twsg1,Tyro3,Vav1,Vav3,Vcam1,Whamm,Wnt3a,Wnt4,Wnt5a,Zfhx3,Zfp703","_rn_":"7"},{"1":"bone remodeling","2":"83","3":"31","4":"14.32","5":"2.50e-04","6":"Cd38,Csf1r,Cthrc1,Ctss,Efna2,Egfr,Enpp1,Epha2,Gja1,Gsk3b,Idua,Inpp4b,Inpp5d,Itgav,Itgb3,Lgr4,Lrrk1,Mc4r,Nox4,P3h4,Pdk4,Pth1r,Ptn,Rac2,S1pr1,Src,Syk,Tcirg1,Tmem64,Tnfrsf11a,Wnt16","_rn_":"14"},{"1":"blood vessel endothelial cell migration","2":"108","3":"32","4":"18.64","5":"1.02e-03","6":"Acvrl1,Alox12,Angpt1,Angpt2,Atp5b,Cdh5,Dll4,Efna1,Epha2,Ephb4,Fgfbp1,Hdac7,Hspb1,Igf1,Itgb1,Klf4,Lemd3,Mef2c,Mmrn2,Nos3,Pdgfb,Pik3c2a,Pik3r3,Plk2,Prcp,Prkd2,Rgcc,Rhoj,Scarb1,Slit2,Tgfb1,Vhl","_rn_":"31"},{"1":"blood vessel remodeling","2":"46","3":"17","4":"7.94","5":"1.11e-03","6":"Acvrl1,Ahr,Angpt2,Bmpr2,Cbs,Dll4,Epas1,Erg,Gja1,Igf1,Itga4,Mef2c,Nos3,Rspo3,Tbx1,Tgfb2,Tgfbr3","_rn_":"33"},{"1":"negative regulation of ossification","2":"36","3":"14","4":"6.21","5":"1.70e-03","6":"Bcl2,Bcor,Dkk1,Ecm1,Enpp1,Kremen1,Kremen2,Mef2c,Rflnb,Smad7,Smurf1,Srgn,Tgfb1,Trpm4","_rn_":"45"},{"1":"positive regulation of MAPK cascade","2":"406","3":"93","4":"70.07","5":"1.88e-03","6":"Ackr3,Ager,Angpt1,Ankrd6,Arhgef5,Arl6ip5,C5ar1,Cav2,Cavin3,Ccl24,Ccl3,Ccl6,Ccl9,Ccn2,Cdon,Cflar,Csf1r,Cspg4,Cx3cl1,Dixdc1,Dkk1,Drd4,Dstyk,Dusp19,Dusp22,Dvl3,Edn3,Efna1,Egfr,Epha4,Ffar4,Fgd2,Flt1,Frs2,Gas6,Gpr183,Hmgb1,Icam1,Igf1,Igfbp3,Ighm,Itga1,Itgb3,Kit,Kitl,Klhdc10,Laptm5,Lpar1,Madd,Map3k20,Map3k21,Mapk8ip1,Mapkbp1,Mdfi,Met,Ncf1,Ndst1,Nenf,Ngf,Nod1,Nox4,Nqo2,Nrxn1,Ntf3,P2ry1,P2ry6,Pde8a,Pdgfb,Pik3cg,Pik3r6,Plce1,Prdx2,Prkd2,Psap,Ptprc,Ramp3,Sash1,Sorbs3,Spi1,Src,Syk,Tab1,Tbx1,Tgfb1,Tgfb2,Timp2,Tlr4,Tnfrsf11a,Tnik,Tpd52l1,Trem2,Wnt16,Wnt5a","_rn_":"47"},{"1":"muscle cell fate commitment","2":"13","3":"7","4":"2.24","5":"2.93e-03","6":"AW551984,Dkk1,Mef2c,Tbx1,Tbx2,Tbx3,Wnt3a","_rn_":"70"}],"options":{"columns":{"min":{},"max":[10]},"rows":{"min":[10],"max":[10]},"pages":{}}}
  </script>
</div>


```r
## GO term relationship
geneTon <- shake_topGOtableResult(GO.sacral)
gs_dendro(geneTon,
          n_gs = 50,
          gs_dist_type = "kappa", 
          clust_method = "ward.D2",
          color_leaves_by = "gs_pvalue",
          color_branches_by = "clusters",
          create_plot = TRUE)
```

![](02_somitogenesisDevelopment_files/figure-html/sacral_dend-1.png)<!-- -->

### Chromatin remodelling during developmental progression

There are also a large number of chromatin loci that change across development. Once again, most change consistently between the somite trios.


```r
dars <- read.table(paste0(dir, "ATAC-seq/results/04_DAregions_stages.tsv"), 
                   stringsAsFactors = FALSE)
dars$somiteSpecific <- ifelse(rowSums(dars[,6:8])>0, 1, 0)
c(average = sum(dars$ave), stage_wise = sum(dars$somiteSpecific))
```

```
##    average stage_wise 
##      28404       4609
```

The accessibility dynamics can also be clustered


```r
## count data for DARs
data <- peakCounts[row.names(dars),]
## convert into z-scores
data <- 2^data
data <- t(apply(data, 1, function(x) (x-mean(x))/sd(x) ))

## split genes into clusters
dist <- dist(data)
hc <- hclust(dist)
# data is too large to use dynamic tree cut; cut the tree directly with k=12
cut.dars <- cutree(hc, k=12)
  
## order clusters to follow developmental progression
clusters.dars <- factor(paste0("cluster",cut.dars), 
                        levels=paste0("cluster", c(10,3, 1,4,9,2,5,12, 7,6,8, 11)))
names(clusters.dars) <- row.names(data)
table(clusters.dars)
```

```
## clusters.dars
## cluster10  cluster3  cluster1  cluster4  cluster9  cluster2  cluster5 cluster12 
##       959      3756      4612      1733      1738      6032      3565       280 
##  cluster7  cluster6  cluster8 cluster11 
##      1159      4551      4238       390
```

and the observed dynamics are more diverse than the observed gene expression. Still, most regions peak at a specific stage, but we observe greater diversity in the thoracic somites, with some regions specific to one of the three stages.


```r
## order columns based on somite number
m <- meta.atac
m$somite <- factor(m$somite, levels=c("SIII","SII","SI"))
stopifnot(identical(colnames(data), m$sample)) # make sure the metadata corresponds with count data matrix
order <- order(m$stage, m$somite)

## heatmap annotation
ha  <- HeatmapAnnotation(df = data.frame(stage = paste0("stage", m[order,]$stage), 
                                         somite = m[order,]$somite), 
                         col = list(stage = cols.stage, somite = cols.somite))

Heatmap(data[,order], 
        cluster_columns = FALSE, 
        col=colorRamp2(breaks = c(-2,0,2,4), 
                       colors = c("steelblue","white","indianred3","indianred4")),
        name = "z-score", 
        show_row_names = FALSE, 
        show_column_names = FALSE, 
        top_annotation = ha,
        row_split = clusters.dars,
        cluster_row_slices = FALSE)
```

![](02_somitogenesisDevelopment_files/figure-html/heatmap_dars-1.png)<!-- -->

The regions with highest accessibility in thoracic somites are the most abundant, accounting for over half of all DA regions (54.4%).


```r
## add cluster identity to DARs
dars$cluster <- clusters.dars

# annotate by vertebral fate
dars$fate <- ifelse(dars$cluster %in% paste0("cluster", c(10,3)), "cervical",
                    ifelse(dars$cluster %in% paste0("cluster", c(1,4,9,2,5,12)), "thoracic",
                           ifelse(dars$cluster %in% paste0("cluster", c(7,6)), "lumbar", "sacral")))
dars[dars$cluster=="cluster11",]$fate <- NA
# more fine-grained annotation
dars$fate.sub <- ifelse(dars$cluster %in% paste0("cluster", c(10,3)), "cervical",
                    ifelse(dars$cluster == "cluster1", "cervical-thoracic",
                           ifelse(dars$cluster %in% paste0("cluster", c(4,9,2,5,12)), "thoracic",
                                  ifelse(dars$cluster == "cluster7", "lumbar", 
                                         ifelse(dars$cluster == "cluster6", "lumbar-sacral", "sacral")))))
dars[dars$cluster=="cluster11",]$fate.sub <- NA
table(dars$fate.sub)[c(1,2,6,3:5)]
```

```
## 
##          cervical cervical-thoracic          thoracic            lumbar 
##              4715              4612             13348              1159 
##     lumbar-sacral            sacral 
##              4551              4238
```

In terms of the location of these DA regions in the genome, they are slightly depleted of promoters, with a corresponding increase in distal and intergenic regions. But the difference is not as striking as observed between the somite trios. However, if we split by their accessibility profile, regions open in thoracic somites show a much more pronounced depletion of promoters and increase of distal elements; whereas regions preferentially open in lumbar somites are enriched in promoter elements.


```r
## annotation of regions based on genomic context
peakAnn <- read.table(paste0(dir, "ATAC-seq/results/05_peaks_classAnnotation.tsv"), sep="\t")
dars.ann <- peakAnn[substr(row.names(dars), 4, 30),]
stopifnot(identical(row.names(dars), paste0("chr", row.names(dars.ann))))
dars.ann$cluster <- dars$cluster
dars.ann$fate <- dars$fate
dars.ann$fate.sub <- dars$fate.sub

df <- data.frame(all = round(colSums(peakAnn[,7:11])/nrow(peakAnn)*100, 2),
                 DA = round(colSums(dars.ann[,7:11])/nrow(dars.ann)*100, 2),
                 cervical = colSums(dars.ann[which(dars.ann$fate=="cervical"),7:11])/
                   length(which(dars.ann$fate=="cervical"))*100,
                 thoracic = colSums(dars.ann[which(dars.ann$fate=="thoracic"),7:11])/
                   length(which(dars.ann$fate=="thoracic"))*100,
                 lumbar = colSums(dars.ann[which(dars.ann$fate=="lumbar"),7:11])/
                   length(which(dars.ann$fate=="lumbar"))*100,
                 sacral = colSums(dars.ann[which(dars.ann$fate=="sacral"),7:11])/
                   length(which(dars.ann$fate=="sacral"))*100)
df <- data.frame(class=rep(row.names(df),6),
                 DA = c(rep(FALSE, 5), rep(TRUE, 25)),
                 context = c(rep(c("all", "DA", "cervical", "thoracic", "lumbar", "sacral"), each=5)),
                 prop = c(df$all, df$DA, df$cervical, df$thoracic, df$lumbar, df$sacral))
df$context <- factor(df$context, levels=c("all", "DA", "cervical", "thoracic", "lumbar", "sacral"))

ggplot(df, aes(context, prop, fill=class)) + 
  geom_bar(stat="identity") +
  scale_fill_manual(values = brewer.pal(n=5, "Spectral")) +
  ylab("proportion of all peaks") +
  th
```

![](02_somitogenesisDevelopment_files/figure-html/genomic_context-1.png)<!-- -->

In contrast with the observations with the somite trios, a large number of DA regions directly overlap DE genes: 12824 peaks (38.85%) overlap 4693 genes (43.89%).


```r
## find overlaps between DA peaks and DE genes
da <- GRanges(seqnames = dars$seqnames, 
              ranges = IRanges(dars$start, dars$end),
              cluster = dars$cluster,
              fate = dars$fate,
              fate.sub = dars$fate.sub)
gene_ann <- read.table(paste0(dir, "RNA-seq/data/Mus_musculus.GRCm38.96.ann"), row.names = 1)
de <- gene_ann[row.names(degs),]
de$fate <- degs[row.names(de),]$fate
de <- GRanges(seqnames = paste0("chr", de$V3),
              ranges = IRanges(de$V4, de$V5),
              strand = de$V6,
              gene = de$V2,
              fate = de$fate)
# length(unique(subjectHits(findOverlaps(de, da))))/nrow(dars)
# 12824 (38.85%)
# length(unique(queryHits(findOverlaps(de, da))))/nrow(degs)
# 4693 (43.89%)
```


#### Overrepresented biological processes {.tabset}

GREAT analysis takes the genes in the proximity of the DA regions, which are putatively regulated by these peaks, and performs enrichment analysis. 56% of the DE genes lie in close proximity to the set of DA regions. Thus, the enriched GO terms are similar to those observed based on the DE gene analysis. 


```r
## background is all open chromatin regions
background <- GRanges(seqnames = paste0("chr", peakAnn$seqnames),
                      ranges = IRanges(peakAnn$start, peakAnn$end))

## run GREAT
job = submitGreatJob(gr = da,
                     bg = background,
                     species = "mm10")

# retrieve which genes are associated with each peak
great.peak_geneAnn <- plotRegionGeneAssociationGraphs(job, plot=FALSE)
# length(intersect(degs$gene, unique(great.peak_geneAnn$gene)))
# 6023 (56.34%)  # DE genes that are the 'closest' gene to a DA peak

## enrichment results
go.dars <- getEnrichmentTables(job, ontology = "GO Biological Process")
pheno.dars <- getEnrichmentTables(job, category = "Phenotype")

## GO biological process
go.dars[[1]][c(1,6,7,13:15,40,59), c(2,6,11:9,13)]
```

<div data-pagedtable="false">
  <script data-pagedtable-source type="application/json">
{"columns":[{"label":[""],"name":["_rn_"],"type":[""],"align":["left"]},{"label":["name"],"name":[1],"type":["chr"],"align":["left"]},{"label":["Hyper_Fold_Enrichment"],"name":[2],"type":["dbl"],"align":["right"]},{"label":["Total_Genes_Annotated"],"name":[3],"type":["dbl"],"align":["right"]},{"label":["Hyper_Background_Gene_Hits"],"name":[4],"type":["dbl"],"align":["right"]},{"label":["Hyper_Foreground_Gene_Hits"],"name":[5],"type":["dbl"],"align":["right"]},{"label":["Hyper_Adjp_BH"],"name":[6],"type":["dbl"],"align":["right"]}],"data":[{"1":"pattern specification process","2":"1.254823","3":"440","4":"437","5":"391","6":"4.281859e-46","_rn_":"1"},{"1":"regionalization","2":"1.252569","3":"344","4":"342","5":"310","6":"1.873102e-37","_rn_":"6"},{"1":"epithelium development","2":"1.147534","3":"1041","4":"1008","5":"835","6":"4.883623e-37","_rn_":"7"},{"1":"skeletal system development","2":"1.197881","3":"462","4":"461","5":"409","6":"1.585783e-32","_rn_":"13"},{"1":"blood vessel development","2":"1.191819","3":"502","4":"498","5":"430","6":"2.574307e-31","_rn_":"14"},{"1":"vasculature development","2":"1.187865","3":"531","4":"527","5":"451","6":"2.904496e-31","_rn_":"15"},{"1":"cell differentiation","2":"1.059532","3":"3498","4":"3370","5":"2694","6":"5.044873e-23","_rn_":"40"},{"1":"muscle tissue development","2":"1.190195","3":"329","4":"327","5":"282","6":"1.896923e-19","_rn_":"59"}],"options":{"columns":{"min":{},"max":[10]},"rows":{"min":[10],"max":[10]},"pages":{}}}
  </script>
</div>

As observed with the somite trios, many of the significant terms using DE genes are also found when using the DA regions, and these terms tend to have much lower p-values compared to those that are only identified with DE genes.


```r
GO.all$p.value_great <- go.dars[[1]][match(GO.all$GO.ID, go.dars[[1]]$ID),]$Hyper_Raw_PValue
GO.all$sig_great <- GO.all$p.value_great<0.05

ggplot(GO.all[!is.na(GO.all$p.value_great),], aes(sig_great, p.value_classic, colour=sig_great)) +
  geom_boxplot() +
  scale_color_manual(values=c('TRUE'="red", 'FALSE'="black")) +
  xlab("term is significant with GREAT") +
  ylab("p-value of enrichment with DE genes") +
  th + theme(legend.position = "none")
```

![](02_somitogenesisDevelopment_files/figure-html/compere_degs-1.png)<!-- -->

However, we can also look at enrichment of mouse phenotype data based on single KO mice. Most of the highly significant enriched phenotypes relate to skeleton abnormalities, validating that we are picking up relevant regions. There is also a very strong enrichment for embryonic-lethal genes, which is expected for important developmental regulators.


```r
## mouse phenotype
pheno.dars[['Mouse Phenotype Single KO']][c(1:3,10,11:20,33), c(1,2,6,11:9,13)]
```

<div data-pagedtable="false">
  <script data-pagedtable-source type="application/json">
{"columns":[{"label":[""],"name":["_rn_"],"type":[""],"align":["left"]},{"label":["ID"],"name":[1],"type":["chr"],"align":["left"]},{"label":["name"],"name":[2],"type":["chr"],"align":["left"]},{"label":["Hyper_Fold_Enrichment"],"name":[3],"type":["dbl"],"align":["right"]},{"label":["Total_Genes_Annotated"],"name":[4],"type":["dbl"],"align":["right"]},{"label":["Hyper_Background_Gene_Hits"],"name":[5],"type":["dbl"],"align":["right"]},{"label":["Hyper_Foreground_Gene_Hits"],"name":[6],"type":["dbl"],"align":["right"]},{"label":["Hyper_Adjp_BH"],"name":[7],"type":["dbl"],"align":["right"]}],"data":[{"1":"MP:0010832","2":"lethality during fetal growth through weaning","3":"1.099021","4":"1760","5":"1748","6":"1452","7":"2.684791e-31","_rn_":"1"},{"1":"MP:0002081","2":"perinatal lethality","3":"1.119604","4":"983","5":"978","6":"810","7":"4.082378e-24","_rn_":"2"},{"1":"MP:0002114","2":"abnormal axial skeleton morphology","3":"1.110581","4":"1157","5":"1150","6":"954","7":"4.077667e-21","_rn_":"3"},{"1":"MP:0004625","2":"abnormal rib joint","3":"1.409417","4":"73","5":"73","6":"65","7":"4.064631e-19","_rn_":"10"},{"1":"MP:0009250","2":"abnormal appendicular skeleton morphology","3":"1.145105","4":"606","5":"601","6":"513","7":"6.987291e-19","_rn_":"11"},{"1":"MP:0002932","2":"abnormal joint morphology","3":"1.207325","4":"276","5":"276","6":"238","7":"1.883921e-18","_rn_":"12"},{"1":"MP:0008148","2":"abnormal sternocostal joint morphology","3":"1.442099","4":"59","5":"59","6":"55","7":"2.298604e-18","_rn_":"13"},{"1":"MP:0005508","2":"abnormal skeleton morphology","3":"1.080918","4":"1832","5":"1815","6":"1454","7":"2.876237e-18","_rn_":"14"},{"1":"MP:0000157","2":"abnormal sternum morphology","3":"1.245194","4":"167","5":"167","6":"147","7":"3.640513e-18","_rn_":"15"},{"1":"MP:0000137","2":"abnormal vertebrae morphology","3":"1.162887","4":"509","5":"506","6":"413","7":"4.930706e-18","_rn_":"16"},{"1":"MP:0000267","2":"abnormal heart development","3":"1.183360","4":"299","5":"299","6":"254","7":"1.182803e-17","_rn_":"17"},{"1":"MP:0005390","2":"skeleton phenotype","3":"1.077184","4":"1928","5":"1908","6":"1515","7":"1.835904e-17","_rn_":"18"},{"1":"MP:0000459","2":"abnormal presacral vertebrae morphology","3":"1.216134","4":"280","5":"279","6":"232","7":"7.422608e-17","_rn_":"19"},{"1":"MP:0004703","2":"abnormal vertebral column morphology","3":"1.131615","4":"686","5":"680","6":"558","7":"2.601176e-16","_rn_":"20"},{"1":"MP:0004624","2":"abnormal thoracic cage morphology","3":"1.150259","4":"432","5":"432","6":"357","7":"9.528605e-14","_rn_":"33"}],"options":{"columns":{"min":{},"max":[10]},"rows":{"min":[10],"max":[10]},"pages":{}}}
  </script>
</div>

We can repeat this analysis using only the regions with preferential accessibility in somites from one fate.

##### Cervical

The 9327 regions with higher accessibility in **cervical** somites again show similar terms for those with highest enrichments: 


```r
## run GREAT
job.cervical = submitGreatJob(gr = da[grepl("cervical", da$fate.sub)],
                              bg = background,
                              species = "mm10")
## enrichment results
go.dars.cervical <- getEnrichmentTables(job.cervical, ontology = "GO Biological Process")
pheno.dars.cervical <- getEnrichmentTables(job.cervical, category = "Phenotype")

## GO biological process
go.dars.cervical[[1]][c(1,2,7,16,32,66,67,77,99), c(1,2,6,11:9,13)]
```

<div data-pagedtable="false">
  <script data-pagedtable-source type="application/json">
{"columns":[{"label":[""],"name":["_rn_"],"type":[""],"align":["left"]},{"label":["ID"],"name":[1],"type":["chr"],"align":["left"]},{"label":["name"],"name":[2],"type":["chr"],"align":["left"]},{"label":["Hyper_Fold_Enrichment"],"name":[3],"type":["dbl"],"align":["right"]},{"label":["Total_Genes_Annotated"],"name":[4],"type":["dbl"],"align":["right"]},{"label":["Hyper_Background_Gene_Hits"],"name":[5],"type":["dbl"],"align":["right"]},{"label":["Hyper_Foreground_Gene_Hits"],"name":[6],"type":["dbl"],"align":["right"]},{"label":["Hyper_Adjp_BH"],"name":[7],"type":["dbl"],"align":["right"]}],"data":[{"1":"GO:0003002","2":"regionalization","3":"1.451247","4":"344","5":"342","6":"243","7":"5.016346e-26","_rn_":"1"},{"1":"GO:0007389","2":"pattern specification process","3":"1.393841","4":"440","5":"437","6":"294","7":"9.812918e-25","_rn_":"2"},{"1":"GO:0060429","2":"epithelium development","3":"1.198765","4":"1041","5":"1008","6":"572","7":"3.281730e-15","_rn_":"7"},{"1":"GO:0002009","2":"morphogenesis of an epithelium","3":"1.264010","4":"498","5":"496","6":"295","7":"9.435542e-14","_rn_":"16"},{"1":"GO:0009952","2":"anterior/posterior pattern specification","3":"1.367619","4":"219","5":"219","6":"149","7":"1.083166e-10","_rn_":"32"},{"1":"GO:0035282","2":"segmentation","3":"1.464296","4":"99","5":"99","6":"70","7":"1.871746e-08","_rn_":"66"},{"1":"GO:0030154","2":"cell differentiation","3":"1.074391","4":"3498","5":"3370","6":"1673","7":"1.938688e-08","_rn_":"67"},{"1":"GO:0001501","2":"skeletal system development","3":"1.201484","4":"462","5":"461","6":"294","7":"6.278988e-08","_rn_":"77"},{"1":"GO:0030855","2":"epithelial cell differentiation","3":"1.197797","4":"543","5":"511","6":"283","7":"3.418710e-07","_rn_":"99"}],"options":{"columns":{"min":{},"max":[10]},"rows":{"min":[10],"max":[10]},"pages":{}}}
  </script>
</div>

And when looking at the effect of KOs, these regions are still enriched for phenotypes related to skeleton abnormalities, including some specific to cervical vertebrae.


```r
## mouse phenotype
pheno.dars.cervical[['Mouse Phenotype Single KO']][c(19,20,24,30,51,69,77,90,102,109,112), c(1,2,6,11:9,13)]
```

<div data-pagedtable="false">
  <script data-pagedtable-source type="application/json">
{"columns":[{"label":[""],"name":["_rn_"],"type":[""],"align":["left"]},{"label":["ID"],"name":[1],"type":["chr"],"align":["left"]},{"label":["name"],"name":[2],"type":["chr"],"align":["left"]},{"label":["Hyper_Fold_Enrichment"],"name":[3],"type":["dbl"],"align":["right"]},{"label":["Total_Genes_Annotated"],"name":[4],"type":["dbl"],"align":["right"]},{"label":["Hyper_Background_Gene_Hits"],"name":[5],"type":["dbl"],"align":["right"]},{"label":["Hyper_Foreground_Gene_Hits"],"name":[6],"type":["dbl"],"align":["right"]},{"label":["Hyper_Adjp_BH"],"name":[7],"type":["dbl"],"align":["right"]}],"data":[{"1":"MP:0000459","2":"abnormal presacral vertebrae morphology","3":"1.322138","4":"280","5":"279","6":"156","7":"3.021614e-08","_rn_":"19"},{"1":"MP:0002084","2":"abnormal developmental patterning","3":"1.258582","4":"401","5":"401","6":"206","7":"4.749867e-08","_rn_":"20"},{"1":"MP:0002114","2":"abnormal axial skeleton morphology","3":"1.140260","4":"1157","5":"1150","6":"643","7":"6.400114e-08","_rn_":"24"},{"1":"MP:0004624","2":"abnormal thoracic cage morphology","3":"1.228940","4":"432","5":"432","6":"247","7":"4.015794e-07","_rn_":"30"},{"1":"MP:0001688","2":"abnormal somite development","3":"1.319238","4":"177","5":"177","6":"107","7":"4.802009e-06","_rn_":"51"},{"1":"MP:0002823","2":"abnormal rib development","3":"2.086345","4":"17","5":"17","6":"12","7":"9.051436e-06","_rn_":"69"},{"1":"MP:0004703","2":"abnormal vertebral column morphology","3":"1.156566","4":"686","5":"680","6":"375","7":"1.858508e-05","_rn_":"77"},{"1":"MP:0000137","2":"abnormal vertebrae morphology","3":"1.179428","4":"509","5":"506","6":"279","7":"3.207384e-05","_rn_":"90"},{"1":"MP:0003048","2":"abnormal cervical vertebrae morphology","3":"1.346201","4":"126","5":"126","6":"83","7":"5.995657e-05","_rn_":"102"},{"1":"MP:0004620","2":"cervical vertebral fusion","3":"1.558204","4":"40","5":"40","6":"27","7":"7.374018e-05","_rn_":"109"},{"1":"MP:0005508","2":"abnormal skeleton morphology","3":"1.084025","4":"1832","5":"1815","6":"942","7":"8.077651e-05","_rn_":"112"}],"options":{"columns":{"min":{},"max":[10]},"rows":{"min":[10],"max":[10]},"pages":{}}}
  </script>
</div>

##### Thoracic

The 17960 regions with higher accessibility in **thoracic** somites show a similar behaviour, recapitulating the results from DE genes. 


```r
## run GREAT
job.thoracic = submitGreatJob(gr = da[grepl("thoracic", da$fate.sub)],
                              bg = background,
                              species = "mm10")
```

```
## 
  |                                                                            
  |                                                                      |   0%
  |                                                                            
  |===                                                                   |   4%
  |                                                                            
  |======                                                                |   9%
  |                                                                            
  |=========                                                             |  13%
  |                                                                            
  |============                                                          |  18%
  |                                                                            
  |================                                                      |  22%
  |                                                                            
  |===================                                                   |  27%
  |                                                                            
  |======================                                                |  31%
  |                                                                            
  |=========================                                             |  35%
  |                                                                            
  |============================                                          |  40%
  |                                                                            
  |===============================                                       |  44%
  |                                                                            
  |==================================                                    |  49%
  |                                                                            
  |=====================================                                 |  53%
  |                                                                            
  |========================================                              |  58%
  |                                                                            
  |===========================================                           |  62%
  |                                                                            
  |===============================================                       |  66%
  |                                                                            
  |==================================================                    |  71%
  |                                                                            
  |=====================================================                 |  75%
  |                                                                            
  |========================================================              |  80%
  |                                                                            
  |===========================================================           |  84%
  |                                                                            
  |==============================================================        |  89%
  |                                                                            
  |=================================================================     |  93%
  |                                                                            
  |====================================================================  |  97%
```

```r
## enrichment results
go.dars.thoracic <- getEnrichmentTables(job.thoracic, ontology = "GO Biological Process")
pheno.dars.thoracic <- getEnrichmentTables(job.thoracic, category = "Phenotype")

## GO biological process
go.dars.thoracic[[1]][c(8,13,16:21,27,50:53,58:61,65,85,89:90), c(1,2,6,11:9,13)]
```

<div data-pagedtable="false">
  <script data-pagedtable-source type="application/json">
{"columns":[{"label":[""],"name":["_rn_"],"type":[""],"align":["left"]},{"label":["ID"],"name":[1],"type":["chr"],"align":["left"]},{"label":["name"],"name":[2],"type":["chr"],"align":["left"]},{"label":["Hyper_Fold_Enrichment"],"name":[3],"type":["dbl"],"align":["right"]},{"label":["Total_Genes_Annotated"],"name":[4],"type":["dbl"],"align":["right"]},{"label":["Hyper_Background_Gene_Hits"],"name":[5],"type":["dbl"],"align":["right"]},{"label":["Hyper_Foreground_Gene_Hits"],"name":[6],"type":["dbl"],"align":["right"]},{"label":["Hyper_Adjp_BH"],"name":[7],"type":["dbl"],"align":["right"]}],"data":[{"1":"GO:0007389","2":"pattern specification process","3":"1.245216","4":"440","5":"437","6":"312","7":"1.553112e-20","_rn_":"8"},{"1":"GO:1901861","2":"regulation of muscle tissue development","3":"1.404016","4":"134","5":"134","6":"96","7":"1.170762e-18","_rn_":"13"},{"1":"GO:0060429","2":"epithelium development","3":"1.149134","4":"1041","5":"1008","6":"664","7":"3.510753e-18","_rn_":"16"},{"1":"GO:0001568","2":"blood vessel development","3":"1.206966","4":"502","5":"498","6":"341","7":"2.338179e-17","_rn_":"17"},{"1":"GO:0016202","2":"regulation of striated muscle tissue development","3":"1.392225","4":"131","5":"131","6":"94","7":"2.338179e-17","_rn_":"18"},{"1":"GO:0060485","2":"mesenchyme development","3":"1.297957","4":"194","5":"194","6":"151","7":"5.427805e-17","_rn_":"19"},{"1":"GO:0001944","2":"vasculature development","3":"1.200124","4":"531","5":"527","6":"353","7":"6.427231e-17","_rn_":"20"},{"1":"GO:0048634","2":"regulation of muscle organ development","3":"1.381767","4":"135","5":"135","6":"96","7":"1.038625e-16","_rn_":"21"},{"1":"GO:0048762","2":"mesenchymal cell differentiation","3":"1.349658","4":"134","5":"134","6":"105","7":"4.669797e-16","_rn_":"27"},{"1":"GO:0001501","2":"skeletal system development","3":"1.175188","4":"462","5":"461","6":"321","7":"2.204370e-12","_rn_":"50"},{"1":"GO:0003002","2":"regionalization","3":"1.207517","4":"344","5":"342","6":"249","7":"2.204443e-12","_rn_":"51"},{"1":"GO:0045595","2":"regulation of cell differentiation","3":"1.095847","4":"1670","5":"1588","6":"998","7":"2.475362e-12","_rn_":"52"},{"1":"GO:1901862","2":"negative regulation of muscle tissue development","3":"1.511605","4":"55","5":"55","6":"37","7":"3.109811e-12","_rn_":"53"},{"1":"GO:0045843","2":"negative regulation of striated muscle tissue development","3":"1.511755","4":"52","5":"52","6":"35","7":"8.903899e-12","_rn_":"58"},{"1":"GO:0055026","2":"negative regulation of cardiac muscle tissue development","3":"1.587399","4":"31","5":"31","6":"23","7":"1.088219e-11","_rn_":"59"},{"1":"GO:0048635","2":"negative regulation of muscle organ development","3":"1.505798","4":"54","5":"54","6":"36","7":"1.206367e-11","_rn_":"60"},{"1":"GO:0045844","2":"positive regulation of striated muscle tissue development","3":"1.449984","4":"65","5":"65","6":"51","7":"1.344191e-11","_rn_":"61"},{"1":"GO:0045597","2":"positive regulation of cell differentiation","3":"1.125996","4":"899","5":"880","6":"576","7":"2.513000e-11","_rn_":"65"},{"1":"GO:0048706","2":"embryonic skeletal system development","3":"1.290577","4":"133","5":"133","6":"103","7":"2.107441e-10","_rn_":"85"},{"1":"GO:0060537","2":"muscle tissue development","3":"1.199106","4":"329","5":"327","6":"234","7":"2.708770e-10","_rn_":"89"},{"1":"GO:0035282","2":"segmentation","3":"1.348408","4":"99","5":"99","6":"71","7":"2.845078e-10","_rn_":"90"}],"options":{"columns":{"min":{},"max":[10]},"rows":{"min":[10],"max":[10]},"pages":{}}}
  </script>
</div>

And similar KO phenotypes, although with some related to more specific derivatives, such as myotome or cartilage. This perhaps reflects the higher statistical power that comes with such a large number of regions.


```r
## mouse phenotype
pheno.dars.thoracic[['Mouse Phenotype Single KO']][c(2,13,14,24,29,39,41,48,56,75,88,95,109,123,124), c(1,2,6,11:9,13)]
```

<div data-pagedtable="false">
  <script data-pagedtable-source type="application/json">
{"columns":[{"label":[""],"name":["_rn_"],"type":[""],"align":["left"]},{"label":["ID"],"name":[1],"type":["chr"],"align":["left"]},{"label":["name"],"name":[2],"type":["chr"],"align":["left"]},{"label":["Hyper_Fold_Enrichment"],"name":[3],"type":["dbl"],"align":["right"]},{"label":["Total_Genes_Annotated"],"name":[4],"type":["dbl"],"align":["right"]},{"label":["Hyper_Background_Gene_Hits"],"name":[5],"type":["dbl"],"align":["right"]},{"label":["Hyper_Foreground_Gene_Hits"],"name":[6],"type":["dbl"],"align":["right"]},{"label":["Hyper_Adjp_BH"],"name":[7],"type":["dbl"],"align":["right"]}],"data":[{"1":"MP:0008272","2":"abnormal endochondral bone ossification","3":"1.519414","4":"93","5":"92","6":"66","7":"1.184216e-17","_rn_":"2"},{"1":"MP:0004915","2":"abnormal Reichert's cartilage morphology","3":"2.807876","4":"3","5":"3","6":"3","7":"4.894990e-13","_rn_":"13"},{"1":"MP:0009250","2":"abnormal appendicular skeleton morphology","3":"1.174149","4":"606","5":"601","6":"417","7":"5.648246e-13","_rn_":"14"},{"1":"MP:0008277","2":"abnormal sternum ossification","3":"1.664482","4":"37","5":"37","6":"27","7":"1.357011e-11","_rn_":"24"},{"1":"MP:0004625","2":"abnormal rib joint","3":"1.456314","4":"73","5":"73","6":"56","7":"2.407011e-11","_rn_":"29"},{"1":"MP:0008148","2":"abnormal sternocostal joint morphology","3":"1.484992","4":"59","5":"59","6":"47","7":"1.372677e-10","_rn_":"39"},{"1":"MP:0000157","2":"abnormal sternum morphology","3":"1.265413","4":"167","5":"167","6":"124","7":"2.603622e-10","_rn_":"41"},{"1":"MP:0003939","2":"abnormal myotome morphology","3":"2.242654","4":"8","5":"8","6":"8","7":"4.296220e-10","_rn_":"48"},{"1":"MP:0008275","2":"failure of endochondral bone ossification","3":"2.831472","4":"3","5":"3","6":"3","7":"9.977152e-10","_rn_":"56"},{"1":"MP:0005508","2":"abnormal skeleton morphology","3":"1.082351","4":"1832","5":"1815","6":"1136","7":"1.703677e-09","_rn_":"75"},{"1":"MP:0002114","2":"abnormal axial skeleton morphology","3":"1.101233","4":"1157","5":"1150","6":"739","7":"4.164810e-09","_rn_":"88"},{"1":"MP:0005390","2":"skeleton phenotype","3":"1.077829","4":"1928","5":"1908","6":"1185","7":"5.841143e-09","_rn_":"95"},{"1":"MP:0000137","2":"abnormal vertebrae morphology","3":"1.157109","4":"509","5":"506","6":"309","7":"1.682977e-08","_rn_":"109"},{"1":"MP:0003419","2":"delayed endochondral bone ossification","3":"1.469905","4":"45","5":"45","6":"34","7":"3.042241e-08","_rn_":"123"},{"1":"MP:0005369","2":"muscle phenotype","3":"1.100279","4":"1139","5":"1128","6":"725","7":"3.406332e-08","_rn_":"124"}],"options":{"columns":{"min":{},"max":[10]},"rows":{"min":[10],"max":[10]},"pages":{}}}
  </script>
</div>

##### Lumbar

The 5710 regions with higher accessibility in **lumbar** somites show less significant enrichments, but the most prominent relate to angiogenesis and differentiation of derivatives, in agreement with the DE results.


```r
## run GREAT
job.lumbar = submitGreatJob(gr = da[grepl("lumbar", da$fate.sub)],
                              bg = background,
                              species = "mm10")
```

```
## 
  |                                                                            
  |                                                                      |   0%
  |                                                                            
  |===                                                                   |   5%
  |                                                                            
  |=======                                                               |   9%
  |                                                                            
  |==========                                                            |  14%
  |                                                                            
  |=============                                                         |  19%
  |                                                                            
  |================                                                      |  23%
  |                                                                            
  |====================                                                  |  28%
  |                                                                            
  |=======================                                               |  33%
  |                                                                            
  |==========================                                            |  38%
  |                                                                            
  |==============================                                        |  42%
  |                                                                            
  |=================================                                     |  47%
  |                                                                            
  |====================================                                  |  52%
  |                                                                            
  |=======================================                               |  56%
  |                                                                            
  |===========================================                           |  61%
  |                                                                            
  |==============================================                        |  66%
  |                                                                            
  |=================================================                     |  70%
  |                                                                            
  |=====================================================                 |  75%
  |                                                                            
  |========================================================              |  80%
  |                                                                            
  |===========================================================           |  85%
  |                                                                            
  |==============================================================        |  89%
  |                                                                            
  |==================================================================    |  94%
  |                                                                            
  |===================================================================== |  99%
```

```r
## enrichment results
go.dars.lumbar <- getEnrichmentTables(job.lumbar, ontology = "GO Biological Process")
pheno.dars.lumbar <- getEnrichmentTables(job.lumbar, category = "Phenotype")

## GO biological process
go.dars.lumbar[[1]][c(2,6,15,17,33,34,45,51,65,75), c(1,2,6,11:9,13)]
```

<div data-pagedtable="false">
  <script data-pagedtable-source type="application/json">
{"columns":[{"label":[""],"name":["_rn_"],"type":[""],"align":["left"]},{"label":["ID"],"name":[1],"type":["chr"],"align":["left"]},{"label":["name"],"name":[2],"type":["chr"],"align":["left"]},{"label":["Hyper_Fold_Enrichment"],"name":[3],"type":["dbl"],"align":["right"]},{"label":["Total_Genes_Annotated"],"name":[4],"type":["dbl"],"align":["right"]},{"label":["Hyper_Background_Gene_Hits"],"name":[5],"type":["dbl"],"align":["right"]},{"label":["Hyper_Foreground_Gene_Hits"],"name":[6],"type":["dbl"],"align":["right"]},{"label":["Hyper_Adjp_BH"],"name":[7],"type":["dbl"],"align":["right"]}],"data":[{"1":"GO:0045746","2":"negative regulation of Notch signaling pathway","3":"2.122161","4":"33","5":"33","6":"21","7":"8.239062e-05","_rn_":"2"},{"1":"GO:0097084","2":"vascular smooth muscle cell development","3":"2.607880","4":"9","5":"9","6":"7","7":"5.368925e-04","_rn_":"6"},{"1":"GO:0030512","2":"negative regulation of transforming growth factor beta receptor signaling pathway","3":"1.673284","4":"56","5":"56","6":"26","7":"4.940597e-03","_rn_":"15"},{"1":"GO:0008593","2":"regulation of Notch signaling pathway","3":"1.489767","4":"77","5":"76","6":"47","7":"5.641071e-03","_rn_":"17"},{"1":"GO:0043588","2":"skin development","3":"1.300150","4":"246","5":"216","6":"86","7":"1.191564e-02","_rn_":"33"},{"1":"GO:0035886","2":"vascular smooth muscle cell differentiation","3":"1.734369","4":"23","5":"23","6":"16","7":"1.191564e-02","_rn_":"34"},{"1":"GO:0051145","2":"smooth muscle cell differentiation","3":"1.473824","4":"46","5":"46","6":"33","7":"1.461555e-02","_rn_":"45"},{"1":"GO:0010830","2":"regulation of myotube differentiation","3":"1.609733","4":"56","5":"54","6":"28","7":"1.539439e-02","_rn_":"51"},{"1":"GO:0042693","2":"muscle cell fate commitment","3":"2.018692","4":"11","5":"11","6":"9","7":"2.249740e-02","_rn_":"65"},{"1":"GO:0060215","2":"primitive hemopoiesis","3":"2.344946","4":"7","5":"7","6":"5","7":"2.893183e-02","_rn_":"75"}],"options":{"columns":{"min":{},"max":[10]},"rows":{"min":[10],"max":[10]},"pages":{}}}
  </script>
</div>

Interestingly, the KO phenotypes are not as clearly enriched for skeleton related phenotypes as with the other sets of DA regions. Instead, the strongest enrichments relate to erythrocytes, again consistent with the results from the DE genes.


```r
## mouse phenotype
pheno.dars.lumbar[['Mouse Phenotype Single KO']][c(4,7:10,12,15,17,23,26:30,50,52,57,59), c(1,2,6,11:9,13)]
```

<div data-pagedtable="false">
  <script data-pagedtable-source type="application/json">
{"columns":[{"label":[""],"name":["_rn_"],"type":[""],"align":["left"]},{"label":["ID"],"name":[1],"type":["chr"],"align":["left"]},{"label":["name"],"name":[2],"type":["chr"],"align":["left"]},{"label":["Hyper_Fold_Enrichment"],"name":[3],"type":["dbl"],"align":["right"]},{"label":["Total_Genes_Annotated"],"name":[4],"type":["dbl"],"align":["right"]},{"label":["Hyper_Background_Gene_Hits"],"name":[5],"type":["dbl"],"align":["right"]},{"label":["Hyper_Foreground_Gene_Hits"],"name":[6],"type":["dbl"],"align":["right"]},{"label":["Hyper_Adjp_BH"],"name":[7],"type":["dbl"],"align":["right"]}],"data":[{"1":"MP:0003396","2":"abnormal embryonic hematopoiesis","3":"1.597634","4":"54","5":"54","6":"33","7":"0.005281491","_rn_":"4"},{"1":"MP:0013657","2":"abnormal blood cell morphology","3":"1.123536","4":"1567","5":"1531","6":"561","7":"0.009195859","_rn_":"7"},{"1":"MP:0002447","2":"abnormal erythrocyte morphology","3":"1.188314","4":"712","5":"701","6":"263","7":"0.009195859","_rn_":"8"},{"1":"MP:0014060","2":"abnormal platelet alpha-granule morphology","3":"2.896385","4":"9","5":"9","6":"5","7":"0.009195859","_rn_":"9"},{"1":"MP:0013659","2":"abnormal erythroid lineage cell morphology","3":"1.177998","4":"768","5":"756","6":"283","7":"0.009195859","_rn_":"10"},{"1":"MP:0004612","2":"fusion of vertebral bodies","3":"2.471868","4":"11","5":"11","6":"8","7":"0.010077357","_rn_":"12"},{"1":"MP:0012283","2":"decreased sternebra number","3":"3.785048","4":"5","5":"5","6":"5","7":"0.012728834","_rn_":"15"},{"1":"MP:0004229","2":"abnormal embryonic erythropoiesis","3":"1.714995","4":"29","5":"29","6":"17","7":"0.014905147","_rn_":"17"},{"1":"MP:0000061","2":"fragile skeleton","3":"2.153562","4":"25","5":"25","6":"12","7":"0.016363757","_rn_":"23"},{"1":"MP:0013661","2":"abnormal myeloid cell number","3":"1.211837","4":"442","5":"435","6":"160","7":"0.021579232","_rn_":"26"},{"1":"MP:0001586","2":"abnormal erythrocyte cell number","3":"1.273851","4":"241","5":"238","6":"86","7":"0.022433818","_rn_":"27"},{"1":"MP:0004230","2":"abnormal embryonic erythrocyte morphology","3":"1.636778","4":"35","5":"35","6":"18","7":"0.022433818","_rn_":"28"},{"1":"MP:0009931","2":"abnormal skin appearance","3":"1.223734","4":"386","5":"382","6":"161","7":"0.022433818","_rn_":"29"},{"1":"MP:0002655","2":"abnormal keratinocyte morphology","3":"1.667539","4":"48","5":"47","6":"25","7":"0.022433818","_rn_":"30"},{"1":"MP:0000141","2":"abnormal vertebral body morphology","3":"1.452620","4":"55","5":"55","6":"29","7":"0.040323003","_rn_":"50"},{"1":"MP:0002258","2":"abnormal cricoid cartilage morphology","3":"1.855895","4":"15","5":"15","6":"9","7":"0.044698759","_rn_":"52"},{"1":"MP:0002429","2":"abnormal blood cell morphology/development","3":"1.130588","4":"850","5":"831","6":"334","7":"0.049885363","_rn_":"57"},{"1":"MP:0003794","2":"delayed somite formation","3":"2.320759","4":"13","5":"13","6":"6","7":"0.053092712","_rn_":"59"}],"options":{"columns":{"min":{},"max":[10]},"rows":{"min":[10],"max":[10]},"pages":{}}}
  </script>
</div>

##### Sacral

The 8789 regions with higher accessibility in **sacral** somites are very strongly enriched for processes related to the skeletal system and differentiation of its components.


```r
## run GREAT
job.sacral = submitGreatJob(gr = da[grepl("sacral", da$fate.sub)],
                              bg = background,
                              species = "mm10")
```

```
## 
  |                                                                            
  |                                                                      |   0%
  |                                                                            
  |===                                                                   |   4%
  |                                                                            
  |=====                                                                 |   8%
  |                                                                            
  |========                                                              |  12%
  |                                                                            
  |===========                                                           |  16%
  |                                                                            
  |==============                                                        |  19%
  |                                                                            
  |================                                                      |  23%
  |                                                                            
  |===================                                                   |  27%
  |                                                                            
  |======================                                                |  31%
  |                                                                            
  |=========================                                             |  35%
  |                                                                            
  |===========================                                           |  39%
  |                                                                            
  |==============================                                        |  43%
  |                                                                            
  |=================================                                     |  47%
  |                                                                            
  |===================================                                   |  51%
  |                                                                            
  |======================================                                |  54%
  |                                                                            
  |=========================================                             |  58%
  |                                                                            
  |============================================                          |  62%
  |                                                                            
  |==============================================                        |  66%
  |                                                                            
  |=================================================                     |  70%
  |                                                                            
  |====================================================                  |  74%
  |                                                                            
  |======================================================                |  78%
  |                                                                            
  |=========================================================             |  82%
  |                                                                            
  |============================================================          |  86%
  |                                                                            
  |===============================================================       |  90%
  |                                                                            
  |=================================================================     |  93%
  |                                                                            
  |====================================================================  |  97%
```

```r
## enrichment results
go.dars.sacral <- getEnrichmentTables(job.sacral, ontology = "GO Biological Process")
pheno.dars.sacral <- getEnrichmentTables(job.sacral, category = "Phenotype")

## GO biological process
go.dars.sacral[[1]][c(1,2,5:8,12,16,26,31,33,38,40,61,72,74,81,82,84,89,108,113,119,120), c(1,2,6,11:9,13)]
```

<div data-pagedtable="false">
  <script data-pagedtable-source type="application/json">
{"columns":[{"label":[""],"name":["_rn_"],"type":[""],"align":["left"]},{"label":["ID"],"name":[1],"type":["chr"],"align":["left"]},{"label":["name"],"name":[2],"type":["chr"],"align":["left"]},{"label":["Hyper_Fold_Enrichment"],"name":[3],"type":["dbl"],"align":["right"]},{"label":["Total_Genes_Annotated"],"name":[4],"type":["dbl"],"align":["right"]},{"label":["Hyper_Background_Gene_Hits"],"name":[5],"type":["dbl"],"align":["right"]},{"label":["Hyper_Foreground_Gene_Hits"],"name":[6],"type":["dbl"],"align":["right"]},{"label":["Hyper_Adjp_BH"],"name":[7],"type":["dbl"],"align":["right"]}],"data":[{"1":"GO:0001501","2":"skeletal system development","3":"1.296260","4":"462","5":"461","6":"294","7":"2.694976e-14","_rn_":"1"},{"1":"GO:0060272","2":"embryonic skeletal joint morphogenesis","3":"2.839320","4":"12","5":"12","6":"9","7":"7.909624e-14","_rn_":"2"},{"1":"GO:0048706","2":"embryonic skeletal system development","3":"1.493064","4":"133","5":"133","6":"91","7":"5.717414e-12","_rn_":"5"},{"1":"GO:0048704","2":"embryonic skeletal system morphogenesis","3":"1.582416","4":"100","5":"100","6":"70","7":"1.425090e-11","_rn_":"6"},{"1":"GO:0048705","2":"skeletal system morphogenesis","3":"1.363692","4":"235","5":"234","6":"160","7":"5.902977e-11","_rn_":"7"},{"1":"GO:0072498","2":"embryonic skeletal joint development","3":"2.473779","4":"14","5":"14","6":"9","7":"6.667095e-11","_rn_":"8"},{"1":"GO:0060351","2":"cartilage development involved in endochondral bone morphogenesis","3":"1.995513","4":"30","5":"30","6":"22","7":"1.224541e-09","_rn_":"12"},{"1":"GO:0051216","2":"cartilage development","3":"1.387301","4":"163","5":"163","6":"110","7":"4.000808e-09","_rn_":"16"},{"1":"GO:0030154","2":"cell differentiation","3":"1.077952","4":"3498","5":"3370","6":"1629","7":"2.627054e-08","_rn_":"26"},{"1":"GO:0060350","2":"endochondral bone morphogenesis","3":"1.661358","4":"56","5":"56","6":"39","7":"5.790231e-08","_rn_":"31"},{"1":"GO:0045595","2":"regulation of cell differentiation","3":"1.115715","4":"1670","5":"1588","6":"836","7":"7.651085e-08","_rn_":"33"},{"1":"GO:0060348","2":"bone development","3":"1.343806","4":"179","5":"178","6":"114","7":"1.883096e-07","_rn_":"38"},{"1":"GO:0060394","2":"negative regulation of pathway-restricted SMAD protein phosphorylation","3":"2.291591","4":"12","5":"12","6":"9","7":"2.279778e-07","_rn_":"40"},{"1":"GO:0045746","2":"negative regulation of Notch signaling pathway","3":"1.872585","4":"33","5":"33","6":"24","7":"2.587335e-06","_rn_":"61"},{"1":"GO:0030512","2":"negative regulation of transforming growth factor beta receptor signaling pathway","3":"1.650284","4":"56","5":"56","6":"32","7":"8.025166e-06","_rn_":"72"},{"1":"GO:0045599","2":"negative regulation of fat cell differentiation","3":"1.522855","4":"50","5":"50","6":"34","7":"9.767144e-06","_rn_":"74"},{"1":"GO:0045165","2":"cell fate commitment","3":"1.227933","4":"247","5":"246","6":"161","7":"1.170928e-05","_rn_":"81"},{"1":"GO:0009952","2":"anterior/posterior pattern specification","3":"1.273341","4":"219","5":"219","6":"139","7":"1.180105e-05","_rn_":"82"},{"1":"GO:0001944","2":"vasculature development","3":"1.170683","4":"531","5":"527","6":"314","7":"1.332227e-05","_rn_":"84"},{"1":"GO:0009954","2":"proximal/distal pattern formation","3":"1.575135","4":"35","5":"35","6":"27","7":"1.708564e-05","_rn_":"89"},{"1":"GO:0061061","2":"muscle structure development","3":"1.176420","4":"473","5":"469","6":"259","7":"4.073608e-05","_rn_":"108"},{"1":"GO:0001568","2":"blood vessel development","3":"1.164107","4":"502","5":"498","6":"298","7":"4.649899e-05","_rn_":"113"},{"1":"GO:0051145","2":"smooth muscle cell differentiation","3":"1.471080","4":"46","5":"46","6":"34","7":"6.278040e-05","_rn_":"119"},{"1":"GO:0032330","2":"regulation of chondrocyte differentiation","3":"1.454796","4":"47","5":"47","6":"34","7":"6.396316e-05","_rn_":"120"}],"options":{"columns":{"min":{},"max":[10]},"rows":{"min":[10],"max":[10]},"pages":{}}}
  </script>
</div>

Similarly strong enrichments for a variety of skeletal phenotypes are observed with gene KOs.


```r
## mouse phenotype
pheno.dars.sacral[['Mouse Phenotype Single KO']][c(1,2,10,12,13,19,26,30,37,40,44,46,48,49,55,59,64,66,68,70,71,74,75,77,81,85,90,93,95), c(1,2,6,11:9,13)]
```

<div data-pagedtable="false">
  <script data-pagedtable-source type="application/json">
{"columns":[{"label":[""],"name":["_rn_"],"type":[""],"align":["left"]},{"label":["ID"],"name":[1],"type":["chr"],"align":["left"]},{"label":["name"],"name":[2],"type":["chr"],"align":["left"]},{"label":["Hyper_Fold_Enrichment"],"name":[3],"type":["dbl"],"align":["right"]},{"label":["Total_Genes_Annotated"],"name":[4],"type":["dbl"],"align":["right"]},{"label":["Hyper_Background_Gene_Hits"],"name":[5],"type":["dbl"],"align":["right"]},{"label":["Hyper_Foreground_Gene_Hits"],"name":[6],"type":["dbl"],"align":["right"]},{"label":["Hyper_Adjp_BH"],"name":[7],"type":["dbl"],"align":["right"]}],"data":[{"1":"MP:0004650","2":"increased lumbar vertebrae number","3":"5.532871","4":"5","5":"5","6":"3","7":"1.232076e-09","_rn_":"1"},{"1":"MP:0003036","2":"vertebral transformation","3":"1.745516","4":"89","5":"89","6":"57","7":"1.232076e-09","_rn_":"2"},{"1":"MP:0004622","2":"sacral vertebral fusion","3":"2.329159","4":"12","5":"12","6":"8","7":"3.277429e-07","_rn_":"10"},{"1":"MP:0004617","2":"sacral vertebral transformation","3":"2.555487","4":"14","5":"14","6":"11","7":"4.211951e-07","_rn_":"12"},{"1":"MP:0002932","2":"abnormal joint morphology","3":"1.281805","4":"276","5":"276","6":"164","7":"7.376143e-07","_rn_":"13"},{"1":"MP:0003050","2":"abnormal sacral vertebrae morphology","3":"1.507661","4":"110","5":"108","6":"49","7":"2.510701e-06","_rn_":"19"},{"1":"MP:0002114","2":"abnormal axial skeleton morphology","3":"1.124512","4":"1157","5":"1150","6":"586","7":"1.101209e-05","_rn_":"26"},{"1":"MP:0003049","2":"abnormal lumbar vertebrae morphology","3":"1.393157","4":"155","5":"154","6":"76","7":"1.686678e-05","_rn_":"30"},{"1":"MP:0005508","2":"abnormal skeleton morphology","3":"1.093775","4":"1832","5":"1815","6":"864","7":"3.408417e-05","_rn_":"37"},{"1":"MP:0005390","2":"skeleton phenotype","3":"1.090436","4":"1928","5":"1908","6":"903","7":"4.517374e-05","_rn_":"40"},{"1":"MP:0000137","2":"abnormal vertebrae morphology","3":"1.185985","4":"509","5":"506","6":"263","7":"5.985034e-05","_rn_":"44"},{"1":"MP:0004616","2":"lumbar vertebral transformation","3":"1.919262","4":"36","5":"36","6":"23","7":"6.298569e-05","_rn_":"46"},{"1":"MP:0002620","2":"abnormal monocyte morphology","3":"1.342986","4":"194","5":"189","6":"94","7":"6.298569e-05","_rn_":"48"},{"1":"MP:0008271","2":"abnormal bone ossification","3":"1.224273","4":"341","5":"337","6":"170","7":"7.826946e-05","_rn_":"49"},{"1":"MP:0004609","2":"vertebral fusion","3":"1.326683","4":"130","5":"130","6":"78","7":"1.111476e-04","_rn_":"55"},{"1":"MP:0001533","2":"abnormal skeleton physiology","3":"1.178156","4":"538","5":"530","6":"250","7":"1.949704e-04","_rn_":"59"},{"1":"MP:0000060","2":"delayed bone ossification","3":"1.348624","4":"110","5":"110","6":"74","7":"2.484970e-04","_rn_":"64"},{"1":"MP:0003396","2":"abnormal embryonic hematopoiesis","3":"1.474462","4":"54","5":"54","6":"36","7":"2.494922e-04","_rn_":"66"},{"1":"MP:0004612","2":"fusion of vertebral bodies","3":"2.236807","4":"11","5":"11","6":"8","7":"3.827439e-04","_rn_":"68"},{"1":"MP:0004618","2":"thoracic vertebral transformation","3":"1.661145","4":"47","5":"47","6":"29","7":"3.830270e-04","_rn_":"70"},{"1":"MP:0009250","2":"abnormal appendicular skeleton morphology","3":"1.145293","4":"606","5":"601","6":"315","7":"4.080512e-04","_rn_":"71"},{"1":"MP:0004703","2":"abnormal vertebral column morphology","3":"1.141370","4":"686","5":"680","6":"343","7":"4.558055e-04","_rn_":"74"},{"1":"MP:0012283","2":"decreased sternebra number","3":"3.337287","4":"5","5":"5","6":"5","7":"4.585982e-04","_rn_":"75"},{"1":"MP:0002447","2":"abnormal erythrocyte morphology","3":"1.156021","4":"712","5":"701","6":"311","7":"4.720346e-04","_rn_":"77"},{"1":"MP:0001586","2":"abnormal erythrocyte cell number","3":"1.252834","4":"241","5":"238","6":"108","7":"5.143511e-04","_rn_":"81"},{"1":"MP:0002896","2":"abnormal bone mineralization","3":"1.313818","4":"165","5":"162","6":"79","7":"5.880081e-04","_rn_":"85"},{"1":"MP:0008148","2":"abnormal sternocostal joint morphology","3":"1.436910","4":"59","5":"59","6":"42","7":"7.236217e-04","_rn_":"90"},{"1":"MP:0013657","2":"abnormal blood cell morphology","3":"1.098804","4":"1567","5":"1531","6":"662","7":"7.885951e-04","_rn_":"93"},{"1":"MP:0000154","2":"rib fusion","3":"1.414067","4":"70","5":"70","6":"43","7":"1.009222e-03","_rn_":"95"}],"options":{"columns":{"min":{},"max":[10]},"rows":{"min":[10],"max":[10]},"pages":{}}}
  </script>
</div>

####

---

We can next check the enrichment results for phenotypes specific to each type of vertebrae. We observe, for example, that *cervical vertebral transformation* phenotypes show a clear enrichment in the regions most accessible in the cervical somites, but not with any of the other three groups of DA loci. Similar results are obtained with each fate except lumbar regions: the enrichment is significant with the corresponding somites, and often not with regions that peak at other fates. As observed above, lumbar regions do not show enrichment for skeletal phenotypes.


```r
getPvals <- function(id=NULL){
c(all = pheno.dars[["Mouse Phenotype Single KO"]][
  pheno.dars[["Mouse Phenotype Single KO"]]$name == id,]$Hyper_Adjp_BH,
  cervical = pheno.dars.cervical[["Mouse Phenotype Single KO"]][
  pheno.dars.cervical[["Mouse Phenotype Single KO"]]$name == id,]$Hyper_Adjp_BH,
  thoracic = pheno.dars.thoracic[["Mouse Phenotype Single KO"]][
  pheno.dars.thoracic[["Mouse Phenotype Single KO"]]$name == id,]$Hyper_Adjp_BH,
  lumbar = pheno.dars.lumbar[["Mouse Phenotype Single KO"]][
  pheno.dars.lumbar[["Mouse Phenotype Single KO"]]$name == id,]$Hyper_Adjp_BH,
  sacral = pheno.dars.sacral[["Mouse Phenotype Single KO"]][
  pheno.dars.sacral[["Mouse Phenotype Single KO"]]$name == id,]$Hyper_Adjp_BH)
}

tmp <- rbind("abnormal cervical vertebrae morphology" = getPvals(id = "abnormal cervical vertebrae morphology"),
             "cervical vertebral fusion" = getPvals(id = "cervical vertebral fusion"),
             "cervical vertebral transformation" = getPvals(id = "cervical vertebral transformation"),

             "thoracic vertebral transformation" = getPvals(id = "thoracic vertebral transformation"),
             "abnormal thoracic cage morphology" = getPvals(id = "abnormal thoracic cage morphology"),
             "abnormal thoracic vertebrae morphology" = getPvals(id = "abnormal thoracic vertebrae morphology"),
             "abnormal thoracic cavity morphology" = getPvals(id = "abnormal thoracic cavity morphology"),
             "thoracic vertebral fusion" = getPvals(id = "thoracic vertebral fusion"),
             
             "abnormal rib joint" = getPvals(id = "abnormal rib joint"),
             "abnormal rib morphology" = getPvals(id = "abnormal rib morphology"),
             "abnormal rib development" = getPvals(id = "abnormal rib development"),
             "rib fusion" = getPvals(id = "rib fusion"),
             "increased rib number" = getPvals(id = "increased rib number"),
             "decreased rib number" = getPvals(id = "decreased rib number"),
             
             "abnormal lumbar vertebrae morphology" = getPvals(id = "abnormal lumbar vertebrae morphology"),
             "increased lumbar vertebrae number" = getPvals(id = "increased lumbar vertebrae number"),
             "decreased lumbar vertebrae number" = getPvals(id = "decreased lumbar vertebrae number"),
             "lumbar vertebral transformation" = getPvals(id = "lumbar vertebral transformation"),
             "lumbar vertebral fusion" = getPvals(id = "lumbar vertebral fusion"),
             
             "abnormal presacral vertebrae morphology" = getPvals(id = "abnormal presacral vertebrae morphology"),
             "decreased presacral vertebrae number" = getPvals(id = "decreased presacral vertebrae number"),
             "increased presacral vertebrae number" = getPvals(id = "increased presacral vertebrae number"),
             
             "abnormal sacral vertebrae morphology" = getPvals(id = "abnormal sacral vertebrae morphology"),
             "sacral vertebral fusion" = getPvals(id = "sacral vertebral fusion"),
             "increased sacral vertebrae number" = getPvals(id = "increased sacral vertebrae number"),
             "absent sacral vertebrae" = getPvals(id = "absent sacral vertebrae"),
             "sacral vertebral transformation" = getPvals(id = "sacral vertebral transformation"),
             "small sacral vertebrae" = getPvals(id = "small sacral vertebrae"))

h1 <- Heatmap(t(-log10(tmp[3:1,-1])), 
        col = colorRamp2(breaks = c(0,6),
                         colors = c("white","darkolivegreen4")),
        name = "-log10(FDR)",
        cluster_columns = FALSE, cluster_rows = FALSE)

h2 <- Heatmap(t(-log10(tmp[c(9:11,5,4,6),-1])), 
        col = colorRamp2(breaks = c(0,15),
                         colors = c("white","cornflowerblue")),
        name = "-log10(FDR)",
        cluster_columns = FALSE, cluster_rows = FALSE)

h3 <- Heatmap(t(-log10(tmp[c(16,15,18,17),-1])), 
        col = colorRamp2(breaks = c(0,15),
                         colors = c("white","darkorange")),
        name = "-log10(FDR)",
        cluster_columns = FALSE, cluster_rows = FALSE)

h4 <- Heatmap(t(-log10(tmp[c(27,24,23),-1])),
        col = colorRamp2(breaks = c(0,15),
                         colors = c("white","indianred1")),
        name = "-log10(FDR)",
        cluster_columns = FALSE, cluster_rows = FALSE)

## plot side by side
grid.newpage()
pushViewport(viewport(layout = grid.layout(nr = 1, nc = 4)))
pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 1))
draw(h1, newpage = FALSE)
upViewport()
pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 2))
draw(h2, newpage = FALSE)
upViewport()
pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 3))
draw(h3, newpage = FALSE)
upViewport()
pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 4))
draw(h4, newpage = FALSE)
upViewport()
```

![](02_somitogenesisDevelopment_files/figure-html/enriched_by_fate-1.png)<!-- -->

```r
## key terms
# Heatmap(-log10(tmp[c(1, 5,10, 15, 23),-1]),
#         col = brewer.pal(n=8, "Blues"),
#         name = "-log10(FDR)",
#         cluster_columns = FALSE, cluster_rows = FALSE)
```

These results demonstrate that the regions we have classified as preferentially open in one particular vertebral fate contain regulatory elements close to genes that, when knocked out, result in phenotypes affecting such vertebrae. 


#### Motif enrichment analysis

To understand which transcription factors are mediating some of the changes in expression, we search for overrepresentd binding motifs between the DA regions and an equivalent set of static open chromatin loci.


```r
## retrieve the FASTA sequences of the DA regions
# all DA
write.table(paste0("samtools faidx ", dir, "REFERENCE/Mus_musculus.GRCm38.dna.primary_assembly.fa ",
                   substr(row.names(dars), 4, 50), " >> ", dir,
                   "RNA+ATAC/results/02_DAregions_stage.fasta"), 
            paste0(dir, "RNA+ATAC/scripts/dataProcessing/02_getDNAseq_DAregions_stage.sh"), 
            quote = FALSE, row.names = FALSE, col.names = FALSE)

## similar set of nonDA to use as a control instead of shuffled seqs
remove <- which(row.names(peakAnn) %in% substr(row.names(dars),4,50))
peaks_nonDA <- peakAnn[-remove,1:3]
peaks_nonDA <- GRanges(peaks_nonDA$seqnames, IRanges(peaks_nonDA$start, peaks_nonDA$end),
                       name=row.names(peaks_nonDA))
## should have the same length distribution
l <- dars$end-dars$start+1 ## DA regions lengths
freqs <- table(cut(l,6)) ## define 10 intervals and count number of seqs in each
intervals <- gsub("[(](.+),(.+)[]]","\\1-\\2",names(freqs))

## randomly sample nonDA regions from all peaks that have the same length distribution
sel <- c()
for(i in 1:length(intervals)){
  limits <- as.numeric(unlist(strsplit(intervals[i], "-")))
  sel <- c(sel, sample(peaks_nonDA[width(peaks_nonDA) > limits[1] & width(peaks_nonDA) <= limits[2],],
                       freqs[i]))
}
control <- do.call('c', sel)
# plot(density(l))
# lines(density(width(control)), col="red")

write.table(paste0("samtools faidx ", dir, "REFERENCE/Mus_musculus.GRCm38.dna.primary_assembly.fa ",
                   control$name, " >> ", dir,
                   "RNA+ATAC/results/02_nonDA_controlRegions_stage.fasta"), 
            paste0(dir, "RNA+ATAC/scripts/dataProcessing/02_getDNAseq_nonDA_controlRegions_stage.sh"), 
            quote = FALSE, row.names = FALSE, col.names = FALSE)

## submit to AME (http://meme-suite.org/tools/ame) to identify enriched TF motifs
## ran against the human and mouse HOCOMOCOv11_full motif databases
## all other parameters default
## downloaded results and saved in dir/RNA+ATAC/results/02_AME_DAregions_stage.tsv

## do the same for each fate
# cervical
write.table(paste0("samtools faidx ", dir, "REFERENCE/Mus_musculus.GRCm38.dna.primary_assembly.fa ",
                   substr(row.names(dars[which(dars$fate=="cervical"),]), 4, 50), " >> ", dir,
                   "RNA+ATAC/results/02_DAregions_stage_cervical.fasta"), 
            paste0(dir, "RNA+ATAC/scripts/dataProcessing/02_getDNAseq_DAregions_stage_cervical.sh"), 
            quote = FALSE, row.names = FALSE, col.names = FALSE)
# thoracic
write.table(paste0("samtools faidx ", dir, "REFERENCE/Mus_musculus.GRCm38.dna.primary_assembly.fa ",
                   substr(row.names(dars[which(dars$fate=="thoracic"),]), 4, 50), " >> ", dir,
                   "RNA+ATAC/results/02_DAregions_stage_thoracic.fasta"), 
            paste0(dir, "RNA+ATAC/scripts/dataProcessing/02_getDNAseq_DAregions_stage_thoracic.sh"), 
            quote = FALSE, row.names = FALSE, col.names = FALSE)
# lumbar
write.table(paste0("samtools faidx ", dir, "REFERENCE/Mus_musculus.GRCm38.dna.primary_assembly.fa ",
                   substr(row.names(dars[which(dars$fate=="lumbar"),]), 4, 50), " >> ", dir,
                   "RNA+ATAC/results/02_DAregions_stage_lumbar.fasta"), 
            paste0(dir, "RNA+ATAC/scripts/dataProcessing/02_getDNAseq_DAregions_stage_lumbar.sh"), 
            quote = FALSE, row.names = FALSE, col.names = FALSE)
## regions tend to be larger, so the background might not be perfectly matched but it looks close enough
# sacral
write.table(paste0("samtools faidx ", dir, "REFERENCE/Mus_musculus.GRCm38.dna.primary_assembly.fa ",
                   substr(row.names(dars[which(dars$fate=="sacral"),]), 4, 50), " >> ", dir,
                   "RNA+ATAC/results/02_DAregions_stage_sacral.fasta"), 
            paste0(dir, "RNA+ATAC/scripts/dataProcessing/02_getDNAseq_DAregions_stage_sacral.sh"), 
            quote = FALSE, row.names = FALSE, col.names = FALSE)
## length distribution for all others is well matched with controls
```

First, we compute the enrichment for all 33013 DA regions, irrespective of their accessibility patterns. Not surprisingly, a large number of TFs are significantly overrepresented.


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
motifs_da <- read.table(paste0(dir, "RNA+ATAC/results/02_AME_DAregions_stage.tsv"),
                        sep = "\t", header = TRUE, stringsAsFactors = FALSE)[,-c(2,4)]

## remove TFs that are not expressed
not_expr <- which(motif_ann[motif_ann$Model %in% motifs_da$motif_ID,]$Transcription.factor %notin% geneCounts$gene)
not_expr <- motif_ann[motif_ann$Model %in% motifs_da$motif_ID,][not_expr,]$Model

motifs_da <- motifs_da[motifs_da$motif_ID %notin% not_expr,]

## annotate TF gene
motifs_da$TF <- motif_ann[match(motifs_da$motif_ID, motif_ann$Model),]$Transcription.factor
motifs_da$family <- motif_ann[match(motifs_da$motif_ID, motif_ann$Model),]$TF.family
motifs_da$subfamily <- motif_ann[match(motifs_da$motif_ID, motif_ann$Model),]$TF.subfamily

length(unique(motifs_da$TF)) # 268
```

```
## [1] 268
```

```r
# table(motifs_da$subfamily)
```

Many of the enriched TFs are similar to those observed with the DA regions across somite trios; strong representation of homeodomain, forkhead, Sox-related, and nuclear receptor families, among many others.


```r
motifs_da$label <- ifelse(motifs_da$rank < 35, motifs_da$TF, "")
ggplot(motifs_da, aes(log2(X.TP/X.FP), -log10(adj_p.value), label=label)) + 
  geom_point() + 
  xlab("log2 fold-change (DA / nonDA)") +
  ylab("-log10(adj p-value)") +
  geom_text_repel(size = 3, max.overlaps = 200) +
  th
```

![](02_somitogenesisDevelopment_files/figure-html/motifs_enriched-1.png)<!-- -->

However, we can also perform the enrichment analysis for regions with preferential accessibility in one fate, to pinpoint any changes in motif usage.


```r
read_AME <- function(file=NULL){
  motifs <- read.table(file, sep = "\t", header = TRUE, stringsAsFactors = FALSE)[,-c(2,4)]
  # remove non-expressed TFs
  motifs <- motifs[motifs$motif_ID %notin% not_expr,]
  # annotate TFs
  motifs$TF <- motif_ann[match(motifs$motif_ID, motif_ann$Model),]$Transcription.factor
  # a few motifs are missing from the annotation - add
  if(sum(grepl("EVX2_MOUSE.H11MO.0.A", motifs$motif_ID))>0){ motifs[motifs$motif_ID == "EVX2_MOUSE.H11MO.0.A",]$TF <- "Evx2" }
  if(sum(grepl("EVX1_MOUSE.H11MO.0.C", motifs$motif_ID))>0){ motifs[motifs$motif_ID == "EVX1_MOUSE.H11MO.0.C",]$TF <- "Evx1" }
  if(sum(grepl("CDX4_MOUSE.H11MO.0.A", motifs$motif_ID))>0){ motifs[motifs$motif_ID == "CDX4_MOUSE.H11MO.0.A",]$TF <- "Cdx4" }
  
  motifs$family <- motif_ann[match(motifs$motif_ID, motif_ann$Model),]$TF.family
  motifs$subfamily <- motif_ann[match(motifs$motif_ID, motif_ann$Model),]$TF.subfamily
  ## label top hits
  motifs$label <- ifelse(motifs$rank < 25, motifs_da$TF, "")
  return(motifs)
}

# cervical
motifs_da.cervical <- read_AME(file = paste0(dir, "RNA+ATAC/results/02_AME_DAregions_stage_cervical.tsv"))
# thoracic
motifs_da.thoracic <- read_AME(file = paste0(dir, "RNA+ATAC/results/02_AME_DAregions_stage_thoracic.tsv"))
# lumbar
motifs_da.lumbar <- read_AME(file = paste0(dir, "RNA+ATAC/results/02_AME_DAregions_stage_lumbar.tsv"))
# sacral
motifs_da.sacral <- read_AME(file = paste0(dir, "RNA+ATAC/results/02_AME_DAregions_stage_sacral.tsv"))
```

The regions with highest accessibility in the lumbar somites show **very** different enrichments to all the rest, with over two thirds of the significant TFs not identified with the other sets of DA regions. A large number of TFs are shared by the other three fates. The cervical and thoracic somites share many motifs, but thoracic regions also have many unique enrichments.


```r
all <- unique(c(motifs_da$motif_ID, motifs_da.cervical$motif_ID, motifs_da.thoracic$motif_ID,
                motifs_da.lumbar$motif_ID, motifs_da.sacral$motif_ID))
compare <- data.frame(motif=all,
                      tf=motif_ann[match(all, motif_ann$Model),]$Transcription.factor)
compare$all <- ifelse(compare$motif %in% motifs_da$motif_ID, 1, 0)
compare$cervical <- ifelse(compare$motif %in% motifs_da.cervical$motif_ID, 1, 0)
compare$thoracic <- ifelse(compare$motif %in% motifs_da.thoracic$motif_ID, 1, 0)
compare$lumbar <- ifelse(compare$motif %in% motifs_da.lumbar$motif_ID, 1, 0)
compare$sacral <- ifelse(compare$motif %in% motifs_da.sacral$motif_ID, 1, 0)
upset(compare[,4:7])
```

![](02_somitogenesisDevelopment_files/figure-html/motif_each_upset-1.png)<!-- -->

Although it is complex to directly compare the results from the different enrichment analyses given the differences in power that stem from the large variation in the number of regions used for each test, the difference in the profile of lumbar somites is definitely a salient feature. This is perhaps explained, to some extent, by the observation that many of the lumbar regions occur in promoters and not as many in distal elements, which is in sharp contrast to the rest of the stages.


### Transcription factor dynamics

To get an insight into potential key regulators of the different processes happening along development, we next focus on the expression of transcription factors specifically. A large number of TFs are differentially expressed between stages, and we can observe similar dynamics as with the rest of the transcriptome, with different groups of regulators showing maximum expression in particular groups of somites.


```r
## list of all mouse TFs (1636)
tfs <- read.table(paste0(dir, "RNA+ATAC/data/Mus_musculus_TF_list.txt"), sep="\t", header = TRUE)
## restrict to those expressed (1310)
tfs <- tfs[tfs$Ensembl %in% row.names(geneCounts),]
## and DE (838; 64% of all expressed)
tfs <- tfs[tfs$Ensembl %in% row.names(degs),]

## expression
data <- geneCounts[tfs$Ensembl,-1]
data <- 2^data-1
data <- t(apply(data, 1, function(x) (x-mean(x))/sd(x) ))

## split by k means
set.seed(942)
clusters.tfs <- kmeans(data, centers=7)
clusters.tfs <- factor(paste0("cluster", clusters.tfs$cluster), 
                       levels = paste0("cluster",c(2,4,7,5,3,1,6)))
names(clusters.tfs) <- row.names(data)
# table(clusters.tfs)

## order columns based on somite number
m <- meta.rna
m$somite <- factor(m$somite, levels=c("SIII","SII","SI"))
stopifnot(identical(colnames(data), m$sample)) # make sure the metadata corresponds with count data matrix
order <- order(m$stage, m$somite)

## heatmap annotation
ha  <- HeatmapAnnotation(df = data.frame(stage = paste0("stage", m[order,]$stage), 
                                         somite = m[order,]$somite), 
                         col = list(stage = cols.stage, somite = cols.somite))

Heatmap(data[,order], 
        cluster_columns = FALSE, 
        col=colorRamp2(breaks = c(-4,-2,0,2,3), 
                       colors = c("steelblue4","steelblue","white","indianred3","indianred4")),
        name = "z-score", 
        show_row_names = FALSE, 
        show_column_names = FALSE, 
        top_annotation = ha, 
        row_split = clusters.tfs,
        cluster_row_slices = FALSE)
```

![](02_somitogenesisDevelopment_files/figure-html/tf_expr-1.png)<!-- -->

When looking at the families of TFs represented in each cluster it stands out that `cluster3` almost exclusively contains TFs from the C2H2 zinc finger family; although this is the most abundant TF family in the mammalian genome and is thus represented in all clusters, the number of these TFs in `cluster3` is much higher than in all other expression patterns.


```r
tfs$cluster <- clusters.tfs[match(tfs$Ensembl, names(clusters.tfs))]

par(mfrow=c(3,3))
for(i in 1:7){
  plot(table(tfs[which(tfs$cluster == paste0("cluster",i)),]$Family),
       ylab = "number of TFs", main = paste0("cluster",i), las=2, cex=0.75)
}
```

![](02_somitogenesisDevelopment_files/figure-html/families-1.png)<!-- -->

We can also observe that some TF families are predominantly found in a particular cluster. These might be important in regulating some of the processes observed at the corresponding stages. For example, GATA factors are involved in blood development and we see them expressed in lumbar somites, consistent with the wider enrichment results seen above.


```r
tmp <- tfs[tfs$Family %in% names(table(tfs$Family)[table(tfs$Family) > 5]),]
Heatmap(prop.table(table(tmp$cluster, tmp$Family),2),
        col = brewer.pal(n=9, "Greys"),
        cluster_rows = FALSE)
```

![](02_somitogenesisDevelopment_files/figure-html/unnamed-chunk-1-1.png)<!-- -->




```r
## save results
write.table(degs, paste0(dir, "RNA+ATAC/results/02_DEgenes_summary_stage_fate.tsv"), 
            quote = FALSE, sep = "\t")
write.table(GO.all, paste0(dir, "RNA+ATAC/results/02_GOenrichment_DEGs_stages.tsv"), 
            quote = FALSE, sep = "\t", row.names = FALSE)
write.table(GO.cervical, paste0(dir, "RNA+ATAC/results/02_GOenrichment_DEGs_stages_cervical.tsv"), 
            quote = FALSE, sep = "\t", row.names = FALSE)
write.table(GO.thoracic, paste0(dir, "RNA+ATAC/results/02_GOenrichment_DEGs_stages_thoracic.tsv"), 
            quote = FALSE, sep = "\t", row.names = FALSE)
write.table(GO.lumbar, paste0(dir, "RNA+ATAC/results/02_GOenrichment_DEGs_stages_lumbar.tsv"), 
            quote = FALSE, sep = "\t", row.names = FALSE)
write.table(GO.sacral, paste0(dir, "RNA+ATAC/results/02_GOenrichment_DEGs_stages_sacral.tsv"), 
            quote = FALSE, sep = "\t", row.names = FALSE)
write.table(tfs, paste0(dir, "RNA+ATAC/results/02_DE_TFs_clusters.tsv"), 
            quote = FALSE, sep = "\t", row.names = FALSE)

write.table(dars, paste0(dir, "RNA+ATAC/results/02_DAregions_summary_stage_fate.tsv"), 
            quote = FALSE, sep = "\t")
write.table(go.dars, paste0(dir, "RNA+ATAC/results/02_GOenrichment_DARs_stages.tsv"),
            quote = FALSE, sep = "\t", row.names = FALSE)
saveRDS(pheno.dars, paste0(dir, "RNA+ATAC/results/02_mousePhenotype_enrichment_DARs_stages.Rds"))
saveRDS(pheno.dars.cervical, paste0(dir, "RNA+ATAC/results/02_mousePhenotype_enrichment_DARs_stages_cervical.Rds"))
saveRDS(pheno.dars.thoracic, paste0(dir, "RNA+ATAC/results/02_mousePhenotype_enrichment_DARs_stages_thoracic.Rds"))
saveRDS(pheno.dars.lumbar, paste0(dir, "RNA+ATAC/results/02_mousePhenotype_enrichment_DARs_stages_lumbar.Rds"))
saveRDS(pheno.dars.sacral, paste0(dir, "RNA+ATAC/results/02_mousePhenotype_enrichment_DARs_stages_sacral.Rds"))
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
##  [1] org.Mm.eg.db_3.13.0         inlmisc_0.5.5              
##  [3] wordcloud2_0.2.1            chromVAR_1.14.0            
##  [5] UpSetR_1.4.0                circlize_0.4.14            
##  [7] ComplexHeatmap_2.8.0        dynamicTreeCut_1.63-1      
##  [9] RColorBrewer_1.1-3          ggrepel_0.9.1              
## [11] ggpubr_0.4.0                ggplot2_3.3.5              
## [13] rGREAT_1.24.0               DESeq2_1.32.0              
## [15] SummarizedExperiment_1.22.0 MatrixGenerics_1.4.3       
## [17] matrixStats_0.62.0          DEFormats_1.20.0           
## [19] GeneTonic_1.4.1             pcaExplorer_2.18.0         
## [21] topGO_2.44.0                SparseM_1.81               
## [23] GO.db_3.13.0                AnnotationDbi_1.54.1       
## [25] Biobase_2.52.0              graph_1.70.0               
## [27] GenomicRanges_1.44.0        GenomeInfoDb_1.28.4        
## [29] IRanges_2.26.0              S4Vectors_0.30.2           
## [31] BiocGenerics_0.38.0        
## 
## loaded via a namespace (and not attached):
##   [1] rappdirs_0.3.3              rtracklayer_1.52.1         
##   [3] AnnotationForge_1.34.1      R.methodsS3_1.8.1          
##   [5] pkgmaker_0.32.2             tidyr_1.2.0                
##   [7] bit64_4.0.5                 knitr_1.38                 
##   [9] R.utils_2.11.0              DelayedArray_0.18.0        
##  [11] data.table_1.14.2           KEGGREST_1.32.0            
##  [13] TFBSTools_1.30.0            RCurl_1.98-1.6             
##  [15] doParallel_1.0.17           generics_0.1.2             
##  [17] cowplot_1.1.1               terra_1.5-21               
##  [19] RSQLite_2.2.12              bit_4.0.4                  
##  [21] tzdb_0.3.0                  webshot_0.5.3              
##  [23] xml2_1.3.3                  httpuv_1.6.5               
##  [25] assertthat_0.2.1            DirichletMultinomial_1.34.0
##  [27] viridis_0.6.2               xfun_0.30                  
##  [29] hms_1.1.1                   jquerylib_0.1.4            
##  [31] evaluate_0.15               promises_1.2.0.1           
##  [33] TSP_1.2-0                   fansi_1.0.3                
##  [35] restfulr_0.0.13             progress_1.2.2             
##  [37] caTools_1.18.2              dendextend_1.15.2          
##  [39] dbplyr_2.1.1                Rgraphviz_2.36.0           
##  [41] igraph_1.3.0                DBI_1.1.2                  
##  [43] geneplotter_1.70.0          htmlwidgets_1.5.4          
##  [45] purrr_0.3.4                 ellipsis_0.3.2             
##  [47] crosstalk_1.2.0             dplyr_1.0.8                
##  [49] backports_1.4.1             annotate_1.70.0            
##  [51] gridBase_0.4-7              biomaRt_2.48.3             
##  [53] vctrs_0.4.1                 Cairo_1.5-15               
##  [55] abind_1.4-5                 cachem_1.0.6               
##  [57] withr_2.5.0                 ggforce_0.3.3              
##  [59] BSgenome_1.60.0             rgdal_1.5-30               
##  [61] checkmate_2.0.0             GenomicAlignments_1.28.0   
##  [63] prettyunits_1.1.1           cluster_2.1.3              
##  [65] backbone_2.0.3              lazyeval_0.2.2             
##  [67] seqLogo_1.58.0              crayon_1.5.1               
##  [69] genefilter_1.74.1           labeling_0.4.2             
##  [71] edgeR_3.34.1                pkgconfig_2.0.3            
##  [73] tweenr_1.0.2                seriation_1.3.5            
##  [75] rlang_1.0.2                 lifecycle_1.0.1            
##  [77] miniUI_0.1.1.1              colourpicker_1.1.1         
##  [79] registry_0.5-1              filelock_1.0.2             
##  [81] BiocFileCache_2.0.0         GOstats_2.58.0             
##  [83] polyclip_1.10-0             rngtools_1.5.2             
##  [85] Matrix_1.4-1                raster_3.5-15              
##  [87] carData_3.0-5               base64enc_0.1-3            
##  [89] GlobalOptions_0.1.2         pheatmap_1.0.12            
##  [91] png_0.1-7                   viridisLite_0.4.0          
##  [93] rjson_0.2.21                bitops_1.0-7               
##  [95] shinydashboard_0.7.2        R.oo_1.24.0                
##  [97] visNetwork_2.1.0            Biostrings_2.60.2          
##  [99] blob_1.2.3                  shape_1.4.6                
## [101] rintrojs_0.3.0              stringr_1.4.0              
## [103] readr_2.1.2                 rstatix_0.7.0              
## [105] shinyAce_0.4.1              ggsignif_0.6.3             
## [107] CNEr_1.28.0                 scales_1.2.0               
## [109] memoise_2.0.1               GSEABase_1.54.0            
## [111] magrittr_2.0.3              plyr_1.8.7                 
## [113] zlibbioc_1.38.0             threejs_0.3.3              
## [115] compiler_4.1.0              BiocIO_1.2.0               
## [117] clue_0.3-60                 Rsamtools_2.8.0            
## [119] cli_3.2.0                   XVector_0.32.0             
## [121] Category_2.58.0             MASS_7.3-56                
## [123] tidyselect_1.1.2            stringi_1.7.6              
## [125] shinyBS_0.61.1              highr_0.9                  
## [127] yaml_2.3.5                  locfit_1.5-9.5             
## [129] sass_0.4.1                  tools_4.1.0                
## [131] rstudioapi_0.13             TFMPvalue_0.0.8            
## [133] foreach_1.5.2               gridExtra_2.3              
## [135] farver_2.1.0                digest_0.6.29              
## [137] pracma_2.3.8                shiny_1.7.1                
## [139] Rcpp_1.0.8.3                car_3.0-12                 
## [141] broom_0.8.0                 later_1.3.0                
## [143] shinyWidgets_0.6.4          httr_1.4.2                 
## [145] colorspace_2.0-3            XML_3.99-0.9               
## [147] splines_4.1.0               RBGL_1.68.0                
## [149] expm_0.999-6                sp_1.4-6                   
## [151] plotly_4.10.0               xtable_1.8-4               
## [153] jsonlite_1.8.0              poweRlaw_0.70.6            
## [155] heatmaply_1.3.0             R6_2.5.1                   
## [157] pillar_1.7.0                htmltools_0.5.2            
## [159] mime_0.12                   NMF_0.24.0                 
## [161] glue_1.6.2                  fastmap_1.1.0              
## [163] DT_0.22                     BiocParallel_1.26.2        
## [165] bs4Dash_2.0.3               codetools_0.2-18           
## [167] utf8_1.2.2                  lattice_0.20-45            
## [169] bslib_0.3.1                 tibble_3.1.6               
## [171] curl_4.3.2                  gtools_3.9.2               
## [173] magick_2.7.3                survival_3.3-1             
## [175] limma_3.48.3                rmarkdown_2.13             
## [177] munsell_0.5.0               GetoptLong_1.0.5           
## [179] GenomeInfoDbData_1.2.6      iterators_1.0.14           
## [181] reshape2_1.4.4              gtable_0.3.0               
## [183] shinycssloaders_1.0.0
```

