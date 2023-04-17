### Code to reproduce all figures
library(SummarizedExperiment)
library(ComplexHeatmap)
library(circlize)
library(GeneTonic)
library(chromVAR)
library(csaw)
library(edgeR)
library(DESeq2)
library(limma)
library(zoo)
library(UpSetR)
library(ggplot2)
library(ggpubr)
library(ggrepel)
library(ggrastr)
library(RColorBrewer)
library(colorspace)
library(inlmisc)

dir <- "/Users/ibarra01/OneDrive - CRUK Cambridge Institute/github/somitogenesis2020/"
out <- "/Users/ibarra01/OneDrive - CRUK Cambridge Institute/github/somitogenesis2020/Figures/plot_test/"

palette(brewer.pal(n=12, "Set3")[-c(1:2)])

th <- theme_bw() + theme(axis.text.x = element_text(size=12), 
                         axis.title.x = element_text(size=12), 
                         axis.text.y = element_text(size=12), 
                         axis.title.y = element_text(size=12), 
                         panel.grid.major = element_blank(), 
                         panel.grid.minor = element_blank(), 
                         axis.line = element_line(colour = "black"), 
                         panel.border = element_blank(), 
                         plot.title = element_text(face="bold", hjust = 0.5))

## cluster colours
cols.stage <- c("darkolivegreen4", brewer.pal(n=9, "Blues")[c(4,6,7)], "darkorange", "indianred1")
names(cols.stage) <- paste0("stage", c(8,18,21,25,27,35))
cols.somite <- hcl.colors(n=5, palette = "PurpOr")[3:1]
names(cols.somite) <- paste0("S", c("I","II","III"))

# Figure 1 ###############

## C: HOX gene expression heatmap ======
# normalised, batch-corrected data
geneCounts <- read.table(paste0(dir, "RNA-seq/data/geneCounts.NORM_batchCorrected_14PCs.tsv"),
                         check.names = FALSE)
# metadata
meta.rna <- read.table(paste0(dir, "RNA-seq/data/metadata_RNAseq.tsv"), 
                       stringsAsFactors = FALSE, header = TRUE)
meta.rna$stage <- factor(paste0("stage", meta.rna$stage), levels = paste0("stage", c(8,18,21,25,27,35)))
meta.rna$somite <- factor(meta.rna$somite, levels=paste0("S", c("III", "II", "I")))
# keep only good quality samples
meta.rna <- meta.rna[meta.rna$sample %in% colnames(geneCounts),]
stopifnot(identical(colnames(geneCounts)[-1], meta.rna$sample))
## remove samples from stage27 that are likely mis-staged
meta.rna <- meta.rna[meta.rna$QC == 1 & meta.rna$wrongStage == 0,]
geneCounts <- geneCounts[,which(colnames(geneCounts) %in% c("gene", meta.rna$sample)),]


# HOX genes
hox <- geneCounts[grep("^Hox",geneCounts$gene),]
hox <- hox[grep("s", hox$gene, invert=TRUE),] # remove antisense and non-coding genes
row.names(hox) <- hox$gene
hox$gene <- NULL

# visualise with z-scores
hox.cnt <- t(apply(hox, 1, function(x) (x-mean(x))/sd(x) ))

# paralogous group
hox.num <- as.numeric(substr(row.names(hox.cnt), 5,6))
names(hox.num) <- row.names(hox.cnt)

# order
order.col <- order(meta.rna$stage)
order.row <- order(hox.num)

# sample annotation
ha.sample <- HeatmapAnnotation(df = data.frame(stage = meta.rna[order.col,]$stage), 
                               col = list(stage = cols.stage))
stages <- factor(meta.rna[order.col,]$stage, levels=paste0("stage",c(8,18,21,25,27,35)))

pdf(paste0(out, "Figure1C_Hox_heatmap.pdf"), useDingbats = FALSE, width = 5, height = 5)
draw(
  Heatmap(hox.cnt[order.row, order.col], 
          cluster_columns = FALSE, 
          cluster_rows = FALSE,
          col=colorRamp2(breaks = c(-3,0,3), colors = c("steelblue","aliceblue","indianred3")), 
          name = "z-score", 
          show_row_names = FALSE, 
          show_column_names = FALSE, 
          top_annotation = ha.sample,
          row_title = NULL,
          column_title = NULL,
          column_split = stages, 
          cluster_column_slices = FALSE,
          row_split = hox.num[order.row],
          border_gp = gpar(col = "grey20", lwd = 0.5), 
          heatmap_legend_param = list(legend_direction = "horizontal")),
  heatmap_legend_side = "bottom")
dev.off()

## D: HOX gene expression PCA ======
# PCA
pca <- prcomp(t(hox))
# variance
eigs <- pca$sdev^2
prop.var <- eigs / sum(eigs) * 100
# data to plot
df <- as.data.frame(pca$x)

pdf(paste0(out, "Figure1D_Hox_PCA.pdf"), useDingbats = FALSE, width = 5.5, height = 4.5)
ggplot(df, aes(PC1, PC2, colour=meta.rna$stage, shape=meta.rna$somite)) + 
  geom_point(cex=2.25) + 
  xlab(paste("PC1 -", round(prop.var[1],2), "%")) + 
  ylab(paste("PC2 -", round(prop.var[2],2), "%")) + 
  labs(colour="stage", shape="somite") + 
  scale_color_manual(values = cols.stage) + 
  th
dev.off()


## E: peak annotation ======
peaks.classes <- read.table(paste0(dir, "ATAC-seq/results/05_peaks_classAnnotation.tsv"), sep="\t")
peaks.classes <- colSums(peaks.classes[,7:12])
df <- data.frame(class = names(peaks.classes),
                 number = peaks.classes, 
                 prop = peaks.classes/sum(peaks.classes)*100)
# relevel to get desired order
df$class <- factor(df$class, levels = rev(c("promoter", "exonic", "intronic", "proximal", "distal", "intergenic")))

pdf(paste0(out, "Figure1E_peakAnnotation.pdf"), useDingbats = FALSE, width = 6, height = 2)
ggplot(df, aes(x = 1, y = prop, fill=class, label=round(prop))) + 
  geom_bar(position="stack", stat="identity") +
  scale_fill_manual(values = c(brewer.pal(n=11, "PiYG")[7:10], "gold2", "darkgoldenrod2")) +
  geom_text(size = 3, position = position_stack(vjust = 0.5)) +
  xlab("") + 
  ylab("proportion of peaks") +
  coord_flip() +
  th + theme(axis.text.y = element_blank(),
             axis.ticks.y = element_blank(),
             legend.position = "bottom")
dev.off()



# Figure 2 ###############

## B: MA pplot SI vs SIII ======
somiteIIIvsI.all <- read.table(paste0(dir, "RNA-seq/results/04_DEresults_somiteTrios_somiteIvsIII.tsv"))
# reorder to plot DE on top
somiteIIIvsI.all <- somiteIIIvsI.all[order(somiteIIIvsI.all$DE),]
somiteIIIvsI.all <- somiteIIIvsI.all[order(somiteIIIvsI.all$FDR),]
somiteIIIvsI.all <- rbind(somiteIIIvsI.all[somiteIIIvsI.all$DE == "NS",], 
                          somiteIIIvsI.all[somiteIIIvsI.all$DE != "NS",])
# label largest changes after excluding very lowly expressed genes
tmp <- somiteIIIvsI.all[order(somiteIIIvsI.all$logFC),]
labels <- tmp[tmp$logCPM>2 & tmp$DE != "NS",]$genes[1:5]
tmp <- somiteIIIvsI.all[order(-somiteIIIvsI.all$logFC),]
labels <- c(labels, tmp[tmp$logCPM>2 & tmp$DE != "NS",]$genes[1:5])
labels <- c(labels, "Col2a1")
somiteIIIvsI.all$label <- ifelse(somiteIIIvsI.all$genes %in% labels, somiteIIIvsI.all$genes, "")

# remove overlapping points for easier figure editing
keep <- iSEE::subsetPointsByGrid(somiteIIIvsI.all$logFC, -log10(somiteIIIvsI.all$FDR), resolution = 300)

pdf(paste0(out, "Figure2B_volcanoPlot.pdf"), useDingbats = FALSE, width = 5, height = 5)
ggplot(somiteIIIvsI.all, aes(-logFC, -log10(FDR), col=as.factor(DE), label=label)) + 
  geom_point(aes(size=logCPM), alpha=0.75) +
  scale_size(range = c(0.5, 3)) +
  scale_color_manual(values = c(down="#CC78AF", NS="grey", up="#583C88")) +
  geom_hline(yintercept = -log10(0.05), lty=2, colour="grey40") +
  geom_text_repel(max.overlaps = 50, size = 2, nudge_x = 1, nudge_y = 0.25, show.legend = FALSE) +
  xlab(expression('log'[2]*' fold-change')) + 
  ylab(expression('-log'[10]*' adjusted p-value')) + 
  xlim(c(-4.25,4.5)) +
  th + theme(legend.position = "none",
             panel.background = element_rect(fill = "transparent", colour = NA),  
             plot.background = element_rect(fill = "transparent", colour = NA))
dev.off()

## C: examples ======
# plot by somite number
df <- meta.rna
# differentiate between the two somites 25 from different stages
df$group <- factor(paste(df$somiteNumber, df$somite, sep="|"))
df$group <- factor(df$group, levels = levels(df$group)[c(16:18, 1:15)])

plots <- list()
df$expr <- as.numeric(geneCounts[geneCounts$gene=="Tbx22",-1])
plots[[1]] <- ggplot(df, aes(group, expr, fill=somite)) +
  geom_boxplot(alpha=0.5) +
  scale_fill_manual(values = cols.somite) +
  geom_vline(xintercept = seq(3.5,15.5,3), lty=3, col="grey") +
  ggtitle("Tbx22") +
  xlab("somite number") +
  ylab(expression('log'[2]*' CPM')) + 
  scale_x_discrete(labels=c(6:8,16:21,23:25,25:27,33:35)) +
  th + theme(plot.title = element_text(face="italic", hjust = 0.5),
             legend.position = "bottom")
df$expr <- as.numeric(geneCounts[geneCounts$gene=="Col2a1",-1])
plots[[2]] <- ggplot(df, aes(group, expr, fill=somite)) +
  geom_boxplot(alpha=0.5) +
  scale_fill_manual(values = cols.somite) +
  geom_vline(xintercept = seq(3.5,15.5,3), lty=3, col="grey") +
  ggtitle("Col2a1") +
  xlab("somite number") +
  ylab(expression('log'[2]*' CPM')) + 
  scale_x_discrete(labels=c(6:8,16:21,23:25,25:27,33:35)) +
  th + theme(plot.title = element_text(face="italic", hjust = 0.5),
             legend.position = "bottom")

pdf(paste0(out, "Figure2C_geneExpr_examples.pdf"), useDingbats = FALSE, width = 5, height = 5)
ggarrange(plotlist = plots, nrow=2, align="hv", common.legend = TRUE)
dev.off()

## D: GO enrichment ======
load(paste0(dir, "RNA+ATAC/results/01_geneTonic_objects.RData"))
# define terms to highlight
terms <- c("GO:0009952", "GO:0003002", "GO:0001756", "GO:0030198", "GO:0048706",
           "GO:0007155", "GO:0051216", "GO:2000050", "GO:0048745", "GO:0070373",
           "GO:0050673", "GO:0030857", "GO:0035023", "GO:0007219", "GO:0030336", "GO:0030335")
GO.somite_res$highlight <- GO.somite_res$gs_id %in% terms

# annotate which terms are also significant from DA regions
GO.somite <- read.table(paste0(dir, "RNA+ATAC/results/01_GO_enrichment_DEGs_somiteTrios.tsv"), sep="\t", header = TRUE)
GO.somite_res$shared <- GO.somite[match(GO.somite_res$gs_id, GO.somite$GO.ID),]$sig_great

# scatter plot
p <- ggplot(GO.somite_res[order(GO.somite_res$shared),], aes(-z_score, -log10(gs_pvalue),
                                                        size = DE_count, label=gs_description)) +
  geom_point(shape=21, aes(colour=shared, fill = -aggr_score)) + 
  scale_color_manual(values = c("grey90", "grey20")) +
  scale_fill_gradientn(colours = brewer.pal(n=11, "RdYlBu"),
                       limits=c(-1,1), labels = c("< -1", seq(-0.5,0.5,0.5), "> 1")) +
  xlab("z score") +
  ylab("-log10 p-value") +
  geom_vline(xintercept = 0, lty=2, colour="grey70") +
  scale_size(range=c(0,10), breaks=c(30,50,200,500)) +
  th 
pdf(paste0(out, "Figure2D_GOresults.pdf"), useDingbats = FALSE, width = 7, height = 3.5)
p
dev.off()

# add labels to plot
p <- p + geom_point(aes(colour=highlight)) + 
  geom_text_repel(data = GO.somite_res[GO.somite_res$gs_id %in% terms, ],
                         colour = "grey50",
                         min.segment.length = unit(0, 'lines'),
                         nudge_y = 0.3,
                         size=3, max.overlaps = 200)
pdf(paste0(out, "Figure2D_GOresults_labels.pdf"), useDingbats = FALSE, width = 7, height = 3.5)
p
dev.off()


## E: DA genomic locations ======
peakAnn <- read.table(paste0(dir, "ATAC-seq/results/05_peaks_classAnnotation.tsv"), sep="\t")
dars <- read.table(paste0(dir, "ATAC-seq/results/04_DAregions_somiteTrios.tsv"), 
                   stringsAsFactors = FALSE)
dars.ann <- peakAnn[substr(row.names(dars), 4, 30),]

df <- data.frame(all = round(colSums(peakAnn[,7:12])/nrow(peakAnn)*100, 2),
                 DA = round(colSums(dars.ann[,7:12])/nrow(dars.ann)*100, 2))
df <- data.frame(class=rep(row.names(df),2),
                 DA = c(rep(FALSE, 6), rep(TRUE, 6)),
                 prop = c(df$all, df$DA))
df$class <- factor(df$class, levels=rev(c("promoter", "exonic", "intronic", "proximal", "distal", "intergenic")))

pdf(paste0(out, "Figure2E_DA_locations.pdf"), useDingbats = FALSE, width = 3, height = 5)
ggplot(df, aes(DA, prop, fill=class, label=round(prop))) + 
  geom_bar(stat="identity", colour="black", width=0.5) +
  scale_fill_manual(values = c(brewer.pal(n=11, "PiYG")[7:10], "gold2", "darkgoldenrod2")) +
  geom_text(size = 3, position = position_stack(vjust = 0.5)) +
  ylab("proportion of all peaks") +
  th
dev.off()


## F: TFBS enrichment ======
motif_dyn_pct <- read.table(paste0(dir, "RNA+ATAC/results/01_motifEnrichment_somiteTrios.tsv"))
motif_dyn_pct$TF <- toupper(motif_dyn_pct$TF)

pdf(paste0(out, "Figure2F_motif_enrichment.pdf"), useDingbats = FALSE, width = 6, height = 2.5)
ggplot(motif_dyn_pct, aes(-z, da/nonDa)) +
  geom_point(pch=21, colour="grey40") +
  geom_point(data=motif_dyn_pct[abs(motif_dyn_pct$z)>0.3,], aes(fill=z>0),
             pch=21, alpha=0.5, colour="grey40") +
  scale_fill_manual(values = unname(cols.somite[c(1,3)])) +
  geom_vline(xintercept = 0, lty=2, colour="grey80") +
  geom_text_repel(data=motif_dyn_pct[motif_dyn_pct$da/motif_dyn_pct$nonDa>1.6 | 
                                       abs(motif_dyn_pct$z)>0.3,], 
                  aes(label=TF), colour="grey45", size=2) +
  xlab("z score") +
  ylab("peaks with TFBS\nDA / nonDA") +
  # xlim(-0.9, 0.5) +
  ylim(c(1,2)) +
  th + theme(legend.position = "none")
dev.off()


## G: motif activity ======
dev <- readRDS(paste0(dir, "ATAC-seq/results/06_motif_deviationScores.Rds"))
df <- as.data.frame(colData(dev))
df$somite <- factor(df$somite, levels = c("SIII", "SII", "SI"))

plots <- list()
for(tf in c("Hoxb3", "Hoxa7", "Hoxd13", "Cdx2", "Evx2", 
            "Meis1", "Rxra", "Nr2c1", "Nr2f2", "Zbtb12")){
  i <- grep(tf, row.names(dev))
  df$z <- assay(dev, 'z')[i,]
  
  plots[[tf]] <- ggplot(df, aes(somite, z, fill=somite)) +
    geom_boxplot(alpha=0.5) +
    scale_fill_manual(values = cols.somite) +
    geom_hline(yintercept = 0, lty=2, colour="grey") +
    ggtitle(tf) +
    xlab("") +
    th + theme(plot.title = element_text(face="italic", hjust = 0.5))
}

pdf(paste0(out, "Figure2G_motif_activity.pdf"), useDingbats = FALSE, width = 10, height = 5)
ggarrange(plotlist = plots, ncol=5, nrow=2, legend = "none", align = "hv")
dev.off()


# Figure 3 ###############

## B: DE genes heatmap ======
data.stage.degs <- readRDS(paste0(dir, "RNA+ATAC/results/02_stage_DEGs_data_for_heatmap.Rds"))
degs <- read.table(paste0(dir, "RNA+ATAC/results/02_DEgenes_summary_stage_fate.tsv"), header=TRUE)
stopifnot(identical(row.names(degs), row.names(data.stage.degs)))

# annotate somite and stage
stopifnot(identical(colnames(data.stage.degs), meta.rna$sample)) # make sure the metadata corresponds with count data matrix
order <- order(meta.rna$stage, meta.rna$somite)
ha  <- HeatmapAnnotation(df = data.frame(stage = meta.rna[order,]$stage, 
                                         somite = meta.rna[order,]$somite), 
                         col = list(stage = cols.stage, somite = cols.somite),
                         simple_anno_size = unit(0.25, "cm"),
                         show_annotation_name = FALSE,
                         annotation_legend_param = list(stage = list(nrow=2),
                                                        somite = list(nrow = 1)))
# group by clusters
clusters <- factor(degs$cluster, levels = paste0("cluster",c(8,2,6,7, 4,5, 3, 11,12,9,1, 10)))
fates <- c("cervical", rep("thoracic",3), "lumbar", "sacral")
names(fates) <- paste0("stage",c(8,18,21,25,27,35))
fates <- factor(fates[meta.rna[order,]$stage], levels=c("cervical", "thoracic", "lumbar", "sacral"))

# label genes of interest
ids <- row.names(geneCounts)[geneCounts$gene %in% c("Rdh10", "Aldh1a2", "Rara", "Rxra", # RA signalling
                                                    paste0("Bmp", c(1,3,4,7)), paste0("Smad", c(3,5)),
                                                    "Sox9", "Mef2c", "Gata2", "Tal1", "Ets1", "Tbx20")]
ha.ids <- rowAnnotation(genes = anno_mark(at = which(row.names(data.stage.degs) %in% ids), 
                                          labels = geneCounts[row.names(data.stage.degs[which(row.names(data.stage.degs) %in% ids),]),1]))

# mark cholesterol pathway genes
GO.thoracic <- read.table(paste0(dir, "RNA+ATAC/results/02_GOenrichment_DEGs_stages_thoracic.tsv"), header = TRUE, sep="\t")
genes <- unlist(strsplit(GO.thoracic[GO.thoracic$GO.ID=="GO:0006695",]$genes, ",", fixed=TRUE))
ids <- row.names(geneCounts[geneCounts$gene %in% genes,])

cols <- c("red", "white")
names(cols) <- c("cholesterol", "other")
ha.chol <-  rowAnnotation(df = data.frame(cholesterol = ifelse(row.names(data.stage.degs) %in% ids, "cholesterol", "other")),
                          col = list(cholesterol = cols),
                          show_annotation_name = FALSE,
                          show_legend = FALSE)

pdf(paste0(out, "Figure3B_DEGs_heatmap.pdf"), useDingbats = FALSE, width = 6, height = 8)
draw(
  Heatmap(data.stage.degs[,order], 
          cluster_columns = FALSE, 
          col=colorRamp2(breaks = c(-3,-2,-log2(1.5),0,log2(1.5),2,3), 
                         colors = c("steelblue4","steelblue",rep("white",3),"indianred3","indianred4")),
          show_row_names = FALSE, 
          show_column_names = FALSE, 
          show_row_dend = FALSE,
          row_title = NULL,
          column_title = NULL,
          left_annotation = ha.chol,
          top_annotation = ha, 
          column_split = fates,
          row_split = clusters,
          cluster_row_slices = FALSE,
          right_annotation = ha.ids,
          border_gp = gpar(col = "grey20", lwd = 0.5), 
          heatmap_legend_param = list(title = "z-score", direction = "horizontal"),
          use_raster = TRUE),
  heatmap_legend_side = "bottom",
  annotation_legend_side = "bottom",
  merge_legend = TRUE)
dev.off()

## C: GO enrichment ======
GO.all <- read.table(paste0(dir, "RNA+ATAC/results/02_GOenrichment_DEGs_stages.tsv"), header = TRUE, sep="\t")
GO.cervical <- read.table(paste0(dir, "RNA+ATAC/results/02_GOenrichment_DEGs_stages_cervical.tsv"), header = TRUE, sep="\t")
GO.thoracic <- read.table(paste0(dir, "RNA+ATAC/results/02_GOenrichment_DEGs_stages_thoracic.tsv"), header = TRUE, sep="\t")
GO.lumbar <- read.table(paste0(dir, "RNA+ATAC/results/02_GOenrichment_DEGs_stages_lumbar.tsv"), header = TRUE, sep="\t",
                        quote='', as.is = TRUE)
GO.sacral <- read.table(paste0(dir, "RNA+ATAC/results/02_GOenrichment_DEGs_stages_sacral.tsv"), header = TRUE, sep="\t")



pvals <- data.frame(term = unique(c(GO.cervical[c(1,11,30),]$Term,
                                    GO.thoracic[c(1,2,7,14,45),]$Term,
                                    GO.lumbar[c(1,5,21,22),]$Term,
                                    GO.sacral[c(1,33,70,177),]$Term)))
pvals$cervical <- GO.cervical[match(pvals$term, GO.cervical$Term),]$p.value_elim
pvals$thoracic <- GO.thoracic[match(pvals$term, GO.thoracic$Term),]$p.value_elim
pvals$lumbar <- GO.lumbar[match(pvals$term, GO.lumbar$Term),]$p.value_elim
pvals$sacral <- GO.sacral[match(pvals$term, GO.sacral$Term),]$p.value_elim
pvals[is.na(pvals)] <- 1

sizes <- data.frame(term = unique(c(GO.cervical[c(1,11,30),]$Term,
                                    GO.thoracic[c(1,2,7,14,45),]$Term,
                                    GO.lumbar[c(1,5,21,22),]$Term,
                                    GO.sacral[c(1,33,70,177),]$Term)))
sizes$cervical <- GO.cervical[match(sizes$term, GO.cervical$Term),]$Significant
sizes$thoracic <- GO.thoracic[match(sizes$term, GO.thoracic$Term),]$Significant
sizes$lumbar <- GO.lumbar[match(sizes$term, GO.lumbar$Term),]$Significant
sizes$sacral <- GO.sacral[match(sizes$term, GO.sacral$Term),]$Significant
sizes[is.na(sizes)] <- 0

hc <- hclust(dist(pvals[,-1]))
levels <- pvals$term[hc$order]

df <- reshape2::melt(pvals)
tmp <- reshape2::melt(sizes)
df$n_sig <- tmp$value

df$term <- factor(df$term, levels = levels[c(3:6,2:1,7:9,11:10,15:16,14:12)])

pdf(paste0(out, "Figure3C_GOenrichment.pdf"), useDingbats = FALSE, width = 7, height = 6)
ggplot(df, aes(term, variable, colour=-log10(value), size=log10(n_sig))) +
  geom_point() +
  scale_color_gradientn(colours = brewer.pal(n=9, "Purples")[-1]) +
  xlab("") +
  ylab("") +
  labs(colour = "-log10 p-value",
       size = "log10 # DEGs") +
  scale_y_discrete(position = "right") +
  th + theme(axis.text.x = element_text(size=10, angle=90, hjust=1, vjust=0.5),
             axis.text.y = element_text(size=10),
             legend.position = "bottom",
             legend.justification='left')
dev.off()

## D: RA signalling ======
# gene expression
df <- meta.rna
df$somite <- factor(df$somite, levels = c("SI", "SII", "SIII"))
df$expr <- as.numeric(geneCounts[geneCounts$gene=="Rara",-1])

pdf(paste0(out, "Figure3D_Rara.pdf"), useDingbats = FALSE, width = 3, height = 4)
ggplot(df, aes(stage, expr, fill=stage)) +
  geom_boxplot(alpha=0.75) +
  scale_fill_manual(values = cols.stage) +
  ggtitle("Rara") +
  xlab("") +
  ylab(expression('log'[2]*' CPM')) + 
  th + theme(axis.text.x = element_blank(),
             axis.ticks.x = element_blank(),
             plot.title = element_text(face="italic", hjust = 0.5),
             legend.position = "none")
dev.off()

## E: RA signalling ======
# TF activity
df <- as.data.frame(colData(dev))
df$somite <- factor(df$somite, levels = c("SI", "SII", "SIII"))
plots <- list()
for(tf in c("Rara", "Rxra")){
  i <- grep(tf, row.names(dev))
  df$z <- assay(dev, 'z')[i,]
  
  plots[[tf]] <- ggplot(df, aes(stage, z, fill=stage)) +
    geom_boxplot(alpha=0.75) +
    scale_fill_manual(values = cols.stage) +
    geom_hline(yintercept = 0, lty=2, colour="grey") +
    ggtitle(toupper(tf)) +
    xlab("") +
    th + theme(axis.text.x = element_blank(),
               axis.ticks.x = element_blank(),
               plot.title = element_text(hjust = 0.5),
               legend.position = "none")
}

pdf(paste0(out, "Figure3E_RAreceptors_dev_scores.pdf"), useDingbats = FALSE, width = 6, height = 4)
ggarrange(plotlist = plots, ncol=2, legend = "none", align="hv")
dev.off()


## F: cholesterol ======
genes <- unlist(strsplit(GO.thoracic[GO.thoracic$GO.ID=="GO:0006695",]$genes, ",", fixed=TRUE))
# highlight the genes in this GO, except for one outlier that clusters separately in the heatmap
data <- geneCounts[geneCounts$gene %in% setdiff(genes, "Srebf1"),]
row.names(data) <- data$gene
data$gene <- NULL
data <- t(apply(data, 1, function(x) (x-mean(x))/sd(x) ))

stopifnot(identical(colnames(data), meta.rna$sample))
order <- order(meta.rna$stage, meta.rna$somite)

## heatmap annotation
ha  <- HeatmapAnnotation(df = data.frame(stage = meta.rna[order,]$stage),
                         col = list(stage = cols.stage),
                         simple_anno_size = unit(0.25, "cm"))

pdf(paste0(out, "Figure3F_cholesterol_biosynthesis.pdf"), useDingbats = FALSE, width = 10, height = 2)
Heatmap(data[,order],
        cluster_columns = FALSE,
        cluster_rows = FALSE,
        col = colorRamp2(breaks = c(-3,-2,0,2,3), 
                         colors = c("steelblue4","steelblue","white","indianred3","indianred4")),
        row_names_gp = gpar(fontsize=8, fontface="italic"),
        column_names_gp = gpar(fontsize=10),
        column_split = fates,
        border_gp = gpar(col = "grey20", lwd = 0.5), 
        name = "z-score",
        show_row_names = TRUE,
        show_column_names = FALSE,
        top_annotation = ha)
dev.off()


## G: FRiP ======
# ATAC-seq metadata
meta.atac <- read.table(paste0(dir, "ATAC-seq/data/metadata_ATACseq_GQ.tsv"), stringsAsFactors = FALSE, header = TRUE)
meta.atac <- meta.atac[meta.atac$QCpass==1,]
meta.atac$stage <- factor(paste0("stage", meta.atac$stage), levels = paste0("stage", c(8,18,21,25,27,35)))
meta.atac$somite <- factor(meta.atac$somite, levels=paste0("S", c("III", "II", "I")))
meta.atac$FRiP <- meta.atac$readsInPeakSet/meta.atac$goodQuality

pdf(paste0(out, "Figure3G_FRiP.pdf"), useDingbats = FALSE, width = 3, height = 4)
ggplot(meta.atac, aes(stage, FRiP, fill=stage)) +
  geom_boxplot(alpha=0.75) +
  scale_fill_manual(values = cols.stage) +
  xlab("") +
  ylab("fraction of reads in peaks") +
  th + theme(axis.text.x = element_blank(),
             axis.ticks.x = element_blank(),
             legend.position = "none")
dev.off()

## H: FOXC1 ======
# TF activity
df <- as.data.frame(colData(dev))
i <- grep("Foxc1", row.names(dev))
df$z <- assay(dev, 'z')[i,]
  
pdf(paste0(out, "Figure3H_Foxc1_dev_scores.pdf"), useDingbats = FALSE, width = 3, height = 4)
ggplot(df, aes(stage, z, fill=stage)) +
  geom_boxplot(alpha=0.75) +
  scale_fill_manual(values = cols.stage) +
  geom_hline(yintercept = 0, lty=2, colour="grey") +
  ggtitle(toupper("Foxc1")) +
  xlab("") +
  th + theme(axis.text.x = element_blank(),
             axis.ticks.x = element_blank(),
             plot.title = element_text(hjust = 0.5),
             legend.position = "none")
dev.off()

## I: Myod1 ======
df <- as.data.frame(colData(dev))
i <- grep("Myod", row.names(dev))
df$z <- assay(dev, 'z')[i,]

pdf(paste0(out, "Figure3I_Myod1_devScore.pdf"), useDingbats = FALSE, width = 3, height = 4)
ggplot(df, aes(stage, z, fill=stage)) +
  geom_boxplot(alpha=0.75) +
  scale_fill_manual(values = cols.stage) +
  geom_hline(yintercept = 0, lty=2, colour="grey") +
  ggtitle("MYOD1") +
  xlab("") +
  th + theme(axis.text.x = element_blank(),
             axis.ticks.x = element_blank(),
             plot.title = element_text(hjust = 0.5),
             legend.position = "none")
dev.off()

## J-K: Smad5 ======
# gene expression
df <- meta.rna
df$expr <- as.numeric(geneCounts[geneCounts$gene=="Smad5",-1])

plots <- list()
plots[[1]] <- ggplot(df, aes(stage, expr, fill=stage)) +
  geom_boxplot(alpha=0.75) +
  scale_fill_manual(values = cols.stage) +
  ggtitle("Smad5") +
  xlab("") +
  ylab(expression('log'[2]*' CPM')) + 
  th + theme(axis.text.x = element_blank(),
             axis.ticks.x = element_blank(),
             plot.title = element_text(face="italic", hjust = 0.5),
             legend.position = "none")

# TF activity
df <- as.data.frame(colData(dev))
i <- grep("Smad5", row.names(dev))
df$z <- assay(dev, 'z')[i,]
  
plots[[2]] <- ggplot(df, aes(stage, z, fill=stage)) +
  geom_boxplot(alpha=0.75) +
  scale_fill_manual(values = cols.stage) +
  geom_hline(yintercept = 0, lty=2, colour="grey") +
  ggtitle("SMAD5") +
  xlab("") +
  th + theme(axis.text.x = element_blank(),
             axis.ticks.x = element_blank(),
             plot.title = element_text(hjust = 0.5),
             legend.position = "none")

pdf(paste0(out, "Figure3J-K_Smad5_expr_devScore.pdf"), useDingbats = FALSE, width = 6, height = 4)
ggarrange(plotlist = plots, legend = "none", align="hv")
dev.off()


## L: Mef2c ======
# gene expression
df <- meta.rna
df$expr <- as.numeric(geneCounts[geneCounts$gene=="Mef2c",-1])

pdf(paste0(out, "Figure3L_Mef2c_expr.pdf"), useDingbats = FALSE, width = 3, height = 4)
ggplot(df, aes(stage, expr, fill=stage)) +
  geom_boxplot(alpha=0.75) +
  scale_fill_manual(values = cols.stage) +
  ggtitle("Mef2c") +
  xlab("") +
  ylab(expression('log'[2]*' CPM')) + 
  th + theme(axis.text.x = element_blank(),
             axis.ticks.x = element_blank(),
             plot.title = element_text(face="italic", hjust = 0.5),
             legend.position = "none")
dev.off()


## M: KO enrichment ======
pheno.dars.cervical <- readRDS(paste0(dir, "RNA+ATAC/results/02_mousePhenotype_enrichment_DARs_stages_cervical.Rds"))
pheno.dars.thoracic <- readRDS(paste0(dir, "RNA+ATAC/results/02_mousePhenotype_enrichment_DARs_stages_thoracic.Rds"))
pheno.dars.lumbar <- readRDS(paste0(dir, "RNA+ATAC/results/02_mousePhenotype_enrichment_DARs_stages_lumbar.Rds"))
pheno.dars.sacral <- readRDS(paste0(dir, "RNA+ATAC/results/02_mousePhenotype_enrichment_DARs_stages_sacral.Rds"))

getPvals <- function(id=NULL){
  c(cervical = pheno.dars.cervical[["Mouse Phenotype Single KO"]][pheno.dars.cervical[["Mouse Phenotype Single KO"]]$name == id,]$Hyper_Adjp_BH,
    thoracic = pheno.dars.thoracic[["Mouse Phenotype Single KO"]][pheno.dars.thoracic[["Mouse Phenotype Single KO"]]$name == id,]$Hyper_Adjp_BH,
    lumbar = pheno.dars.lumbar[["Mouse Phenotype Single KO"]][pheno.dars.lumbar[["Mouse Phenotype Single KO"]]$name == id,]$Hyper_Adjp_BH,
    sacral = pheno.dars.sacral[["Mouse Phenotype Single KO"]][pheno.dars.sacral[["Mouse Phenotype Single KO"]]$name == id,]$Hyper_Adjp_BH)
}
tmp <- rbind("abnormal cervical vertebrae morphology" = getPvals(id = "abnormal cervical vertebrae morphology"),
             "abnormal thoracic cage morphology" = getPvals(id = "abnormal thoracic cage morphology"),
             "abnormal rib morphology" = getPvals(id = "abnormal rib morphology"),
             "abnormal lumbar vertebrae morphology" = getPvals(id = "abnormal lumbar vertebrae morphology"),
             "abnormal sacral vertebrae morphology" = getPvals(id = "abnormal sacral vertebrae morphology"))

df <- reshape2::melt(tmp)
df$Var1 <- factor(df$Var1, levels = rev(row.names(tmp)))
df <- df[order(df$Var1),]
df$x <- c(rev(seq(0.25,1,0.25)), rev(seq(2.25,3,0.25)), rev(seq(4.25,5,0.25)), rev(seq(6.25,7,0.25)), rev(seq(8.25,9,0.25)))

pdf(paste0(out, "Figure3M_KOenrichment.pdf"), useDingbats = FALSE, width = 6, height = 2)
ggplot(df, aes(x, -log10(value), fill=Var2)) +
  geom_point(aes(size=-log10(value)), pch=21, colour="grey20") +
  geom_segment(aes(x=x, xend=x, y=0, yend=-log10(value)), colour="grey20") +
  scale_fill_manual(values = c("darkolivegreen4", "cornflowerblue", "darkorange", "indianred1")) +
  geom_hline(yintercept = -log10(0.05), lty=2, colour="grey") +
  scale_x_continuous(breaks = seq(0.5,8.5,2), labels = levels(df$Var1), position = "top") +
  scale_y_reverse() +
  coord_flip() +
  xlab("") +
  ylab("-log10 p-value") +
  th + theme(legend.position = "none")
dev.off()



# Figure 4 ###############

## A: enhancer catalogues ======
cisCor <- readRDS(paste0(dir, "RNA+ATAC/results/03_gene_peak_correlations.Rds"))
links <- read.table(paste0(dir, "RNA+ATAC/results/03_gene_peak_links.tsv"), header = TRUE)
encode <- read.table(paste0(dir, "ATAC-seq/results/05_peaks_overlapExternalData.tsv"), sep="\t")
row.names(encode) <- paste0("chr", row.names(encode))
encode$encode <- ifelse(encode$pELS | encode$dELS, "enhancer",
                        ifelse(encode$PLS, "promoter",
                               ifelse(encode$CTCF_bound, "CTCFonly",
                                      ifelse(encode$H3K4me3, "H3K4me3only", "not_in_encode"))))

df <- data.frame(type = rep(c("enhancer_cCRE", "not_in_encode", "FANTOM5"), each=7),
                 prop = c(sum(encode$encode=="enhancer")/nrow(encode),
                          sum(encode[unique(cisCor$PeakRanges),]$encode=="enhancer")/nrow(encode[unique(cisCor$PeakRanges),]),
                          sum(links[links$rObs>0.3,]$encode=="enhancer")/nrow(links[links$rObs>0.3,]),
                          sum(links[links$rObs>0.4,]$encode=="enhancer")/nrow(links[links$rObs>0.4,]),
                          sum(links[links$rObs>0.5,]$encode=="enhancer")/nrow(links[links$rObs>0.5,]),
                          sum(links[links$rObs>0.6,]$encode=="enhancer")/nrow(links[links$rObs>0.6,]),
                          sum(links[links$rObs>0.7,]$encode=="enhancer")/nrow(links[links$rObs>0.7,]),
                          sum(encode$encode=="not_in_encode")/nrow(encode),
                          sum(encode[unique(cisCor$PeakRanges),]$encode=="not_in_encode")/nrow(encode[unique(cisCor$PeakRanges),]),
                          sum(links[links$rObs>0.3,]$encode=="not_in_encode")/nrow(links[links$rObs>0.3,]),
                          sum(links[links$rObs>0.4,]$encode=="not_in_encode")/nrow(links[links$rObs>0.4,]),
                          sum(links[links$rObs>0.5,]$encode=="not_in_encode")/nrow(links[links$rObs>0.5,]),
                          sum(links[links$rObs>0.6,]$encode=="not_in_encode")/nrow(links[links$rObs>0.6,]),
                          sum(links[links$rObs>0.7,]$encode=="not_in_encode")/nrow(links[links$rObs>0.7,]),
                          sum(encode$FANTOM)/nrow(encode),
                          sum(encode[unique(cisCor$PeakRanges),]$FANTOM)/nrow(encode[unique(cisCor$PeakRanges),]),
                          sum(links[links$rObs>0.3,]$fantom)/nrow(links[links$rObs>0.3,]),
                          sum(links[links$rObs>0.4,]$fantom)/nrow(links[links$rObs>0.4,]),
                          sum(links[links$rObs>0.5,]$fantom)/nrow(links[links$rObs>0.5,]),
                          sum(links[links$rObs>0.6,]$fantom)/nrow(links[links$rObs>0.6,]),
                          sum(links[links$rObs>0.7,]$fantom)/nrow(links[links$rObs>0.7,])))

pdf(paste0(out, "Figure4A_enhancer_overlap.pdf"), useDingbats = FALSE, width = 3, height = 4.5)
ggplot(df) +
  geom_point(aes(rep(c(2.5,2.75,seq(3,5,0.5)),3), prop, colour=type), size=2, alpha=rep(c(0.6,0.6, rep(1,5)),3)) +
  geom_line(data=df[-c(1:2,8:9,15:16),], aes(rep(seq(3,5,0.5),3), prop, colour=type), lty=3) +
  scale_color_manual(values = c("indianred3", "darkorange", "grey60")) +
  xlab("") +
  ylab("fraction of peaks") +
  ylim(c(0,1)) +
  scale_x_continuous(breaks=c(2.5,2.75,seq(3,5,0.5)),
                     labels=c("all peaks", "peaks near DEGs", paste0("linked peaks > 0.", 3:7))) +
  labs(colour="") +
  th + theme(axis.text.x = element_text(angle=90, hjust=1, vjust = 0.5),
             panel.grid.major.y = element_line(size=0.5, colour="grey90"),
             panel.grid.minor.y = element_line(size=0.25, colour="grey95"),
             legend.position = "bottom")
dev.off()

## B: highly-connected genes ======
n_links <- as.data.frame(table(links$Gene))
n_links$lab <- ifelse(n_links$Freq>14, as.character(n_links$Var1), "")

pdf(paste0(out, "Figure4B_highly_regulated_genes.pdf"), useDingbats = FALSE, width = 4, height = 4.5)
ggplot(n_links[order(n_links$Freq),], aes(Freq, 1:nrow(n_links), colour=Freq>5, label=lab)) +
  geom_point(alpha=0.75) +
  # geom_text_repel(size=3, angle=-45, min.segment.length = 1,
  #                 nudge_y = -400, nudge_x = 0.5, direction = "y") +
  scale_color_manual(values = c("grey", "steelblue4")) +
  xlab("number of linked peaks") +
  ylab("DE genes") +
  th + theme(axis.text.y = element_blank(),
             axis.ticks.y = element_blank(),
             legend.position = "none")
dev.off()

## C: regulation scores ======
fig.d <- readRDS(paste0(dir, "RNA+ATAC/results/03_FigR_network.Rds"))

# remove overlapping points for easier figure editing
keep <- iSEE::subsetPointsByGrid(fig.d$Corr.log10P, fig.d$Enrichment.log10P, resolution = 300)

tmp <- fig.d[keep,]
pdf(paste0(out, "Figure4C_regulation_scores.pdf"), useDingbats = FALSE, width = 8, height = 4.5)
ggplot(tmp, aes(Corr.log10P, Enrichment.log10P, fill=Score)) +
  geom_point_rast(data = tmp[tmp$Motif == "Nr6a1" & tmp$Enrichment.log10P>0,], size=2.5, pch=24, colour="grey30") +
  geom_point_rast(data = tmp[tmp$Motif != "Nr6a1",], size=1.5, pch=21, colour=rgb(1,1,1,0)) +
  scale_fill_gradientn(colours = diverge_hcl(n=7, palette = "Tropic"),
                        limits=c(-2,2),
                        oob = scales::squish,
                        breaks=scales::breaks_pretty(n=3)) +
  geom_vline(xintercept = 0, lty=2, colour="grey") +
  xlab(expression('(signed) -log'[10]*' p-value correlation')) +
  ylab(expression('-log'[10]*' p-value TF motif enrichment')) +
  labs(shape = "Nr6a1 enrichment") +
  th
dev.off()


## D: regulation scores heatmap ======
# get cluster membership
clusters <- read.table(paste0(dir, "RNA+ATAC/results/03_regulons_clusters.tsv"), header = TRUE)
clusters$cluster <- factor(clusters$cluster, levels=paste0("cluster", c(3,1,2,4)))

# filter interactions
score.cut <- 1.25
DORCsToKeep <- fig.d %>% dplyr::filter(abs(Score) >= score.cut) %>% dplyr::pull(DORC) %>% unique()
TFsToKeep <- fig.d %>% dplyr::filter(abs(Score) >= score.cut) %>% dplyr::pull(Motif) %>% unique()
net.d <- fig.d %>% dplyr::filter(DORC %in% DORCsToKeep & Motif %in% TFsToKeep) %>%
  reshape2::dcast(DORC ~ Motif) %>%
  tibble::column_to_rownames("DORC") %>% as.matrix()

# heatmap
h <- Heatmap(net.d,
             column_km = 2,
             row_split = clusters[order(clusters$gene),]$cluster,
             cluster_row_slices = FALSE)
h <- draw(h)

# reorder rows to have Hox at the bottom
tmp <- row_order(h)
tmp[['cluster4']] <- tmp[['cluster4']][c(70:72,1:46,59:69,47:58)] # reorder rows respecting dendrogram
tmp <- as.data.frame(unlist(tmp))
colnames(tmp) <- "idx"
tmp$cluster <- substr(row.names(tmp), 1,8)
tmp$cluster <- factor(tmp$cluster, levels=paste0("cluster", c(3,1,2,4)))
tmp$gene <- row.names(net.d)[tmp$idx]

# label genes of interest
tgfb <- c("Trp53","Eid2","Cited2","Bmp2","Skil","Tgfbr2","Tgfb3","Tgfb1","Creb1","Col1a2","Col3a1","Spi1","Ccl2","Ltbp1","Appl1","Tab1","Ptprk","Mstn","Eng","Fut8","Chst11","Zyx","Map3k7","Map3k1","Smad1","Itga8","Acvr1","Npnt","Onecut2","Dcp1a","Flcn","Mir145b","Acvr1c","Gcnt2","Parp1","Bmpr1b","Fermt2","Usp9x","Zmiz1","Smad3","Cav2","Gdnf","Appl2","Jun","Zfyve9","Smad9","Onecut1","Nlk","Dusp22","Arrb2","Cldn5","Usp15","Hpgd","Acvrl1","Glg1","Ltbp4","Lpxn","Nrros","Smad7","Itgb5","Itgb6","Pxn","Nodal","Amhr2","Cdh5","Pml","Tgfbr3","Bmp8a","Twsg1","Foxh1","Mtmr4","Stat3","Fos","Ptk2","Wfikkn2","Itgb8","Bmp8b","Smad2","Lrrc32","Sirt1","Bambi","Tgfbr3l","Mir145a","Mir143","Smad5","Cited1","Ltbp3","Smad6","Src","Tgfb2","Tgfbr1","Smad4","Adam9","Hipk2")
bmp <- c("Wnt5a","Wnt3a","Wnt1","Mir329","Spg20","Bmp6","Bmp5","Bmp4","Bmp2","Megf8","Tgfb3","Tgfb1","Comp","Bmp15","Dlx3","Sostdc1","Gdf15","Bmp7","Gdf9","Atf2","Bmper","Scx","Acvr2a","Dlx5","Mapk3","Smad1","Acvr2b","Acvr1","Chrd","Pparg","Smurf1","Abl1","Ecsit","Dsg4","Ror2","Egr1","Bmpr1b","Usp9x","Ext1","Etv2","Smad3","Tmem100","Pdcd4","Lefty1","Rgma","Bmncr","Smad9","Sost","Fkbp8","Zfp128","Tmprss6","Zcchc12","Tcf7l2","Mir23a","Mir210","Mir16-1","Mir147","Sfrp1","Usp15","Sfrp2","Acvrl1","Bmpr1a","Slc39a5","Bmp10","Gdf2","Nog","Msx2","Msx1","Smad7","Rgmb","Bmpr2","Fam83g","Hjv","Notch2","Nodal","Slc33a1","Chrdl1","Ddx5","Tgfbr3","Bmp8a","Twsg1","Vsir","Gdf1","Fst","Lefty2","Hivep1","Lef1","Nanog","Hes5","Hes1","Bmp8b","Runx2","Smad2","Mir125a","Mir122","Mirlet7f-2","Mirlet7f-1","Mirlet7d","Mirlet7b","Mirlet7a-1","Smad5","Smad6","Smpd3","Tgfb2","Smad4","Hfe","Id1","Gdf7","Gdf6","Gdf5","Gdf3")
ids <- c("Cdx1","Cdx2","Hey1","Lef1", "Tbx3","Foxp4","Sall4", "Cdh5","Dlk1","Epha2","Id3", "Bmp7","Smad5","Wnt5a",
         intersect(row.names(net.d), c(tgfb,bmp)))
ha.ids <- rowAnnotation(genes = anno_mark(at = tmp[match(ids, tmp$gene),]$idx, 
                                          labels = ids))

# label TFS of interest
ids <- c("Cdx1", "Gbx2", "Sall4", "Gli1","Tcf3", "Rara", "Lef1", "Nr6a1", "Ebf1", "Sox6", "Sox9", "Mef2c", "Snai2","Zfp637","Hoxa9","Foxp1",
         "Smad5","Sox11","Pax9")
ha.tfs <- columnAnnotation(TFs = anno_mark(at = match(ids, colnames(net.d)), labels = toupper(ids),
                                           which = "column", side = "bottom"))

pdf(paste0(out, "Figure4D_reg_heatmap.pdf"), useDingbats = FALSE, width = 6, height = 8)
draw(
  Heatmap(net.d,
          row_order = tmp$idx,
          col = colorRamp2(breaks = c(-2,-1,0,1,2),
                           colors = diverge_hcl(n=5, palette = "Tropic")),
          show_row_names = FALSE, 
          show_column_names = FALSE,
          column_names_gp = gpar(fontsize=8),
          column_km = 2,
          row_split = clusters[order(clusters$gene),]$cluster,
          cluster_row_slices = FALSE,
          row_title = NULL,
          column_title = NULL,
          right_annotation = ha.ids,
          bottom_annotation = ha.tfs,
          border_gp = gpar(col = "grey20", lwd = 0.5), 
          heatmap_legend_param = list(title = "regulation score", direction = "horizontal"),
          use_raster = TRUE),
  heatmap_legend_side = "bottom")
dev.off()

## E: expression heatmap ======
# gene counts for the same regulated genes
data <- geneCounts[match(tmp$gene, geneCounts$gene),]
row.names(data) <- data$gene
data$gene <- NULL
data <- t(apply(data, 1, function(x) (x-mean(x))/sd(x) ))

stopifnot(identical(colnames(data), meta.rna$sample)) # make sure the metadata corresponds with count data matrix
order <- order(meta.rna$stage, meta.rna$somite)
ha  <- HeatmapAnnotation(df = data.frame(stage = meta.rna[order,]$stage, 
                                         somite = meta.rna[order,]$somite), 
                         col = list(stage = cols.stage, somite = cols.somite),
                         simple_anno_size = unit(0.25, "cm"),
                         show_annotation_name = FALSE,
                         annotation_legend_param = list(stage = list(nrow=2),
                                                        somite = list(nrow = 1)))
# group by clusters
fates <- c("cervical", rep("thoracic",3), "lumbar", "sacral")
names(fates) <- paste0("stage",c(8,18,21,25,27,35))
fates <- factor(fates[meta.rna[order,]$stage], levels=c("cervical", "thoracic", "lumbar", "sacral"))

pdf(paste0(out, "Figure4E_expr_heatmap.pdf"), useDingbats = FALSE, width = 5, height = 8)
draw(
  Heatmap(data[,order], 
        cluster_columns = FALSE, 
        cluster_rows = FALSE,
        col = colorRamp2(breaks = c(-4,-2,0,2,4),
                   colors = c("steelblue4","steelblue","white","indianred3","indianred4")),
        show_row_names = FALSE, 
        show_column_names = FALSE,
        show_row_dend = FALSE,
        row_title = NULL,
        column_title = NULL,
        top_annotation = ha,
        column_split = fates,
        row_split = tmp$cluster,
        cluster_row_slices = FALSE,
        border_gp = gpar(col = "grey20", lwd = 0.5), 
        heatmap_legend_param = list(title = "z-score", direction = "horizontal"),
        use_raster = TRUE),
  heatmap_legend_side = "bottom",
  annotation_legend_side = "bottom",
  merge_legend = TRUE)
dev.off()



## F: Sall4 ======
# gene expression
df <- meta.rna
df$expr <- as.numeric(geneCounts[geneCounts$gene=="Sall4",-1])

pdf(paste0(out, "Figure4F_Sall4_expr.pdf"), useDingbats = FALSE, width = 3.5, height = 4)
ggplot(df, aes(stage, expr, fill=somite)) +
  geom_boxplot(alpha=0.5) +
  scale_fill_manual(values = cols.somite) +
  ggtitle("Sall4") +
  xlab("") +
  ylab(expression('log'[2]*' CPM')) + 
  th + theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5),
             plot.title = element_text(face="italic", hjust = 0.5),
             legend.position = "none")
dev.off()

## G: expression heatmap ======
# Nr6a1 - Gdf11
df <- data.frame(stage = meta.rna$stage,
                 Nr6a1 = as.numeric(geneCounts[geneCounts$gene=="Nr6a1",-1]),
                 Gdf11 = as.numeric(geneCounts[geneCounts$gene=="Gdf11",-1]))

pdf(paste0(out, "Figure4G_Nr6a2_Gdf11.pdf"), useDingbats = FALSE, width = 4, height = 4)
ggplot(df, aes(Gdf11, Nr6a1, colour=stage)) +
  geom_point(size=3, alpha=0.85) +
  scale_color_manual(values = cols.stage) +
  # annotation_logticks() +
  xlab(expression('log'[2]*' CPM Gdf11')) +
  ylab(expression('log'[2]*' CPM Nr6a1')) +
  th + theme(legend.position = "none")
dev.off()








# Figure S1 ###############

## Left-right DE analysis ======
res <- read.table(paste0(dir, "RNA-seq/results/00_DEresults_side_allSomites.tsv"))

pdf(paste0(out, "FigureS1.pdf"), useDingbats = FALSE, width = 5, height = 5)
ggplot(res, aes(logCPM, logFC, colour=FDR < 0.05)) +
  geom_point() +
  scale_color_manual(values = c("black", "indianred")) +
  geom_hline(yintercept = 0, lty=2, colour="red") +
  xlab(expression('log'[2]*' mean expression')) +
  ylab(expression('log'[2]*' fold-change')) +
  th + theme(legend.position = "none")
dev.off()


# Figure S2 ###############

## A: PCA normalised data ======
data <- read.table(paste0(dir, "RNA-seq/data/geneCounts.RAW.tsv"), check.names = FALSE)
dataNorm <- read.table(paste0(dir, "RNA-seq/data/geneCounts.NORM_logCPM.tsv"), check.names = FALSE)
meta <- read.table(paste0(dir, "RNA-seq/data/metadata_RNAseq.tsv"), stringsAsFactors = FALSE, header = TRUE)
meta <- meta[meta$sample %in% colnames(dataNorm),]

# variance-stabilisation
tmp <- data[row.names(dataNorm), colnames(dataNorm)]
data.vst <- vst(as.matrix(tmp[,-1]))
vars <- rowVars(data.vst)
names(vars) <- row.names(dataNorm)

# PCA
pca <- prcomp(t(dataNorm[names(vars[order(vars, decreasing = TRUE)])[1:1000],-1]))

# plot
df <- as.data.frame(pca$x)
df$stage <- paste0("stage", meta[match(row.names(df), meta$sample),'stage'])
df$date <- meta[match(row.names(df), meta$sample),'date']

plots <- list()
plots[[1]] <- ggplot(df, aes(PC1, PC2, colour=stage)) + 
  geom_point() + 
  scale_color_manual(values = cols.stage) +
  labs(col="stage") + 
  ggtitle("stage") + 
  th + theme(legend.position = "bottom")
plots[[2]] <- ggplot(df, aes(PC1, PC2, colour=date)) + 
  geom_point() + 
  labs(col="date") + 
  ggtitle("date") + 
  th + theme(legend.position = "bottom")

pdf(paste0(out, "FigureS2A.pdf"), useDingbats = FALSE, width = 9, height = 5)
ggarrange(plotlist = plots, ncol=2)
dev.off()

## B: PCA batch-corrected data ======
dataNorm.corr <- read.table(paste0(dir, "RNA-seq/data/geneCounts.NORM_batchCorrected_14PCs.tsv"), check.names = FALSE)
pcs <- read.table(paste0(dir, "RNA-seq/results/02_pcs_residuals.tab"))

# variance-stabilisation
tmp <- data[row.names(dataNorm.corr), colnames(dataNorm.corr)]
data.vst <- vst(as.matrix(tmp[,-1]))
data.vst <- removeBatchEffect(data.vst, design = model.matrix(~0+group, meta[meta$QC==1,]), covariates = pcs[,1:14])
vars <- rowVars(data.vst)
names(vars) <- row.names(dataNorm)

# PCA
pca <- prcomp(t(dataNorm.corr[names(vars[order(vars, decreasing = TRUE)])[1:1000],-1]))

df <- as.data.frame(pca$x)
df$stage <- paste0("stage", meta[match(row.names(df), meta$sample),'stage'])
df$date <- meta[match(row.names(df), meta$sample),'date']

plots <- list()
plots[[1]] <- ggplot(df, aes(PC1, PC2, colour=stage)) + 
  geom_point() + 
  scale_color_manual(values = cols.stage) +
  labs(col="stage") + 
  ggtitle("stage") + 
  th + theme(legend.position = "bottom")
plots[[2]] <- ggplot(df, aes(PC1, PC2, colour=date)) + 
  geom_point() + 
  labs(col="date") + 
  ggtitle("date") + 
  th + theme(legend.position = "bottom")

pdf(paste0(out, "FigureS2B.pdf"), useDingbats = FALSE, width = 9, height = 5)
ggarrange(plotlist = plots, ncol=2)
dev.off()


# Figure S3 ###############
cols <- GetColors(34)[c(25,24,15:13)]

## A: insert size distribution ======
diagnostics <- readRDS(paste0(dir, "ATAC-seq/results/01_insertSizeFrequencies.Rds"))
samples <- c("e23_SII-2.noDUPs.GQ.bam", "e16_SI-2.noDUPs.GQ.bam",
             "e26_SII-2.noDUPs.GQ.bam", "e2_SII-2.noDUPs.GQ.bam",
             "e13_SIII-2.noDUPs.GQ.bam")

pdf(paste0(out, "FigureS3A.pdf"), useDingbats = FALSE, width = 12, height = 3)
par(mfrow=c(1,5))
for(i in 1:length(samples)){
  tmp <- diagnostics[[samples[i]]]
  tmp <- unlist(sapply(1:nrow(tmp), function(i) rep(tmp[i,2], tmp[i,3])))
  plot(density(tmp), type="l", lwd=5, xlim=c(0,1e3), col=cols[i],
       xlab("insert size"), main="")
}
dev.off()

## B: TSS enrichment ======
tss <- readRDS(paste0(dir, "ATAC-seq/results/02_TSSinsertionCounts.Rds"))

# normalise to background
tss.norm <- t(do.call("cbind", lapply(tss, function(x) colMeans(x)/mean(colMeans(x[,c(1:100,1901:2001)])))))

samples <- gsub(".noDUPs.GQ.bam", "", samples)
pdf(paste0(out, "FigureS3B.pdf"), useDingbats = FALSE, width = 12, height = 3)
par(mfrow=c(1,5))
for(i in 1:length(samples)){
  plot(rollmean(tss.norm[samples[i],], k=25), type="l", lwd=5, 
       xlab="", ylab="", 
       col=cols[i], ylim=c(0,10), axes=FALSE)
  box(bty="l"); axis(1, at=c(0,1000,2000), labels = c("-1kb","TSS","1kb")); axis(2,las=2)
  abline(h=5, lty=2, lwd=2)
}
dev.off()

## C: number of peaks ======
meta.atac.all <- read.table(paste0(dir, "ATAC-seq/data/metadata_ATACseq.tsv"), stringsAsFactors = FALSE, header = TRUE)

pdf(paste0(out, "FigureS3C.pdf"), useDingbats = FALSE, width = 4, height = 4)
ggplot(meta.atac.all, aes(as.factor(insSizeDist), nPeaks/1e3)) +
  geom_violin() +
  geom_jitter(width=0.1, colour="grey40") +
  geom_jitter(data=meta.atac.all[meta.atac.all$sample %in% samples,], 
              aes(as.factor(insSizeDist), nPeaks/1e3, colour=as.factor(insSizeDist))) +
  scale_color_manual(values = cols) +
  geom_hline(yintercept = 15, colour="grey", lty=2) +
  xlab("insert size score") +
  ylab("number of peaks x 1000") +
  th + theme(legend.position = "none")
dev.off()

## D: FRiP ======
pdf(paste0(out, "FigureS3D.pdf"), useDingbats = FALSE, width = 4, height = 4)
ggplot(meta.atac.all, aes(as.factor(insSizeDist), readsInPeaks/goodQuality*100)) +
  geom_violin() +
  geom_jitter(width=0.1, colour="grey40") +
  geom_jitter(data=meta.atac.all[meta.atac.all$sample %in% samples,], 
              aes(as.factor(insSizeDist), readsInPeaks/goodQuality*100, colour=as.factor(insSizeDist))) +
  scale_color_manual(values = cols) +
  geom_hline(yintercept = 3, colour="grey", lty=2) +
  xlab("insert size score") +
  ylab("number of peaks x 1000") +
  th + theme(legend.position = "none")
dev.off()

## E: QC pass ======
qc <- data.frame(nuclosome = ifelse(as.numeric(meta.atac.all$insSizeDist)>=2, 1, 0), 
                 nPeaks = ifelse(meta.atac.all$nPeaks > 15000, 1, 0), 
                 frip = ifelse(meta.atac.all$readsInPeaks/meta.atac.all$goodQuality*100 >= 3, 1, 0), 
                 tss = ifelse(meta.atac.all$TSSscore > 4, 1, 0))
row.names(qc) <- meta.atac.all$sample

pdf(paste0(out, "FigureS3E.pdf"), useDingbats = FALSE, width = 4, height = 4)
upset(qc, mainbar.y.label = "number of samples", sets.x.label = "number of samples\nthat pass", text.scale=1.25)
dev.off()

## F: DNA size ======
pdf(paste0(out, "FigureS3F.pdf"), useDingbats = FALSE, width = 4, height = 4)
ggplot(meta.atac.all, aes(as.factor(insSizeDist), size/1e3, fill=as.factor(insSizeDist))) +
  geom_boxplot() +
  scale_fill_manual(values = cols) +
  xlab("insert size score") +
  ylab("DNA size (kb)") +
  th + theme(legend.position = "none")
dev.off()

## G: Seq depth ======
pdf(paste0(out, "FigureS3G.pdf"), useDingbats = FALSE, width = 4, height = 4)
ggplot(meta.atac.all, aes(as.factor(insSizeDist), librarySize/1e6, fill=as.factor(insSizeDist))) +
  geom_boxplot() +
  scale_fill_manual(values = cols) +
  xlab("insert size score") +
  ylab("million fragments") +
  th + theme(legend.position = "none")
dev.off()



# Figure S4 ###############

## A: MA plot background ======
background <- readRDS(paste0(dir, "ATAC-seq/results/03_backgroundCounts_10kbBins.Rds"))
adj.counts <- cpm(asDGEList(background), log=TRUE)

pdf(paste0(out, "FigureS4A.pdf"), useDingbats = FALSE, width = 8, height = 4)
par(mfrow=c(1, 2))
for (i in c(11,40)){
  cur.x <- adj.counts[,1]
  cur.y <- adj.counts[,i]
  smoothScatter(x=(cur.x+cur.y)/2, y=cur.x-cur.y, 
                main=paste("1 vs", i), 
                xlab="A", ylab="M",
                ylim=c(-4,4))
  abline(h=0, col="red", lty=2, lwd=2)
}
dev.off()

## B: MA plot window + scaled norm ======
filtered.data <- readRDS(file=paste0(dir, "ATAC-seq/results/03_windowCounts_filteredWindows.Rds"))

# average counts
abval <- aveLogCPM(asDGEList(filtered.data))
o <- order(abval)

# raw counts
adjc <- log2(assay(filtered.data)+0.5)

# normalised counts - scaling factors
filtered.data <- normFactors(filtered.data, se.out = TRUE)
re.adjc <- cpm(asDGEList(filtered.data), log=TRUE)

pdf(paste0(out, "FigureS4B.pdf"), useDingbats = FALSE, width = 8, height = 8)
par(mfrow=c(2,2))
for(i in c(11,40)){
  mval <- adjc[,1]-adjc[,i]
  fit <- loessFit(x=abval, y=mval)
  smoothScatter(abval, mval, 
                main=paste("Raw 1 vs",i), 
                ylab="M", xlab="Average logCPM", 
                ylim=c(-5,5))
  lines(abval[o], fit$fitted[o], col="red")
  abline(h=0)
}
for(i in c(11,40)){
  mval <- re.adjc[,1]-re.adjc[,i]
  fit <- loessFit(x=abval, y=mval)
  smoothScatter(abval, mval, 
                main=paste("Normalised 1 vs",i), 
                ylab="M", xlab="Average logCPM", 
                ylim=c(-5,5))
  lines(abval[o], fit$fitted[o], col="red")
  abline(h=0)
}
dev.off()

## C: trended norm ======
re.adjc <- adjc - assay(filtered.data, "offset")/log(2)

pdf(paste0(out, "FigureS4C.pdf"), useDingbats = FALSE, width = 8, height = 4)
par(mfrow=c(1,2))
for(i in c(11,40)){
  mval <- re.adjc[,1]-re.adjc[,i]
  fit <- loessFit(x=abval, y=mval)
  smoothScatter(abval, mval, 
                main=paste("Normalised 1 vs",i), 
                ylab="M", xlab="Average logCPM", 
                ylim=c(-5,6.5))
  lines(abval[o], fit$fitted[o], col="red")
  abline(h=0)
}
dev.off()

## D: PCA ======

# variable regions
vars <- rowVars(as.matrix(re.adjc))
tmp <- re.adjc[order(vars, decreasing=TRUE)[1:5000],]

# PCA
pca <- prcomp(t(tmp))
df <- cbind(pca$x, meta.atac)

pdf(paste0(out, "FigureS4D.pdf"), useDingbats = FALSE, width = 8, height = 5)
plots <- list()
plots[[1]] <- ggplot(df, aes(PC1, PC2, color=stage)) + 
  geom_point() + 
  scale_color_manual(values = cols.stage) +
  labs(colour="stage") +
  th + theme(legend.position = "bottom",
             axis.text.x = element_blank(),
             axis.text.y = element_blank(),
             axis.ticks = element_blank())
plots[[2]] <- ggplot(df, aes(PC1, PC2, colour=readsInPeakSet/goodQuality*100)) + 
  geom_point() + 
  scale_color_gradientn(colors=rev(brewer.pal(n=10,"RdYlBu"))) + 
  labs(colour="FRiP") +
  th + theme(legend.position = "bottom",
             axis.text.x = element_blank(),
             axis.text.y = element_blank(),
             axis.ticks = element_blank())
ggarrange(plotlist = plots, align = "hv")
dev.off()

## D: Batch correction ======
pcs <- read.table(paste0(dir, "ATAC-seq/results/03_pcs_residuals.tab"))

# regress out technical noise
meta.atac$group <- factor(paste(meta.atac$stage, meta.atac$somite, sep="."))
design <- model.matrix(~0+group, meta.atac)
colnames(design) <- paste0("stage", levels(meta.atac$group))
norm.counts.corr.pca <- removeBatchEffect(re.adjc, design=design, covariates = pcs[,1:18])

# variable regions
vars <- rowVars(norm.counts.corr.pca)
tmp <- norm.counts.corr.pca[order(vars, decreasing=TRUE)[1:5000],]

# PCA
pca <- prcomp(t(tmp))
df <- cbind(pca$x, meta.atac)

pdf(paste0(out, "FigureS4E.pdf"), useDingbats = FALSE, width = 8, height = 5)
plots <- list()
plots[[1]] <- ggplot(df, aes(PC1, PC2, color=stage)) + 
  geom_point() + 
  scale_color_manual(values = cols.stage) +
  labs(colour="stage") +
  th + theme(legend.position = "bottom",
             axis.text.x = element_blank(),
             axis.text.y = element_blank(),
             axis.ticks = element_blank())
plots[[2]] <- ggplot(df, aes(PC1, PC2, colour=readsInPeakSet/goodQuality*100)) + 
  geom_point() + 
  scale_color_gradientn(colors=rev(brewer.pal(n=10,"RdYlBu"))) + 
  labs(colour="FRiP") +
  th + theme(legend.position = "bottom",
             axis.text.x = element_blank(),
             axis.text.y = element_blank(),
             axis.ticks = element_blank())
ggarrange(plotlist = plots, align = "hv")
dev.off()



# Figure S5 ###############

## A: DA accessibility patterns ======
peakCounts <- readRDS(paste0(dir, "ATAC-seq/results/04_peakCounts_csawMerged.NORM.batchCorrected_18PCs.Rds"))
dars <- read.table(paste0(dir, "RNA+ATAC/results/02_DAregions_summary_stage_fate.tsv"))

data <- peakCounts[row.names(dars),]
data <- 2^data
data <- t(apply(data, 1, function(x) (x-mean(x))/sd(x) ))

## order clusters to follow developmental progression
clusters.dars <- factor(dars$cluster, levels=paste0("cluster", c(10,3, 1,4,9,2,5,12, 7,6,8, 11)))
names(clusters.dars) <- row.names(data)

stopifnot(identical(colnames(data), meta.atac$sample)) # make sure the metadata corresponds with count data matrix
order <- order(meta.atac$stage, meta.atac$somite)

## heatmap annotation
ha  <- HeatmapAnnotation(df = data.frame(stage = meta.atac[order,]$stage, 
                                         somite = meta.atac[order,]$somite), 
                         col = list(stage = cols.stage, somite = cols.somite))
fates <- c("cervical", rep("thoracic",3), "lumbar", "sacral")
names(fates) <- paste0("stage",c(8,18,21,25,27,35))
fates <- factor(fates[meta.atac[order,]$stage], levels=c("cervical", "thoracic", "lumbar", "sacral"))

pdf(paste0(out, "FigureS5A.pdf"), useDingbats = FALSE, width = 6, height = 8)
draw(Heatmap(data[,order], 
             cluster_columns = FALSE, 
             col=colorRamp2(breaks = c(-2,0,2,4), 
                            colors = c("steelblue","white","indianred3","indianred4")),
             name = "z-score", 
             show_row_names = FALSE, 
             show_column_names = FALSE, 
             show_row_dend = FALSE,
             row_title = NULL,
             column_title = NULL,
             top_annotation = ha,
             column_split = fates,
             row_split = clusters.dars,
             cluster_row_slices = FALSE,
             border_gp = gpar(col = "grey20", lwd = 0.5), 
             heatmap_legend_param = list(title = "z-score", direction = "horizontal"),
             use_raster = TRUE),
     heatmap_legend_side = "bottom",
     annotation_legend_side = "bottom",
     merge_legend = TRUE)
dev.off()

## B: stage DARs genomic location ======
dars.ann <- peakAnn[substr(row.names(dars), 4, 30),]

df <- data.frame(all = round(colSums(peakAnn[,7:12])/nrow(peakAnn)*100, 2),
                 DA = round(colSums(dars.ann[,7:12])/nrow(dars.ann)*100, 2))
df <- data.frame(class=rep(row.names(df),2),
                 DA = c(rep(FALSE, 6), rep(TRUE, 6)),
                 prop = c(df$all, df$DA))
df$class <- factor(df$class, levels=rev(c("promoter", "exonic", "intronic", "proximal", "distal", "intergenic")))

pdf(paste0(out, "FigureS5B.pdf"), useDingbats = FALSE, width = 3, height = 5)
ggplot(df, aes(DA, prop, fill=class, label=round(prop))) + 
  geom_bar(stat="identity", colour="black", width=0.5) +
  scale_fill_manual(values = c(brewer.pal(n=11, "PiYG")[7:10], "gold2", "darkgoldenrod2")) +
  geom_text(size = 3, position = position_stack(vjust = 0.5)) +
  ylab("proportion of all peaks") +
  th
dev.off()

## C: Myf5 ======
df <- as.data.frame(colData(dev))
i <- grep("Myf5", row.names(dev))
df$z <- assay(dev, 'z')[i,]

pdf(paste0(out, "FigureS5C.pdf"), useDingbats = FALSE, width = 3, height = 4)
ggplot(df, aes(stage, z, fill=stage)) +
  geom_boxplot(alpha=0.75) +
  scale_fill_manual(values = cols.stage) +
  geom_hline(yintercept = 0, lty=2, colour="grey") +
  ggtitle("MYF5") +
  xlab("") +
  th + theme(axis.text.x = element_blank(),
             axis.ticks.x = element_blank(),
             plot.title = element_text(hjust = 0.5),
             legend.position = "none")
dev.off()

# ## D: Sox9 ======
df <- meta.rna
df$expr <- as.numeric(geneCounts[geneCounts$gene=="Sox9",-1])

pdf(paste0(out, "FigureS5D.pdf"), useDingbats = FALSE, width = 3, height = 4)
ggplot(df, aes(stage, expr, fill=stage)) +
  geom_boxplot(alpha=0.75) +
  scale_fill_manual(values = cols.stage) +
  ggtitle("Sox9") +
  xlab("") +
  ylab(expression('log'[2]*' CPM')) +
  th + theme(axis.text.x = element_blank(),
             axis.ticks.x = element_blank(),
             plot.title = element_text(face="italic", hjust = 0.5),
             legend.position = "none")
dev.off()