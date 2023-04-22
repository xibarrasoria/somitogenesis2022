# helper.R
load("data/data.RData")
library(Gviz)
library(TxDb.Mmusculus.UCSC.mm10.ensGene)
txdb <- TxDb.Mmusculus.UCSC.mm10.ensGene
ensembl <- biomaRt::useEnsembl(biomart = "genes", dataset = "mmusculus_gene_ensembl", version = 96)


## colour palette
cols.stage <- c("darkolivegreen4", brewer.pal(n=9, "Blues")[c(4,6,7)], "darkorange", "indianred1")
names(cols.stage) <- paste0("stage", c(8,18,21,25,27,35))
cols.somite <- hcl.colors(n=5, palette = "PurpOr")[3:1]
names(cols.somite) <- paste0("S", c("I","II","III"))

## add vertical lines
vline <- function(x = 0, color = "grey") {
  list(
    type = "line", 
    y0 = 0, 
    y1 = 1, 
    yref = "paper",
    x0 = x, 
    x1 = x, 
    line = list(color = color, dash="dot", width=0.5)
  )
}

## list of all genes in dataset
genesAvailable <- as.character(geneCounts$gene)
genesAvailable <- genesAvailable[order(nchar(genesAvailable))]

## list of available TF with chromVar results
tfsAvailable <- as.character(chromvar$gene)
tfsAvailable <- tfsAvailable[order(nchar(tfsAvailable))]

## samples shared across data types
shared_samples <- intersect(meta.rna$sample, meta.atac$sample)

#### Gene expression
## plot expression across samples
boxplotExpr <- function(group_by='stage', colour_by='somite', gene="Hoxa1"){
  ## get data
  df <- data.frame(sample = meta.rna$sample,
                   group = meta.rna[,group_by],
                   expr = as.numeric(geneCounts[geneCounts$gene==gene,-1][1,]),
                   col = meta.rna[,colour_by])
  
  ## set colour palette
  if(colour_by == 'somite'){
    cols <- cols.somite
  }else{
    cols <- cols.stage
  }
  
  ## plot
  p <- plot_ly(data = df,
          type = "box",
          x = ~group, y = ~expr,
          color = ~col, colors = cols,
          width = 700, height = 450
  ) %>% layout( 
      images = list(
        list(
          source =  base64enc::dataURI(file = "images/somiteTrios.png"),
          x = 1.05, y = 0.1, 
          sizex = 0.85, sizey = 0.85,
          xref = "paper", yref = "paper", 
          xanchor = "left", 
          yanchor = "bottom"
        ),
        list(
          source =  base64enc::dataURI(file = "images/somiteStages.png"),
          x = 0, y = -0.1, 
          sizex = 1.05, sizey = 1.05,
          xref = "paper", yref = "paper", 
          xanchor = "bottom", 
          yanchor = "center"
          )
      ),
      margin = list(t = 50, b = 95, r = 120)
    ) 
  if(!identical(df$group, df$col)){ 
    p <- p %>%
      layout(boxmode = "group",
      autosize  = FALSE,
      title = gene,
      xaxis = c(list(title = ""), showline= T),
      yaxis = c(list(title = "log2 CPM", showline=T, 
                     # range=c(floor(min(df$expr)), ceiling(max(df$expr))),
                     hoverformat = ".2f")),
      hoverinfo = 'text',
      text="",
      showlegend = FALSE
      )
  }else{
    p <- p %>%
      layout(autosize  = FALSE,
             title = gene,
             xaxis = c(list(title = ""), showline= T),
             yaxis = c(list(title = "log2 CPM", showline=T, 
                            # range=c(floor(min(df$expr)), ceiling(max(df$expr))),
                            hoverformat = ".2f")),
             hoverinfo = 'text',
             text="",
             showlegend = FALSE
      )
  }
  
  if(group_by == "somite"){
    p <- p %>% 
      layout(shapes = list(vline(0.5), vline(1.5)))
  }
  if(group_by == "stage"){
    p <- p %>% 
      layout(shapes = list(vline(0.5), vline(1.5), vline(2.5), vline(3.5), vline(4.5)))
  }
  
  return(p)
}

## summary of DE results
# trios
summaryDEtrios <- function(gene="Hoxa1"){
  if(gene != ""){
    ## average
    all <- DEtriosAll[DEtriosAll$genes==gene,][1,]
    
    ## per-stage
    stage <- list()
    for(s in paste0("stage", c(8,18,21,25,27,35))){
      stage[[s]] <- DEtriosStage[[s]][DEtriosStage[[s]]$genes==gene,][1,]
    }
    
    ## plot
    plots <- list()
    rng <- range(c(all[,2:4], unlist(lapply(stage, function(x) x[2:4])) ))
    
    # avg
    df <- data.frame(var1 = factor(c("SI", "SI", "SII"), levels = c("SI", "SII")), 
                     var2 = factor(c("SII", "SIII", "SIII"), levels = c("SII", "SIII")), 
                     fc = as.numeric(all[,2:4]), 
                     fdr = all$FDR)
    df <- df[order(abs(df$fc)),]
    plots[[1]] <- ggplot(df, aes(var1, var2, fill = fc)) + 
      geom_tile(color = ifelse(df$fdr < 0.05 & abs(df$fc) > log2(1.5), "black", "grey"), 
                size = ifelse(df$fdr < 0.05 & abs(df$fc) > log2(1.5), 0.5, 0.1)) + 
      scale_fill_gradient2(low = "steelblue", high = "indianred3", mid = "white", midpoint = 0, limit = rng, name="log2 FC") + 
      xlab("average") + 
      ylab("") +
      coord_fixed() + 
      theme_minimal() + 
      theme(panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank(), 
            axis.text = element_text(family='sans', face=1, size=10), 
            axis.title.x = element_text(family='sans', face=1, size=12))
    
    # stage
    i=2
    for(s in paste0("stage", c(8,18,21,25,27,35))){
      df <- data.frame(var1 = factor(c("SI", "SI", "SII"), levels = c("SI", "SII")), 
                       var2 = factor(c("SII", "SIII", "SIII"), levels = c("SII", "SIII")), 
                       fc = c(stage[[s]][,2], stage[[s]][,3], stage[[s]][,4]), 
                       fdr = stage[[s]]$FDR)
      df <- df[order(abs(df$fc)),]
      plots[[i]] <- ggplot(df, aes(var1, var2, fill = fc)) + 
        geom_tile(color = ifelse(df$fdr < 0.05 & abs(df$fc) > log2(1.5), "black", "grey"), 
                  size = ifelse(df$fdr < 0.05 & abs(df$fc) > log2(1.5), 0.5, 0.1)) + 
        scale_fill_gradient2(low = "steelblue", high = "indianred3", mid = "white", midpoint = 0, limit = rng, name="log2 FC") + 
        xlab(s) + 
        ylab("") + 
        coord_fixed() + 
        theme_minimal() + 
        theme(panel.grid.major = element_blank(), 
              panel.grid.minor = element_blank(), 
              axis.text = element_text(family='sans', face=1, size=10), 
              axis.title.x = element_text(family='sans', face=1, size=12))
      i <- i+1
    }
    p <- ggarrange(plotlist = plots, ncol = 7, nrow = 1, legend = "none")
    return(p)
  }
}

# stages
summaryDEstages <- function(gene="Hoxa1"){
  if(gene != ""){
    ## average
    all <- DEstageAll[DEstageAll$genes==gene,][1,]
    
    ## per-somite
    somite <- list()
    for(s in paste0("somite", c("I","II","III"))){
      somite[[s]] <- DEstageSomite[[s]][DEstageSomite[[s]]$genes==gene,][1,]
    }
    
    ## plot
    plots <- list()
    rng <- range(c(all[2:16], somite[["somiteI"]][2:16], somite[["somiteII"]][2:16], somite[["somiteIII"]][2:16]))
    
    # avg
    stages <- c(8,18,21,25,27,35)
    df <- data.frame(var1 = factor(rep(c(8,18,21,25,27), times=c(5:1)), levels = stages),
                     var2 = factor(c(stages[2:6], stages[3:6], stages[4:6], stages[5:6], stages[6]), levels = stages),
                     fc = as.numeric(all[,2:16]),
                     fdr = all$FDR)
    df <- df[order(abs(df$fc)),]
    plots[[1]] <- ggplot(df, aes(var1, var2, fill = fc)) + 
      geom_tile(color = ifelse(df$fdr < 0.05 & abs(df$fc) > log2(1.5), "black", "grey"), 
                size = ifelse(df$fdr < 0.05 & abs(df$fc) > log2(1.5), 0.5, 0.1)) + 
      scale_fill_gradient2(low = "steelblue", high = "indianred3", mid = "white", midpoint = 0, limit = rng, name="log2 FC") + 
      xlab("average") + 
      ylab("") +
      coord_fixed() + 
      theme_minimal() + 
      theme(panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank(), 
            axis.text = element_text(family='sans', face=1, size=10), 
            axis.title.x = element_text(family='sans', face=1, size=12))
    
    # somite
    i=2
    for(s in paste0("somite", c("III","II","I"))){
      df <- data.frame(var1=factor(rep(c(8,18,21,25,27), times=c(5:1)), levels = stages),
                       var2 = factor(c(stages[2:6], stages[3:6], stages[4:6], stages[5:6], stages[6]), levels = stages),
                       fc=as.numeric(somite[[s]][,2:16]),
                       fdr=somite[[s]]$FDR)
      df <- df[order(abs(df$fc)),]
      plots[[i]] <- ggplot(df, aes(var1, var2, fill = fc)) + 
        geom_tile(color = ifelse(df$fdr < 0.05 & abs(df$fc) > log2(1.5), "black", "grey"), 
                  size = ifelse(df$fdr < 0.05 & abs(df$fc) > log2(1.5), 0.5, 0.1)) + 
        scale_fill_gradient2(low = "steelblue", high = "indianred3", mid = "white", midpoint = 0, limit = rng, name="log2 FC") + 
        xlab(s) + 
        ylab("") + 
        coord_fixed() + 
        theme_minimal() + 
        theme(panel.grid.major = element_blank(), 
              panel.grid.minor = element_blank(), 
              axis.text = element_text(family='sans', face=1, size=10), 
              axis.title.x = element_text(family='sans', face=1, size=12))
      i <- i+1
    }
    p <- ggarrange(plotlist = plots, ncol = 4, nrow = 1, legend = "none")
    return(p)
  }
}



#### Peak accessibility
## show peaks in a genomic region with gene models
plotPeaks <- function(type="gene", gene=NULL, region=NULL, flank=5){
  if(type == "gene"){
    ## get gene's genomic location
    ann <- gene_ann[gene_ann$gene==gene,]
    chr <- paste0("chr", ann$chr)
    range <- GRanges(chr,
                      IRanges(ann$start - flank*1e3,
                              ann$end + flank*1e3),
                      gene = ann$gene)
  }else{
    ## convert to GRanges
    chr <- unlist(strsplit(region, "[:-]"))[1]
    range <- GRanges(chr,
                      IRanges(as.numeric(unlist(strsplit(region, "[:-]"))[2]) - flank*1e3,
                              as.numeric(unlist(strsplit(region, "[:-]"))[3]) + flank*1e3))
  }
  ## find overlapping peaks
  overlaps <- peaks[subjectHits(findOverlaps(range, peaks))]

  atrack <- AnnotationTrack(overlaps,
                            id = paste0("peak_", 1:length(overlaps)),
                            name="peaks",
                            stacking="dense",
                            fill = "#66808080")

  itrack <- IdeogramTrack(genome = "mm10", chromosome = chr)
  gtrack <- GenomeAxisTrack()
  biomTrack <- BiomartGeneRegionTrack(biomart = ensembl,
                                      chromosome = chr, start = start(range), end = end(range),
                                      transcriptAnnotation = "symbol",
                                      stacking = "squish",
                                      name = "ENSEMBL")

  p <- plotTracks(list(itrack, gtrack, atrack, biomTrack),
                  sizes = c(0.1,0.15,0.15,0.6),
                  from = start(range), to = end(range),
                  featureAnnotation = "id", fontcolor.feature = "grey20", fontsize.feature = 8,
                  extend.right = 1e3, extend.left = 1e3)
  return(p)

}

## plot accessibility across samples
boxplotAcc <- function(type="gene", group_by='stage', colour_by='somite', peak=NULL){
  ## select relevant peak
  peak <- unlist(strsplit(peak, " "))[2]
  
  ## get data
  df <- data.frame(sample = meta.atac$sample,
                   group = meta.atac[,group_by],
                   acc = peakCounts[peak,],
                   col = meta.atac[,colour_by])

  ## set colour palette
  if(colour_by == 'somite'){
    cols <- cols.somite
  }else{
    cols <- cols.stage
  }

  ## plot
  p <- plot_ly(data = df,
               type = "box",
               x = ~group, y = ~acc,
               color = ~col, colors = cols,
               width = 630, height = 405
  ) %>% layout( 
      images = list(
        list(
          source =  base64enc::dataURI(file = "images/somiteTrios.png"),
          x = 1.05, y = 0.1, 
          sizex = 0.85, sizey = 0.85,
          xref = "paper", yref = "paper", 
          xanchor = "left", 
          yanchor = "bottom"
        ),
        list(
          source =  base64enc::dataURI(file = "images/somiteStages.png"),
          x = 0, y = -0.1, 
          sizex = 1.05, sizey = 1.05,
          xref = "paper", yref = "paper", 
          xanchor = "bottom", 
          yanchor = "center"
        )
      ),
      margin = list(t = 50, b = 95, r = 120)
    ) 
  if(!identical(df$group, df$col)){
    p <- p %>%
      layout(boxmode = "group",
             autosize  = FALSE,
             title = peak,
             xaxis = c(list(title = ""), showline= T),
             yaxis = c(list(title = "log2 CPM", showline=T,
                            hoverformat = ".2f")),
             hoverinfo = 'text',
             text="",
             showlegend = FALSE
      ) 
  }else{
    p <- p %>%
      layout(autosize  = FALSE,
             title = peak,
             xaxis = c(list(title = ""), showline= T),
             yaxis = c(list(title = "log2 CPM", showline=T,
                            hoverformat = ".2f")),
             hoverinfo = 'text',
             text="",
             showlegend = FALSE
      ) 
  }

  if(group_by == "somite"){
    p <- p %>%
      layout(shapes = list(vline(0.5), vline(1.5)))
  }
  if(group_by == "stage"){
    p <- p %>%
      layout(shapes = list(vline(0.5), vline(1.5), vline(2.5), vline(3.5), vline(4.5)))
  }

  return(p)
}



#### Chromatin deviation scores
## plot TF activity across samples
boxplotChromvar <- function(group_by='stage', colour_by='somite', gene="Hoxa1"){
  ## get data
  df <- data.frame(sample = meta.atac$sample,
                   group = meta.atac[,group_by],
                   dev = as.numeric(chromvar[chromvar$gene==gene,-1]),
                   col = meta.atac[,colour_by])
  
  ## set colour palette
  if(colour_by == 'somite'){
    cols <- cols.somite
  }else{
    cols <- cols.stage
  }
  
  ## plot
  p <- plot_ly(data = df,
               type = "box",
               x = ~group, y = ~dev,
               color = ~col, colors = cols,
               width = 700, height = 450
  ) %>% layout( 
      images = list(
        list(
          source =  base64enc::dataURI(file = "images/somiteTrios.png"),
          x = 1.05, y = 0.1, 
          sizex = 0.85, sizey = 0.85,
          xref = "paper", yref = "paper", 
          xanchor = "left", 
          yanchor = "bottom"
        ),
        list(
          source =  base64enc::dataURI(file = "images/somiteStages.png"),
          x = 0, y = -0.1, 
          sizex = 1.05, sizey = 1.05,
          xref = "paper", yref = "paper", 
          xanchor = "bottom", 
          yanchor = "center"
        )
      ),
      margin = list(t = 50, b = 95, r = 120)
    ) 
  if(!identical(df$group, df$col)){
    p <- p %>%
      layout(boxmode = "group",
             autosize  = FALSE,
             title = gene,
             xaxis = c(list(title = ""), showline= T),
             yaxis = c(list(title = "z-score", showline=T, 
                            range=c(floor(min(df$dev)), ceiling(max(df$dev))),
                            hoverformat = ".2f")),
             hoverinfo = 'text',
             text="",
             showlegend = FALSE
      ) 
  }else{
    p <- p %>%
      layout(autosize  = FALSE,
             title = gene,
             xaxis = c(list(title = ""), showline= T),
             yaxis = c(list(title = "z-score", showline=T, 
                            range=c(floor(min(df$dev)), ceiling(max(df$dev))),
                            hoverformat = ".2f")),
             hoverinfo = 'text',
             text="",
             showlegend = FALSE
      ) 
  }
  
  if(group_by == "somite"){
    p <- p %>% 
      layout(shapes = list(vline(0.5), vline(1.5)))
  }
  if(group_by == "stage"){
    p <- p %>% 
      layout(shapes = list(vline(0.5), vline(1.5), vline(2.5), vline(3.5), vline(4.5)))
  }
  
  return(p)
}



#### DE results - trios
## display DE results table
printDEtriosTable <- function(level="average", contrast="stage8"){
  if(level == "average"){
    tbl <- data.frame(gene = DEtriosAll$genes,
                      logCPM = round(DEtriosAll$logCPM, 2),
                      logFC = round(DEtriosAll$maxFC, 2),
                      FDR = DEtriosAll$FDR,
                      in_avg = as.logical(DEtrios_summary[row.names(DEtriosAll),'ave']),
                      row.names = row.names(DEtriosAll))
    tbl <- tbl[order(tbl$FDR),]
    tbl$FDR <- as.numeric(format(tbl$FDR, digits=4, width=5))
    tbl <- tbl[row.names(tbl) %in% row.names(DEtrios_summary),]
  }else{
    tbl <- data.frame(gene = DEtriosStage[[contrast]]$genes,
                      logCPM = round(DEtriosStage[[contrast]]$logCPM, 2),
                      logFC = round(DEtriosStage[[contrast]]$maxFC, 2),
                      FDR = DEtriosStage[[contrast]]$FDR,
                      in_avg = as.logical(DEtrios_summary[row.names(DEtriosStage[[contrast]]),'ave']),
                      row.names = row.names(DEtriosStage[[contrast]]))
    tbl <- tbl[order(tbl$FDR),]
    tbl$FDR <- as.numeric(format(tbl$FDR, digits=4, width=5))
    tbl <- tbl[row.names(tbl) %in% row.names(DEtrios_summary),]
  }
  return(tbl)
}

## plot expression across samples
boxplotExpr_DEtrios <- function(level="average", contrast="stage8", group_by='somite', colour_by='somite', sel=1){
  ## compute table to get gene
  if(level == "average"){
    tbl <- data.frame(gene = DEtriosAll$genes,
                      logCPM = round(DEtriosAll$logCPM, 2),
                      logFC = round(DEtriosAll$maxFC, 2),
                      FDR = DEtriosAll$FDR,
                      row.names = row.names(DEtriosAll))
    tbl <- tbl[order(tbl$FDR),]
    tbl <- tbl[row.names(tbl) %in% row.names(DEtrios_summary),]
  }else{
    tbl <- data.frame(gene = DEtriosStage[[contrast]]$genes,
                      logCPM = round(DEtriosStage[[contrast]]$logCPM, 2),
                      logFC = round(DEtriosStage[[contrast]]$maxFC, 2),
                      FDR = DEtriosStage[[contrast]]$FDR,
                      row.names = row.names(DEtriosStage[[contrast]]))
    tbl <- tbl[order(tbl$FDR),]
    tbl <- tbl[row.names(tbl) %in% row.names(DEtrios_summary),]
  }
  
  ## gene to plot
  gene <- tbl[sel,'gene']
  ## get data
  df <- data.frame(sample = meta.rna$sample,
                   group = meta.rna[,group_by],
                   expr = as.numeric(geneCounts[geneCounts$gene==gene,-1][1,]),
                   col = meta.rna[,colour_by])
  
  ## set colour palette
  if(colour_by == 'somite'){
    cols <- cols.somite
  }else{
    cols <- cols.stage
  }
  
  ## plot
  p <- plot_ly(data = df,
               type = "box",
               x = ~group, y = ~expr,
               color = ~col, colors = cols,
               width = 700, height = 450
  ) %>% layout( 
    images = list(
      list(
        source =  base64enc::dataURI(file = "images/somiteTrios.png"),
        x = 1.05, y = 0.1, 
        sizex = 0.85, sizey = 0.85,
        xref = "paper", yref = "paper", 
        xanchor = "left", 
        yanchor = "bottom"
      ),
      list(
        source =  base64enc::dataURI(file = "images/somiteStages.png"),
        x = 0, y = -0.1, 
        sizex = 1.05, sizey = 1.05,
        xref = "paper", yref = "paper", 
        xanchor = "bottom", 
        yanchor = "center"
      )
    ),
    margin = list(t = 50, b = 95, r = 120)
  ) 
  if(!identical(df$group, df$col)){ 
    p <- p %>%
      layout(boxmode = "group",
             autosize  = FALSE,
             title = gene,
             xaxis = c(list(title = ""), showline= T),
             yaxis = c(list(title = "log2 CPM", showline=T, 
                            # range=c(floor(min(df$expr)), ceiling(max(df$expr))),
                            hoverformat = ".2f")),
             hoverinfo = 'text',
             text="",
             showlegend = FALSE
      )
  }else{
    p <- p %>%
      layout(autosize  = FALSE,
             title = gene,
             xaxis = c(list(title = ""), showline= T),
             yaxis = c(list(title = "log2 CPM", showline=T, 
                            # range=c(floor(min(df$expr)), ceiling(max(df$expr))),
                            hoverformat = ".2f")),
             hoverinfo = 'text',
             text="",
             showlegend = FALSE
      )
  }
  
  if(group_by == "somite"){
    p <- p %>% 
      layout(shapes = list(vline(0.5), vline(1.5)))
  }
  if(group_by == "stage"){
    p <- p %>% 
      layout(shapes = list(vline(0.5), vline(1.5), vline(2.5), vline(3.5), vline(4.5)))
  }
  
  return(p)
}

## summary of DE results
# trios
summaryDEtrios_tableTrios <- function(level="average", contrast="stage8", sel=1){
  ## compute table to get gene
  if(level == "average"){
    tbl <- data.frame(gene = DEtriosAll$genes,
                      logCPM = round(DEtriosAll$logCPM, 2),
                      logFC = round(DEtriosAll$maxFC, 2),
                      FDR = DEtriosAll$FDR,
                      row.names = row.names(DEtriosAll))
    tbl <- tbl[order(tbl$FDR),]
    tbl <- tbl[row.names(tbl) %in% row.names(DEtrios_summary),]
  }else{
    tbl <- data.frame(gene = DEtriosStage[[contrast]]$genes,
                      logCPM = round(DEtriosStage[[contrast]]$logCPM, 2),
                      logFC = round(DEtriosStage[[contrast]]$maxFC, 2),
                      FDR = DEtriosStage[[contrast]]$FDR,
                      row.names = row.names(DEtriosStage[[contrast]]))
    tbl <- tbl[order(tbl$FDR),]
    tbl <- tbl[row.names(tbl) %in% row.names(DEtrios_summary),]
  }
  
  ## gene to plot
  gene <- tbl[sel,'gene']
  
  if(gene != ""){
    ## average
    all <- DEtriosAll[DEtriosAll$genes==gene,]

    ## per-stage
    stage <- list()
    for(s in paste0("stage", c(8,18,21,25,27,35))){
      stage[[s]] <- DEtriosStage[[s]][DEtriosStage[[s]]$genes==gene,]
    }

    ## plot
    plots <- list()
    rng <- range(c(all[,2:4], unlist(lapply(stage, function(x) x[2:4])) ))

    # avg
    df <- data.frame(var1 = factor(c("SI", "SI", "SII"), levels = c("SI", "SII")),
                     var2 = factor(c("SII", "SIII", "SIII"), levels = c("SII", "SIII")),
                     fc = as.numeric(all[,2:4]),
                     fdr = all$FDR)
    df <- df[order(abs(df$fc)),]
    plots[[1]] <- ggplot(df, aes(var1, var2, fill = fc)) +
      geom_tile(color = ifelse(df$fdr < 0.05 & abs(df$fc) > log2(1.5), "black", "grey"),
                size = ifelse(df$fdr < 0.05 & abs(df$fc) > log2(1.5), 0.5, 0.1)) +
      scale_fill_gradient2(low = "steelblue", high = "indianred3", mid = "white", midpoint = 0, limit = rng, name="log2 FC") +
      xlab("average") +
      ylab("") +
      coord_fixed() +
      theme_minimal() +
      theme(panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            axis.text = element_text(family='sans', face=1, size=10),
            axis.title.x = element_text(family='sans', face=1, size=12))

    # stage
    i=2
    for(s in paste0("stage", c(8,18,21,25,27,35))){
      df <- data.frame(var1 = factor(c("SI", "SI", "SII"), levels = c("SI", "SII")),
                       var2 = factor(c("SII", "SIII", "SIII"), levels = c("SII", "SIII")),
                       fc = c(stage[[s]][,2], stage[[s]][,3], stage[[s]][,4]),
                       fdr = stage[[s]]$FDR)
      df <- df[order(abs(df$fc)),]
      plots[[i]] <- ggplot(df, aes(var1, var2, fill = fc)) +
        geom_tile(color = ifelse(df$fdr < 0.05 & abs(df$fc) > log2(1.5), "black", "grey"),
                  size = ifelse(df$fdr < 0.05 & abs(df$fc) > log2(1.5), 0.5, 0.1)) +
        scale_fill_gradient2(low = "steelblue", high = "indianred3", mid = "white", midpoint = 0, limit = rng, name="log2 FC") +
        xlab(s) +
        ylab("") +
        coord_fixed() +
        theme_minimal() +
        theme(panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              axis.text = element_text(family='sans', face=1, size=10),
              axis.title.x = element_text(family='sans', face=1, size=12))
      i <- i+1
    }
    p <- ggarrange(plotlist = plots, ncol = 7, nrow = 1, legend = "none")
    return(p)
  }
}

# stages
summaryDEstages_tableTrios <- function(level="average", contrast="stage8", sel=1){
  ## compute table to get gene
  if(level == "average"){
    tbl <- data.frame(gene = DEtriosAll$genes,
                      logCPM = round(DEtriosAll$logCPM, 2),
                      logFC = round(DEtriosAll$maxFC, 2),
                      FDR = DEtriosAll$FDR,
                      row.names = row.names(DEtriosAll))
    tbl <- tbl[order(tbl$FDR),]
    tbl <- tbl[row.names(tbl) %in% row.names(DEtrios_summary),]
  }else{
    tbl <- data.frame(gene = DEtriosStage[[contrast]]$genes,
                      logCPM = round(DEtriosStage[[contrast]]$logCPM, 2),
                      logFC = round(DEtriosStage[[contrast]]$maxFC, 2),
                      FDR = DEtriosStage[[contrast]]$FDR,
                      row.names = row.names(DEtriosStage[[contrast]]))
    tbl <- tbl[order(tbl$FDR),]
    tbl <- tbl[row.names(tbl) %in% row.names(DEtrios_summary),]
  }
  
  ## gene to plot
  gene <- tbl[sel,'gene']
  
  if(gene != ""){
    ## average
    all <- DEstageAll[DEstageAll$genes==gene,]

    ## per-somite
    somite <- list()
    for(s in paste0("somite", c("I","II","III"))){
      somite[[s]] <- DEstageSomite[[s]][DEstageSomite[[s]]$genes==gene,]
    }

    ## plot
    plots <- list()
    rng <- range(c(all[2:16], somite[["somiteI"]][2:16], somite[["somiteII"]][2:16], somite[["somiteIII"]][2:16]))

    # avg
    stages <- c(8,18,21,25,27,35)
    df <- data.frame(var1 = factor(rep(c(8,18,21,25,27), times=c(5:1)), levels = stages),
                     var2 = factor(c(stages[2:6], stages[3:6], stages[4:6], stages[5:6], stages[6]), levels = stages),
                     fc = as.numeric(all[,2:16]),
                     fdr = all$FDR)
    df <- df[order(abs(df$fc)),]
    plots[[1]] <- ggplot(df, aes(var1, var2, fill = fc)) +
      geom_tile(color = ifelse(df$fdr < 0.05 & abs(df$fc) > log2(1.5), "black", "grey"),
                size = ifelse(df$fdr < 0.05 & abs(df$fc) > log2(1.5), 0.5, 0.1)) +
      scale_fill_gradient2(low = "steelblue", high = "indianred3", mid = "white", midpoint = 0, limit = rng, name="log2 FC") +
      xlab("average") +
      ylab("") +
      coord_fixed() +
      theme_minimal() +
      theme(panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            axis.text = element_text(family='sans', face=1, size=10),
            axis.title.x = element_text(family='sans', face=1, size=12))

    # somite
    i=2
    for(s in paste0("somite", c("III","II","I"))){
      df <- data.frame(var1=factor(rep(c(8,18,21,25,27), times=c(5:1)), levels = stages),
                       var2 = factor(c(stages[2:6], stages[3:6], stages[4:6], stages[5:6], stages[6]), levels = stages),
                       fc=as.numeric(somite[[s]][,2:16]),
                       fdr=somite[[s]]$FDR)
      df <- df[order(abs(df$fc)),]
      plots[[i]] <- ggplot(df, aes(var1, var2, fill = fc)) +
        geom_tile(color = ifelse(df$fdr < 0.05 & abs(df$fc) > log2(1.5), "black", "grey"),
                  size = ifelse(df$fdr < 0.05 & abs(df$fc) > log2(1.5), 0.5, 0.1)) +
        scale_fill_gradient2(low = "steelblue", high = "indianred3", mid = "white", midpoint = 0, limit = rng, name="log2 FC") +
        xlab(s) +
        ylab("") +
        coord_fixed() +
        theme_minimal() +
        theme(panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              axis.text = element_text(family='sans', face=1, size=10),
              axis.title.x = element_text(family='sans', face=1, size=12))
      i <- i+1
    }
    p <- ggarrange(plotlist = plots, ncol = 4, nrow = 1, legend = "none")
    return(p)
  }
}



#### DE results - stages
## display DE results table
printDEstagesTable <- function(level="average", contrast="somiteI"){
  if(level == "average"){
    tbl <- data.frame(gene = DEstageAll$genes,
                      logCPM = round(DEstageAll$logCPM, 2),
                      logFC = round(DEstageAll$maxFC, 2),
                      FDR = DEstageAll$FDR,
                      in_avg = as.logical(DEstage_summary[row.names(DEstageAll),'average']),
                      row.names = row.names(DEstageAll))
    tbl <- tbl[order(tbl$FDR),]
    tbl$FDR <- as.numeric(format(tbl$FDR, digits=4, width=5))
    tbl <- tbl[row.names(tbl) %in% row.names(DEstage_summary),]
  }else{
    tbl <- data.frame(gene = DEstageSomite[[contrast]]$genes,
                      logCPM = round(DEstageSomite[[contrast]]$logCPM, 2),
                      logFC = round(DEstageSomite[[contrast]]$maxFC, 2),
                      FDR = DEstageSomite[[contrast]]$FDR,
                      in_avg = as.logical(DEstage_summary[row.names(DEstageSomite[[contrast]]),'average']),
                      row.names = row.names(DEstageSomite[[contrast]]))
    tbl <- tbl[order(tbl$FDR),]
    tbl$FDR <- as.numeric(format(tbl$FDR, digits=4, width=5))
    tbl <- tbl[row.names(tbl) %in% row.names(DEstage_summary),]
  }
  return(tbl)
}

## plot expression across samples
boxplotExpr_DEstages <- function(level="average", contrast="somiteI", group_by='stage', colour_by='stage', sel=1){
  ## compute table to get gene
  if(level == "average"){
    tbl <- data.frame(gene = DEstageAll$genes,
                      logCPM = round(DEstageAll$logCPM, 2),
                      logFC = round(DEstageAll$maxFC, 2),
                      FDR = DEstageAll$FDR,
                      row.names = row.names(DEstageAll))
    tbl <- tbl[order(tbl$FDR),]
    tbl <- tbl[row.names(tbl) %in% row.names(DEstage_summary),]
  }else{
    tbl <- data.frame(gene = DEstageSomite[[contrast]]$genes,
                      logCPM = round(DEstageSomite[[contrast]]$logCPM, 2),
                      logFC = round(DEstageSomite[[contrast]]$maxFC, 2),
                      FDR = DEstageSomite[[contrast]]$FDR,
                      row.names = row.names(DEstageSomite[[contrast]]))
    tbl <- tbl[order(tbl$FDR),]
    tbl <- tbl[row.names(tbl) %in% row.names(DEstage_summary),]
  }
  
  ## gene to plot
  gene <- tbl[sel,'gene']
  ## get data
  df <- data.frame(sample = meta.rna$sample,
                   group = meta.rna[,group_by],
                   expr = as.numeric(geneCounts[geneCounts$gene==gene,-1][1,]),
                   col = meta.rna[,colour_by])
  
  ## set colour palette
  if(colour_by == 'somite'){
    cols <- cols.somite
  }else{
    cols <- cols.stage
  }
  
  ## plot
  p <- plot_ly(data = df,
               type = "box",
               x = ~group, y = ~expr,
               color = ~col, colors = cols,
               width = 700, height = 450
  ) %>% layout( 
    images = list(
      list(
        source =  base64enc::dataURI(file = "images/somiteTrios.png"),
        x = 1.05, y = 0.1, 
        sizex = 0.85, sizey = 0.85,
        xref = "paper", yref = "paper", 
        xanchor = "left", 
        yanchor = "bottom"
      ),
      list(
        source =  base64enc::dataURI(file = "images/somiteStages.png"),
        x = 0, y = -0.1, 
        sizex = 1.05, sizey = 1.05,
        xref = "paper", yref = "paper", 
        xanchor = "bottom", 
        yanchor = "center"
      )
    ),
    margin = list(t = 50, b = 95, r = 120)
  ) 
  if(!identical(df$group, df$col)){ 
    p <- p %>%
      layout(boxmode = "group",
             autosize  = FALSE,
             title = gene,
             xaxis = c(list(title = ""), showline= T),
             yaxis = c(list(title = "log2 CPM", showline=T, 
                            # range=c(floor(min(df$expr)), ceiling(max(df$expr))),
                            hoverformat = ".2f")),
             hoverinfo = 'text',
             text="",
             showlegend = FALSE
      )
  }else{
    p <- p %>%
      layout(autosize  = FALSE,
             title = gene,
             xaxis = c(list(title = ""), showline= T),
             yaxis = c(list(title = "log2 CPM", showline=T, 
                            # range=c(floor(min(df$expr)), ceiling(max(df$expr))),
                            hoverformat = ".2f")),
             hoverinfo = 'text',
             text="",
             showlegend = FALSE
      )
  }
  
  if(group_by == "somite"){
    p <- p %>% 
      layout(shapes = list(vline(0.5), vline(1.5)))
  }
  if(group_by == "stage"){
    p <- p %>% 
      layout(shapes = list(vline(0.5), vline(1.5), vline(2.5), vline(3.5), vline(4.5)))
  }
  
  return(p)
}

## summary of DE results
# trios
summaryDEtrios_tableStages <- function(level="average", contrast="somiteI", sel=1){
  ## compute table to get gene
  if(level == "average"){
    tbl <- data.frame(gene = DEstageAll$genes,
                      logCPM = round(DEstageAll$logCPM, 2),
                      logFC = round(DEstageAll$maxFC, 2),
                      FDR = DEstageAll$FDR,
                      row.names = row.names(DEstageAll))
    tbl <- tbl[order(tbl$FDR),]
    tbl <- tbl[row.names(tbl) %in% row.names(DEstage_summary),]
  }else{
    tbl <- data.frame(gene = DEstageSomite[[contrast]]$genes,
                      logCPM = round(DEstageSomite[[contrast]]$logCPM, 2),
                      logFC = round(DEstageSomite[[contrast]]$maxFC, 2),
                      FDR = DEstageSomite[[contrast]]$FDR,
                      row.names = row.names(DEstageSomite[[contrast]]))
    tbl <- tbl[order(tbl$FDR),]
    tbl <- tbl[row.names(tbl) %in% row.names(DEstage_summary),]
  }
  
  ## gene to plot
  gene <- tbl[sel,'gene']
  
  if(gene != ""){
    ## average
    all <- DEtriosAll[DEtriosAll$genes==gene,]
    
    ## per-stage
    stage <- list()
    for(s in paste0("stage", c(8,18,21,25,27,35))){
      stage[[s]] <- DEtriosStage[[s]][DEtriosStage[[s]]$genes==gene,]
    }
    
    ## plot
    plots <- list()
    rng <- range(c(all[,2:4], unlist(lapply(stage, function(x) x[2:4])) ))
    
    # avg
    df <- data.frame(var1 = factor(c("SI", "SI", "SII"), levels = c("SI", "SII")),
                     var2 = factor(c("SII", "SIII", "SIII"), levels = c("SII", "SIII")),
                     fc = as.numeric(all[,2:4]),
                     fdr = all$FDR)
    df <- df[order(abs(df$fc)),]
    plots[[1]] <- ggplot(df, aes(var1, var2, fill = fc)) +
      geom_tile(color = ifelse(df$fdr < 0.05 & abs(df$fc) > log2(1.5), "black", "grey"),
                size = ifelse(df$fdr < 0.05 & abs(df$fc) > log2(1.5), 0.5, 0.1)) +
      scale_fill_gradient2(low = "steelblue", high = "indianred3", mid = "white", midpoint = 0, limit = rng, name="log2 FC") +
      xlab("average") +
      ylab("") +
      coord_fixed() +
      theme_minimal() +
      theme(panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            axis.text = element_text(family='sans', face=1, size=10),
            axis.title.x = element_text(family='sans', face=1, size=12))
    
    # stage
    i=2
    for(s in paste0("stage", c(8,18,21,25,27,35))){
      df <- data.frame(var1 = factor(c("SI", "SI", "SII"), levels = c("SI", "SII")),
                       var2 = factor(c("SII", "SIII", "SIII"), levels = c("SII", "SIII")),
                       fc = c(stage[[s]][,2], stage[[s]][,3], stage[[s]][,4]),
                       fdr = stage[[s]]$FDR)
      df <- df[order(abs(df$fc)),]
      plots[[i]] <- ggplot(df, aes(var1, var2, fill = fc)) +
        geom_tile(color = ifelse(df$fdr < 0.05 & abs(df$fc) > log2(1.5), "black", "grey"),
                  size = ifelse(df$fdr < 0.05 & abs(df$fc) > log2(1.5), 0.5, 0.1)) +
        scale_fill_gradient2(low = "steelblue", high = "indianred3", mid = "white", midpoint = 0, limit = rng, name="log2 FC") +
        xlab(s) +
        ylab("") +
        coord_fixed() +
        theme_minimal() +
        theme(panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              axis.text = element_text(family='sans', face=1, size=10),
              axis.title.x = element_text(family='sans', face=1, size=12))
      i <- i+1
    }
    p <- ggarrange(plotlist = plots, ncol = 7, nrow = 1, legend = "none")
    return(p)
  }
}

# stages
summaryDEstages_tableStages <- function(level="average", contrast="somiteI", sel=1){
  ## compute table to get gene
  if(level == "average"){
    tbl <- data.frame(gene = DEstageAll$genes,
                      logCPM = round(DEstageAll$logCPM, 2),
                      logFC = round(DEstageAll$maxFC, 2),
                      FDR = DEstageAll$FDR,
                      row.names = row.names(DEstageAll))
    tbl <- tbl[order(tbl$FDR),]
    tbl <- tbl[row.names(tbl) %in% row.names(DEstage_summary),]
  }else{
    tbl <- data.frame(gene = DEstageSomite[[contrast]]$genes,
                      logCPM = round(DEstageSomite[[contrast]]$logCPM, 2),
                      logFC = round(DEstageSomite[[contrast]]$maxFC, 2),
                      FDR = DEstageSomite[[contrast]]$FDR,
                      row.names = row.names(DEstageSomite[[contrast]]))
    tbl <- tbl[order(tbl$FDR),]
    tbl <- tbl[row.names(tbl) %in% row.names(DEstage_summary),]
  }
  
  ## gene to plot
  gene <- tbl[sel,'gene']
  
  if(gene != ""){
    ## average
    all <- DEstageAll[DEstageAll$genes==gene,]
    
    ## per-somite
    somite <- list()
    for(s in paste0("somite", c("I","II","III"))){
      somite[[s]] <- DEstageSomite[[s]][DEstageSomite[[s]]$genes==gene,]
    }
    
    ## plot
    plots <- list()
    rng <- range(c(all[2:16], somite[["somiteI"]][2:16], somite[["somiteII"]][2:16], somite[["somiteIII"]][2:16]))
    
    # avg
    stages <- c(8,18,21,25,27,35)
    df <- data.frame(var1 = factor(rep(c(8,18,21,25,27), times=c(5:1)), levels = stages),
                     var2 = factor(c(stages[2:6], stages[3:6], stages[4:6], stages[5:6], stages[6]), levels = stages),
                     fc = as.numeric(all[,2:16]),
                     fdr = all$FDR)
    df <- df[order(abs(df$fc)),]
    plots[[1]] <- ggplot(df, aes(var1, var2, fill = fc)) +
      geom_tile(color = ifelse(df$fdr < 0.05 & abs(df$fc) > log2(1.5), "black", "grey"),
                size = ifelse(df$fdr < 0.05 & abs(df$fc) > log2(1.5), 0.5, 0.1)) +
      scale_fill_gradient2(low = "steelblue", high = "indianred3", mid = "white", midpoint = 0, limit = rng, name="log2 FC") +
      xlab("average") +
      ylab("") +
      coord_fixed() +
      theme_minimal() +
      theme(panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            axis.text = element_text(family='sans', face=1, size=10),
            axis.title.x = element_text(family='sans', face=1, size=12))
    
    # somite
    i=2
    for(s in paste0("somite", c("III","II","I"))){
      df <- data.frame(var1=factor(rep(c(8,18,21,25,27), times=c(5:1)), levels = stages),
                       var2 = factor(c(stages[2:6], stages[3:6], stages[4:6], stages[5:6], stages[6]), levels = stages),
                       fc=as.numeric(somite[[s]][,2:16]),
                       fdr=somite[[s]]$FDR)
      df <- df[order(abs(df$fc)),]
      plots[[i]] <- ggplot(df, aes(var1, var2, fill = fc)) +
        geom_tile(color = ifelse(df$fdr < 0.05 & abs(df$fc) > log2(1.5), "black", "grey"),
                  size = ifelse(df$fdr < 0.05 & abs(df$fc) > log2(1.5), 0.5, 0.1)) +
        scale_fill_gradient2(low = "steelblue", high = "indianred3", mid = "white", midpoint = 0, limit = rng, name="log2 FC") +
        xlab(s) +
        ylab("") +
        coord_fixed() +
        theme_minimal() +
        theme(panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              axis.text = element_text(family='sans', face=1, size=10),
              axis.title.x = element_text(family='sans', face=1, size=12))
      i <- i+1
    }
    p <- ggarrange(plotlist = plots, ncol = 4, nrow = 1, legend = "none")
    return(p)
  }
}



#### peak-gene correlation
scatter_peak_gene <- function(sel=1, colour_by="stage"){
  ## set colour palette
  if(colour_by == 'somite'){
    cols <- cols.somite
  }else{
    cols <- cols.stage
  }
  
  ## get data
  link <- links[sel,]
  df <- data.frame(peak = peakCounts[link$peak, shared_samples],
                   gene = as.numeric(geneCounts[geneCounts$gene == link$gene, shared_samples]),
                   stage = meta.rna[match(shared_samples, meta.rna$sample),'stage'],
                   somite = meta.rna[match(shared_samples, meta.rna$sample),'somite'],
                   col = meta.rna[match(shared_samples, meta.rna$sample),colour_by])
  
  # linear regression
  fit <- lm(gene ~ peak, data = df)
  
  p <- plot_ly(data = df,
               type = "scatter",
               x = ~peak, y = ~gene,
               size=2,
               color = ~col, colors = cols,
               hoverinfo = 'text',
               text = ~paste('stage:', stage, '\nsomite: ', somite),
               width = 550, height = 350
  ) %>%
    add_lines(x = ~peak, y = fitted(fit), line = list(width = 2, color="grey"))  %>%
    layout(autosize  = FALSE,
           title = paste(link$peak, link$gene, sep=" -- "),
           xaxis = c(list(title = "peak accessibility"), showline= T, zeroline = FALSE),
           yaxis = c(list(title = "gene expression", showline=T, zeroline = FALSE)),
           showlegend = FALSE
    )
  
  return(p)
}





##################################################
### prepare environment
# dir <- "/Users/xi629080/Documents/somitogenesis2022/"
# dir <- "/Users/ibarra01/OneDrive - CRUK Cambridge Institute/github/somitogenesis2020/"

# ## load data
# # metadata
# meta.rna <- read.table(paste0(dir, "RNA-seq/data/metadata_RNAseq.tsv"),
#                        stringsAsFactors = FALSE, header = TRUE)
# meta.rna <- meta.rna[meta.rna$QC == 1 & meta.rna$wrongStage == 0,]
# meta.rna$stage <- factor(paste0("stage", meta.rna$stage), levels = paste0("stage", c(8,18,21,25,27,35)))
# meta.rna$somite <- factor(meta.rna$somite, levels=paste0("S", c("III", "II", "I")))
# meta.rna$somiteNumber <- as.factor(meta.rna$somiteNumber)
# 
# meta.atac <- read.table(paste0(dir, "ATAC-seq/data/metadata_ATACseq_GQ.tsv"),
#                         stringsAsFactors = FALSE, header = TRUE)
# meta.atac <- meta.atac[meta.atac$QCpass==1,]
# meta.atac$somiteNumber <- ifelse(meta.atac$somite == "SI", meta.atac$stage,
#                                  ifelse(meta.atac$somite == "SII", meta.atac$stage-1, meta.atac$stage-2))
# meta.atac$stage <- factor(paste0("stage", meta.atac$stage), levels = paste0("stage", c(8,18,21,25,27,35)))
# meta.atac$somite <- factor(meta.atac$somite, levels=paste0("S", c("III", "II", "I")))
# meta.atac$somiteNumber <- as.factor(meta.atac$somiteNumber)
# 
# ## expression data
# geneCounts <- read.table(paste0(dir, "RNA-seq/data/geneCounts.NORM_batchCorrected_14PCs.tsv"),
#                          check.names = FALSE, stringsAsFactors = FALSE)
# geneCounts <- geneCounts[,which(colnames(geneCounts) %in% c("gene", meta.rna$sample)),]
# stopifnot(identical(colnames(geneCounts)[-1], meta.rna$sample))
# 
# ## accessibility data
# peakCounts <- readRDS(paste0(dir, "ATAC-seq/results/04_peakCounts_csawMerged.NORM.batchCorrected_18PCs.Rds"))
# stopifnot(identical(colnames(peakCounts), meta.atac$sample))
# 
# ## peak ranges
# tmp <- unlist(lapply(strsplit(row.names(peakCounts), ":"), '[[', 2))
# peaks <- GRanges(seqnames = unlist(lapply(strsplit(row.names(peakCounts), ":"), '[[', 1)),
#                  IRanges(as.numeric(unlist(lapply(strsplit(tmp, "-"), '[[', 1))),
#                          as.numeric(unlist(lapply(strsplit(tmp, "-"), '[[', 2)))),
#                  name = row.names(peakCounts))
# 
# ## gene annotation
# gene_ann <- read.table("RNA-seq/data/Mus_musculus.GRCm38.96.ann", row.names = 1)
# colnames(gene_ann) <- c("gene", "chr", "start", "end", "strand", "biotype")
# 
# 
# ## DE results
# # somite trios
# tmp <- list()
# for(contrast in paste0("somite", c("IvsII", "IIvsIII", "IvsIII"))){
#   tmp[[contrast]] <- read.table(paste0(dir, "RNA-seq/results/04_DEresults_somiteTrios_", contrast, ".tsv"), stringsAsFactors = FALSE)
#   tmp[[contrast]] <- tmp[[contrast]][order(tmp[[contrast]]$genes),]
# }
# # summarise in single table mathing per-stage structure
# DEtriosAll <- data.frame(genes = tmp[[1]]$genes,
#                          logFC.IIvsI = tmp[['somiteIvsII']]$logFC,
#                          logFC.IIIvsI = tmp[['somiteIvsIII']]$logFC,
#                          logFC.IIIvsII = tmp[['somiteIIvsIII']]$logFC,
#                          logCPM = tmp[[1]]$logCPM,
#                          F = NA,
#                          PValue = NA,
#                          FDR = pmin(tmp[['somiteIvsII']]$FDR, tmp[['somiteIvsIII']]$FDR, tmp[['somiteIIvsIII']]$FDR))
# row.names(DEtriosAll) <- row.names(tmp[[1]])
# DEtriosAll$maxFC <- apply(DEtriosAll[,2:4], 1, function(x) x[which.max(abs(x))])
# DEtriosAll$DE <- ifelse(DEtriosAll$FDR < 0.05 & abs(DEtriosAll$maxFC) > log2(1.5), 1, 0)
# DEtriosAll <- DEtriosAll[order(DEtriosAll$FDR),]
# 
# DEtriosStage <- list()
# for(contrast in paste0("stage", c(8,18,21,25,27,35))){
#   DEtriosStage[[contrast]] <- read.table(paste0(dir, "RNA-seq/results/04_DEresults_somiteTrios_", contrast, ".tsv"), stringsAsFactors = FALSE)
# }
# 
# DEtrios_summary <- read.table("RNA-seq/results/04_DEresults_summary_somiteTrios.tsv")
# DEtrios_summary <- DEtrios_summary[rowSums(DEtrios_summary)>0,]
# 
# # stages
# DEstageAll <- read.table(paste0(dir, "RNA-seq/results/04_DEresults_stage_average.tsv"), stringsAsFactors = FALSE)
# 
# DEstageSomite <- list()
# for(contrast in paste0("somite", c("I", "II", "III"))){
#   DEstageSomite[[contrast]] <- read.table(paste0(dir, "RNA-seq/results/04_DEresults_stage_",contrast,".tsv"), stringsAsFactors = FALSE)
# }
# 
# DEstage_summary <- read.table("RNA-seq/results/04_DEresults_summary_stage.tsv")
# DEstage_summary <- DEstage_summary[rowSums(DEstage_summary)>0,]
# 
# 
# ## chromatin deviation
# dev <- readRDS(paste0(dir, "ATAC-seq/results/06_motif_deviationScores.Rds"))
# chromvar <- as.data.frame(cbind(gene = rowData(dev)$name, assay(dev, 'z')))
# 
# 
# ## peak-gene correlations
# links <- read.table("RNA+ATAC/results/03_gene_peak_links.tsv", header = TRUE)[,c(2:5,7,9:10)]
# colnames(links) <- c("peak", "gene", "pearson_r", "p_value", "DE", "DA", "distance")
# links$pearson_r <- round(links$pearson_r, 3)
# links$p_value <- round(links$p_value  , 4)
# links <- links[,c("peak", "gene", "distance", "pearson_r", "p_value", "DE", "DA"), ]
# links <- links[order(links$pearson_r, decreasing = TRUE),]
# fig.d <- readRDS("RNA+ATAC/results/03_FigR_network.Rds")
# 
# save(meta.rna, meta.atac, geneCounts, peakCounts, gene_ann, peaks,
#      DEtriosAll, DEtriosStage, DEtrios_summary,
#      DEstageAll, DEstageSomite, DEstage_summary,
#      chromvar, links,
#      file=paste0(dir, "shiny/data/data.RData"))



