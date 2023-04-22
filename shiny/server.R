# server

shinyServer(
  function(input, output, session) {
    
    #### gene expression
    updateSelectizeInput(session = session, inputId = 'gene', choices = c(Choose = '', genesAvailable), 
                         selected = "Hoxa1", server = TRUE)
    
    output$plotGeneExpr <- renderPlotly({
      boxplotExpr(group_by=input$group, colour_by=input$colour, gene=input$gene)
    })
    
    output$summaryDEtrios <- renderPlot({
      summaryDEtrios(gene=input$gene)
    }, height=90, width=620)

    output$summaryDEstages <- renderPlot({
      summaryDEstages(gene=input$gene)
    }, height=120, width=655)

    output$legend <- renderImage({
      outfile <- tempfile(fileext='.png')
      p <- readPNG("images/summaryDE_legend_horiz.png")
      writePNG(p, target=outfile)
      list(src = outfile, contentType = 'image/png', width = 255, height = 45)
    })

    
    #### peak accessibility
    updateSelectizeInput(session = session, inputId = 'geneATAC', choices = c(Choose = '', genesAvailable), 
                         selected = "Hoxa1", server = TRUE)
    
    output$tracks <- renderPlot({
      if(input$inputATAC == "gene"){
        validate(
          need(input$geneATAC %in% geneCounts$gene, 'Please type or select a gene in the dropdown menu.')
        )
      }
      if(input$inputATAC == "region"){
        region <- unlist(strsplit(input$regionATAC, "[:-]"))
        validate(
          need(as.numeric(region[3]) > as.numeric(region[2]), 'End coordinate must be greater than start.')
        )
      }

      plotPeaks(type=input$inputATAC, gene=input$geneATAC, region=input$regionATAC, flank=input$flank)
    }, width = 1100, height = 260)

    
    output$peaks_in_range <- renderUI({
      if(input$inputATAC == "gene"){
        validate(
          need(input$geneATAC %in% geneCounts$gene, '')
        )
        ## get gene's genomic location
        ann <- gene_ann[gene_ann$gene==input$geneATAC,]
        range <- GRanges(paste0("chr", ann$chr),
                          IRanges(ann$start - input$flank*1e3,
                                  ann$end + input$flank*1e3),
                          gene = ann$gene)
      }else{
        region <- unlist(strsplit(input$regionATAC, "[:-]"))
        validate(
          need(region[3] > region[2], 'End coordinate must be greater than start.')
        )
        ## convert to GRanges
        range <- GRanges(region[1],
                         IRanges(as.numeric(region[2]) - input$flank*1e3,
                                 as.numeric(region[3]) + input$flank*1e3))
      }
      ## find overlapping peaks
      overlaps <- peaks[subjectHits(findOverlaps(range, peaks))]$name
      
      if(length(overlaps)>0){
        choiceList = paste0("Peak_", 1:length(overlaps), ": ", overlaps)
      }else{
        choiceList = "No peaks in range"
      }
      
      selectInput(inputId = 'peaks_in_range_list', label = "", 
                  choices = c(Choose = '', choiceList),
                  selectize = TRUE, selected = 1)
    })
    
    output$plotPeakAcc <- renderPlotly({
      if(input$inputATAC == "gene"){
        validate(
          need(input$geneATAC %in% geneCounts$gene, '')
        )
        ## get gene's genomic location
        ann <- gene_ann[gene_ann$gene==input$geneATAC,]
        range <- GRanges(paste0("chr", ann$chr),
                         IRanges(ann$start - input$flank*1e3,
                                 ann$end + input$flank*1e3),
                         gene = ann$gene)
      }else{
        region <- unlist(strsplit(input$regionATAC, "[:-]"))
        validate(
          need(as.numeric(region[3]) > as.numeric(region[2]), 'End coordinate must be greater than start.')
        )
        ## convert to GRanges
        range <- GRanges(region[1],
                         IRanges(as.numeric(region[2]) - input$flank*1e3,
                                 as.numeric(region[3]) + input$flank*1e3))
      }
      ## find overlapping peaks
      overlaps <- peaks[subjectHits(findOverlaps(range, peaks))]$name
      
      validate(
        need(input$geneATAC %in% geneCounts$gene && length(overlaps)>0, 'Please select a different gene/region to plot accessibility levels.')
      )
      
      boxplotAcc(type=input$inputATAC,
                 group_by=input$groupATAC, colour_by=input$colourATAC, 
                 peak=input$peaks_in_range_list)
    })
    
    
    #### TF activity
    updateSelectizeInput(session = session, inputId = 'tf', choices = c(Choose = '', tfsAvailable), 
                         selected = "Hoxa1", server = TRUE)
    
    output$plotTFactivty <- renderPlotly({
      validate(
        need(input$tf %in% chromvar$gene, 'Please select a TF from the dropdown menu.')
      )
      boxplotChromvar(group_by=input$groupTF, colour_by=input$colourTF, gene=input$tf)
    })
    
    
    #### DE results - trios
    output$DEtriosTable <- DT::renderDataTable(
      datatable( printDEtriosTable(level=input$DEtriosLevel, contrast=input$contrastTrios), 
                 colnames = c("gene", "log2_CPM", "log2_FC", "FDR", "shared_across_stages"), rownames = FALSE,
                 class="compact", options=list(pageLength = 6)
      ), 
      server=TRUE
    )
    
    output$plotGeneExpr_DEtrios <- renderPlotly({
      validate(
        if(input$DEtriosLevel == "average"){
          need(input$DEtriosTable_row_last_clicked <= nrow(DEtrios_summary), '')
        }else{
          need(input$DEtriosTable_row_last_clicked <= nrow(DEtrios_summary), '')
        }
      )
      sel = input$DEtriosTable_row_last_clicked
      if(length(sel)){
        boxplotExpr_DEtrios(level=input$DEtriosLevel, contrast=input$contrastTrios, sel=sel,
                             group_by=input$groupTrios, colour_by=input$colourTrios)
      }
    })
    
    output$summaryDEtrios_tableTrios <- renderPlot({
      validate(
        if(input$DEtriosLevel == "average"){
          need(input$DEtriosTable_row_last_clicked <= nrow(DEtrios_summary), '')
        }else{
          need(input$DEtriosTable_row_last_clicked <= nrow(DEtrios_summary), '')
        }
      )
      sel = input$DEtriosTable_row_last_clicked
      summaryDEtrios_tableTrios(level=input$DEtriosLevel, contrast=input$contrastTrios, sel=sel)
    }, height=90, width=620)
    
    output$summaryDEstages_tableTrios <- renderPlot({
      validate(
        if(input$DEtriosLevel == "average"){
          need(input$DEtriosTable_row_last_clicked <= nrow(DEtrios_summary), '')
        }else{
          need(input$DEtriosTable_row_last_clicked <= nrow(DEtrios_summary), '')
        }
      )
      sel = input$DEtriosTable_row_last_clicked
      summaryDEstages_tableTrios(level=input$DEtriosLevel, contrast=input$contrastTrios, sel=sel)
    }, height=120, width=655)
    
    output$legend_trios <- renderImage({
      outfile <- tempfile(fileext='.png')
      p <- readPNG("images/summaryDE_legend_horiz.png")
      writePNG(p, target=outfile)
      list(src = outfile, contentType = 'image/png', width = 255, height = 45)
    })
    
    
    #### DE results - stages
    output$DEstagesTable <- DT::renderDataTable(
      datatable( printDEstagesTable(level=input$DEstagesLevel, contrast=input$contrastStages), 
                 colnames = c("gene", "log2_CPM", "log2_FC", "FDR", "shared_across_somites"), rownames = FALSE,
                 class="compact", options=list(pageLength = 6)
      ), 
      server=TRUE
    )
    
    output$plotGeneExpr_DEstages <- renderPlotly({
      validate(
        if(input$DEstagesLevel == "average"){
          need(input$DEstagesTable_row_last_clicked <= nrow(DEstage_summary), '')
        }else{
          need(input$DEstagesTable_row_last_clicked <= nrow(DEstage_summary), '')
        }
      )
      sel = input$DEstagesTable_row_last_clicked
      if(length(sel)){
        boxplotExpr_DEstages(level=input$DEstagesLevel, contrast=input$contrastStages, sel=sel,
                             group_by=input$groupStages, colour_by=input$colourStages)
      }
    })
    
    output$summaryDEtrios_tableStages <- renderPlot({
      validate(
        if(input$DEstagesLevel == "average"){
          need(input$DEstagesTable_row_last_clicked <= nrow(DEstage_summary), '')
        }else{
          need(input$DEstagesTable_row_last_clicked <= nrow(DEstage_summary), '')
        }
      )
      sel = input$DEstagesTable_row_last_clicked
      summaryDEtrios_tableStages(level=input$DEstagesLevel, contrast=input$contrastStages, sel=sel)
    }, height=90, width=620)
    
    output$summaryDEstages_tableStages <- renderPlot({
      validate(
        if(input$DEstagesLevel == "average"){
          need(input$DEstagesTable_row_last_clicked <= nrow(DEstage_summary), '')
        }else{
          need(input$DEstagesTable_row_last_clicked <= nrow(DEstage_summary), '')
        }
      )
      sel = input$DEstagesTable_row_last_clicked
      summaryDEstages_tableStages(level=input$DEstagesLevel, contrast=input$contrastStages, sel=sel)
    }, height=120, width=655)
    
    output$legend_stages <- renderImage({
      outfile <- tempfile(fileext='.png')
      p <- readPNG("images/summaryDE_legend_horiz.png")
      writePNG(p, target=outfile)
      list(src = outfile, contentType = 'image/png', width = 255, height = 45)
    })    
    
    
    
    #### peak-gene correlation
    output$linksTable <- DT::renderDataTable(
      datatable( links, 
                 rownames = FALSE,
                 class="compact", options=list(pageLength = 10)
      ), 
      server=TRUE
    )
    
    output$plotPeak_vs_Gene <- renderPlotly({
      sel = input$linksTable_row_last_clicked
      if(length(sel)){
        scatter_peak_gene(sel=sel, colour_by=input$colourLinks)
      }
    })
    
  }
)



