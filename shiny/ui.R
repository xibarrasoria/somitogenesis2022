# ui
library(shiny)
library(shinycssloaders)
library(DT)
library(plotly)
library(ggplot2)
library(RColorBrewer)
library(ggpubr)
library(png)
library(GenomicRanges)
library(SummarizedExperiment)
source("helper.R")

shinyUI(fluidPage(
  titlePanel(h4("Molecular pofiling of mouse somites")),
  tabsetPanel(
    tabPanel("About",
             fluidRow(
               column(12,
                      br(),
                      HTML( "<p style='font-size:16px'>Website accompanying the paper 
                      <b><i>A transcriptional and regulatory map of mouse somite maturation,</b></i>
                      by Ibarra-Soria, Thierion et al., <i>Developmental Cell</i>, 2023. </p>" ),
                      hr(),
                      HTML( "<p>Raw and processed data are available in the ArrayExpress repository and can be accessed through the 
                      <a href='https://www.ebi.ac.uk/biostudies/'>BioStudies database</a> under accession numbers 
                      <a href='https://www.ebi.ac.uk/biostudies/arrayexpress/studies/E-MTAB-12511'>E-MTAB-12511</a> for the RNA-seq dataset
                      <a href='https://www.ebi.ac.uk/biostudies/arrayexpress/studies/E-MTAB-12539'>E-MTAB-12539</a> for the ATAC-seq dataset. 
                      Both FASTQ files and count tables are available.</p>
                      <p>Relevant metadata is available with the paper in Table S1.</p>" ),
                      HTML( "<p>Code associated with this study is available in 
                            <a href='https://github.com/xibarrasoria/somitogenesis2022'>GitHub</a>.</p>" ),
                      hr()
               )
             ),
             fluidRow(
               column(10,
                      HTML( "<p style='font-size:16px'><b>How to use the website </b></p>" ),
                      HTML( "<p style='font-size:15px'>
                            <li><b>For all plots:</b>
                            <ul class='square'>
                            <li>It is possible to select the colour to switch between somite trios and developmental stages.</li>
                            <li>For boxplots, it is also possible to select how to group samples for each box. The combination of both
                            group and colour options allow visualising dynamics driven by somite maturation or developmental progression.</li>
                            <li>All values are log2 CPM (counts-per-million) of batch-corrected data, except for the TF activity tab.</li>
                            <li>On hover, a menu is shown on the top right corner. The first icon (camera) allows downloading the plot.</li>
                            </ul></li>
                            <li><b>For all tables:</b>
                            <ul class='square'>
                            <li>Clicking on a row will display the corresponding gene expression plot.</li>
                            <li>Clicking on a column header sorts the table by ascending values of that column. Clicking again, sorts in
                            descending order.</li>
                            </ul></li>
                            </p>"),
                      hr(),
                      HTML( "<p style='font-size:15px'><b>Tab description</b></p>
                            <li><b>Gene expression: </b>type or select a gene to plot its expression across samples.
                            <ul class='square'>
                            <li>Below, a summary representation of results from differential expression analyses between
                            the trios and across stages. The grids correspond to each pairwise comparison and the colour indicates
                            the fold-change.</li>
                            <li>Outlined tiles indicate that the difference is statistically significant.</li>
                            <li>Results for average and per-stage/per-somite results are shown.</li>
                            <li>The colour scale is shared across all <i>Trios</i> and <i>Stages</i> plots.</li>
                            </ul></li>" ),
                      HTML( "<li><b>Peak accessibility: </b>input a gene or a genomic region (chr:start-end) to
                            display a schematic showing the peaks found in the region and their genomic context. 
                            Use the drop-down menu to select a peak and plot accessibility levels across samples.</li>" ),
                      br(),
                      HTML( "<li><b>TF activity: </b>type or select a transcription factor to plot the chromatin deviation scores
                            computed by chromVar.</li>" ),
                      br(),
                      HTML( "<li><b>DE results tabs: </b>explore the results from differential expression analyses between somite 
                            <b>trios</b> or across <b>stages</b>. On the left-side bar, select whether to show results from testing all 
                            samples together (average) or the per-somite/per-stage results. Click on a row to plot gene expression.
                            <ul class='square'>
                            <li>The column <i>shared_across_stages/somites</i> indicates whether the gene was significant when considering
                            average expression across stages/somites, i.e. whether the gene shows consistent DE results across stages/somites.
                            <p>When exploring the results of stage/somite-specific analyses it can be useful to look at those that were not 
                            significant in the average test. To do this, click on the <i>shared_across_stages/somites</i> column header to bring 
                            genes with 'false' to the top</p>
                            </ul></li>" ),
                      HTML( "<li><b>Peak-gene correlation: </b>explore the peak-to-gene linking results from FigR. If the gene (peak) was 
                            significantly DE, the <i>DE</i> (<i>DA</i>) column indicates the vertebral fate where the element shows highest 
                            activity. Click on a row to plot the peak accessibility versus gene expression levels across samples. </li>" ),
                      br(),
                      br()
                      
               )
             )
    ),
    tabPanel("Gene expression",
      sidebarLayout(
        sidebarPanel(width=2,
                     fluidRow(
                       selectizeInput("gene", label = h5("Gene of interest"), choices = NULL, 
                                      options = list(placeholder = 'type a gene name'))
                     ),
                     fluidRow(
                       radioButtons("group", "group by:", 
                                    choices = list("somite"='somite', "stage"='stage', "somiteNumber"='somiteNumber'), 
                                    selected = 'stage')
                       ),
                     fluidRow(
                       radioButtons("colour", "colour by:", 
                                    choices = list("somite"='somite', "stage"='stage'),
                                    selected = 'somite')
                       )
        ),
        mainPanel(
          fluidRow(
            column(8,
                   div(style = "height:20px"),
                   plotlyOutput("plotGeneExpr")
            )
          ),
          fluidRow(
            div(style = "height:40px"),
            h5("Differential expression summary"),
            column(8,
                   h6("Trios"),
                   plotOutput("summaryDEtrios")
            )
          ),
          fluidRow(
            column(9,
                   div(style = "margin-top:-300px"),
                   h6("Stages"),
                   div(style = "margin-left:-1em",
                       plotOutput("summaryDEstages")
                   ),
                   div(style = "margin-top:-270px"),
                   div(style = "margin-left:15em",
                       imageOutput("legend")
                   )
            )
          )
        )
      )
    ),
    
    
    tabPanel("Chromatin accessibility",
      sidebarLayout(
        sidebarPanel(width=2,
                     fluidRow(
                       radioButtons("inputATAC", "search by:", 
                                    choices = list("gene"='gene', "region"='region'), 
                                    selected = 'gene'),
                       conditionalPanel(
                         condition = "input.inputATAC == 'gene'",
                         selectizeInput("geneATAC", label = h5("Gene of interest"), choices = NULL, 
                                        options = list(placeholder = 'type a gene name'))
                         ),
                       conditionalPanel(
                         condition = "input.inputATAC == 'region'",
                         textInput("regionATAC", label = h5("Genomic location"), 
                                   placeholder = "chr:start-end", value = "chr6:52150590-52163317")
                         )
                       ),
                     fluidRow(
                       numericInput("flank", "Flanking range (kb):", 5, min = 0, max = 250)
                       ),
                     fluidRow(
                       radioButtons("groupATAC", "group by:", 
                                    choices = list("somite"='somite', "stage"='stage', "somiteNumber"='somiteNumber'), 
                                    selected = 'stage')
                       ),
                     fluidRow(
                       radioButtons("colourATAC", "colour by:", 
                                    choices = list("somite"='somite', "stage"='stage'),
                                    selected = 'somite')
                       )
               ),
        mainPanel(
          fluidRow(
            withSpinner( plotOutput("tracks"), type = 6, size = 0.5),
            div(style = "margin-top:-130px"),
            uiOutput("peaks_in_range"),
            conditionalPanel(
              condition = "input.peaks_in_range_list != ''",
              plotlyOutput("plotPeakAcc")
            )
          )
        )
      )
    ),
    
    
    tabPanel("TF activity",
      sidebarLayout(
        sidebarPanel(width=2,
                     fluidRow(
                       selectizeInput("tf", label = h5("TF of interest"), choices = NULL, 
                                      options = list(placeholder = 'type a gene name'))
                       ),
                     fluidRow(
                       radioButtons("groupTF", "group by:", 
                                    choices = list("somite"='somite', "stage"='stage', "somiteNumber"='somiteNumber'), 
                                    selected = 'stage')
                     ),
                     fluidRow(
                       radioButtons("colourTF", "colour by:", 
                                    choices = list("somite"='somite', "stage"='stage'),
                                    selected = 'somite')
                     )
        ),
        mainPanel(
          fluidRow(
            column(8,
                   div(style = "height:20px"),
                   plotlyOutput("plotTFactivty")
            )
          ),
          fluidRow(
            column(10,
                   div(style = "height:40px"),
                   br(),
                   hr(),
                   HTML( "<p>The boxplot shows the z-score of the deviation scores computed by chromVAR, based on the 
                   <b>accessibility levels of all regions in the genome with a predicted binding site for the TF of intrest</b>.</p>
                   <p>Positive (negative) values indicate greater (lower) accessiblity compared to background chromatin regions.</p>" ),
                   hr()
            )
          )
        )
      )
    ),
    
    
    tabPanel("DE results - trios",
      sidebarLayout(
        sidebarPanel(width=2,
                     fluidRow(
                       radioButtons("DEtriosLevel", "DE results by:", 
                                    choices = list("average"='average', "stage"='stage'), 
                                    selected = 'average')
                     ),
                     conditionalPanel(
                       condition = "input.DEtriosLevel == 'stage'",
                       selectInput("contrastTrios", "Stage:", 
                                   choices = list("stage8"='stage8', "stage18"='stage18', "stage21"='stage21', 
                                                  "stage25"='stage25', "stage27"='stage27', "stage35"='stage35'), 
                                   selected = 'stage8')
                     ),
                     fluidRow(
                       radioButtons("groupTrios", "plot by:", 
                                    choices = list("somite"='somite', "stage"='stage'), 
                                    selected = 'stage')
                     ),
                     fluidRow(
                       radioButtons("colourTrios", "colour by:", 
                                    choices = list("somite"='somite', "stage"='stage'),
                                    selected = 'somite')
                     )
        ),
        mainPanel(
          fluidRow(
                   div(style = "height:10px"),
                   h6('Click on a row to plot gene expression'),
                   DT::dataTableOutput("DEtriosTable")
          ),
          conditionalPanel(
            condition = "input.DEtriosTable_row_last_clicked",
            fluidRow(
              column(8,
                     div(style = "height:20px"),
                     plotlyOutput("plotGeneExpr_DEtrios")
              )
            ),
            fluidRow(
              div(style = "height:40px"),
              h5("Differential expression summary"),
              column(8,
                     h6("Trios"),
                     plotOutput("summaryDEtrios_tableTrios")
              )
            ),
            fluidRow(
              column(9,
                     div(style = "margin-top:-300px"),
                     h6("Stages"),
                     div(style = "margin-left:-1em",
                         plotOutput("summaryDEstages_tableTrios")
                     ),
                     div(style = "margin-top:-270px"),
                     div(style = "margin-left:15em",
                         imageOutput("legend_trios")
                    )
              )
            )
          )
        )
      )
    ),
    
    
    tabPanel("DE results - stages",
      sidebarLayout(
        sidebarPanel(width=2,
                     fluidRow(
                       radioButtons("DEstagesLevel", "DE results by:", 
                                    choices = list("average"='average', "somite"='somite'), 
                                    selected = 'average')
                     ),
                     conditionalPanel(
                       condition = "input.DEstagesLevel == 'somite'",
                       selectInput("contrastStages", "Somite:", 
                                   choices = list("somiteI"='somiteI', "somiteII"='somiteII', "somiteIII"='somiteIII'), 
                                   selected = 'somiteI')
                     ),
                     fluidRow(
                       radioButtons("groupStages", "plot by:", 
                                    choices = list("somite"='somite', "stage"='stage'), 
                                    selected = 'stage')
                     ),
                     fluidRow(
                       radioButtons("colourStages", "colour by:", 
                                    choices = list("somite"='somite', "stage"='stage'),
                                    selected = 'stage')
                     )
        ),
        mainPanel(
          fluidRow(
            div(style = "height:10px"),
            h6('Click on a row to plot gene expression'),
            DT::dataTableOutput("DEstagesTable")
          ),
          conditionalPanel(
            condition = "input.DEstagesTable_row_last_clicked",
            fluidRow(
              column(8,
                     div(style = "height:20px"),
                     plotlyOutput("plotGeneExpr_DEstages")
              )
            ),
            fluidRow(
              div(style = "height:40px"),
              h5("Differential expression summary"),
              column(8,
                     h6("Trios"),
                     plotOutput("summaryDEtrios_tableStages")
              )
            ),
            fluidRow(
              column(9,
                     div(style = "margin-top:-300px"),
                     h6("Stages"),
                     div(style = "margin-left:-1em",
                         plotOutput("summaryDEstages_tableStages")
                     ),
                     div(style = "margin-top:-270px"),
                     div(style = "margin-left:15em",
                         imageOutput("legend_stages")
                    )
              )
            )
          )
        )
      )
    ),
    
    
    tabPanel("Peak-gene correlation",
      sidebarLayout(
        sidebarPanel(width=2,
                     fluidRow(
                       radioButtons("colourLinks", "colour by:", 
                                    choices = list("somite"='somite', "stage"='stage'),
                                    selected = 'stage')
                     )
        ),
        mainPanel(
          fluidRow(
            div(style = "height:20px"),
            DT::dataTableOutput("linksTable")
          ),
          fluidRow(
            div(style = "height:30px"),
            plotlyOutput("plotPeak_vs_Gene")
          )
        )
        )
      
    )
    
    
  )
))

