require(shiny)
require(shinyjs)
require(shinythemes)
require(DT)
require(plotly)

panel_width = 3
result_with = 10

appCSS <- ".mandatory_star { color: blue; }"

# Define UI for application that draws a histogram
shinyUI(
  fluidPage(theme = shinytheme("simplex"),
    navbarPage("GSECA",

             tabPanel("Analysis",

                      # titlePanel("Settings"),
                      # h4("Gene Set Enrichment Class Analysis"),
                      # h4("Settings"),
                      # img(src='GSECA.banner.png', align = "top", width="390", height="120"),
                      shinyjs::useShinyjs(),
                      shinyjs::inlineCSS(appCSS),

                      sidebarLayout(
                        sidebarPanel(
                                 fileInput('exp_matrix', 'Expression Matrix (tab-separated file)*', accept = c('text/tab-separated-values', c('.tsv','.gz') ))
                                 , fileInput('sample_labels', 'Sample type labels (tab-separated file)*', accept = c('text/tab-separated-values', c('.tsv','.gz') ))
                                 , selectInput("symbol", 'Gene ID', choices=c( "HUGO symbol"= "symbol"
                                                                               ,"ENSEMBL Gene Id"="ensembl_gene_id"
                                 ), selected="symbol")
                                 # , selectInput("nClass", 'Number of ECs (default = 7)', choices=c(2:10), selected=7)
                                 , sliderInput("nClass", 'Number of ECs (default = 7)', min=2, max=10, value = 7)
                                 , textInput("analysis", "Analysis Title:", value = "GSECA")
                                 , selectInput("gs_dataset", "Select MsigDB GeneSet:",choices = list('Hallmark (H)'='h.all.v6.0.symbols.gmt',
                                                                                                     'Positional (C1)'='c1.all.v6.0.symbols.gmt',
                                                                                                     'Curated (C2)'='c2.all.v6.0.symbols.gmt',
                                                                                                     'Chemical and genetic pertubations (C2)'='c2.cgp.v6.0.symbols.gmt',
                                                                                                     'Biocarta (C2)'='c2.cp.biocarta.v6.0.symbols.gmt',
                                                                                                     'KEGG (C2)'='c2.cp.kegg.v6.0.symbols.gmt',
                                                                                                     'Reactome (C2)'='c2.cp.reactome.v6.0.symbols.gmt',
                                                                                                     'All canonical pathways (C2)'='c2.cp.v6.0.symbols.gmt',
                                                                                                     'Motif (C3)'='c3.all.v6.0.symbols.gmt',
                                                                                                     'MiRNA (C3)'='c3.mir.v6.0.symbols.gmt',
                                                                                                     'Transcription factors (C3)'='c3.tft.v6.0.symbols.gmt',
                                                                                                     'Computational (C4)'='c4.all.v6.0.symbols.gmt',
                                                                                                     'Cancer gene neighborhoods (C4)'='c4.cgn.v6.0.symbols.gmt',
                                                                                                     'Cancer Modules (C4)'='c4.cm.v6.0.symbols.gmt',
                                                                                                     "GO (C5)"='c5.all.v6.0.symbols.gmt',
                                                                                                     "GO BP (C5)"='c5.bp.v6.0.symbols.gmt',
                                                                                                     "GO CC (C5)"='c5.cc.v6.0.symbols.gmt',
                                                                                                     "GO MF (C5)"='c5.mf.v6.0.symbols.gmt',
                                                                                                     'Oncogenic signatures (C6)'='c6.all.v6.0.symbols.gmt',
                                                                                                     'Immunological signatures (C7)'='c7.all.v6.0.symbols.gmt')
                                               ,selected = 'c2.cp.kegg.v6.0.symbols.gmt')

                                , selectInput("stat.test", 'Statistical Test', choices=c( "Fisher's Exact Test"='fisher'
                                                                                        ,"Chi-Square Test"='chisq'
                                                                                        ), selected='fisher')
                                
                                , selectInput("correction", 'Multiple test correction', choices=c( "Bonferroni"="bonferroni"
                                                                                          ,"FDR"='fdr'), selected='fdr')
                                , textInput("p.adj", "Adjusted p-value cutoff",value = 0.05)

                                # , textInput("pw_sim", "Power - Number of simulation:",value = 100)
                                #  CPUs=3
                                #
                                 , checkboxInput("customGS", label = tags$b("Use Custom Gene Set"), value = FALSE)
                                 , conditionalPanel(condition = "input.customGS == true",fileInput('gene_set', label = NULL,accept = c('Gene Matrix Transposed','.gmt'))),
                                 # radioButtons("genome","Genome for Monte Carlo",choices=c("HG19","HG38"),selected="HG19",inline = TRUE),
                                 selectInput("empirical", 'Empirical P-Value', choices=c("True","False"), selected="False"),
                                 selectInput("bootstrapping", 'Bootstrapping Sample Size', choices=c("True","False"), selected="False"),
                                 numericInput("nsim", "Number of random sampling (Monte Carlo or Bootstrapping)", value=1000),
                                 p("Mandatory fields are marked with *"),
                                 actionButton("submit", "Run GSECA", class = "btn-primary"),
                                 downloadButton('downloadData', 'Download Results'),
                                 p("")
                                 # p("Download the example cohorts A and B:"),
                                 # downloadButton('downloadExample', 'Download Example Datasets')
                               ),

                        # MAIN PANEL
                        mainPanel(
                          img(src='GSECA.banner.png', align = "top", width="450", height="130"),
                          tabsetPanel(
                              tabPanel('Results',   DT::dataTableOutput('gseca_results'))
                            , tabPanel("EnrichmentClass Map",    plotOutput("gseca_heatmap", width = "100%"))
                            # , tabPanel("Score",      plotOutput("gseca_score",height = "1000px"))
                          )
                         )
                       )
             ),

             tabPanel("Help",
                      titlePanel("Help page"),
                      includeMarkdown("../README.md")
             )

    )
  )
)
