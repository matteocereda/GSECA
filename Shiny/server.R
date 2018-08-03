source("../Scripts/config.R")

shinyInput <- function(FUN, len, id, ...) {
  inputs <- character(len)
  for (i in seq_len(len)) {
    inputs[i] <- as.character(FUN(paste0(id, i), ...))
  }
  inputs
}

options(shiny.maxRequestSize = 30*1024^2)
options(scipen=3)

shinyServer(function(input, output) {
      
  current = reactiveValues(res =NULL)
      
      RUN_GSECA = eventReactive(input$submit, {
              setwd("../")
              M             = read.delim(input$exp_matrix$datapath, header=T, stringsAsFactors = F)
              L             = read.delim(input$sample_labels$datapath, header=T, stringsAsFactors = F)[,1]
              symbol        = input$symbol
              nClass        = input$nClass
              gene.set.path = ifelse(!input$customGS,get.gene.set(input$gs_dataset),input$gene_set$datapath)
              geneset       = read.gmt.file(gene.set.path)
              s.test        = "fisher" #input$stat.test
              correction    = input$correction
              p_adj_th      = input$p.adj
              empirical     = input$empirical=="True"
              bootstrapping = input$bootstrapping=="True"
              analysis      = input$analysis
              nsim          = input$nsim
              outdir        = "../Results"
              
              cpus = detectCores()
              
              if (!is.null(cpus)) 
                {
                  cpus = cpus - 1
                } else {
                  cpus = 3
                }
              GSECA_executor( M          = M
                            , L          = L
                            , symbol     = symbol
                            , geneset    = geneset
                            , s.test     = s.test
                            , correction = correction 
                            , p_adj_th   = p_adj_th
                            , nClass     = nClass
                            , analysis   = analysis
                            , outdir     = outdir
                            , N.CORES    = cpus
                            , EMPIRICAL  = empirical
                            , BOOTSTRP   = bootstrapping
                            , nsim       = nsim )
        
                                                
      })
            
      # TAB results
      output$gseca_results   = DT::renderDataTable({
        results = RUN_GSECA()
        # print(head(results))
        current$res   = results[['gseca']][,]
        current$heat  = results[['ECmap']]
        # current$score =  results[['score']]
        # colnames(results) = toupper(colnames(results))
        DT::datatable(current$res, options = list(orderClasses = TRUE, pageLength = 20),selection= 'single',rownames = FALSE)
        # %>%
          # formatRound(columns=c('PERC','or','pw'), digits=2) %>% formatSignif(columns=c('pv','p.adj'), digits=3)
      })

      # TAB grapihics

      output$gseca_heatmap = renderPlot({
        # browser()
        # p= ggplotly(current$heat)
        # plot_ly(p)
        # grid.draw(current$heat)
        grid.newpage()
        pushViewport(viewport(layout=grid.layout(nrow=4, ncol=6
                                                 , heights=c(.05, .8, .1, .05)
                                                 , widths= c( .025, .3,.2,.2,.2, .025 ))))
        
        # grid.text("GSECA", vp = viewport(layout.pos.row = 1, layout.pos.col = 1:6), gp=gpar(fontface='bold'), just="center")
        
        pushViewport(viewport(layout.pos.row = 2, layout.pos.col = 2:5))
        
        # sink(file='tmp.txt')
        draw(current$heat, newpage = FALSE)
        
        upViewport()
        # sink()
        
        cm = ColorMapping(name = "FDR", col_fun = colorRamp2(c(0, 0.05, 0.1), c("red", "yellow", "blue")))
        pushViewport(viewport(layout.pos.row = 3, layout.pos.col = 2))
        color_mapping_legend(cm, legend_direction='horizontal')
        upViewport()
        
        # m = round(limits[2]*100,0)
        # legLabels <- as.character(c(0.01, round(m/2,1) , m))
        # boxSize <- unit(0.15, "inches")
        
        pushViewport(viewport(layout.pos.row = 3, layout.pos.col = 3:5))
        grid.text(expression(bold(paste(Delta, "(%)"))), x = unit(0.1, "native"), y =  unit(0.5, "native")+ unit(.9, "lines"), gp=gpar(fontface='bold') )

        upViewport()
      
        grid.text("GSECA v.1 2018", vp = viewport(layout.pos.row = 4, layout.pos.col = 1:4), just="center")

      })
  })
