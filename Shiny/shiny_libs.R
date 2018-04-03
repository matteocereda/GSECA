library(DT)

get.gene.set = function(ix) {
  paste0("../gene_sets/msigdb/",ix)
}

get_mixture_expr_class_parqu=function(data
                                      , normalization = "FPKM" # FPKM. RSEM, TPM, RPKM
                                      , ne_value=0.01
                                      , method="H_quantiles" # means, H_quantiles, quantiles
                                      , lvls=E.LEVELS
                                      , nc=2
                                      , set.cores = "manual"
){

  suppressPackageStartupMessages(library(plyr))
  suppressPackageStartupMessages(library(mixtools))
  suppressPackageStartupMessages(library(parallel))
  # suppressPackageStartupMessages(library(snow))

  rna=data
  rna$log_value=NA
  rna$log_value[which(rna$value>ne_value)]=log2(rna$value[which(rna$value>ne_value)])

  l=dlply(rna, ~Tumor_Sample_Barcode)

  # nclust = ncl
  # print(nclust)
  if(set.cores == "auto"){
    ncor <- detectCores()
    print(paste0("Parallelizing, ",ncor," cores detected ... "))
  } else if (set.cores == "manual") {
    ncor <- nc
    print(paste0("Parallelizing, ",ncor," cores ... "))
  }

  clus <- makeCluster(ncor)
  clusterEvalQ(clus, library(mixtools))
  clusterExport(clus, c("l","GMM"),envir=environment())
  print("Getting thresholds ...")


  S1= -4
  S2= 4

  if (normalization=="FPKM" & ne_value==0.1){

    S1= -4
    S2= 4

  } else if (normalization=="FPKM" & ne_value==0.01){

    S1= -5
    S2= 2.8

  } else if(normalization=="RSEM"){
    S1= 0
    S2= 8

  }else if(normalization=="TPM"){


  }else if(normalization=="RPKM"){


  }

  l = parLapply(clus, l, function(x) tryCatch(GMM(x, s1=S1, s2=S2), error = function(x) return(NA)) )

  stopCluster(clus)

  print(paste0("Samples with errors: ", which(sapply(l, length)==1)))

  rna = do.call(rbind, l)

  print("Setting expression classes ...")
  rna$expr_class=NA

  if(method=="means" & length(lvls)==4){

    rna$expr_class[which(rna$log_value > rna$mean2)] = lvls[1]
    rna$expr_class[which(rna$log_value <= rna$mean2 & rna$log_value > rna$mean1 )] = lvls[2]
    rna$expr_class[which(rna$log_value <= rna$mean1)] = lvls[3]
    rna$expr_class[which(rna$value <= ne_value)] = lvls[4]

  } else if(method=="H_quantiles" & length(lvls)==4){

    rna$expr_class[which(rna$log_value > rna$H_th2)] = lvls[1]
    rna$expr_class[which(rna$log_value <= rna$H_th2 & rna$log_value > rna$H_th1 )] = lvls[2]
    rna$expr_class[which(rna$log_value <= rna$H_th1)] = lvls[3]
    rna$expr_class[which(rna$value <= ne_value)] = lvls[4]

  } else if(method=="quantiles" & length(lvls)==5 ){

    rna$expr_class[which(rna$log_value > rna$H_th2)] = lvls[1]
    rna$expr_class[which(rna$log_value <= rna$H_th2 & rna$log_value > rna$H_th1 )] = lvls[2]
    rna$expr_class[which(rna$log_value <= rna$H_th1 & rna$log_value > rna$L_th2)] = lvls[3]
    rna$expr_class[which(rna$log_value <= rna$L_th2 & rna$log_value > rna$L_th1)] = lvls[4]
    rna$expr_class[which(rna$log_value <= rna$L_th1 | rna$value <= ne_value)] = lvls[5]

  }else{
    stop(message("get_mixture_expr_class_parqu: method and classes do not match"))
  }
  return(rna)
}

get_summary_expr_class_parallel = function(  pl
                                             , r
                                             , lvls   = E.LEVELS
                                             , cls    = S.LABELS
                                             , method = 'fisher' # 'chisq', 'wilcox'
                                             , PW_SIM = 100
                                             , TAIL = "two.sided"
                                             , POWER = 0.9,
                                             nc=3
                                             , set.cores = "manual"
){
  suppressPackageStartupMessages(library(parallel, verbose=F))
  # retaining only genes in pathways
  print("Expression Classes:")
  print(lvls)
  print("Cohorts:")
  print(cls)

  total_genes = unique(unlist(pl))

  rna = subset(r, symbol%in%total_genes)

  pll = mapply( function(x,y) data.frame("gene"=x, "gene_set"=y), pl, names(pl), SIMPLIFY = F)

  if(set.cores == "auto"){
    ncor <- detectCores()
    print(paste0("Parallelizing, ",ncor," cores detected ... "))
  } else if (set.cores == "manual") {
    ncor <- nc
    print(paste0("Parallelizing, ",ncor," cores ... "))
  }

  clus <- makeCluster(ncor)
  clusterEvalQ(clus, library(plyr) )
  clusterEvalQ(clus, library(statmod) )
  clusterEvalQ(clus, library(pwr) )

  if(method=='fisher'){
    clusterExport(cl=clus, c("rna",'lvls','cls', 'method',  "TAIL", "PW_SIM", "FISHER"), envir = environment())
  }else if( method=='chisq'){
    clusterExport(cl=clus, c("rna",'lvls','cls', 'method',  "TAIL", "PW_SIM", "CHISQ"), envir = environment())
  }

    print("SEA ...")

  res = parLapply(clus
                  , pll
                  , get_summary_expr_class
                  , r=rna
                  , lvls   = E.LEVELS
                  , cls    = S.LABELS
                  , method = method
                  , TAIL= TAIL
                  , PW_SIM = PW_SIM
                  , POWER = POWER
  )

  stopCluster(clus)
  res
}

heatmap_GSECA_shiny = function( gseca, filename=NULL, p_adj=0.1, pow=0.9, sc = NULL ){

  gseca$gene_set =fix_names(gseca$gene_set)

  score = unique(gseca[,c('gene_set','score')])
  score = score[order(score[,2],decreasing = T),]
  score$gene_set = fix_names(score$gene_set)
  
  cat("Trhesholds: \n")
  cat(paste0("P-value ",p_adj,"\n"))
  cat(paste0("Power ",pow,"\n"))

  gseca$gene_set = factor(gseca$gene_set, levels=rev(score$gene_set))
  score$gene_set = factor(score$gene_set, levels=rev(score$gene_set))

  sel = unique(subset(gseca, p.adj<=p_adj & pw>=pow )$gene_set)

  if(!is.null(sc)){
    qt = quantile(score$score, seq(0,1,0.05))
    sel_score = score$gene_set[which(score$score>=qt[paste0(sc,'%')])]
    sel = intersect(sel, sel_score)
  }
  toplot = subset(gseca, gene_set%in%sel)[,c('gene_set','class','type','PERC','p.adj','direction')]

  toplot = ddply(toplot, .(gene_set,class), summarise,
                 delta=PERC[type=='CASE']-PERC[type=='CNTR']
                 ,p.adj=p.adj[type=='CASE']
                 ,direction=direction[type=='CASE']
  )

  toplot$p.adj[which(toplot$p.adj>p_adj)]=NA

   library(RColorBrewer)

  myPalette <-
    colorRampPalette(c("#00007F", "blue",
                       "#007FFF", "cyan",
                       "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))

  myPalette <-colorRampPalette(brewer.pal(9, "RdYlBu"))

  require(scales)
  require(gtable)


  ## To format heatmap y.ticks with appropriate width (5 chars),
  ## to align with gg_rows y.tics
  ytickform <- function(x){
    lab <- sprintf("%05s",x)
  }

  heat =ggplot(toplot, aes(x=class,y=gene_set, fill=p.adj)) +
    geom_tile(fill = "white",colour = "black") + #
    theme_bw(20) +
    coord_equal() +
    theme(
      legend.position = 'left',
      panel.grid.major = element_blank(), panel.grid.minor = element_blank(),  panel.border = element_blank(),  panel.background = element_blank(),
      axis.text.x = element_text(angle=90, vjust=0.5,hjust=1)) + xlab("") +ylab("") +
    geom_point(data=subset(toplot, delta>=0), aes(size=delta),shape=24,colour="black",na.rm=TRUE,show.legend = T) +
    geom_point(data=subset(toplot, delta<0), aes(size=abs(delta)),shape=25,colour="black",na.rm=TRUE,show.legend = T) +
    scale_fill_gradientn(colours = myPalette(5),na.value="white")+
    scale_size_continuous(range = c(1,5)
                          #,limits=c(0,6),breaks=c(0,1,2), labels=c("0","1","2")
    )

  barscore =ggplot(subset(score, gene_set%in%unique(toplot$gene_set)), aes(y=score,x=gene_set))+
    geom_bar(stat='identity',fill='lightgrey', col='black',width = .9)+theme_bw()+coord_flip()+
    # theme(axis.text.y = element_blank())+
    xlab('')

  return(list('heat'=heat,'barscore'=barscore))

  # g <- ggplotGrob(gg_hm )
  # g <- gtable_add_cols(g, unit(3,"cm"))
  # g <- gtable_add_grob(g, ggplotGrob(gg_rows),
  #                      t = 1, l=ncol(g), b=nrow(g)-1, r=ncol(g))

   # grid.draw(g)

}


GSECA_executor = function(M, L, symbol, geneset, s.test, tail.test, pw_sim, pw, correction, p_adj_th, norm,
                 analysis, outdir){

  print.logo()
  analysis=create.output.folder(analysis, outdir)
  print("preparing dataset ...")
  expr = get_expression_dataset(M, L, symbol)


  # 02. MIXTURE MODEL ====
  print("Mixture modelling ...")
  
  expr <- get_mixture_expr_class(expr
                                 , normalization = norm
                                 , ne_value=0.01
                                 , nc = 3
  )
  
  # # 04. SEA ====
  print("Running SEA ...")

  # sea = get_summary_expr_class_parallel(pl=geneset, r=expr
  #                                         , method = s.test
  #                                         , PW_SIM = pw_sim
  #                                         , TAIL = tail.test
  #                                         , POWER = pw
  # )
  
  gcr  = gene_class_representation(expr$expr_class
                                   , pl = geneset
                                   , method = s.test
                                   , TAIL = tail.test
                                   , PW_SIM = pw_sim
  )
  
  print("DONE")

  # 05. GSECA ====
  print("Running GSECA ...")
  # gseca = global_summary_expr_class(sea, geneset
  #                                   , method = s.test
  #                                   , PW_SIM = pw_sim
  #                                   , TAIL = tail.test
  #                                   , POWER = pw
  # )
  # browser()
  # load("../Results/Job18.2017-06-21.11_45_08.GSECA/results.Rdata")

  gseca = GSECA(  pl=geneset
                  , gcr
                  , method = s.test # fisher, chisq, wilcoxon
                  , PW_SIM = pw_sim
                  , TAIL = tail.test
                  , POWER = pw
                  , correction = correction)
  
  
  # gseca_plot = heatmap_GSECA_shiny(gseca
  #                                  , p_adj=p_adj_th
  #                                  , pow=pw
  #                                  )
  
  # PSUMLOG = 5e-5
  PSUMLOG = 0.01
  # PADJ    = 0.1
  # PW      = 0
  hm = heatmap_GSECA_sumlog(gseca
                            , filename=paste0(analysis,"/heatmap.pdf")
                            , p_adj = p_adj_th
                            , pow = pw
                            , psumlog=PSUMLOG)

  print("Saving files...")
  save(expr, gcr, gseca, file=paste0(analysis,"/results.Rdata") )
  # pdf(file=paste0(analysis,"/heatmap.pdf"),width = unit(8.27, "inches"), height = unit(11.69, "inches"))
  # print(gseca_plot$heat)
  # print(gseca_plot$barscore)
  # dev.off()
  write.csv(gseca, file=paste0(analysis,"/results.csv"), row.names = F)

  res=list('gseca'=gseca,
           # 'heat'=gseca_plot$heat,
           # 'score'=gseca_plot$barscore,
           'heat'= hm,
           'analysis'=analysis)
 return(res)
}
