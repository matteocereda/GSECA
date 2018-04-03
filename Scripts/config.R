# REQUIREMENTS ==========

options(stringsAsFactors = F)

pkgs = c( 'ggrepel'
         ,'grid'
         ,'gridExtra'
         ,'statmod'
         ,'R.utils'
         ,'reshape2'
         ,'mixtools'
         # ,'pwr' # ???????
         ,'scales'
         ,'gtable'
         ,'plyr'
         ,'parallel'
         ,'doParallel'
         ,'foreach'
         ,'snowfall'
         ,'rlecuyer'
         ,'RColorBrewer'
         ,'metap'
         ,'ComplexHeatmap'
         , 'circlize'
         )

not_installed = which( ! pkgs %in% rownames(installed.packages()) )
if(length(not_installed)>0){
  for(i in not_installed) install.packages(pkgs[i])
}

for(i in pkgs) suppressWarnings(suppressPackageStartupMessages(library(i, character.only =T)))

# ROUTINES ===========

print.logo = function() {
cat("██████╗ ███████╗███████╗ ██████╗ █████╗\n")
cat("██╔════╝ ██╔════╝██╔════╝██╔════╝██╔══██╗\n")
cat("██║  ███╗███████╗█████╗  ██║     ███████║\n")
cat("██║   ██║╚════██║██╔══╝  ██║     ██╔══██║\n")
cat("╚██████╔╝███████║███████╗╚██████╗██║  ██║\n")
cat("╚═════╝ ╚══════╝╚══════╝ ╚═════╝╚═╝  ╚═╝\n")
}


lapply_pb = function(X, FUN, ...){
  env <- environment()
  pb_Total <- length(X)
  counter <- 0
  pb <- txtProgressBar(min = 0, max = pb_Total, style = 3)

  # wrapper around FUN
  wrapper <- function(...){
    curVal <- get("counter", envir = env)
    assign("counter", curVal +1 ,envir=env)
    setTxtProgressBar(get("pb", envir=env), curVal +1)
    FUN(...)
  }
  res <- lapply(X, wrapper, ...)
  close(pb)
  res
}

mapply_pb = function(FUN, X, Y,  ...){
  env <- environment()
  pb_Total <- length(X)
  counter <- 0
  pb <- txtProgressBar(min = 0, max = pb_Total, style = 3)

  # wrapper around FUN
  wrapper <- function(...){
    curVal <- get("counter", envir = env)
    assign("counter", curVal +1 ,envir=env)
    setTxtProgressBar(get("pb", envir=env), curVal +1)
    FUN(...)
  }
  res <- mapply(wrapper, X, Y, ...)
  close(pb)
  res
}

ids_from_list_names = function(x,y){

  x = mapply(function(a,b){
    a$id = b
    return(a)
  },x,y, SIMPLIFY = F)

  return(x)
}

fix_names = function(p) {
  library(R.utils)
  return(capitalize(gsub(pattern = "_", replacement = " ",fixed = T,x = tolower(gsub(pattern = "KEGG_",replacement = "",fixed = T,x = as.character(as.character(p)))))))
}

read.gmt.file = function(pathMsigDbFile) {
  inputFile <- pathMsigDbFile
  con  <- file(inputFile, open = "r")

  c = 1
  pathway.list <- vector(mode="list",length=0)
  while (length(oneLine <- readLines(con, n = 1, warn = FALSE)) > 0)
  {
    myVector <- do.call("rbind",strsplit(oneLine, "\t"))
    t = vector(mode="list",length=1)
    t[[1]] = myVector[3:length(myVector)]
    names(t) = myVector[1]
    pathway.list = c(pathway.list,t)
    c = c+1
  }

  close(con)
  return(pathway.list)
}

create.output.folder = function (s, odir) {
  if(!dir.exists(odir)) dir.create(odir)

  job     = paste0('Job',length(grep(".Job", list.dirs(path = odir, recursive =F)))+1)
  timer    = format(Sys.time(), "%Y-%m-%d.%H_%M_%S")

  analisys = paste(job,timer,s, sep=".")
  analisys = paste0(odir,"/", analisys )

  dir.create( analisys )
  # for(i in s$sample ) dir.create( paste0( analisys, "/", i))
  return(analisys)
}

# STATS ===========

get_padj = function(x,m){
  x$p.adj=p.adjust(x$pv, method=m) #n = length(x$gene_set)
  x
}

summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE, conf.interval=.95, .drop=TRUE) {
  library(plyr)
  # New version of length which can handle NA's: if na.rm==T, don't count them
  length2 <- function (x, na.rm=FALSE) {
    if (na.rm) sum(!is.na(x))
    else       length(x)
  }
  # This does the summary. For each group's data frame, return a vector with
  # N, mean, and sd
  datac <- ddply(data, groupvars, .drop=.drop,
                 .fun = function(xx, col) {
                   c(N    = length2(xx[[col]], na.rm=na.rm),
                     mean   = mean   (xx[[col]], na.rm=na.rm),
                     median = median   (xx[[col]], na.rm=na.rm),
                     sd   = sd     (xx[[col]], na.rm=na.rm)
                   )
                 },
                 measurevar
  )
  # Rename the "mean" column
  # datac <- rename(datac, c("mean" = measurevar))
  datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean
  # Confidence interval multiplier for standard error
  # Calculate t-statistic for confidence interval:
  # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
  ciMult <- qt(conf.interval/2 + .5, datac$N-1)
  datac$ci <- datac$se * ciMult
  return(datac)
}


# bootTest <- rep(FALSE, 1000) #store results of thousand simulations
# n <- 12000  # guess sample size - we could iterate/search across values, I just changed until I got roughly 80% power
# p1 <- 0.1   # frequency of mutations
# p2 <- 0.01  # frequency of 'associated' mutations
# effectSize <- 2 #hypothesised
# for (i in seq(along=bootTest)) { #run many simulations
#   r1 <- runif(n)<p1 # randomly generate TRUEs with correct frequency
#   r2 <- runif(n)<ifelse(r1, p2*effectSize, p2) # create correlated TRUEs
#   bootTest[i] <- chisq.test(table(r1, r2))$p.value < 0.05 # store whether true-positive is found
# }
# mean(bootTest) # will be the proportion of true positives in total positives = power

power.chisq.test=function (p1, p2, n1, n2, alpha = 0.05, nsim = 100){
  y1 <- rbinom(nsim, size = n1, prob = p1)
  y2 <- rbinom(nsim, size = n2, prob = p2)
  y <- cbind(y1, n1 - y1, y2, n2 - y2)
  # p.value <- rep(0, nsim)
  p.value <- apply(y, 1, function(x)  chisq.test(matrix(x, 2, 2))$p.value )
  # for (i in 1:nsim) p.value[i] <- chisq.test(matrix(y[i, ], 2, 2))$p.value
  mean(p.value < alpha)
}

power.fisher.test=function (p1, p2, n1, n2, alpha = 0.05, nsim = 100, alternative='two.sided'){
  y1 <- rbinom(nsim, size = n1, prob = p1)
  y2 <- rbinom(nsim, size = n2, prob = p2)
  y <- cbind(y1, n1 - y1, y2, n2 - y2)
  # p.value <- rep(0, nsim)
  p.value <- apply(y, 1, function(x)  fisher.test(matrix(x, 2, 2), alternative=alternative)$p.value )
  # for (i in 1:nsim) p.value[i] <- chisq.test(matrix(y[i, ], 2, 2))$p.value
  mean(p.value < alpha)
}

FISHER = function (u, TAIL= TAIL, PW_SIM = PW_SIM){
  u = unlist(u)

  a = as.numeric(u["CASE.N.class"]); b = as.numeric(u['CASE.N.rest'])
  c = as.numeric(u["CNTR.N.class"]); d = as.numeric(u['CNTR.N.rest'])

  if( (b==0 & c==0) | (a==0 & d==0) ){ # ??? --- controllare con google -----
    return(  c(
      'pv' = 0,
      'or' = Inf,
      'pw' = 1
    ))
  }

  m = matrix( as.numeric(u[c("CASE.N.class","CASE.N.rest","CNTR.N.class", "CNTR.N.rest")]), nc=2, byrow = T)
  f = fisher.test(m, alternative = TAIL )  #, alt="greater"

  pw = power.fisher.test( a/(a+b), c/(c+d), (a+b),(c+d), alpha=.05, nsim= PW_SIM, alternative=TAIL)

  c(
    'pv' = f$p,
    'or' = f$estimate,
    'pw' = pw
  )
}

CHISQ = function (u){
  u = unlist(u)
  a = as.numeric(u["CASE.N.class"])
  b = as.numeric(u['CASE.N.rest'])
  c = as.numeric(u["CNTR.N.class"])
  d = as.numeric(u['CNTR.N.rest'])
  if( (b==0 & c==0) | (a==0 & d==0) ){ # ??? --- controllare con google -----
    return(  c(
      'pv' = 0,
      'or' = Inf,
      'pw' = 1
    ))
  }

  m = matrix( as.numeric(u[c("CASE.N.class","CASE.N.rest","CNTR.N.class", "CNTR.N.rest")]), nc=2, byrow = T)
  f = suppressWarnings(chisq.test( m ))

  or = (a*d)/(c*b)

  pw = power.chisq.test( a/(a+b), c/(c+d), (a+b),(c+d), alpha=.05, nsim= PW_SIM)
  
  c(
    'pv' = f$p.value,
    'or' = or,
    'pw' = pw
  )
}

# PREPARE INPUT ======

get_expression_dataset = function(M, L, gene_symbol="ensembl_gene_id"){
  if(is.null(M) | is.null(L)){
    stop(message('Please upload expression matrix and sample lables'))
  }else{

    if(length(L)!=(ncol(M)-1)) stop(message("Number of labels not equal to number of samples in the matrix"))
    require(reshape2)

    colnames(M)[1] = gene_symbol
    names(L)=colnames(M)[2:ncol(M)]
    expr= melt(M,id.vars = gene_symbol)
    colnames(expr) = c(gene_symbol, 'Tumor_Sample_Barcode', 'value')
    # expr[,1] = toupper(expr[,1])

    load("annotation/gencode.v26.annotation.Rdata")

    # define biomart object
    # require(biomaRt)
    # mart <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl", host="www.ensembl.org")
    # # query biomart
    # gene_symb <- getBM(attributes = c("ensembl_gene_id"
    #                                   #,'hgnc_symbol'
    #                                   ,'external_gene_name'
    #                                   ,'entrezgene'
    #                                   ,'gene_biotype')
    #                    # , filters = 'ensembl_gene_id'
    #                    # , values = expr$ensembl_gene_id
    #                    , mart = mart)

    if(gene_symbol=="ensembl_gene_id"){
      expr$symbol    = genes[match(expr$ensembl_gene_id, genes$gene_id), c('gene_name')]
      expr$gene_type = genes[match(expr$ensembl_gene_id, genes$gene_id), c('gene_type')]
      # expr = cbind(expr, gene_symb[match(expr$ensembl_gene_id, gene_symb$ensembl_gene_id), which(colnames(gene_symb)!="ensembl_gene_id")])
    }else if(gene_symbol=='symbol'){
      # gene_symb <- getBM(attributes = c("ensembl_gene_id"
      #                                   #,'hgnc_symbol'
      #                                   ,'external_gene_name'
      #                                   ,'entrezgene'
      #                                   ,'gene_biotype')
      #                    , filters = 'external_gene_name'
      #                    , values = expr$symbol
      #                    , mart = mart)
      expr$ensembl_gene_id = genes[match(expr$symbol, genes$gene_name), c('gene_id')]
      expr$gene_type       = genes[match(expr$symbol, genes$gene_name), c('gene_type')]
      # expr = cbind(expr, gene_symb[match(expr$symbol, gene_symb$external_gene_name), which(colnames(gene_symb)!="external_gene_name")])
    }

    expr$type=L[expr$Tumor_Sample_Barcode]
    
    expr$ensembl_gene_id  = as.character(expr$ensembl_gene_id)
    expr$symbol           = as.character(expr$symbol)
    
    return(expr)
  }
}


# MIXTURE MODEL =============

GMM = function(x
               , mean_lower_boundary =(-5) # lower boundary for 1st mean
               , mean_upper_boundary =2.8 # upper boundary for 2nd mean
               , my_perc  = c('20%','80%') # quantiles of expression distribution to start with
               , iterations = 50 # max number of iterations to reach correct results
               , nm_k = 2
               , nm_fast = TRUE
               , nm_maxit = 10000
               , nm_epsilon =  0.001
               ){

  na = subset(x, is.na(log_value) )
  x  = subset(x, !is.na(log_value) )

  perc = quantile(x$log_value, seq(0,1,.05))[ my_perc ]

  #initialization of internal variables
  normMix     = NA
  infloop     = 0
  threshold_1 = 1
  threshold_2 = 0
  mean_1      = (-10)
  mean_2      = 10

  while( (threshold_1>threshold_2) | (threshold_1>0) | ( mean_1<mean_lower_boundary ) | (mean_2>mean_upper_boundary) ){

    normMix = normalmixEM(  x$log_value
                          , mu      = perc
                          , k       = nm_k
                          , fast    = nm_fast
                          , maxit   = nm_maxit
                          , epsilon = nm_epsilon
                          )

    mean_1 = min(normMix$mu)
    mean_2 = max(normMix$mu)

    threshold_1 = mean_1+normMix$sigma[which.min(normMix$mu)]
    # ----- ? why did we comment this ? ------
    threshold_2 = mean_2 #-normMix$sigma[which.max(normMix$mu)]

    infloop=infloop+1 # to avoid infinite loops

    cat( 't1:\t',threshold_1 ,'t2:\t',threshold_2
        ,'m1:\t' ,mean_1 ,'m2:\t',mean_2
        ,'infloop:\t',infloop,'\t')

    if(infloop > iterations){
      break
    }
  }

  # na$comp.1 = NA  # na$comp.2 = NA  # x$comp.1 =normMix$posterior[,1]  # x$comp.2 =normMix$posterior[,2]

  x = as.data.frame(rbind(x,na))
  x = unrowname(x)

  # x$perc1   = perc[my_perc[1]]
  # x$perc2   = perc[my_perc[2]]
  x$mean1   = min(normMix$mu)
  x$mean2   = max(normMix$mu)
  x$st1     = normMix$sigma[which.min(normMix$mu)]
  x$st2     = normMix$sigma[which.max(normMix$mu)]
  x$lambda1 = normMix$lambda[which.min(normMix$mu)]
  x$lambda2 = normMix$lambda[which.max(normMix$mu)]

  qu1       = qnorm(c(0.25,0.75), x$mean1, x$st1)
  qu2       = qnorm(c(0.25,0.75), x$mean2, x$st2)
  x$L_th1   = qu1[1]
  x$L_th2   = qu1[2]
  x$H_th1   = qu2[1]
  x$H_th2   = qu2[2]

  return(x)
}

GMM_posterior = function(x
               , mean_lower_boundary =(-5) # lower boundary for 1st mean
               , mean_upper_boundary =2.8 # upper boundary for 2nd mean
               , my_perc  = c('20%','80%') # quantiles of expression distribution to start with
               , iterations = 50 # max number of iterations to reach correct results
               , nm_k = 2
               , nm_fast = TRUE
               , nm_maxit = 10000
               , nm_epsilon =  0.001
){

  na = subset(x, is.na(log_value) | !is.finite(log_value) )
  x  = subset(x, !is.na(log_value) & is.finite(log_value) )

  perc = quantile(x$log_value, seq(0,1,.05))[ my_perc ]

  #initialization of internal variables
  normMix     = NA
  infloop     = 0
  threshold_1 = 1
  threshold_2 = 0
  mean_1      = (-10)
  mean_2      = 10
  
  while( (threshold_1>threshold_2) | (threshold_1>0) |
         ( mean_1<mean_lower_boundary ) | (mean_2>mean_upper_boundary) ){

    normMix = normalmixEM(  x$log_value
                            , mu      = perc
                            , k       = nm_k
                            , fast    = nm_fast
                            , maxit   = nm_maxit
                            , epsilon = nm_epsilon
    )

    mean_1 = min(normMix$mu)
    mean_2 = max(normMix$mu)

    threshold_1 = mean_1+normMix$sigma[which.min(normMix$mu)]
    # ----- ? why did we comment this ? ------
    threshold_2 = mean_2 #-normMix$sigma[which.max(normMix$mu)]

    infloop=infloop+1 # to avoid infinite loops

    cat( 't1:\t',threshold_1 ,'t2:\t',threshold_2
         ,'m1:\t' ,mean_1 ,'m2:\t',mean_2
         ,'infloop:\t',infloop,'\t')

    if(infloop > iterations){
      break
    }
  }

  na$comp.1 = NA
  na$comp.2 = NA
  x$comp.1  = normMix$posterior[,1]
  x$comp.2  = normMix$posterior[,2]

  x = as.data.frame(rbind(x,na))
  x = unrowname(x)

  x$mean1   = min(normMix$mu)
  x$mean2   = max(normMix$mu)
  x$st1     = normMix$sigma[which.min(normMix$mu)]
  x$st2     = normMix$sigma[which.max(normMix$mu)]
  x$lambda1 = normMix$lambda[which.min(normMix$mu)]
  x$lambda2 = normMix$lambda[which.max(normMix$mu)]

  return(x)
}

set_expression_class = function(x
                                ,lev
                                ,thresholds=c('L_th1', 'L_th2','H_th1', 'H_th2')
                                ,ne_value=0.01
){

  th = unique(x[,thresholds])

  if(th[1]<log2(ne_value)) th[1] = ne_value

  classes=cut(x$log_value, breaks = c(-Inf, th, Inf), labels = (lev))

  if(sum(is.na(classes))>0) classes[is.na(classes)]='NE'
  x$expr_class = classes
  return(x)
}

set_expression_class_on_posterior = function(x
                                , posterior_cutoff=.9
){

  x$expr_class = "ME"

  i.na = which(is.na(x$log_value))
  if(length(i.na)>0) x$expr_class [ i.na ] = "NE"

  x$expr_class [ which(x$comp.1>=posterior_cutoff) ] = 'LE'
  x$expr_class [ which(x$comp.2>=posterior_cutoff) ] = 'HE'

  qt = quantile(subset(x, expr_class=='HE')$log_value)

  x$expr_class [ which(x$expr_class=='HE' & x$log_value<=qt["25%"] ) ] = 'HE1'
  x$expr_class [ which(x$expr_class=='HE' & x$log_value>=qt["25%"] & x$log_value<=qt["50%"]) ] = 'HE2'
  x$expr_class [ which(x$expr_class=='HE' & x$log_value>=qt["50%"] & x$log_value<=qt["75%"]) ] = 'HE3'
  x$expr_class [ which(x$expr_class=='HE' & x$log_value>=qt["75%"] ) ] = 'HE4'

  x$expr_class = factor(x$expr_class, levels = c('NE','LE','ME','HE1','HE2','HE3','HE4'))
  return(x)
}

get_mixture_expr_class=function(rna
                                      , normalization = "FPKM" # FPKM. RSEM, TPM, RPKM
                                      , ne_value=0.01
                                      , ncor=2
){

  message("FMM Normalization:\t", normalization, "\tNE threshold:\t", ne_value, "\tParallelizing ",ncor," cores ...")

  rna = rna[,c("Tumor_Sample_Barcode", "type","ensembl_gene_id","symbol","gene_type","value")]

  rna$log_value=NA
  rna$log_value[which(rna$value>ne_value)] = log2(rna$value[which(rna$value>ne_value)])

  message("FMM Getting thresholds ...")
  custnorm=c("FPKM",'RSEM','TPM','RPKM')
  S1=-4;S2=4
  if (normalization==custnorm[1]){
    if(ne_value==0.1){
      S1=-4;S2=4
    }else if(ne_value==0.01){
      S1=-5;S2=2.8
    }
  } else if(normalization==custnorm[2]){
    S1=0;S2=8
  }

  l = dlply(rna, ~Tumor_Sample_Barcode)

  if(.Platform$OS.type=="unix"){
    # on unix-like architectures
    l = mclapply(l, function(x, y, z){
                   tryCatch( GMM_posterior(x, mean_lower_boundary=y, mean_upper_boundary=z)
                             , error = function(x) return(NA))
                 },y=S1, z=S2, mc.silent = T, mc.cores = ncor)
  } else {
    # on other OS, e.g. Windows
    sfInit(parallel=TRUE, cpus=ncor, type="SOCK")
    sfLibrary(mixtools)
    sfLibrary(plyr)
    sfExport(list = c('S1','S2','GMM_posterior'))
    sfClusterSetupRNG()
    l = sfLapply(l, function(x, y, z){
                   tryCatch( GMM_posterior(x, mean_lower_boundary=y, mean_upper_boundary=z)
                             , error = function(x) return(NA))
                 },y=S1, z=S2)
    sfStop()
  }
  stats = lapply(l, function(x) x[1,c("Tumor_Sample_Barcode",'type','mean1','mean2','st1','st2','lambda1','lambda2')])
  stats = do.call(rbind,stats)
  stats = unrowname(stats)

  message("FMM Setting expression classes ...")
  l = lapply(l, function(x) x[,!colnames(x)%in%c('mean1','mean2','st1','st2','lambda1','lambda2')])
  l = lapply_pb(l, set_expression_class_on_posterior, posterior_cutoff=.9)

  rna = do.call("rbind", l)
  rna = unrowname(rna)

  return(list( 'expr_class'=rna
              ,'mixture_model'=stats))
}


# GENE CLASS REPRESENTATION ===========

# Representation of genes per expression classs

gene_class_representation=function(  r
                                     , pl
                                     , method = 'fisher'
                                     , PW_SIM = 100
                                     , TAIL = "two.sided"
                                     , POWER = 0.9
                                     , correction = 'fdr'
                                     , cutoff = 0.01
){
  require(plyr)
  require(metap)
  message("GRC # gene sets:\t", length(pl),"\n")

  total_genes = unique(unlist(pl))

  cn = c("ensembl_gene_id","symbol","type", "expr_class")

  samples = unique(r[,c("Tumor_Sample_Barcode","type")])
  samples = table(samples$type)

  rna = subset(r, symbol%in%total_genes )[,cn]
  rna = unrowname(rna)

  t = with(rna, table(ensembl_gene_id, expr_class,type ))
  
  case = as.data.frame(t[,,"CASE"]);
  case$N.not_in_class = samples['CASE'] - case$Freq
  colnames(case)[3:4] = c("CASE.N.class","CASE.N.rest")

  cntr = as.data.frame(t[,,"CNTR"]);
  cntr$N.not_in_class = samples['CNTR'] - cntr$Freq
  colnames(cntr)[3:4] = c("CNTR.N.class","CNTR.N.rest")

  ix = paste0(case$ensembl_gene_id,".",case$expr_class)
  iy = paste0(cntr$ensembl_gene_id,".",cntr$expr_class)

  cnts=as.data.frame(cbind(case, cntr[match(ix, iy),3:4]))

  rm(case, cntr)

  cnts$symbol = rna$symbol[ match(cnts$ensembl_gene_id, rna$ensembl_gene_id) ]
  cnts$CASE.P = with(cnts,CASE.N.class/(CASE.N.class+CASE.N.rest))
  cnts$CNTR.P = with(cnts,CNTR.N.class/(CNTR.N.class+CNTR.N.rest))
  cnts$delta  = with(cnts, CASE.P - CNTR.P  )
  cnts$direction = NA
  cnts$direction[which(cnts$delta<0)] = "D"
  cnts$direction[which(cnts$delta>0)] = "E"

  invariant = subset(cnts, is.na(direction))
  cnts      = subset(cnts, !is.na(direction))

  message(paste0("GCR Testing genes with ",method," ..."))

  if(method=="fisher"){
    stats = t(apply(cnts, 1, FISHER, TAIL=TAIL, PW_SIM=PW_SIM ))
  }else if(method=="chisq"){
    stats = t(apply(cnts, 1, CHISQ))
  }

  # controllare che i risultatiti non cambino ----
  # setto i pvalue == 0 pari al minimo valore misurato
  min_pv = min(stats[ which(stats[,1]>0) ,1 ], na.rm = T)

  if(min_pv>0.005) min_pv = 0.005

  stats[which(stats[,1]==0),1] = min_pv

  invariant$pv = 1
  invariant$or = 0
  invariant$pw = 0

  cnts$pv = stats[,1]
  cnts$or = stats[,2]
  cnts$pw = stats[,3]

  cnts = rbind(cnts, invariant)

  # forse la correzione bisogna farla non per classi ma per geneset ----
  message("GCR multiple testing correction ...")
  cnts2 = ddply( cnts[,c("symbol","pv")], .(symbol), get_padj, m=correction,  .progress = 'text')
  cnts$p.adj = cnts2$p.adj[match(cnts$symbol,cnts2$symbol)]
  rm(cnts2)
  

  # Fisher's method method --------
  message("GCR scoring ...")
  
  cnts2 = ddply(cnts[,c("symbol","pv")], .(symbol), mutate
               , sumlog   = sumlog(pv)$p
               # , .drop = F
               , .progress = 'text'
  )
  cnts$sumlog = cnts2$sumlog[match(cnts$symbol,cnts2$symbol)]
  rm(cnts2)
  
  # Ranking ----
  to.ord = unique(cnts$sumlog)
  ord    = 1:length(to.ord)
  n.ord  = sort(to.ord, decreasing = F)
  rank   = ord[ match( to.ord , n.ord ) ]
  tmp    = data.frame(to.ord, rank)
  cnts$rank = tmp$rank[match(cnts$sumlog,tmp$to.ord)]
  rm(tmp, ord, to.ord, n.ord, rank)  
  
  cnts = cnts[order(cnts$rank, cnts$expr_class),]

  cnts$sig = FALSE
  if( !is.na(POWER)){
    cnts$sig = cnts$p.adj<=cutoff & cnts$pw>=POWER
  }else{
    cnts$sig = cnts$p.adj<=cutoff
  }

  cnts = cnts[, c( "ensembl_gene_id","symbol", "expr_class",
                   'CASE.N.class','CASE.N.rest','CNTR.N.class','CNTR.N.rest',
                   'CASE.P','CNTR.P','delta',
                   'direction','or','pv',"p.adj",'pw', 'sig'
                   ,'rank'
                   ,"sumlog"
                   )]
  return(cnts)
}


# GSECA ==============

GSECA_rank = function(pv){
  to.ord = unique(pv)
  ord    = 1:length(to.ord)
  n.ord = sort(to.ord, decreasing = F)
  rank = ord[ match( pv , n.ord ) ]
  return(rank)
}


GSECA_runner = function(gene_set,
                        genes
                        , method='fisher'
                        , TAIL= 'two.sided'
                        , PW_SIM = 100
                 ) {
  require(plyr)

  genes = subset(genes, symbol%in%gene_set)

  px = ddply(genes, .(expr_class), summarise
             , CASE.N.class = sum(CASE.N.class)
             , CASE.N.rest  = sum(CASE.N.rest)
             , CNTR.N.class = sum(CNTR.N.class)
             , CNTR.N.rest  = sum(CNTR.N.rest)
             , .drop=F
  )

  px$CASE.P = with(px,CASE.N.class/(CASE.N.class+CASE.N.rest))
  px$CNTR.P = with(px,CNTR.N.class/(CNTR.N.class+CNTR.N.rest))

  px$delta  = with(px, CASE.P - CNTR.P  )

  # px$delta  = with(px, (CASE.P - CNTR.P)/(CASE.P + CNTR.P)  )*100
  
  px$direction = NA
  px$direction[which(px$delta<0)] = "D"
  px$direction[which(px$delta>0)] = "E"

  if(method=='fisher'){
    stats = t(apply(px, 1, FISHER, TAIL=TAIL, PW_SIM=PW_SIM ))
  }else if(method=='chisq'){
    stats = t(apply(cnts, 1, CHISQ))
  }

 # else if(method=='wilcox'){
  # }

  px$pv = stats[,1]
  px$or = stats[,2]
  px$pw = stats[,3]

  return(px)
}

GSECA_core <- function(gcr
                       , pl
                       , method='fisher'
                       , TAIL= 'two.sided'
                       , PW_SIM=100
                       , correction='fdr'
                       ) {
  
  px = lapply_pb(pl, GSECA_runner
                 , genes = gcr
                 , method=method
                 , TAIL= TAIL
                 , PW_SIM = PW_SIM
  )
  
  
  px = mapply(function(x,y){x$gene_set=y;return(x)}, px, names(px), SIMPLIFY = F)
  px = as.data.frame(do.call('rbind',px))
  px = unrowname(px)
  
  # --- correggo per gene_set ----
  px = ddply(px, .(gene_set), get_padj, m=correction)
  
  message("GSECA scoring ...")
  
  # Fisher's method ----
  
  px = ddply(px, .(gene_set), mutate
             , sumlog   = sumlog(pv)$p
  )
  px[,c(
    "gene_set","expr_class","CASE.N.class","CASE.N.rest","CNTR.N.class","CNTR.N.rest","CASE.P","CNTR.P",
    "delta","direction","pv","or","pw","p.adj"
    ,"sumlog"
  )]
}

GSECA = function(  pl
                  , gcr
                  , method ='fisher' # fisher, chisq, wilcoxon
                  , PW_SIM = 100
                  , TAIL = "two.sided"
                  , POWER = 0.9
                  , correction='fdr'
                  , cutoff = 0.01
                  # , remove_invariant_gene=F
                  ){
# 
#   if(remove_invariant_gene){
#     x = ddply(gcr, .(ensembl_gene_id), summarise, n=length(expr_class[is.na(direction)]))
#     invariant_in_all_classes = x[which(x$n==5),1]
#     gcr = subset(gcr, !ensembl_gene_id%in%invariant_in_all_classes )
#   }

  message("GSECA core ...")

  core = GSECA_core(gcr, pl, method=method, TAIL=TAIL, PW_SIM=PW_SIM, correction=correction)

  # Ranking ----
  core = ddply(core, .(gene_set), mutate
             , rank = GSECA_rank(sumlog)
             )

  core = core[order(core$rank, core$expr_class),]

  core$sig = FALSE
  if( !is.na(POWER)){
    core$sig = core$p.adj<=cutoff & core$pw>=POWER
  }else{
    core$sig = core$p.adj<=cutoff
  }


  core[,c(
   "gene_set","expr_class","CASE.N.class","CASE.N.rest","CNTR.N.class","CNTR.N.rest","CASE.P","CNTR.P",
   "delta","direction","pv","or","pw","p.adj", "sig"
   ,'rank'
   ,"sumlog"
   )]
}




GSECA_get_GS = function(x,y) subset(x, gene_set)

# EXECUTORs =========

# execute_GSECA <- function(x=NULL
#                           , M 
#                           , L
#                           , expr=NULL
#                           , pl
#                           , norm = "FPKM"
#                           , ncor = 2
#                           , s.test = 'fisher'
#                           , tail.test = 'two.sided'
#                           , correction = 'fdr' 
#                           , pw_sim = 1
#                           , pw = 0.9
#                           ){
#   # 0. Prepare Data
#   if(!is.null(expr)){
#     rna = expr
#   } else{
#     rna = get_expression_dataset(M, L, 'symbol')   
#   }
#   if(!is.null(x)) rna = subset(rna, Tumor_Sample_Barcode%in%x)
#   ## GSECA
#   # 1. Mixture Model
#   expr <- get_mixture_expr_class(rna
#                                  , normalization = norm
#                                  , ne_value=0.01
#                                  , nc = ncor
#   )
#   rm(rna)
#   # 2. GCR
#   gcr  = gene_class_representation(expr$expr_class
#                                    , pl
#                                    , method = s.test
#                                    , TAIL = tail.test
#                                    , PW_SIM = pw_sim
#   )
#   rm(expr)
#   # 3. GSECA
#   gseca = GSECA(  pl
#                   , gcr
#                   , method = s.test
#                   , PW_SIM = pw_sim
#                   , TAIL = tail.test
#                   , POWER = pw
#                   , correction = correction)
#   rm(gcr)
#   # 4. Resilts
#   return(gseca)
# }

GSECA_executor_shiny <- function( M
                                 , L
                                 , symbol
                                 , geneset
                                 , s.test = 'fisher'
                                 , tail.test = 'two.sided'
                                 , pw_sim = 1
                                 , pw = 0.9
                                 , correction = 'fdr'
                                 , p_adj_th = 0.1
                                 , norm = "FPKM"
                                 , analysis
                                 , outdir
                                 , N.CORES = 3
                                 , EMPIRICAL = F
                                 , BOOTSTRP = F
                                 , nsim = 10){
  
  print.logo()
  analysis=create.output.folder(analysis, outdir)
  # 01. Prepare Data 
  print("preparing dataset ...")
  expr = get_expression_dataset(M, L, symbol)
  
  # 02. MIXTURE MODEL 
  print("Mixture modelling ...")
  expr <- get_mixture_expr_class(expr
                                 , normalization = norm
                                 , ne_value=0.01
                                 , nc = N.CORES
  )
  
  # 03. GCR 
  print("Gene Class Representation ...")
  gcr  = gene_class_representation(expr$expr_class
                                   , pl = geneset
                                   , method = s.test
                                   , TAIL = tail.test
                                   , PW_SIM = pw_sim
  )
  
  # 04. GSECA 
  print("Running GSECA ...")
  gseca = GSECA(  pl=geneset
                  , gcr
                  , method = s.test # fisher, chisq, wilcoxon
                  , PW_SIM = pw_sim
                  , TAIL = tail.test
                  , POWER = pw
                  , correction = correction
                  )
  
  gseca$p.emp = NA
  gseca$success_rate = NA
  
  if(EMPIRICAL){
    # 05. EMPIRICAL PVALUES
    gseca.emp = GSECA.Bootstrap.empirical( expr$expr_class
                                           , gseca
                                           , pl = geneset
                                           , NSIM = nsim
                                           , nc = N.CORES
                                           , method = s.test
                                           , PW_SIM = pw_sim
                                           , TAIL = tail.test
                                           , POWER = pw
                                           , correction = correction
                                           , cutoff = 0.01
    )
    
    if(!is.null(gseca.emp)){
      gseca$p.emp = gseca.emp$p.emp[match(gseca$gene_set, gseca.emp$gene_set)] 
    }
  }
 
  if(BOOTSTRP){
    # 07. CORRECTION FOR SAMPLE SIZE 
    gseca.boot=GSECA.Bootstrap.sample_size(   expr$expr_class
                                              , gseca
                                              , pl = geneset
                                              , phen = c("CASE","CNTR")
                                              , NSIM = nsim
                                              , nc = N.CORES
                                              # , sig.threshold=0.01
                                              , method = s.test
                                              , PW_SIM = pw_sim
                                              , TAIL = tail.test
                                              , POWER = pw
                                              , correction = correction
                                              , cutoff = 0.01
    )
    
    if(!is.null(gseca.boot)){
      gseca$success_rate = gseca.boot$success_rate[match(gseca$gene_set, 
                                                         gseca.boot$gene_set)]
    }
  }
  
  PSUMLOG = 0.01
  PADJ    = 0.1
  PW      = 0
  PEMP    = 1
  SRATE   = 0.7
  # browser()
  ecmap = GSECA.ECmap( gseca
                       , filename=paste0(analysis,"/ECmap.psumlog_",PSUMLOG,
                                         "_padj_",PADJ,
                                         "_pw_",PW,
                                         "_pemp_",PEMP,
                                         "_srate_",SRATE,".pdf")
                       , p_adj = PADJ
                       , pow = PW
                       , psumlog=PSUMLOG
                       , pemp = PEMP
                       , srate=SRATE)

  # Results 
  print("Saving files...")
  gseca$gene_set <- fix_names(gseca$gene_set)
  save(expr, gcr, gseca, file=paste0(analysis,"/results.Rdata"))
  write.csv(gseca, file=paste0(analysis,"/results.csv"), row.names = F)
  
  res = list( 'gseca' = gseca, 
              'ECmap' = ecmap,
              'analysis' = analysis)
  
  return(res)
}

# SIMULATION ====

execute_GSECA_bootstrap = function(sample,  rna  , pl
                           , type='sample.size'
                           # GRC
                           , GRC.method = 'fisher'
                           , GRC.PW_SIM = 1
                           , GRC.TAIL = "two.sided"
                           , GRC.POWER = 0.9
                           , GRC.correction = 'fdr'
                           , GRC.cutoff = 0.01
                           
                           # GSECA
                           , GSECA.method ='fisher'
                           , GSECA.PW_SIM = 1
                           , GSECA.TAIL = "two.sided"
                           , GSECA.POWER = 0.9
                           , GSECA.correction='fdr'
                           , GSECA.cutoff = 0.01
){
  
  if(type=='p.empirical'){
    
    rna$type = sample$type[ match( as.character(rna$Tumor_Sample_Barcode), sample$barcode)]
  
  }else{
    
    rna = subset(rna, Tumor_Sample_Barcode%in%sample)
    
  }
  
  dd  = gene_class_representation(rna
                                  , pl
                                  , method =   GRC.method
                                  , TAIL =     GRC.TAIL
                                  , PW_SIM =   GRC.PW_SIM
                                  , correction=GRC.correction
                                  , POWER =    GRC.POWER
                                  , cutoff =   GRC.cutoff
  )
  
  res = GSECA(  pl
                , dd
                
                , method =    GSECA.method
                , PW_SIM =    GSECA.PW_SIM
                , TAIL =      GSECA.TAIL
                , POWER =     GSECA.POWER
                , correction= GSECA.correction
                , cutoff =    GSECA.cutoff 
  )
  res
}

GSECA.Bootstrap.empirical <- function( expr
                             , obs.gseca
                             , pl
                             , NSIM=1000
                             , nc=3
                             , method = 'fisher'
                             , PW_SIM = 1
                             , TAIL = "two.sided"
                             , POWER = 0.9
                             , correction = 'fdr'
                             , cutoff = 0.01
                             ){
  rna = expr
  total_genes = unique(unlist(pl))
  rna = subset(rna, symbol %in% total_genes )
  
  boot = unique(rna[,c('Tumor_Sample_Barcode','type')])
  
  samples = 1:NSIM
  samples = lapply(samples, 
                   function(x,y,z) data.frame(barcode=z,type=sample(y))
                   , y=as.character(boot$type), z=as.character(boot$Tumor_Sample_Barcode)  
                   )
  
  if(.Platform$OS.type=="unix"){
      # on unix-like architectures
      l = mclapply(samples, function(x, y, z, ...){
        tryCatch( 
          execute_GSECA_bootstrap(x, rna=y, pl=z, type='p.empirical'
                                  , GRC.method=method, GRC.PW_SIM=PW_SIM
                                  , GRC.TAIL=TAIL,GRC.POWER=POWER
                                  , GRC.correction=correction, GRC.cutoff=cutoff
                                  , GSECA.method=method, GSECA.PW_SIM=PW_SIM
                                  , GSECA.TAIL=TAIL, GSECA.POWER=POWER
                                  , GSECA.correction=correction, GSECA.cutoff=cutoff
                                  )
                  , error = function(x) return(NA)
          )
      }
      , y=rna, z=pl
      , GRC.method=method, GRC.PW_SIM=PW_SIM
      , GRC.TAIL=TAIL,GRC.POWER=POWER
      , GRC.correction=correction, GRC.cutoff=cutoff
      , GSECA.method=method, GSECA.PW_SIM=PW_SIM
      , GSECA.TAIL=TAIL, GSECA.POWER=POWER
      , GSECA.correction=correction, GSECA.cutoff=cutoff
      , mc.silent = T, mc.cores = nc)
      
    } else {
      # on other OS, e.g. Windows
      sfInit(parallel=TRUE, cpus=nc, type="SOCK")
      sfLibrary(mixtools)
      sfLibrary(plyr)
      sfExport(list = c('rna','pl','execute_GSECA_boot_sample_size'))
      sfClusterSetupRNG()
      l = sfLapply(l, function(x, y, z){
        tryCatch( 
          execute_GSECA_bootstrap(x, rna=y, pl=z, type='p.empirical'
                                  , GRC.method=method, GRC.PW_SIM=PW_SIM
                                  , GRC.TAIL=TAIL,GRC.POWER=POWER
                                  , GRC.correction=correction, GRC.cutoff=cutoff
                                  , GSECA.method=method, GSECA.PW_SIM=PW_SIM
                                  , GSECA.TAIL=TAIL, GSECA.POWER=POWER
                                  , GSECA.correction=correction, GSECA.cutoff=cutoff
                                  )
                  , error = function(x) return(NA)
          )
      }
      , y=rna, z=pl
      , GRC.method=method, GRC.PW_SIM=PW_SIM
      , GRC.TAIL=TAIL,GRC.POWER=POWER
      , GRC.correction=correction, GRC.cutoff=cutoff
      , GSECA.method=method, GSECA.PW_SIM=PW_SIM
      , GSECA.TAIL=TAIL, GSECA.POWER=POWER
      , GSECA.correction=correction, GSECA.cutoff=cutoff
      )
      sfStop()
    }
    
    l <- do.call(rbind.data.frame, lapply(l, function(x) unique(x[,c("gene_set", "sumlog")])))
    
    message("GSECA.Bootstrap Success rate ...")
    
    l$obs = obs.gseca$sumlog[match(l$gene_set,obs.gseca$gene_set)]
    
    l <- ddply(l, .(gene_set), summarise
               
               , p.emp = (1+sum(sumlog<obs, na.rm = T))/(1+length(na.omit(sumlog)))
    
               )

  return(l)
}


GSECA.Bootstrap.sample_size <- function( expr
                            , obs.gseca
                            , pl
                            , phen = c("CASE","CNTR")
                            , NSIM=1000
                            , nc=3
                            , cutoff=0.01
                            , method = 'fisher'
                            , PW_SIM = 1
                            , TAIL = "two.sided"
                            , POWER = 0.9
                            , correction = 'fdr'
                            ){
  rna = expr
  total_genes = unique(unlist(pl))
  rna = subset(rna, symbol %in% total_genes )
  
  message("GSECA.Bootstrap Create Phenotype lists ...")
  
  pheno = list()
  pheno[[1]] = as.character(unique(subset(rna, type==phen[1])$Tumor_Sample_Barcode))
  pheno[[2]] = as.character(unique(subset(rna, type==phen[2])$Tumor_Sample_Barcode))
  
  if(length(pheno[[1]])!=length(pheno[[2]])){
    
    ref.idx     = which.max(c(length(pheno[[1]]), length(pheno[[2]])))
    size.idx    = which.min(c(length(pheno[[1]]), length(pheno[[2]])))
    sample.size = length(pheno[[size.idx]])
    
    message("GSECA.Bootstrap Bootstrapping ...")

    registerDoParallel(cores=nc)
    r <- foreach(i=1:NSIM) %dopar% {
      sample_barcode = vector(mode = "numeric", length = 2*sample.size)
      sample_barcode[1:sample.size] = pheno[[size.idx]]
      sample_barcode[(sample.size+1):length(sample_barcode)] = sample(pheno[[ref.idx]]
                                                                      , length(pheno[[size.idx]]))
      return(sample_barcode)
    }

    
    if(.Platform$OS.type=="unix"){
      # on unix-like architectures
      l = mclapply(r, function(x, y, z, ...){
        tryCatch( 
          execute_GSECA_bootstrap(x, rna=y, pl=z, type='sample.size'
                                  , GRC.method=method, GRC.PW_SIM=PW_SIM
                                  , GRC.TAIL=TAIL, GRC.POWER=POWER
                                  , GRC.correction=correction, GRC.cutoff=cutoff
                                  , GSECA.method=method, GSECA.PW_SIM=PW_SIM
                                  , GSECA.TAIL=TAIL, GSECA.POWER=POWER
                                  , GSECA.correction=correction, GSECA.cutoff=cutoff
          )
          , error = function(x) return(NA))
      }
      , y=rna, z=pl
      , GRC.method=method, GRC.PW_SIM=PW_SIM
      , GRC.TAIL=TAIL,GRC.POWER=POWER
      , GRC.correction=correction, GRC.cutoff=cutoff
      , GSECA.method=method, GSECA.PW_SIM=PW_SIM
      , GSECA.TAIL=TAIL, GSECA.POWER=POWER
      , GSECA.correction=correction, GSECA.cutoff=cutoff
      , mc.silent = T, mc.cores = nc)

    } else {
      # on other OS, e.g. Windows
      sfInit(parallel=TRUE, cpus=nc, type="SOCK")
      sfLibrary(mixtools)
      sfLibrary(plyr)
      sfExport(list = c('rna','pl','execute_GSECA_boot_sample_size'))
      sfClusterSetupRNG()
      l = sfLapply(l, function(x, y, z, ...){
        tryCatch( 
          execute_GSECA_bootstrap(x, rna=y, pl=z, type='sample.size'
                                  , GRC.method=method, GRC.PW_SIM=PW_SIM
                                  , GRC.TAIL=TAIL, GRC.POWER=POWER
                                  , GRC.correction=correction, GRC.cutoff=cutoff
                                  , GSECA.method=method, GSECA.PW_SIM=PW_SIM
                                  , GSECA.TAIL=TAIL, GSECA.POWER=POWER
                                  , GSECA.correction=correction, GSECA.cutoff=cutoff
          )
          , error = function(x) return(NA))
      }
      , y=rna, z=pl
      , GRC.method=method, GRC.PW_SIM=PW_SIM
      , GRC.TAIL=TAIL, GRC.POWER=POWER
      , GRC.correction=correction, GRC.cutoff=cutoff
      , GSECA.method=method, GSECA.PW_SIM=PW_SIM
      , GSECA.TAIL=TAIL, GSECA.POWER=POWER
      , GSECA.correction=correction, GSECA.cutoff=cutoff
      )
      sfStop()
    }
    
    l <- do.call(rbind.data.frame, lapply(l, function(x) unique(x[,c("gene_set", "sumlog")])))
    
    message("GSECA.Bootstrap Success rate ...")
    
    l$obs = obs.gseca$sumlog[match(l$gene_set,obs.gseca$gene_set)]
    l$sig.threshold = cutoff
    l <- ddply(l, .(gene_set), summarise
               , p.emp        = (1+sum(sumlog<obs, na.rm = T))/(1+length(na.omit(sumlog)))
               , success_rate = sum(sumlog<sig.threshold, na.rm = T)/length(na.omit(sumlog))
    )
    
  } else{
    print("Same number of samples per phenotype - Skip Sample size correction")
    l=NULL
  }
  
  return(l)
}

# GSECA.Bootstrap.matrix <- function( M
#                              , L
#                              , gene_symbol_type
#                              , pl
#                              , NSIM
#                              , nc=3
#                              , sig.threshold=0.01
#                              , method = 'fisher'
#                              , PW_SIM = 1
#                              , TAIL = "two.sided"
#                              , POWER = 0.9
#                              , correction = 'fdr'
#                              , cutoff = 0.01
#                              
# ){
#   rna = get_expression_dataset(M, L, gene_symbol_type)
# 
#   phen = unique(L)[order(unique(L), decreasing = F)]
#   
#   cat("Create Phenotype lists ...\n")
#   pheno = list()
#   pheno[[1]] = as.character(unique(subset(rna, type==phen[1])$Tumor_Sample_Barcode))
#   pheno[[2]] = as.character(unique(subset(rna, type==phen[2])$Tumor_Sample_Barcode))
#   
#   if(length(pheno[[1]])!=length(pheno[[2]])){
#     
#     ref.idx = which.max(c(length(pheno[[1]]), length(pheno[[2]])))
#     size.idx = which.min(c(length(pheno[[1]]), length(pheno[[2]])))
#     sample.size = length(pheno[[size.idx]])
#     
#     registerDoParallel(cores=nc)
#     cat("Bootstrapping ...\n")
#     r <- foreach(i=1:NSIM) %dopar% {
#       sample_barcode = vector(mode = "numeric", length = 2*sample.size)
#       sample_barcode[1:sample.size] = pheno[[size.idx]]
#       sample_barcode[(sample.size+1):length(sample_barcode)] = sample(pheno[[ref.idx]], length(pheno[[size.idx]]))
#       return(sample_barcode)
#     }
#     
#     cat("Execute GSECA ...\n")
#     l <- mclapply(r, execute_GSECA, expr=rna, pl=pl, ncor=1, mc.cores = nc, mc.silent = T)
#     l <- do.call(rbind.data.frame, lapply(l, function(x) unique(x[,c("gene_set", "sumlog")])))
#     
#     cat("Success rate ...\n")
#     l <- ddply(l, .(gene_set), summarise
#                , success_rate = sum(sumlog<sig.threshold, na.rm = T)/length(na.omit(sumlog))
#     )
#     
#   } else{
#     print("Same number of samples per phenotype - Skip Sample size correction")
#     l=NULL
#   }
#   
#   return(l)
# }

# PLOTTING ============

g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  legend
}

# Expression profile in samples


# Expression profile in samples with
plot_mix_comps <- function(x, mu, sigma, lam) {
  lam * dnorm(x, mu, sigma)
}

plot_mix_density <- function(x, mu1, mu2, sigma1, sigma2, lam1, lam2) {
  lam1 * dnorm(x, mu1, sigma1) + lam2 * dnorm(x, mu2, sigma2)
}

# Plot expression
plot_expression_SE = function(dexpr, pl
                           ,gene_set
                           ,colname='symbol'
                           ,measure="log_value"
                           ,x1=-10
                           ,x2=12
                           ,protein_coding=F
){

  m1=summarySE(dexpr, groupvars = 'type', measurevar = c('mean1'))
  l1=summarySE(dexpr, groupvars = 'type', measurevar = c('lambda1'))
  m2=summarySE(dexpr, groupvars = 'type', measurevar = c('mean2'))
  l2=summarySE(dexpr, groupvars = 'type', measurevar = c('lambda2'))
  s1=summarySE(dexpr, groupvars = 'type', measurevar = c('st1'))
  s2=summarySE(dexpr, groupvars = 'type', measurevar = c('st2'))


  df=rbind(cbind(m1,var='mean1'),cbind(l1,var='lambda1')
           ,cbind(m2,var='mean2'),cbind(l2,var='lambda2')
           , cbind(s1,var='st1'),cbind(s2,var='st2'))
  df=subset(df, type=='CASE')
  ggplot(dexpr,aes(x=log_value))+
    stat_function(geom = "line"
                    , fun = plot_mix_density,
                    args = list(mu1 = subset(df, var=='mean1'), mu2 = subset(df, var=='mean2')
                                , sigma1 =subset(df, var=='st1'), subset(df, var=='st2')
                                , lam1 = subset(df, var=='lamda1'), lam2 = subset(df, var=='lambda2')), color='red'
    )+  xlim(x1,x2)


}

# Plot expression curve for a sample
plot_expression = function(dexpr, pl, my_sample, gene_set
                           ,colname='symbol'
                           ,measure="log_value"
                           ,x1=-10
                           ,x2=12
                           ,protein_coding=F
){
  require(ggrepel)

  dexpr$Tumor_Sample_Barcode = as.character(dexpr$Tumor_Sample_Barcode)
  df = subset(dexpr, Tumor_Sample_Barcode%in%my_sample)

  if(protein_coding) df = subset(df, gene_biotype=='protein_coding')

  if(nrow(df)==0) stop('No sample found')
  df$selected = df[,colname]%in%pl[[gene_set]]

  p = ggplot(df, aes(x=log_value))+
      stat_function(  geom = "line"
                  , fun = plot_mix_density,
                    args = list(mu1 = df$mean1[1], mu2 = df$mean2[1]
                             , sigma1 = df$st1[1], sigma2 = df$st2[1]
                             , lam1 = df$lambda1[1], lam2 = df$lambda2[2]), color='red'
                 )    +
    xlim(x1,x2)


  gg  = ggplot_build(p)$data[[1]]
  gga = unique(ggplot_build(p)$plot[[1]]$type)
  ymax = max(gg$y)
  ymed = ymax/2
  lmax = ifelse( ymax>=0.5, round(ymax), 1)
  lmed = ifelse(lmax==1, 0.5, ifelse( ymed>0.5, round(ymed), 0.5))

  drvrs=subset(df, selected )
  ix=which(is.na(drvrs$log_value))
  if(length(ix)>0) drvrs$log_value[ix] = x1

  plvs=c(1,2); names(plvs)=c("CASE","CNTR")
  plvs = plvs[which(names(plvs)%in%gga)]

  p.drvrs=rbind()
  # for(j in unique(drvrs$type)){

  # dd1 = subset(gg,group == plvs[j])
  dd1=gg
  # y1  =  subset(drvrs,(type)==j )
  y1  =  drvrs
  yy1=xxx1=c()
  for (i in y1[,measure]){
    xxx1 = c(xxx1, dd1$x[which(abs(dd1$x-i)==min(abs(dd1$x-i)))][1])
    yy1  = c(yy1,  dd1$y[which(abs(dd1$x-i)==min(abs(dd1$x-i)))][1])
  }
  if(!is.null(xxx1)){
    y1$cell=xxx1
    y1$y=yy1
    p.drvrs=rbind(p.drvrs,y1)
  }
  # }

  p.drvrs      = as.data.frame(p.drvrs, stringsAsFactors=F)

  p=
    p +
    geom_vline(data=df[1,], aes(xintercept=L_th1), linetype="dashed", color="blue")+
    geom_vline(data=df[1,], aes(xintercept=L_th2), linetype="dashed", color="blue")+
    geom_vline(data=df[1,], aes(xintercept=H_th1), linetype="dashed", color="blue")+
    geom_vline(data=df[1,], aes(xintercept=H_th2), linetype="dashed", color="blue")+
    geom_point(data=p.drvrs, aes(x = cell, y=y), colour="black", show.legend = F)+
    geom_text_repel(data=p.drvrs, aes(x=cell, y=y, label = symbol), colour="black")+
    theme_bw()+theme(legend.position = 'none')+ggtitle(my_sample)
  #        (data=p.drivers, aes(x = cell,y=y,label=gname, family=""),  size=rel(3.5), hjust=rel(.75), vjust=rel(-.25), show.legend = F)
  t=table(subset(df, selected)$expr_class)
  t=as.data.frame(t)
  colnames(t)=c('Expr.Class','Genes')

  p= p+annotation_custom(tableGrob(t,rows = NULL
                                # , theme = ttheme_minimal()
                                ), xmin=x2-3, xmax=x2-1, ymax = ymax/3, ymin=0)
   return(p)
}




plot_mixtures <- function(exval,mix,
                          x_lim =c(-15,15)
) {

  require(ggplot2)
  require(plyr)

  r = ddply(exval, .(expr_class), summarise
            , down = min(log_value, na.rm = T)
            , up = max(log_value, na.rm = T))
  rownames(r)=r[,1]
  rx = unlist(c(r["LE",2:3], r["HE1",2:3], r['HE3',2:3]))

  histo_pal = "#56B4E9"
  comp1_pal = "blue"
  comp2_pal = "red"

  font_family = "sans"
  i_lwd=1
  toplot <-
    ggplot(subset(exval, !is.na(log_value)), aes(x=log_value)) + xlim(x_lim) +
    geom_histogram(aes(y= ..density..), binwidth = 0.5, colour = "black",fill = NA, alpha=0.4) +
    stat_function(geom = "line", fun = plot_mix_comps,
                  args = list(mu = mix$mean1, sigma = mix$st1,lam = mix$lambda1),
                  lwd = i_lwd, linetype = 1, aes(colour = "Less Expressed")) +
    stat_function(geom = "line", fun = plot_mix_comps,
                  args = list(mu = mix$mean2, sigma = mix$st2,lam = mix$lambda2),
                  lwd = i_lwd, linetype = 1, aes(colour = "More Expressed")) +
    # stat_function(geom = "line", fun = plot_mix_density,
    #               args = list(mu1 = mean1, mu2 = mean2, sigma1 = std1, sigma2 = std2,
    #                           lam1 = lambda1, lam2 = lambda2),
    #               lwd = 0.8, alpha=0.8, linetype = 2, aes(colour = "Mixture Density")) +
    scale_colour_manual("Curves", breaks=c("Less Expressed","More Expressed","Mixture Density"),
                        values=c("Less Expressed"=comp1_pal,"More Expressed"=comp2_pal,"Mixture Density"="black"))+
    theme_classic() +

    ylab("Probability Density") +
    xlab("Expression Values") +
    geom_vline(data=data.frame(x = rx), aes(xintercept=x), col='grey')
  toplot
}

ytickform <- function(x){
  lab <- sprintf("%05s",x)
}


heatmap_GCR = function( gcr
                          , filename=NULL
                          # , p_adj=0.1
                          # , pow=0.9
                          , sc = NULL
){

  require(RColorBrewer)
  require(scales)
  require(gtable)



  toplot = subset(gcr,
                  sig
                    # & symbol%in%pl[['KEGG_CYTOKINE_CYTOKINE_RECEPTOR_INTERACTION']]
                  )

  if(nrow(toplot)==0| is.null(toplot)) stop(message("No enriched genes"))

  toplot$p.adj[which(!toplot$sig)]=NA


  myPalette <-colorRampPalette(brewer.pal(9, "RdYlBu"))

  gg_hm =ggplot(toplot, aes(x=expr_class,y=symbol, fill=p.adj)) +
    geom_tile(fill = "white",colour = "black") + #
    theme_bw() +
    coord_equal() +
    theme(
      legend.position = 'left'
      , panel.grid.major = element_blank()
      , panel.grid.minor = element_blank()
      , panel.border = element_blank()
      , panel.background = element_blank()
      , axis.text.x = element_text(angle=90, vjust=0.5,hjust=1)) + xlab("") +ylab("") +
    geom_point(data=subset(toplot, delta>=0), aes(size=delta),shape=24,colour="black",na.rm=TRUE,show.legend = T) +
    geom_point(data=subset(toplot, delta<0), aes(size=abs(delta)),shape=25,colour="black",na.rm=TRUE,show.legend = T) +
    scale_fill_gradientn(colours = myPalette(5),na.value="white")+
    scale_size_continuous(range = c(1,5)
                          #,limits=c(0,6),breaks=c(0,1,2), labels=c("0","1","2")
    )


  pdf(file=filename,width = unit(12.27, "inches"), height = unit(15.69, "inches"))
  grid.draw(g)
  dev.off()


 }

anno_bimodal = function() {
  n = 7
  
  set.seed(30580)
  xseq<-seq(0, 7 ,.01)
  densities<-dnorm(xseq, 2,1)
  densities2<-dnorm(xseq, 5,1)
  
  
  pushViewport(viewport( xscale = c(0,7), 
                         yscale = c(0,1)
  ))
  
  # grid.rect(gp = gpar(fill = "transparent"))
  # grid.lines(xseq, densities
  #            ,  gp = gpar(col = 'blue')
  #            , default.units = "native") 
  grid.polygon(c(xseq,rev(xseq))
               , c(densities, rep(0,length(xseq)))
               , gp = gpar(fill = 'blue', col = 'blue', alpha=0.5)
               , default.units = "native"
  )
  grid.polygon(c(xseq,rev(xseq))
               , c(densities2*2.25, rep(0,length(xseq)))
               , gp = gpar(fill = 'red', col = 'red', alpha=0.5)
               , default.units = "native"
  )
  
  # grid.yaxis(gp = gpar(fontsize = 8))
  upViewport()
}



GSECA.ECmap = function( gseca
                        , filename=NULL
                        , p_adj=0.1
                        , pow=0.9
                        , psumlog=1
                        , pemp=NULL
                        , srate=NULL
                                
){

  require(RColorBrewer)

  message("Thresholds: \nP-value adj ",p_adj,"\nPower ",pow,"\nP-value sumlog ",psumlog,"\n")
  
  sbs = subset(gseca, 
               p.adj<=p_adj 
               & pw>=pow 
               & sumlog<=psumlog)
  
  if(!is.na(pemp) && !is.na(unique(gseca$p.emp))){
    sbs = subset(sbs, p.emp<=pemp)
  }
  
  sr_bar = F
  if(!is.na(srate) && !is.na(unique(gseca$success_rate))){
    sr_bar = T
    sbs = subset(sbs, success_rate>=srate)
  }
  
  sel = unique(sbs$gene_set)
  
  if( is.na(unique(gseca$p.emp)) && is.na(unique(gseca$success_rate))){
    sel <- unique(sbs[order(sbs$sumlog, decreasing = F), "gene_set"])[1:20]
    warning("Bootstrapping controls are missing; top 20 gene sets visualized in EC map.")
  }
  
  if(length(sel)==0| is.null(sel)) stop(message("No enriched gene sets"))

  toplot = subset(gseca, gene_set%in%sel)
  toplot$p.adj[ which(toplot$p.adj>p_adj) ]=NA
  toplot$gene_set = fix_names(toplot$gene_set)

  score = unique(toplot[,c('gene_set','sumlog','p.emp','success_rate')])
  
  gs_order = score[order(score[,2],decreasing = F),'gene_set']

  toplot$gene_set = factor(toplot$gene_set, levels=gs_order)
  score$gene_set = factor(score$gene_set, levels=gs_order)

  ECmap  = dcast(data = toplot, gene_set~expr_class, value.var = 'p.adj' ); rownames(ECmap)  = ECmap[,1];  ECmap  = ECmap[,-1]
  ECmap2 = dcast(data = toplot, gene_set~expr_class, value.var = 'delta' ); rownames(ECmap2) = ECmap2[,1]; ECmap2 = ECmap2[,-1]
  ECmap3 = dcast(data = toplot, gene_set~expr_class, value.var = 'CASE.P' ); rownames(ECmap3) = ECmap3[,1]; ECmap3 = ECmap3[,-1]
  ECmap4 = dcast(data = toplot, gene_set~expr_class, value.var = 'CNTR.P' ); rownames(ECmap4) = ECmap4[,1]; ECmap4 = ECmap4[,-1]
  
  limits = range(abs(ECmap2))
  # for(i in 1:nrow(ECmap2)) for(j in 1:ncol(ECmap2)) ECmap2[i,j] = sqrt(abs(ECmap2[i,j])) * sign(ECmap2[i,j]) * 2
  # mx = apply(ECmap2,1, function(x) max(abs(x)))
  # mx = max(abs(ECmap2))
  # ECmap2 = ECmap2/mx
  # normalize(ECmap2[1,], range =c(-1,1), method="range")

  ECmap = as.matrix(ECmap)
  ECmap2 = as.matrix(ECmap2)
  # ECmap2 =round(ECmap2*100,2)
  
  ha = HeatmapAnnotation(name='Bimodal',col_mean = anno_bimodal, height=unit(1,"cm"))
  
  # col_fun = colorRamp2(c(0, 0.05, 0.1), c("red", "yellow", "blue"))
  
  pycol = c(
    rgb(252,231,62, max=255)
    ,rgb(160,218,69, max=255)
    ,rgb(77,192,111, max=255)
    ,rgb(37,160,134, max=255)
    ,rgb(41,125,140, max=255)
    ,rgb(55,89,137, max=255)
    ,rgb(70,51,124, max=255)
    ,rgb(70,31,91, max=255)
    ,rgb(70,31,104, max=255)
  )
  
  col_fun = colorRamp2(c(0,0.00001,0.0001, 0.05, 0.1),pycol[c(1,3,4,7,9)])
  
  # ECcol = c("grey80"
  #           ,rgb(49,39,129, max=255)
  #           ,rgb(255,238,0, max=255)
  #           ,rgb(254,222,150, max=255)
  #           ,rgb(248,179,67, max=255)
  #           ,rgb(227,152,110, max=255)
  #           ,rgb(192,31,53, max=255))
  ECcol = rep("grey80",7)
  # 
  ECM = Heatmap(ECmap, name = "ECMap"
                , col =  col_fun
                , rect_gp = gpar(col = 'grey80', fill=NA ) #type = "none"
                , cell_fun = function(j, i, x, y, width, height, fill) {
                  hback = ECmap3[i, j]
                  myfill = "white"
                  if(!is.na(ECmap[i, j])) myfill = col_fun(ECmap[i, j])

                   grid.polygon(x=x+unit(outer(c(-.5, .5, .5, -.5), width), 'native'),
                                y=y+unit(outer( c(-.5,-.5,-.5+hback,-.5+hback),height),'native'),
                                gp=gpar(fill=ECcol[j], col = 'black', alpha=0.75, lwd=unit(0.25,'points'))) #ECcol[j]
                   
                   # baseline = -.5+hback

                   grid.lines(x=x+unit(outer(c(-.5, .5), width), 'native'), y =y+unit(outer( c(0,0),height),'native'),gp=gpar(col='grey50', lty='dashed',lwd=unit(0.25,'points'))) #,lty='dashed'
                   grid.lines(x=x+unit(outer(c(-.5, .5), width), 'native'), y =y+unit(outer( c(-.5,-.5),height),'native'),gp=gpar(col='black')) #,lty='dashed'
                   grid.lines(x=x+unit(outer(c(-.5, .5), width), 'native'), y =y+unit(outer( c(.5,.5),height),'native'),gp=gpar(col='black'))
                   if(j==1){
                     grid.lines(x=x+unit(outer(c(-.5,-.5), width), 'native'), y =y+unit(outer( c(-.5,.5),height),'native'),gp=gpar(col='black')) #,lty='dashed'
                   }
                   if(j==length(ECcol)){
                     grid.lines(x=x+unit(outer(c(.5,.5), width), 'native'), y =y+unit(outer( c(-.5,.5),height),'native'),gp=gpar(col='black')) #,lty='dashed'
                   }

                   delta = ECmap2[i, j]*10
                   baseline=0                   
                   b = as.numeric(width)*sqrt(abs(delta))*sign(delta)
                   h = as.numeric(height)*sqrt(abs(delta))*sign(delta)
                   
                   if(myfill!='white'){
                      grid.polygon(  x=x+unit(outer(c(-0.5, 0, .5), b), 'native')
                                , y=y+unit(outer(baseline+c(0, 0.5, 0), h),'native'),
                                id.lengths=rep(3, 1), gp=gpar(fill=myfill,col='black',lwd=unit(0.25,'points')))#ifelse(delta>0,ECcol[j], 'white')
                  }
                  }
                , cluster_rows = FALSE, cluster_columns = FALSE,show_row_names = T, show_column_names = T
                , row_names_side = 'left'
                , show_heatmap_legend = F
                , row_names_gp = gpar(fontsize = 8)
                , column_names_gp = gpar(fontsize = 8)
                , top_annotation = ha
                # , row_names_max_width = unit(2,"cm")
  )
  
  
  if(sr_bar & length(unique(score[match(rownames(ECmap), score$gene_set), "success_rate"]))>1){
    row_ha = rowAnnotation(
      ES = row_anno_barplot(matrix(-10*log10(score[match(rownames(ECmap), score$gene_set), "sumlog"])), 
                            axis = TRUE,  border = T, axis_side = "bottom", 
                            gp = gpar(fill = c("black")), bar_width = 0.6)
      
      # , Pemp = row_anno_barplot(matrix((score[match(rownames(ECmap), score$gene_set), "p.emp"])), axis = TRUE, axis_side = "bottom", gp = gpar(fill = c("black")))
      , SR = row_anno_barplot(matrix(score[match(rownames(ECmap), score$gene_set), "success_rate"]), 
                             axis = TRUE, border = T, axis_side = "bottom", 
                             gp = gpar(fill = c("black")), bar_width = 0.6)
      , gap = unit(c(2, 4), "mm")
      , width = unit(3.5, "cm")
      , show_annotation_name = TRUE
      , annotation_name_side = 'top'
      , annotation_name_rot = 0
      , annotation_name_gp = gpar(fontsize= 8)
    )
  } else {
    row_ha = rowAnnotation(
      ES = row_anno_barplot(matrix(-10*log10(score[match(rownames(ECmap), score$gene_set), "sumlog"])), 
                           axis = TRUE,  border = T, axis_side = "bottom", 
                           gp = gpar(fill = c("black")), bar_width = 0.6)  
      # , Pemp = row_anno_barplot(matrix((score[match(rownames(ECmap), score$gene_set), "p.emp"])), axis = TRUE, axis_side = "bottom", gp = gpar(fill = c("black")))
      , gap = unit(c(2, 4), "mm")
      , width = unit(3.75, "cm")
      , show_annotation_name = TRUE
      , annotation_name_side = 'top'
      , annotation_name_rot = 0
      , annotation_name_gp = gpar(fontsize= 8)
    )
  } 
  
  limits = range(abs(ECmap2)*100)
  ml = matrix(limits, nr=1, dimnames = list("Delta (%)", as.character(round(limits))))
  lgha = HeatmapAnnotation(
     which='column'
    ,text = column_anno_text(as.character(round(limits)),gp = gpar(fontsize=8)) #, rot = 90, offset = 1, just = 'right'
    # annotation_height = unit.c(unit(5, "mm"), unit(5, "mm"))
  )
  
  lgH = Heatmap(ml, name = "Legend"
                , rect_gp = gpar(col = NA, fill=NA)
                , cell_fun = function(j, i, x, y, width, height, fill) {
                  delta = ml[i, j]/10
                  baseline=0        
                  b = as.numeric(width)*sqrt(abs(delta))*sign(delta)
                  h = as.numeric(height)*sqrt(abs(delta))*sign(delta)
                  grid.polygon(  x=x+unit(outer(c(-0.25, 0, .25), b), 'native')
                                 , y=y+unit(outer(baseline+c(-.5, 0, -.5), h),'native'),
                                 id.lengths=rep(3, 1), gp=gpar(fill='white',col='black',lwd=unit(0.25,'points')))#ifelse(delta>0,ECcol[j], 'white')
                }
                , cluster_rows = FALSE, cluster_columns = FALSE,show_row_names = T, show_column_names = F
                , column_names_side = 'bottom'
                , show_heatmap_legend = F
                , width = unit(1, 'cm')
                , row_names_side = 'left'
                , row_names_gp = gpar(fontsize = 8)
                , bottom_annotation = lgha

  )
  
  cm = ColorMapping(name = "FDR", col_fun = colorRamp2(c(0, 0.05, 0.1),pycol[c(1,4,9)]))
  
  if(!is.null(filename)){
    # pdf(file=filename,  width = unit(8.27, "inches"), height = unit(11.69, "inches") )
    pdf(file=filename,  paper = "a4", useDingbats = F)
  
    grid.newpage()
    pushViewport(viewport(layout=grid.layout(  nrow=4, heights = c(.05, .8, .1, .05)
                                             , ncol=6, widths  = c( .05, .1, 3.5, 3.5, .1, .05 ))))
    pushViewport(viewport(layout.pos.row = 2, layout.pos.col = 3:4))
    draw(ECM+row_ha, newpage = FALSE)
    upViewport()

    pushViewport(viewport(layout.pos.row = 3, layout.pos.col = 3))
    color_mapping_legend(cm, legend_direction='horizontal', color_bar='continuous'
                         , labels_gp = gpar(fontsize = 8)
                         , legend_width = unit(2, "cm")
                         , legend_height = unit(0.2, "cm")
                         , at=c(0,0.05,0.1)
                         , title_gp = gpar(fontsize = 8))
    upViewport()
    
    pushViewport(viewport(layout.pos.row = 3, layout.pos.col = 4))
    draw(lgH, newpage = FALSE)
    upViewport()
    
    pushViewport(viewport(layout.pos.row = 4, layout.pos.col = 3:4))
    grid.text("GSECA v.0 2018", just="center", gp=gpar(fontsize=6, col='black'))
    upViewport()
  
     dev.off()
    # return(NULL)
    
  }
  # else{
  
    return(ECM+row_ha)
  
  # }  
}

# CUT CODE
# 
# tmp = read.delim('tmp.txt', h=F)
# system('rm tmp.txt')
# tmp = unique(tmp)
# # 
# b = as.numeric(width)*delta
# h = as.numeric(height)*delta
# 
# w = unique(tmp[,1])
# h = unique(tmp[,2])
# 
# delta =  quantile(abs(ECmap2))
# 
# b0 = w*delta['0%']; b1 = w*delta['50%']; b2 = w*delta['100%']; 
# h0 = h*delta['0%']; h1 = h*delta['50%']; h2 = h*delta['100%']; 
# 
# x0 = unit(0.1, "native"); x1 = unit(0.2, "native"); x2 = unit(0.3, "native")
# yp = unit(0.5, "native") - unit(.5, "lines")
# # 
# # tu = c(-0.5, 0, .5)
# # tl = c(-0.5, 0.5, -0.5)
# 
# tu = c(0, 0.5,1)
# tl = c(-0.5, 0.5, -0.5)
# 
# grid.polygon(x = x0 + unit(outer(tu, b0),'native'),
#              y = yp + unit(outer(tl, h0),'native'),
#                id.lengths=rep(3, 1), gp=gpar(fill=NA), default.units = 'native'
#              )
#   
# grid.polygon(x = x1 + unit(outer(tu, b1),'native'),
#              y = yp + unit(outer(tl, h1),'native'),
#              id.lengths=rep(3, 1), gp=gpar(fill=NA), default.units = 'native' 
# )
# 
# grid.polygon(x = x2 + unit(outer(tu, h2),'native'),
#              y = yp + unit(outer(tl, b2),'native'),
#              id.lengths=rep(3, 1), gp=gpar(fill=NA), default.units = 'native'
# )
# 
# grid.rect(x = , y =  unit(0.5, "native") - unit(.5, "lines") , width = boxSize, height = boxSize, just = "bottom", gp = gpar(fill = NA))
# grid.rect(x = unit(0.2, "native"), y =  unit(0.5, "native") - unit(.5, "lines"), width = boxSize, height = boxSize, just = "bottom", gp = gpar(fill = NA))
# grid.rect(x = unit(0.3, "native"), y =  unit(0.5, "native") - unit(.5, "lines"), width = boxSize, height = boxSize, just = "bottom", gp = gpar(fill = NA))
# 
# grid.text(legLabels[1], x = unit(0.1, "native"), y = unit(0.5, "npc") - unit(1, "lines"), gp = gpar(fontsize = 10))
# grid.text(legLabels[2], x = unit(0.2, "native"), y = unit(0.5, "npc") - unit(1, "lines"), gp = gpar(fontsize = 10))
# grid.text(legLabels[3], x = unit(0.3, "native"), y = unit(0.5, "npc") - unit(1, "lines"), gp = gpar(fontsize = 10))
# 

# leg <- function(legLabels) {
#   nlabels <- length(legLabels)
#   pushViewport(viewport(layout = grid.layout(nlabels, 1)))
#   for (i in 1:nlabels) {
#    pushViewport(viewport(layout.pos.row = i))
#     grid.rect(width = boxSize, height = boxSize, just = "bottom",
#                 + gp = gpar(fill = boxColours[i]))
#     grid.text(legLabels[i], y = unit(0.5, "npc") - unit(1, "lines"))
#     popViewport()
#     }
#   popViewport()
#   upViewport()
# }
# g= Legend(at = c(-2, -1, 0, 1, 2), col_fun = col_fun, title = "main matrix", title_position = "topleft")
# grid.draw(g)




orderGSECA = function(x, PSUMLOG, PEMP, SRATE){
  y = unique(x[,c('gene_set','sumlog', 'success_rate','p.emp')])
  y$method='GSECA'
  colnames(y)[2] = 'adj.P.Val'
  
  y$sig = y$adj.P.Val<=PSUMLOG & y$p.emp<=PEMP & y$success_rate>=SRATE
  
  y=y[order(1-y$sig, y$p.emp, y$adj.P.Val, 1-y$success_rate,  decreasing = F),]
  y$rank = 1:nrow(y)
  
  return(y)
}

