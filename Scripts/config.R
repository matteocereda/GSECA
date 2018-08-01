# REQUIREMENTS ==========
# Set .GlobalEnv options
options(stringsAsFactors = F)

# Install/Load required R packages (dependencies)
pkgs <- c( 'ggrepel'
           ,'grid'
           ,'gridExtra'
           ,'statmod'
           ,'R.utils'
           ,'reshape2'
           ,'mixtools'
           ,'scales'
           ,'gtable'
           ,'plyr'
           ,'parallel'
           ,'doParallel'
           ,'snowfall'
           ,'rlecuyer'
           ,'RColorBrewer'
           ,'metap'
           ,'ComplexHeatmap'
           , 'circlize'
)

not_installed <- which( ! pkgs %in% rownames(installed.packages()) )
if(length(not_installed)>0){
  for(i in not_installed) install.packages(pkgs[i])
}

for(i in pkgs) suppressWarnings(suppressPackageStartupMessages(library(i, character.only =T)))

# ROUTINES ===========
# Utility functions
print.logo <- function() {
  cat("██████╗ ███████╗███████╗ ██████╗ █████╗\n")
  cat("██╔════╝ ██╔════╝██╔════╝██╔════╝██╔══██╗\n")
  cat("██║  ███╗███████╗█████╗  ██║     ███████║\n")
  cat("██║   ██║╚════██║██╔══╝  ██║     ██╔══██║\n")
  cat("╚██████╔╝███████║███████╗╚██████╗██║  ██║\n")
  cat("╚═════╝ ╚══════╝╚══════╝ ╚═════╝╚═╝  ╚═╝\n")
}

lapply_pb <- function(X, FUN, ...){
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

mapply_pb <- function(FUN, X, Y,  ...){
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

ids_from_list_names <- function(x,y){
  
  x <- mapply(function(a,b){
    a$id <- b
    return(a)
  },x,y, SIMPLIFY = F)
  
  return(x)
}

fix_names <- function(p) {
  library(R.utils)
  return(capitalize(gsub(pattern = "_", replacement = " ",fixed = T,x = tolower(gsub(pattern = "KEGG_",replacement = "",fixed = T,x = as.character(as.character(p)))))))
}

read.gmt.file <- function(pathMsigDbFile) {
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

create.output.folder <- function (s, odir) {
  if(!dir.exists(odir)) dir.create(odir)
  
  job      = paste0('Job',length(grep(".Job", list.dirs(path = odir, recursive =F)))+1)
  timer    = format(Sys.time(), "%Y-%m-%d.%H_%M_%S")
  
  analisys = paste(job,timer,s, sep=".")
  analisys = paste0(odir,"/", analisys )
  
  dir.create( analisys )
  # for(i in s$sample ) dir.create( paste0( analisys, "/", i))
  return(analisys)
}

get.gene.set <- function(ix) {
  paste0("../gene_sets/msigdb/",ix)
}

# STATS ===========
# Functions for computing statistics
get_padj <- function(cor) {
  function(x) {
    p.adjust(x, method=cor) 
  }
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

fisher.gonzales <- function (x # matrix
                             , alternative = "two.sided"
                             , or = 1
) {
  PVAL <- NULL
  m <- sum(x[, 1L])
  n <- sum(x[, 2L])
  k <- sum(x[1L, ])
  x <- x[1L, 1L]
  lo <- max(0L, k - n)
  hi <- min(k, m)
  support <- lo:hi
  logdc <- dhyper(support, m, n, k, log = TRUE)
  
  dnhyper <- function(ncp) {
    d <- logdc + log(ncp) * support
    d <- exp(d - max(d))
    d/sum(d)
  }
  
  mnhyper <- function(ncp) {
    if (ncp == 0)
      return(lo)
    if (ncp == Inf)
      return(hi)
    sum(support * dnhyper(ncp))
  }
  # 
  pnhyper <- function(q, ncp = 1, upper.tail = FALSE) {
    if (ncp == 1) {
      return(if (upper.tail) phyper(x - 1, m, n, k,lower.tail = FALSE) else phyper(x, m, n, k))
    }
    if (ncp == 0) {
      return(as.numeric(if (upper.tail) q <= lo else q >=  lo))
    }
    
    if (ncp == Inf) {
      return(as.numeric(if (upper.tail) q <= hi else q >=  hi))
    }
    
    sum(dnhyper(ncp)[if (upper.tail) support >= q else support <= q])
  }
  PVAL <- switch(alternative
                 , less    = pnhyper(x, or)
                 , greater = pnhyper(x, or, upper.tail = TRUE), 
                 two.sided = {
                   relErr <- 1 + 10^(-7)
                   d <- dnhyper(or)
                   sum(d[d <= d[x - lo + 1] * relErr])
                 })
  RVAL <- list(p.value = PVAL)
  mle <- function(x) {
    if (x == lo) 
      return(0)
    if (x == hi) 
      return(Inf)
    mu <- mnhyper(1)
    if (mu > x) 
      uniroot(function(t) mnhyper(t) - x, c(0, 1))$root
    else if (mu < x) 
      1/uniroot(function(t) mnhyper(1/t) - x, c(.Machine$double.eps, 
                                                1))$root
    else 1
  }
  
  ESTIMATE <- c(`odds ratio` = mle(x))
  
  RVAL <- c(RVAL, list(
    estimate = ESTIMATE
  ))
  return(RVAL)
}

FISHER <- function (u, TAIL = TAIL){
  # Compute Fisher's Exact test
  u <- unlist(u)
  
  a <- as.numeric(u["CASE.N.class"]); b <- as.numeric(u['CASE.N.rest'])
  c <- as.numeric(u["CNTR.N.class"]); d <- as.numeric(u['CNTR.N.rest'])
  
  if( (b==0 & c==0) | (a==0 & d==0) ){ 
    return( c(
      'pv' = 0,
      'or' = Inf
    ) )
  }
  
  m <- matrix( as.numeric(u[c("CASE.N.class", "CASE.N.rest", "CNTR.N.class", "CNTR.N.rest")]), nc=2, byrow = T )
  # f <- fisher.test(m, alternative = TAIL)
  f <- fisher.gonzales(m, alternative = TAIL)
  
  
  c(
    'pv' = f$p       ,
    'or' = f$estimate
  )
}
# PREPARE INPUT ======
# Process input files for GSECA analysis. Required inputs are:
# M = Gene Expression matrix
# L = Sample label list
get_expression_dataset <- function(M=NULL, L=NULL, gene_symbol="ensembl_gene_id", annot=NULL){
  
  if(is.null(M) | is.null(L)){
    stop(message('Please upload expression matrix and sample lables'))
  }else{
    
    if(length(L)!=(ncol(M)-1)) stop(message("Number of labels not equal to number of samples in the matrix"))
    require(reshape2)
    
    colnames(M)[1] <- "geneID"
    names(L) <- colnames(M)[2:ncol(M)]
    expr <- melt(M,id.vars = "geneID")
    colnames(expr) <- c('geneID', 'SampleID', 'value')
    expr$geneID    <- as.character(expr$geneID)
    expr$SampleID  <- as.character(expr$SampleID)
    expr$value <- as.numeric(expr$value)
    if (is.null(annot)) {
      # Built-in annotation: GencodeV26
      load("annotation/gencode.v26.annotation.Rdata")
      colnames(genes)[c(6,8)] <- c("ensembl_gene_id", "symbol")
    } else{
      genes <- annot
    }
    
    expr$gene_type <- as.character(genes[match(expr[,'geneID'], genes[,gene_symbol]), c('gene_type')])
    expr$type <- L[expr$SampleID]
    
    return(expr)
  }
}

prepare_pl <- function(pl, expr, annot=NULL) {
  
  if( length(grep("ENS", expr$geneID[1]))!=0 & length(grep("ENS", pl[[1]]))==0 ) {
    if(is.null(annot)){
      load("annotation/gencode.v26.annotation.Rdata") 
    }
    colnames(genes)[c(6,8)] <- c("ensembl_gene_id", "symbol")
    pl <- lapply(pl, function(x) na.omit(genes$ensembl_gene_id[match(x, genes$symbol)]))
  }
  
  return(pl)
}

prepare_expr <- function(pl, expr) {
  
  total_genes    <- unique(unlist(pl))
  rna <- vector(mode = 'list', length = 2)
  rna$expr_class <- lapply(expr$expr_class, 
                           function(x) subset(x, geneID%in%total_genes, select = c("geneID", "expr_class")))
  
  rna$mixture_model <- unique(expr$mixture_model[,c('SampleID','type')])
  
  return(rna)
}

# MIXTURE MODEL =============
# Functions to fit a two-component Gaussian Mixture model to gene expression
# distributions and define the expression classes.
GMM_posterior <- function(x
                          , my_perc  = c('20%','80%') # quantiles of expression distribution to start with
                          , iterations = 100 # max number of iterations to reach correct results
                          , nm_k = 2
                          , nm_fast = TRUE
                          , nm_maxit = 10000
                          , nm_epsilon =  0.001
){
  
  # x <- data.frame("log_value" = x)
  na <- subset(x, is.na(log_value) | !is.finite(log_value) )
  x  <- subset(x, !is.na(log_value) & is.finite(log_value) )
  
  perc = quantile(x$log_value, seq(0,1,.05))[ my_perc ]
  
  #initialization of internal variables
  normMix     = NA
  infloop     = 0
  threshold_1 = 1
  threshold_2 = 0
  mean_1      = (-10)
  mean_2      = 10
  lambda1 = 1
  lambda2 = 0
  st1 = 1
  st2 = 0
  
  while( (threshold_1>threshold_2) | ( lambda2 - lambda1 <= 0.1 ) ) {
    
    normMix <- normalmixEM(  x$log_value
                             , mu      = perc
                             , k       = nm_k
                             , fast    = nm_fast
                             , maxit   = nm_maxit
                             , epsilon = nm_epsilon
    )
    
    mean_1 <- min(normMix$mu)
    mean_2 <- max(normMix$mu)
    lambda1 <- normMix$lambda[which.min(normMix$mu)]
    lambda2 <- normMix$lambda[which.max(normMix$mu)]
    st1     <- normMix$sigma[which.min(normMix$mu)]
    st2     <- normMix$sigma[which.max(normMix$mu)]
    threshold_1 <- mean_1 + normMix$sigma[which.min(normMix$mu)]
    threshold_2 <- mean_2
    
    infloop=infloop+1 # to avoid infinite loops
    
    cat( 't1:\t',threshold_1 ,'t2:\t',threshold_2
         ,'m1:\t' ,mean_1 ,'m2:\t',mean_2
         ,'s1:\t' ,st1 ,'s2:\t',st2
         ,'l1:\t' ,lambda1 ,'l2:\t',lambda2
         ,'infloop:\t',infloop,'\n')
    
    if(infloop > iterations){
      break
    }
  }
  
  na$comp.1 <- NA
  na$comp.2 <- NA
  x$comp.1  <- normMix$posterior[,1]
  x$comp.2  <- normMix$posterior[,2]
  
  x <- as.data.frame(rbind(x,na))
  x <- unrowname(x)
  
  x$mean1   <- min(normMix$mu)
  x$mean2   <- max(normMix$mu)
  x$st1     <- normMix$sigma[which.min(normMix$mu)]
  x$st2     <- normMix$sigma[which.max(normMix$mu)]
  x$lambda1 <- normMix$lambda[which.min(normMix$mu)]
  x$lambda2 <- normMix$lambda[which.max(normMix$mu)]
  
  return(x)
}


assignToComp <- function(x, ne=T, posterior_cutoff=.9) 
{
  
  # The function assign genes to one of the two 
  # components of the FMM, using the posterior
  # probabilities of the model
  message(paste0("posterior_cutoff = ", posterior_cutoff))
  x$expr_class = 'ME'
  
  if(ne) {
    i.na = which(is.na(x$log_value))
    if(length(i.na)>0) x$expr_class [ i.na ] = "NE"
  } 
  
  x$expr_class [ which(x$comp.1>=posterior_cutoff) ] = 'LE'
  x$expr_class [ which(x$comp.2>=posterior_cutoff) ] = 'HE'
  
  return(x)
}

setExpressionClass2 <- function(x, ...) 
{
  # Define 2 expression classes
  message(" -- # 2 classes (LE, HE)")
  x$expr_class <- with(x, ifelse(comp.2 > comp.1, "HE","LE" ) )
  x$expr_class <- factor(x$expr_class, levels=c("LE","HE"))
  
  return(x)
}

setExpressionClass3 <- function(x, ...) 
{
  # Define 3 expression classes
  message(" -- # 3 classes (LE, ME, HE)")
  x <- assignToComp(x, ne=F, ...)
  
  x$expr_class = factor(x$expr_class, levels=c("LE","ME","HE"))
  
  return(x)
}

setExpressionClass4 <- function(x, ...) 
{
  # Define 4 expression classes
  message(" -- # 4 classes (NE, LE, ME, HE)")
  x <- assignToComp(x, ne=T, ...)
  
  x$expr_class = factor(x$expr_class, levels=c("NE","LE","ME","HE"))
  
  return(x)
}

setExpressionClass5 <- function(x, ...) 
{
  # Define 5 expression classes
  message(" -- # 5 classes (NE, LE, ME, HE1, HE2)")
  x <- assignToComp(x, ne=T, ...)
  
  qt = quantile(subset(x, expr_class=='HE')$log_value)
  x$expr_class [ which(x$expr_class=='HE' & x$log_value<qt["50%"]  ) ] = 'HE1'
  x$expr_class [ which(x$expr_class=='HE' & x$log_value>=qt["50%"] ) ] = 'HE2'
  x$expr_class = factor(x$expr_class, levels=c("NE","LE","ME","HE1","HE2"))
  
  return(x)
}

setExpressionClass6 <- function(x, ...) 
{
  # Define 6 expression classes
  message(" -- # 6 classes (NE, LE, ME, HE1, HE2, HE3)")
  x <- assignToComp(x, ne=T, ...)
  
  qt = quantile(subset(x, expr_class=='HE')$log_value, seq(0,1,0.01))
  x$expr_class [ which(x$expr_class=='HE' & x$log_value<qt["33%"]  ) ] = 'HE1'
  x$expr_class [ which(x$expr_class=='HE' & x$log_value>=qt["33%"] & x$log_value<qt["66%"]) ] = 'HE2'
  x$expr_class [ which(x$expr_class=='HE' & x$log_value>=qt["66%"] ) ] = 'HE3'
  x$expr_class = factor(x$expr_class, levels=c("NE","LE","ME","HE1","HE2","HE3"))
  
  return(x)
}

setExpressionClass7 <- function(x, ...)
{
  # Define 7 expression classes
  message(" -- # 7 classes (NE, LE, ME, HE1, HE2, HE3, HE4)")
  x <- assignToComp(x, ne=T, ...)
  
  qt = quantile(subset(x, expr_class=='HE')$log_value)
  
  x$expr_class [ which(x$expr_class=='HE' & x$log_value<=qt["25%"] ) ] = 'HE1'
  x$expr_class [ which(x$expr_class=='HE' & x$log_value>=qt["25%"] & x$log_value<=qt["50%"]) ] = 'HE2'
  x$expr_class [ which(x$expr_class=='HE' & x$log_value>=qt["50%"] & x$log_value<=qt["75%"]) ] = 'HE3'
  x$expr_class [ which(x$expr_class=='HE' & x$log_value>=qt["75%"] ) ] = 'HE4'
  
  x$expr_class = factor(x$expr_class, levels = c('NE','LE','ME','HE1','HE2','HE3','HE4'))
  
  return(x)
}

setExpressionClass8 <- function(x, ...)
{
  # Define 8 expression classes
  message(" -- # 8 classes (NE, LE1, LE2, ME, HE1, HE2, HE3, HE4)")
  x <- assignToComp(x, ne=T, ...)
  
  qt = quantile(subset(x, expr_class=='HE')$log_value)
  
  x$expr_class [ which(x$expr_class=='HE' & x$log_value<=qt["25%"] ) ] = 'HE1'
  x$expr_class [ which(x$expr_class=='HE' & x$log_value>=qt["25%"] & x$log_value<=qt["50%"]) ] = 'HE2'
  x$expr_class [ which(x$expr_class=='HE' & x$log_value>=qt["50%"] & x$log_value<=qt["75%"]) ] = 'HE3'
  x$expr_class [ which(x$expr_class=='HE' & x$log_value>=qt["75%"] ) ] = 'HE4'
  
  qt = quantile(subset(x, expr_class=='LE')$log_value)
  x$expr_class [ which(x$expr_class=='LE' & x$log_value<qt["50%"]  ) ] = 'LE1'
  x$expr_class [ which(x$expr_class=='LE' & x$log_value>=qt["50%"] ) ] = 'LE2'
  
  x$expr_class = factor(x$expr_class, levels = c('NE','LE1',"LE2", 'ME','HE1','HE2','HE3','HE4'))
  
  return(x)
  
}

setExpressionClass9 <- function(x, ...)
{
  # Define 9 expression classes
  message(" -- # 9 classes (NE, LE1, LE2, ME1, ME2, HE1, HE2, HE3, HE4)")
  x <- assignToComp(x, ne=T, ...)
  
  qt = quantile(subset(x, expr_class=='LE')$log_value)
  x$expr_class [ which(x$expr_class=='LE' & x$log_value<qt["50%"]  ) ] = 'LE1'
  x$expr_class [ which(x$expr_class=='LE' & x$log_value>=qt["50%"] ) ] = 'LE2'
  qt = quantile(subset(x, expr_class=='ME')$log_value)
  x$expr_class [ which(x$expr_class=='ME' & x$log_value<qt["50%"]  ) ] = 'ME1'
  x$expr_class [ which(x$expr_class=='ME' & x$log_value>=qt["50%"] ) ] = 'ME2'
  
  x$expr_class = factor(x$expr_class, levels = c('NE','LE1',"LE2", 'ME1', 'ME2', 'HE1','HE2','HE3','HE4'))
  
  return(x)
}

setExpressionClass10 <- function(x, ...)
{
  # Define 10 expression classes
  message(" -- # 10 classes (NE, LE1, LE2, ME1, ME2, HE1, HE2, HE3, HE4, HE5)")
  x <- assignToComp(x, ne=T, ...)
  
  # LE
  qt = quantile(subset(x, expr_class=='LE')$log_value)
  x$expr_class [ which(x$expr_class=='LE' & x$log_value<qt["50%"]  ) ] = 'LE1'
  x$expr_class [ which(x$expr_class=='LE' & x$log_value>=qt["50%"] ) ] = 'LE2'
  
  # HE
  qt = quantile(subset(x, expr_class=='HE')$log_value, seq(0,1,length.out = 6))
  x$expr_class [ which(x$expr_class=='HE' & x$log_value<=qt[2] ) ] = 'HE1'
  x$expr_class [ which(x$expr_class=='HE' & x$log_value>=qt[2] & x$log_value<=qt[3]) ] = 'HE2'
  x$expr_class [ which(x$expr_class=='HE' & x$log_value>=qt[3] & x$log_value<=qt[4]) ] = 'HE3'
  x$expr_class [ which(x$expr_class=='HE' & x$log_value>=qt[4] & x$log_value<=qt[5]) ] = 'HE4'
  x$expr_class [ which(x$expr_class=='HE' & x$log_value>=qt[5] ) ] = 'HE5'
  # ME
  qt = quantile(subset(x, expr_class=='ME')$log_value)
  x$expr_class [ which(x$expr_class=='ME' & x$log_value<qt["50%"]  ) ] = 'ME1'
  x$expr_class [ which(x$expr_class=='ME' & x$log_value>=qt["50%"] ) ] = 'ME2'
  
  
  x$expr_class = factor(x$expr_class, levels = c('NE','LE1',"LE2", 'ME1', "ME2", 'HE1','HE2','HE3','HE4','HE5') )
  
  return(x)
}

set_expression_class_on_posterior <- function(x
                                              , posterior_cutoff=.9
                                              , nClass=7) 
{
  setFunc <- get(grep(nClass, ls(pos = 1), value = T))
  x <- setFunc(x, posterior_cutoff)
  return(x)
}

get_mixture_expr_class <- function(rna
                                   , ne_value=0.01
                                   , nClass=7
                                   , ncor=2){
  
  message(" === Finite Mixture Model === \n", 
          "[*] NE threshold: "                , ne_value, "\n", 
          "[*] Number of Expression Classes: ", nClass  , "\n", 
          " -- Parallelizing "                , ncor    , " cores ...")
  
  samples <- unique(rna[ ,c("SampleID", "type")])
  samples <- unrowname(samples)
  rna <- rna[ ,c("SampleID", "geneID", "value")]
  
  rna$log_value <- NA
  rna$log_value[which(rna$value>ne_value)] <- log2(rna$value[which(rna$value>ne_value)])
  
  message("[*] FMM Getting thresholds ...")
  
  l <- dlply(rna[,c("SampleID","geneID", "log_value")], ~SampleID, 
             function(x) x[,c("geneID","log_value")])
  
  if(.Platform$OS.type=="unix"){
    # on unix-like architectures
    l <- mclapply(l, function(x){
      tryCatch(GMM_posterior(x), error = function(x) return(NA))
    }
    , mc.silent = T, mc.cores = ncor
    )
  } else {
    # on other OS, e.g. Windows
    sfInit(parallel=TRUE, cpus=ncor, type="SOCK")
    sfLibrary(mixtools)
    sfLibrary(plyr)
    sfExport(list = c('GMM_posterior'))
    sfClusterSetupRNG()
    l <- sfLapply(l, function(x){
      tryCatch(GMM_posterior(x), error = function(x) return(NA))
    }
    )
    sfStop()
  }
  
  n_na <-  length(which(is.na(l)))
  if ( n_na != 0 ) {
    warning(paste0(" -- Excluding ", n_na," samples not fitted by the mixture"))
    nm <- names(which(is.na(l)))
    l[which(is.na(l))] <- NULL
    samples <- samples[!samples[, "SampleID"]%in%nm, ]
  }
  
  stats <- lapply(l, function(x) x[1,c('mean1','mean2','st1','st2','lambda1','lambda2')])
  stats <- do.call(rbind,stats)
  stats <- cbind(samples, stats[match(rownames(stats), samples$SampleID),])
  stats <- unrowname(stats)
  
  message("[*] FMM Setting expression classes ...")
  l <- lapply(l, function(x) x[,!colnames(x)%in%c('mean1','mean2','st1','st2','lambda1','lambda2')])
  l <- suppressMessages(lapply_pb(l, set_expression_class_on_posterior, posterior_cutoff=.9, nClass=nClass))
  
  # rna <- do.call("rbind", l)
  
  return(list( 'expr_class'=l
               ,'mixture_model'=stats))
}

# GENE CLASS REPRESENTATION ===========
# Representation of genes per expression classs
gene_class_representation <- function(r, pl
                                      , method = 'fisher'
                                      , TAIL = "two.sided"
                                      , correction = 'fdr'
                                      , cutoff = 0.01
){
  
  require(plyr)
  require(metap)
  
  total_genes <- unique(unlist(pl))
  
  message(" === Gene Class Representation === \n", 
          "[*] # gene sets:\t",   length(pl))
  
  cn <- c("geneID", "expr_class")
  
  samplesID <- unique(r$mixture_model[,c("SampleID","type")])
  typeID <- as.list(t(samplesID[,2])); names(typeID) <- samplesID$SampleID
  samples <- table(samplesID$type)
  rm(samplesID)
  
  rna <- mapply(function(a, b) {
    a <- subset(a, geneID%in%total_genes)[,cn]
    a$type <- b
    return(a)
  }, a=r$expr_class, b=typeID, SIMPLIFY = F)
  
  rna <- do.call(rbind, 
                 lapply(rna, function(x) x[,c('geneID', 'expr_class', 'type')])
  )
  
  t <- with(rna, table(geneID, expr_class,type))
  
  case <- as.data.frame(t[,,"CASE"]);
  case$N.not_in_class <- samples['CASE'] - case$Freq
  colnames(case)[3:4] <- c("CASE.N.class","CASE.N.rest")
  
  cntr <- as.data.frame(t[,,"CNTR"]);
  cntr$N.not_in_class <- samples['CNTR'] - cntr$Freq
  colnames(cntr)[3:4] <- c("CNTR.N.class","CNTR.N.rest")
  
  ix <- paste0(case$geneID,".",case$expr_class)
  iy <- paste0(cntr$geneID,".",cntr$expr_class)
  
  cnts <- as.data.frame(cbind(case, cntr[match(ix, iy),3:4]))
  
  rm(case, cntr)
  
  cnts$CASE.P <- with(cnts,CASE.N.class/(CASE.N.class+CASE.N.rest))
  cnts$CNTR.P <- with(cnts,CNTR.N.class/(CNTR.N.class+CNTR.N.rest))
  cnts$delta  <- with(cnts, CASE.P - CNTR.P  )
  cnts$direction <- NA
  cnts$direction[which(cnts$delta<0)] <- "D"
  cnts$direction[which(cnts$delta>0)] <- "E"
  
  invariant <- subset(cnts, is.na(direction))
  cnts      <- subset(cnts, !is.na(direction))
  
  message(paste0("[*] GCR Testing genes with ",method," ..."))
  
  if(method=="fisher"){
    stats <- t(apply(cnts, 1, FISHER, TAIL=TAIL))
  }
  
  # set pvalue == 0 as minimum measured value
  min_pv = min(stats[ which(stats[,1]>0) ,1 ], na.rm = T)
  
  if(min_pv>0.005) min_pv = 0.005
  
  stats[which(stats[,1]==0),1] = min_pv
  
  invariant$pv <- 1
  invariant$or <- 0
  
  cnts$pv <- stats[,1]
  cnts$or <- stats[,2]
  
  cnts <- rbind(cnts, invariant)
  
  message("[*] GCR multiple testing correction ...")
  
  # Fisher's method
  message("[*] GCR scoring ...")
  
  cnts <- dlply(cnts, .(geneID), here(mutate)
                , p.adj  = get_padj(correction)(pv)
                , sumlog = suppressWarnings(sumlog(pv)$p)
                , .progress = 'text'
  )
  
  cnts <- lapply(cnts, function(x) {
    to.ord = unique(x$sumlog)
    ord    = 1:length(to.ord)
    n.ord  = sort(to.ord, decreasing = F)
    rank   = ord[ match( to.ord , n.ord ) ]
    tmp    = data.frame(to.ord, rank)
    x$rank <- tmp$rank[match(x$sumlog,tmp$to.ord)]
    x <- x[order(x$rank, x$expr_class),]
    x$sig <- FALSE
    
    x$sig = x$p.adj<=cutoff
    
    x <- x[, c( 'geneID', 'expr_class',
                'CASE.N.class','CASE.N.rest',
                'CNTR.N.class','CNTR.N.rest',
                'CASE.P','CNTR.P','delta',
                'direction','or','pv', "p.adj",
                'sig','rank',"sumlog"
    )]
    return(x)
    rm(tmp, ord, to.ord, n.ord, rank)  
  })
  
  cnts <- unrowname(do.call(rbind, cnts))
  return(cnts)
}

# GSECA ==============
# Implementation of the Gene Set Enrichment Class Analysis 
# functions for the execution of GSECA
GSECArunner <- function(gene_set, geneClass
                        , method='fisher'
                        , TAIL= 'two.sided'
) {
  require(plyr)
  
  geneClass <- subset(geneClass, geneID%in%unique(unlist(gene_set)))
  
  if ( length(geneClass[,1]) != 0) {
    # Perform analysis only for
    # gene sets with matched genes
    # in the input matrix
    px <- ddply(geneClass, .(expr_class), summarise
                , CASE.N.class = sum(CASE.N.class)
                , CASE.N.rest  = sum(CASE.N.rest)
                , CNTR.N.class = sum(CNTR.N.class)
                , CNTR.N.rest  = sum(CNTR.N.rest)
                , .drop=F
    )
    
    px$CASE.P <- with(px, CASE.N.class/(CASE.N.class+CASE.N.rest))
    px$CNTR.P <- with(px, CNTR.N.class/(CNTR.N.class+CNTR.N.rest))
    
    px$delta  <- with(px, CASE.P - CNTR.P  )
    
    px$direction <- NA
    px$direction[which(px$delta<0)] <- "D"
    px$direction[which(px$delta>0)] <- "E"
    
    if(method=='fisher'){
      stats <- t(apply(px, 1, FISHER, TAIL=TAIL))
    }
    px$pv = stats[,1]
    px$or = stats[,2]
    
    return(px)
  } else {
    return(NA)
  }
}

GSECA_core <- function(gcr
                       , pl
                       , method='fisher'
                       , TAIL= 'two.sided'
                       , correction='fdr'
) {
  
  px <- mclapply(pl, GSECArunner
                 , geneClass = gcr
                 , method    = method
                 , TAIL      = TAIL
                 , mc.cores = 3, mc.silent = T
  )
  
  n_na <-  length(which(is.na(px)))
  if ( n_na != 0 ) {
    warning(paste0(" -- Excluding ", n_na," gene sets with no genes in input matrix"))
    nm <- names(which(is.na(px)))
    px[which(is.na(px))] <- NULL
  }
  
  px <- mapply(function(x,y){x$gene_set=y;return(x)}, px, names(px), SIMPLIFY = F)
  px <- as.data.frame(do.call('rbind',px))
  px <- unrowname(px)
  
  message("[*] GSECA scoring ...")
  
  # Fisher's method
  px <- dlply(px, .(gene_set), here(mutate)
              , p.adj  = get_padj(correction)(pv)
              , sumlog = suppressWarnings(sumlog(pv)$p)
              , .progress = 'text'
  )
  
  px <- lapply(px, function(x) x[,c("gene_set", "expr_class", "CASE.N.class",
                                    "CASE.N.rest", "CNTR.N.class", "CNTR.N.rest",
                                    "CASE.P", "CNTR.P", "delta", "direction", "pv",
                                    "or","p.adj","sumlog")])
  return(px)
}

GSECA <- function(  pl
                    , gcr
                    , method ='fisher' # fisher
                    , TAIL = "two.sided"
                    , correction='fdr'
                    , cutoff = 0.01
){
  
  message(" === Gene Set Enrichment Class Analysis === ")
  message("[*] GSECA core ...")
  
  core <- GSECA_core(gcr, pl, method=method, TAIL=TAIL, correction=correction)
  
  # Ranking
  core <- lapply(core, function(x) {
    to.ord = unique(x$sumlog)
    ord    = 1:length(to.ord)
    n.ord  = sort(to.ord, decreasing = F)
    rank   = ord[ match( to.ord , n.ord ) ]
    tmp    = data.frame(to.ord, rank)
    x$rank <- tmp$rank[match(x$sumlog,tmp$to.ord)]
    x <- x[order(x$rank, x$expr_class),]
    x$sig <- FALSE
    
    x$sig = x$p.adj<=cutoff
    
    x <- x[,c("gene_set","expr_class","CASE.N.class",
              "CASE.N.rest","CNTR.N.class","CNTR.N.rest",
              "CASE.P","CNTR.P", "delta","direction","pv",
              "or", "p.adj", "sig",'rank',"sumlog")]
    return(x)
    rm(tmp, ord, to.ord, n.ord, rank)  
  })
  
  core <- unrowname(do.call(rbind, core))
  return(core)
}

orderGSECA <- function(x, PSUMLOG, PEMP, SRATE){
  
  # Rank pathways using 3 metrics:
  # AS, SR and emprical P-Value
  
  y = unique(x[,c('gene_set','sumlog', 'success_rate','p.emp')])
  y$method='GSECA'
  colnames(y)[2] = 'adj.P.Val'
  
  y$sig = y$adj.P.Val<=PSUMLOG & y$p.emp<=PEMP & y$success_rate>=SRATE
  
  y=y[order(1-y$sig, y$adj.P.Val,  y$p.emp, 1-y$success_rate,  decreasing = F),]
  y$rank = 1:nrow(y)
  
  return(y)
}

# BOOTSTRAPPING ====
# Implementation of bootstrapping procedures
# [*] Calculation of empirical P value
# [*] Correction for sample size

execute_GSECA_bootstrap <- function(sample,  rna, pl
                                    , type='sample.size'
                                    # GRC
                                    , GRC.method = 'fisher'
                                    , GRC.TAIL = "two.sided"
                                    , GRC.correction = 'fdr'
                                    , GRC.cutoff = 0.01
                                    
                                    # GSECA
                                    , GSECA.method ='fisher'
                                    , GSECA.TAIL = "two.sided"
                                    , GSECA.correction='fdr'
                                    , GSECA.cutoff = 0.01
){
  
  if(type=='p.empirical'){
    
    rna$mixture_model$type <- sample$type[match(as.character(rna$mixture_model$SampleID), sample$barcode)]
    
  }else if(type=='sample.size'){
    
    idx <- match(sample, names(rna$expr_class))
    rna$expr_class[-idx] <- NULL
    rna$mixture_model <- subset(rna$mixture_model, SampleID%in%sample)
    gc()
  }
  rm(sample)
  dd <- gene_class_representation(rna, pl
                                  , method = GRC.method
                                  , TAIL   = GRC.TAIL
                                  , cutoff = GRC.cutoff
                                  , correction = GRC.correction
  )
  rm(rna)
  res <- GSECA(pl, dd
               , method =  GSECA.method
               , TAIL   =  GSECA.TAIL
               , cutoff =  GSECA.cutoff
               , correction = GSECA.correction
  )
  rm(dd)
  rx <- unique(res[,c("gene_set", "sumlog")])
  rm(res)
  return(rx)
}

GSECA.Bootstrap.empirical <- function( rna, obs.gseca, pl
                                       , NSIM=10
                                       , nc=2
                                       , method = 'fisher'
                                       , TAIL = "two.sided"
                                       , correction = 'fdr'
                                       , cutoff = 0.01
){
  
  boot <- unique(rna$mixture_model[,c('SampleID','type')])
  
  samples <- 1:NSIM
  samples <- lapply(samples, 
                    function(x,y,z) data.frame(barcode=z,type=sample(y))
                    , y=as.character(boot$type), z=as.character(boot$SampleID)  
  )
  rm(boot)
  
  functionsToCluster <- as.character(lsf.str(envir = .GlobalEnv))
  dataToCluster      <- ls()
  sfInit(parallel=TRUE, cpus=nc, type="SOCK")
  sfExport(list = c(functionsToCluster, dataToCluster))
  sfLibrary(parallel)
  sfClusterSetupRNG()
  l <- sfLapply(samples, function(x){
    tryCatch(
      execute_GSECA_bootstrap(sample=x, rna=rna, pl=pl, type='p.empirical'
                              , GRC.method=method
                              , GRC.TAIL=TAIL
                              , GRC.correction=correction
                              , GRC.cutoff=cutoff
                              , GSECA.method=method
                              , GSECA.TAIL=TAIL
                              , GSECA.correction=correction
                              , GSECA.cutoff=cutoff
      )
      , error = function(x) return(NA)
    )
  })
  sfStop()
  
  l <- do.call(rbind.data.frame, lapply(l, function(x) unique(x[,c("gene_set", "sumlog")])))
  
  message("[*] Compute P.emp ...")
  
  l$obs <- obs.gseca$sumlog[match(l$gene_set,obs.gseca$gene_set)]
  
  l <- ddply(l, .(gene_set), summarise
             
             , p.emp = (1+sum(sumlog<obs, na.rm = T))/(1+length(na.omit(sumlog)))
             
  )
  
  return(l)
}

GSECA.Bootstrap.sample_size <- function( rna, obs.gseca, pl
                                         , phen = c("CASE","CNTR")
                                         , NSIM=1000
                                         , nc=3
                                         , cutoff=0.01
                                         , method = 'fisher'
                                         , TAIL = "two.sided"
                                         , correction = 'fdr'
){
  
  message("GSECA.Bootstrap Create Phenotype lists ...")
  
  pheno <- list()
  pheno[[1]] <- as.character(unique(subset(rna$mixture_model, type==phen[1])$SampleID))
  pheno[[2]] <- as.character(unique(subset(rna$mixture_model, type==phen[2])$SampleID))
  
  if(length(pheno[[1]])!=length(pheno[[2]])){
    
    ref.idx     <- which.max(c(length(pheno[[1]]), length(pheno[[2]])))
    size.idx    <- which.min(c(length(pheno[[1]]), length(pheno[[2]])))
    sample.size <- length(pheno[[size.idx]])
    
    message("GSECA.Bootstrap Bootstrapping ...")
    
    r <- foreach(i=1:NSIM) %do% {
      sample_barcode = vector(mode = "numeric", length = 2*sample.size)
      sample_barcode[1:sample.size] = pheno[[size.idx]]
      sample_barcode[(sample.size+1):length(sample_barcode)] = sample(pheno[[ref.idx]]
                                                                      , length(pheno[[size.idx]]))
      return(sample_barcode)
    }
    
    functionsToCluster <- as.character(lsf.str(envir = .GlobalEnv))
    dataToCluster      <- ls()
    sfInit(parallel=TRUE, cpus=nc, type="SOCK")
    sfExport(list = c(functionsToCluster, dataToCluster))
    sfLibrary(parallel)
    sfClusterSetupRNG()
    l <- sfLapply(r, function(x){
      tryCatch(
        execute_GSECA_bootstrap(sample=x, rna=rna, pl=pl, type='sample.size'
                                , GRC.method=method
                                , GRC.TAIL=TAIL
                                , GRC.correction=correction
                                , GRC.cutoff=cutoff
                                , GSECA.method=method
                                , GSECA.TAIL=TAIL
                                , GSECA.correction=correction
                                , GSECA.cutoff=cutoff
        )
        , error = function(x) return(NA)
      )
    })
    sfStop()
    
    l <- do.call(rbind.data.frame, lapply(l, function(x) unique(x[,c("gene_set", "sumlog")])))
    
    message("[*] Compute Success Rate ...")
    
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

# EXECUTORs =========
# [*] Main Executor
# Wraps and execute all functions 
# needed for the analysis
GSECA_executor <- function( M
                           , L
                           , symbol
                           , geneset
                           , s.test = 'fisher'
                           , tail.test = 'two.sided'
                           , correction = 'fdr'
                           , p_adj_th = 0.1
                           , nClass = 7
                           , analysis = NULL
                           , outdir = NULL
                           , N.CORES = 3
                           , EMPIRICAL = F
                           , BOOTSTRP = F
                           , nsim = 10
                           , PSUMLOG = 0.01
                           , PADJ    = 0.1
                           , PEMP    = 1
                           , SRATE   = 0.7
){
  
  print.logo()
  
  if( N.CORES > nsim ) N.CORES <- nsim
  
  # 01. Prepare Data 
  expr <- get_expression_dataset(M, L, symbol)
  
  # 02. MIXTURE MODEL
  expr <- get_mixture_expr_class(expr
                                 , nClass   = nClass
                                 , ne_value = 0.01
                                 , nc       = N.CORES
  )
  
  # 03. GCR 
  gcr  <- gene_class_representation(expr
                                   , pl = geneset
                                   , method = s.test
                                   , TAIL = tail.test
  )
  
  # 04. GSECA 
  gseca <- GSECA(  pl=geneset
                  , gcr
                  , method = s.test # fisher
                  , TAIL = tail.test
                  , correction = correction
                  )
  
  gseca$p.emp = NA
  gseca$success_rate = NA
  
  if(EMPIRICAL){
    # 05. EMPIRICAL PVALUES
    message(" === Empirical P-Value estimation === ")
    gseca.emp <- GSECA.Bootstrap.empirical( expr
                                           , gseca
                                           , pl = geneset
                                           , NSIM = nsim
                                           , nc = N.CORES
                                           , method = s.test
                                           , TAIL = tail.test
                                           , correction = correction
                                           , cutoff = 0.01
    ) 
    } else {
      gseca.emp <- NULL
    }
    
    if(!is.null(gseca.emp)){
      gseca$p.emp <- gseca.emp$p.emp[match(gseca$gene_set, gseca.emp$gene_set)] 
    }
  
 
  if(BOOTSTRP){
    # 07. CORRECTION FOR SAMPLE SIZE 
    message(" === Success Rate for sample size === ")
    gseca.boot <- GSECA.Bootstrap.sample_size( expr
                                              , gseca
                                              , pl = geneset
                                              , phen = c("CASE","CNTR")
                                              , NSIM = nsim
                                              , nc = N.CORES
                                              , method = s.test
                                              , TAIL = tail.test
                                              , correction = correction
                                              , cutoff = 0.01
    ) 
    } else {
      gseca.boot <- NULL
    }
    
    if(!is.null(gseca.boot)){
      gseca$success_rate <- gseca.boot$success_rate[match(gseca$gene_set, 
                                                         gseca.boot$gene_set)]
    }
  
  
  
  if( !is.null(analysis) & !is.null(outdir)  ) {
    
    analysis <- create.output.folder(analysis, outdir) 
    
    outfig <- paste0(analysis,"/ECmap.psumlog_",PSUMLOG,
                     "_padj_",PADJ,
                     "_pemp_",PEMP,
                     "_srate_",SRATE,".pdf") 
    
    # Results 
    message(paste0("Saving files in folder: ", analysis))
    gseca$gene_set <- fix_names(gseca$gene_set)
    save(expr, gcr, gseca, file=paste0(analysis,"/results.Rdata"))
    write.csv(gseca, file=paste0(analysis,"/results.csv"), row.names = F)
    
  } else {
    outfig <- NULL
  }
  
  ecmap <- GSECA.ECmap( gseca
                       , filename=outfig
                       , p_adj = PADJ
                       , psumlog=PSUMLOG
                       , pemp = PEMP
                       , srate=SRATE)

  res = list( 'gseca' = gseca, 
              'ECmap' = ecmap,
              'analysis' = analysis)
  
  return(res)
}

# PLOTTING ============
# Visualization of analysis results

g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  legend
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
  
  upViewport()
}



GSECA.ECmap = function( gseca
                        , filename=NULL
                        , p_adj=0.1
                        , psumlog=1
                        , pemp=NULL
                        , srate=NULL
                        , toprank=NULL
                                
){

  require(RColorBrewer)

  message("Thresholds: \nP-value adj ",p_adj,"\nP-value sumlog ",psumlog,"\n")
  
  sbs = subset(gseca, 
               p.adj<=p_adj 
               & sumlog<=psumlog)
  
  if(!is.na(pemp) && !is.na(unique(gseca$p.emp))){
    sbs = subset(sbs, p.emp<=pemp)
  }
  
  sr_bar = F
  if(!is.na(srate) && !is.na(unique(gseca$success_rate))){
    sr_bar = T
    sbs = subset(sbs, success_rate>=srate)
  }
  
  if( is.null(toprank) ) {
    sel = unique(sbs$gene_set)
    
    if( is.na(unique(gseca$p.emp)) && is.na(unique(gseca$success_rate))){
      sel <- unique(sbs[order(sbs$sumlog, decreasing = F), "gene_set"])[1:20]
      warning("Bootstrapping controls are missing; top 20 gene sets visualized in EC map.")
    }
    
    if(length(sel)==0| is.null(sel)) stop(message("No enriched gene sets"))
  } else {
    sel <- orderGSECA(gseca, PSUMLOG = psumlog, PEMP = pemp, SRATE = srate)$gene_set[1:toprank]
  }

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
  ECmap = as.matrix(ECmap)
  ECmap2 = as.matrix(ECmap2)

  ha = HeatmapAnnotation(name='Bimodal',col_mean = anno_bimodal, height=unit(1,"cm"))
  
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
  
  ECcol = rep("grey80",7)
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
                                id.lengths=rep(3, 1), gp=gpar(fill=myfill,col='black',lwd=unit(0.25,'points')))
                  }
                  }
                , cluster_rows = FALSE, cluster_columns = FALSE,show_row_names = T, show_column_names = T
                , row_names_side = 'left'
                , show_heatmap_legend = F
                , row_names_gp = gpar(fontsize = 8)
                , column_names_gp = gpar(fontsize = 8)
                , top_annotation = ha
  )
  
  
  if(sr_bar & length(unique(score[match(rownames(ECmap), score$gene_set), "success_rate"]))>1){
    row_ha = rowAnnotation(
      ES = row_anno_barplot(matrix(-10*log10(score[match(rownames(ECmap), score$gene_set), "sumlog"])), 
                            axis = TRUE,  border = T, axis_side = "bottom", 
                            gp = gpar(fill = c("black")), bar_width = 0.6)
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
    ,text = column_anno_text(as.character(round(limits)),gp = gpar(fontsize=8)) 
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
                                 id.lengths=rep(3, 1), gp=gpar(fill='white',col='black',lwd=unit(0.25,'points')))
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

  }
  return(ECM+row_ha)
}

# TIME CHECKING ===============
getTestStat <- function(test,tp) {
  nsamples <- length(unique(test$SampleID))
  tot <- Reduce(sum, lapply(tp, function(x) x['elapsed']))
  cat(" === Test Log === \n"
      ,"[+] Number of samples:", nsamples, "\n"
      ,"[+] Tot Time:", round(tot,2), " sec \n"
      ,"[+] FMM:"   , round(tp[[1]]['elapsed'], 2), "sec \n"
      ,"[+] GCR:"   , round(tp[[2]]['elapsed'], 2), "sec \n"
      ,"[+] GSECA:" , round(tp[[3]]['elapsed'], 2), "sec \n"
      ,"[+] P.Emp:" , round(tp[[4]]['elapsed'], 2), "sec \n"
      ,"[+] S.Size:", round(tp[[5]]['elapsed'], 2), "sec \n"
  )
}

