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
           ,'circlize'
)

not_installed <- which( ! pkgs %in% rownames(installed.packages()) )

if(any(grepl('ComplexHeatmap', not_installed))) {
  source("http://bioconductor.org/biocLite.R")
  biocLite('ComplexHeatmap')
}

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
  paste0("gene_sets/msigdb/",ix)
}

today=function() format(Sys.time(), "%y%m%d")

# STATS ===========
# Functions for computing statistics
get_padj <- function(cor) {
  function(x) {
    p.adjust(x, method=cor)
  }
}

summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE, conf.interval=.95, .drop=TRUE) {
  library(plyr)
  length2 <- function (x, na.rm=FALSE) {
    if (na.rm) sum(!is.na(x))
    else       length(x)
  }
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
  datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean
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

FISHER <- function (u, TAIL =  "two.sided"){
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
    set.seed(30580)
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



set_expression_class_on_posterior <- function(x
                                              , posterior_cutoff=.9
                                              , nClass=7)
{
  setFunc <- get(grep(nClass, ls(pos = 1), value = T))
  x <- setFunc(x, posterior_cutoff)
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

get_mixture_expr_class <- function(rna
                                   , ne_value=0.01
                                   , nClass=7
                                   , ncor=2){

  message(" === FMM and DD === \n",
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

  message("[*] DD Setting expression classes ...")
  l <- lapply(l, function(x) x[,!colnames(x)%in%c('mean1','mean2','st1','st2','lambda1','lambda2')])
  l <- suppressMessages(lapply_pb(l, set_expression_class_on_posterior, posterior_cutoff=.9, nClass=nClass))

  # rna <- do.call("rbind", l)

  return(list( 'expr_class'=l
               ,'mixture_model'=stats))
}

# GENE CLASS REPRESENTATION  ===========
# Representation of genes per expression classs
gene_class_representation <- function(r, pl, correction = 'fdr' , cutoff = 0.01

){
  require(plyr)
  require(metap)
  total_genes <- unique(unlist(pl))
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

  stats <- t(apply(cnts, 1, FISHER))

  # set pvalue == 0 as minimum measured value
  min_pv = min(stats[ which(stats[,1]>0) ,1 ], na.rm = T)

  if(min_pv>0.005) min_pv = 0.005

  stats[which(stats[,1]==0),1] = min_pv

  invariant$pv <- 1
  invariant$or <- 0

  cnts$pv <- stats[,1]
  cnts$or <- stats[,2]

  cnts <- rbind(cnts, invariant)
  cnts <- dlply(cnts, .(geneID), here(plyr::mutate)
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
GSECArunner <- function(gene_set, geneClass){
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

    stats <- t(apply(px, 1, FISHER))

    px$pv = stats[,1]
    px$or = stats[,2]

    return(px)
  } else {
    return(NA)
  }
}

GSECA_core <- function(gcr, pl, ncores=3, correction='fdr') {

  px <- mclapply(pl, GSECArunner
                 , geneClass = gcr
                 , mc.cores = ncores, mc.silent = T)

  n_na <-  length(which(is.na(px)))
  if ( n_na != 0 ) {
    warning(paste0(" -- Excluding ", n_na," gene sets with no genes in input matrix"))
    # nm <- names(which(is.na(px)))
    px[which(is.na(px))] <- NULL
  }

  px <- mapply(function(x,y){x$gene_set=y;return(x)}, px, names(px), SIMPLIFY = F)
  px <- as.data.frame(do.call('rbind',px))
  px <- unrowname(px)

  message("[*] GSECA scoring ...")

  # Fisher's method
  px <- dlply(px, .(gene_set), here(plyr::mutate)
              , p.adj  = get_padj(correction)(pv)
              , AS = suppressWarnings(sumlog(pv)$p)
              , .progress = 'text'
  )


  px <- lapply(px, function(x) x[,c("gene_set", "expr_class", "CASE.N.class",
                                    "CASE.N.rest", "CNTR.N.class", "CNTR.N.rest",
                                    "CASE.P", "CNTR.P", "delta", "direction", "pv",
                                    "or","p.adj","AS")])
  return(px)
}

GSECA <- function(pl, gcr, correction='fdr', cutoff = 0.01){

  message(" === Gene Set Enrichment Class Analysis === ")
  message("[*] GSECA core ...")

  core <- GSECA_core(gcr, pl, correction = correction)

  # Ranking
  core <- lapply(core, function(x) {
    to.ord = unique(x$AS)
    ord    = 1:length(to.ord)
    n.ord  = sort(to.ord, decreasing = F)
    rank   = ord[ match( to.ord , n.ord ) ]
    tmp    = data.frame(to.ord, rank)
    x$rank <- tmp$rank[match(x$AS,tmp$to.ord)]
    x <- x[order(x$rank, x$expr_class),]
    x$sig <- FALSE

    x$sig = x$p.adj<=cutoff

    x <- x[,c("gene_set","expr_class","CASE.N.class",
              "CASE.N.rest","CNTR.N.class","CNTR.N.rest",
              "CASE.P","CNTR.P", "delta","direction","pv",
              "or", "p.adj", "sig",'rank',"AS")]
    return(x)
    rm(tmp, ord, to.ord, n.ord, rank)
  })

  core <- unrowname(do.call(rbind, core))
  return(core)
}

# BOOTSTRAPPING ====
# Implementation of bootstrapping procedures
# [*] Calculation of empirical P value
# [*] Correction for sample size

execute_GSECA_bootstrap <- function(sample,  rna, pl, type='sample.size', correction="fdr", cutoff = 0.01){

  if(type=='p.empirical'){

    rna$mixture_model$type <- sample$type[match(as.character(rna$mixture_model$SampleID), sample$barcode)]

  }else if(type=='sample.size'){

    idx <- match(sample, names(rna$expr_class))
    rna$expr_class[-idx] <- NULL
    rna$mixture_model <- subset(rna$mixture_model, SampleID%in%sample)
    gc()
  }
  rm(sample)
  dd <- gene_class_representation(rna, pl, correction = correction)
  rm(rna)
  res <- GSECA(pl, dd, correction = correction, cutoff = cutoff )
  rm(dd)
  rx <- unique(res[,c("gene_set", "AS")])
  rm(res)
  return(rx)
}

GSECA.Bootstrap.empirical <- function( rna, obs.gseca, pl
                                       , NSIM=10
                                       , nc=2
                                       , correction = "fdr"
                                       , cutoff = 0.01
){

  boot <- unique(rna$mixture_model[,c('SampleID','type')])
  set.seed(30580)
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
    tryCatch(execute_GSECA_bootstrap(sample=x, rna=rna, pl=pl, type='p.empirical', correction=correction, cutoff=cutoff), error = function(x) return(NA))
  })
  sfStop()

  message("[*] Compute P.emp ...")
  l <- do.call(rbind.data.frame, lapply(l, function(x) unique(x[,c("gene_set", "AS")])))
  l$obs <- obs.gseca$AS[match(l$gene_set,obs.gseca$gene_set)]
  l <- ddply(l, .(gene_set), summarise, p.emp = (1+sum(AS<obs, na.rm = T))/(1+length(na.omit(AS))))
  return(l)
}

GSECA.Bootstrap.sample_size <- function( rna, obs.gseca, pl , phen = c("CASE","CNTR"), NSIM=1000, nc=3 , correction="fdr", cutoff=0.05){

  message("GSECA.Bootstrap Create Phenotype lists ...")

  pheno <- list()
  pheno[[1]] <- as.character(unique(subset(rna$mixture_model, type==phen[1])$SampleID))
  pheno[[2]] <- as.character(unique(subset(rna$mixture_model, type==phen[2])$SampleID))

  if(length(pheno[[1]])!=length(pheno[[2]])){

    ref.idx     <- which.max(c(length(pheno[[1]]), length(pheno[[2]])))
    size.idx    <- which.min(c(length(pheno[[1]]), length(pheno[[2]])))
    sample.size <- length(pheno[[size.idx]])

    message("GSECA.Bootstrap Bootstrapping ...")

    set.seed(30580)
    r <- lapply(1:NSIM, function(i) {
      sample_barcode = vector(mode = "numeric", length = 2*sample.size)
      sample_barcode[1:sample.size] = pheno[[size.idx]]
      sample_barcode[(sample.size+1):length(sample_barcode)] = sample(pheno[[ref.idx]]
                                                                      , length(pheno[[size.idx]]))
      return(sample_barcode)
    })


    functionsToCluster <- as.character(lsf.str(envir = .GlobalEnv))
    dataToCluster      <- ls()
    sfInit(parallel=TRUE, cpus=nc, type="SOCK")
    sfExport(list = c(functionsToCluster, dataToCluster))
    sfLibrary(parallel)
    sfClusterSetupRNG()
    l <- sfLapply(r, function(x){
      tryCatch(execute_GSECA_bootstrap(sample=x, rna=rna, pl=pl, type='sample.size', correction=correction,cutoff=cutoff), error = function(x) return(NA))
    })
    sfStop()

    message("[*] Compute Success Rate ...")
    l <- do.call(rbind.data.frame, lapply(l, function(x) unique(x[,c("gene_set", "AS")])))
    l$obs = obs.gseca$AS[match(l$gene_set,obs.gseca$gene_set)]
    l$sig.threshold = cutoff
    l <- ddply(l, .(gene_set), summarise
               , p.emp        = (1+sum(AS<obs, na.rm = T))/(1+length(na.omit(AS)))
               , success_rate = sum(AS<sig.threshold, na.rm = T)/length(na.omit(AS))
    )
  } else{
    print("Same number of samples per phenotype - Skip Sample size correction")
    l=NULL
  }
  return(l)
}


# UTILS =================

orderGSECA <- function(x, P.AS, P.EMP=NA, S.RATE=NA){
  # Rank pathways using 3 metrics:
  # AS, SR and emprical P-Value
  y = unique(x[,c('gene_set','AS', 'sr','p.emp')])
  y$method='GSECA'
  if(!is.na(P.EMP) & !is.na(S.RATE)){
    y$sig = y$AS<=P.AS & y$p.emp<=P.EMP & y$sr>=S.RATE
    y=y[order(1-y$sig, y$AS,  y$p.emp, 1-y$sr,  decreasing = F),]
    y$rank = 1:nrow(y)
  }else if(!is.na(P.EMP) & is.na(S.RATE) ){
    y$sig = y$AS<=P.AS & y$p.emp<=P.EMP
    y=y[order(1-y$sig, y$AS,  y$p.emp,  decreasing = F),]
    y$rank = 1:nrow(y)
  }else if(is.na(P.EMP) & !is.na(S.RATE)){
    y$sig = y$AS<=P.AS & y$sr>=S.RATE
    y=y[order(1-y$sig, y$AS, 1-y$sr,  decreasing = F),]
    y$rank = 1:nrow(y)
  } else {
    y$sig = y$AS<=P.AS
    y=y[order(1-y$sig, y$AS,  decreasing = F),]
    y$rank = 1:nrow(y)
  }
  return(y)
}

# GSECAordering <- function(x, ASCORE, PEMP=NA, SRATE=NA){
#
#   # Rank pathways using 3 metrics:
#   # AS, SR and emprical P-Value
#   # for GSECA only
#
#   y = unique(x[,c('gene_set','as','sr','p.emp')])
#
#   if(!is.na(PEMP) & !is.na(SRATE)){
#
#     y$sig = y$as<=ASCORE & y$p.emp<=PEMP & y$sr>=SRATE
#
#     y=y[order(1-y$sig, y$as,  y$p.emp, 1-y$sr,  decreasing = F),]
#
#     y$rank = 1:nrow(y)
#
#   }else if(!is.na(PEMP) & is.na(SRATE) ){
#
#     y$sig = y$as<=ASCORE & y$p.emp<=PEMP
#
#     y=y[order(1-y$sig, y$as,  y$p.emp,  decreasing = F),]
#     y$rank = 1:nrow(y)
#
#   }else if(is.na(PEMP) & !is.na(SRATE)){
#
#     y$sig = y$as<=ASCORE & y$sr>=SRATE
#
#     y=y[order(1-y$sig, y$as, 1-y$sr,  decreasing = F),]
#     y$rank = 1:nrow(y)
#
#   } else {
#     y$sig = y$as<=ASCORE
#
#     y=y[order(1-y$sig, y$as, decreasing = F),]
#     y$rank = 1:nrow(y)
#   }
#
#   return(y)
# }

# PLOTTING ============
# Visualization of analysis results

g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  legend
}

ytickform <- function(x){  lab <- sprintf("%05s",x) }

anno_bimodal = AnnotationFunction(
  fun = function(index, k, n) {
    pushViewport(viewport( xscale = c(0,7), yscale = c(0,1)))
    set.seed(30580)
    xseq<-seq(0, 7 ,.01)
    densities<-dnorm(xseq, 2,1)
    densities2<-dnorm(xseq, 5,1)

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
    popViewport()
  }
  , n=7
  , subsetable = TRUE
  , height = unit(1, "cm")
)


# EC mapper
GSECA.ECmap = function( gseca
                        , filename=NULL
                        , p_adj=0.1
                        , AS=1      # Association score threshold
                        , pemp=NULL # Empiricial p-value threshold
                        , SR=NULL   # Success rate threshold
                        , toprank=0

){

  require(RColorBrewer)
  message("Thresholds: \nP-value adj ",p_adj,"\nAssociation score ",AS,"\n")

  pemp  = ifelse(sum(is.na(unique(gseca$p.emp)))==0, pemp, NA)
  SR = ifelse(sum(is.na(unique(gseca$sr)))==0, SR, NA)
  sbs = subset( orderGSECA(subset(gseca,p.adj<=p_adj), AS, pemp, SR ) , sig)

  if(nrow(sbs)==0){
    message("NO altered gene sets")
    return(NULL)
  }

  sel = unique(sbs$gene_set)
  if(length(sel)==0| is.null(sel)) stop(message("No enriched gene sets"))

  if( toprank!=0 ) {
    toprank = ifelse(toprank>length(sel), length(sel), toprank)
    sel = sel[1:toprank]
  }

  toplot = subset(gseca, gene_set%in%sel)
  toplot$p.adj[ which(toplot$p.adj>p_adj) ]=NA
  toplot$gene_set = fix_names(toplot$gene_set)

  score = unique(toplot[,c('gene_set','AS','p.emp','sr')])

  gs_order = score[order(score[,2],decreasing = F),'gene_set']

  toplot$gene_set = factor(toplot$gene_set, levels=gs_order)
  score$gene_set  = factor(score$gene_set, levels=gs_order)

  ECmap  = dcast(data = toplot, gene_set~expr_class, value.var = 'p.adj' ); rownames(ECmap)  = ECmap[,1];  ECmap  = ECmap[,-1]
  ECmap2 = dcast(data = toplot, gene_set~expr_class, value.var = 'delta' ); rownames(ECmap2) = ECmap2[,1]; ECmap2 = ECmap2[,-1]
  ECmap3 = dcast(data = toplot, gene_set~expr_class, value.var = 'CASE.P' ); rownames(ECmap3) = ECmap3[,1]; ECmap3 = ECmap3[,-1]
  # ECmap4 = dcast(data = toplot, gene_set~expr_class, value.var = 'CNTR.P' ); rownames(ECmap4) = ECmap4[,1]; ECmap4 = ECmap4[,-1]

  limits = range(abs(ECmap2))
  ECmap  = as.matrix(ECmap)
  ECmap2 = as.matrix(ECmap2)


  # Annotations
  ha = HeatmapAnnotation(col_mean = anno_bimodal, show_legend = F, show_annotation_name = F)
  if(!is.na(SR) & length(unique(score[match(rownames(ECmap), score$gene_set), "sr"]))>1){
    row_ha = rowAnnotation(
      AS = anno_barplot(matrix(-10*log10(score[match(rownames(ECmap), score$gene_set), "AS"]))
                        , axis = TRUE,  border = T, axis_param = list(side = "bottom")
                        , gp = gpar(fill = c("black")), bar_width = 0.6)
      , SR = anno_barplot(matrix(score[match(rownames(ECmap), score$gene_set), "sr"])
                          , axis = TRUE, border = T, axis_param = list(side = "bottom")
                          , gp = gpar(fill = c("black")), bar_width = 0.6)
      , gap = unit(c(2, 4), "mm")
      , width = unit(2.5, "cm")
      , show_annotation_name = TRUE
      , annotation_name_side = 'top'
      , annotation_name_rot = 0
      , annotation_name_gp = gpar(fontsize= 8)
    )
  } else {
    row_ha = rowAnnotation(
      AS = anno_barplot(matrix(-10*log10(score[match(rownames(ECmap), score$gene_set), "AS"]))
                        , axis = TRUE,  border = T, axis_param = list(side = "bottom")
                        , gp = gpar(fill = c("black")), bar_width = 0.6)
      , gap = unit(c(2, 4), "mm")
      , width = unit(1.25, "cm")
      , show_annotation_name = TRUE
      , annotation_name_side = 'top'
      , annotation_name_rot = 0
      , annotation_name_gp = gpar(fontsize= 8)
    )
  }

  pycol   = c(rgb(252,231,62, max=255),rgb(160,218,69, max=255),rgb(77,192,111, max=255),rgb(37,160,134, max=255),
              rgb(41,125,140, max=255),rgb(55,89,137, max=255),rgb(70,51,124, max=255),rgb(70,31,91, max=255),rgb(70,31,104, max=255))
  col_fun = colorRamp2(c(0,0.00001,0.0001, 0.05, 0.1),pycol[c(1,3,4,7,9)])

  ECcol = rep("grey80",7)
  ECM = Heatmap(ECmap, name = "ECMap"
                , col =  col_fun
                , rect_gp = gpar(col = 'grey80', fill=NA )
                , cluster_rows = FALSE, cluster_columns = FALSE
                , show_row_names = T, show_column_names = T
                , row_names_side = 'left'
                , show_heatmap_legend = T
                , heatmap_legend_param = list(title="FDR",at=c(0, 0.01, 0.05, 0.1)) #, direction = "horizontal"
                , row_names_gp = gpar(fontsize = 8)
                , column_names_gp = gpar(fontsize = 8)
                , top_annotation = ha
                , right_annotation = row_ha
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

  )

  limits = range(abs(ECmap2)*100)
  ml = matrix(limits, nr=1, dimnames = list("Delta (%)", as.character(round(limits))))

  lgH = Heatmap(ml
                , name = "Legend"
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
                , cluster_rows = FALSE, cluster_columns = FALSE
                , show_row_names = F, show_column_names = F
                , column_names_side = 'bottom'
                , show_heatmap_legend = F
                , width = unit(1.5, 'cm')
                , height = unit(1.5, 'cm')
                , row_names_side = 'left'
                , row_names_gp = gpar(fontsize = 8)
                , bottom_annotation =  HeatmapAnnotation(
                  which='column'
                  ,text = anno_text(as.character(round(limits)),gp = gpar(fontsize=8), rot=0, just ='center', location=0.5)
                )
                , top_annotation =  HeatmapAnnotation(
                  which='column'
                  ,text = anno_text(c("Delta","(%)"),gp = gpar(fontsize=8), rot=0, just ='center', location=0.5)
                )


  )

  if(!is.null(filename)){
    pdf(file=filename,  paper = "a4", useDingbats = F)
    grid.newpage()

    l=grid.layout(  nrow=4, heights = c(.05, .8, .1, .05)
                    , ncol=3, widths  = c( .05, .9, .05))
    # grid.show.layout(l)

    marginT <- viewport( layout.pos.row = 1,   layout.pos.col = 1:3, name = "marginT")
    marginL <- viewport( layout.pos.row = 1:4, layout.pos.col = 1,   name = "marginL")
    marginR <- viewport( layout.pos.row = 1:4, layout.pos.col = 3,   name = "marginR")
    marginB <- viewport( layout.pos.row = 4,   layout.pos.col = 1:3, name = "marginB")
    plot   <- viewport(  layout.pos.row = 2,   layout.pos.col = 2,   name = "ECMAP")
    legend <- viewport(  layout.pos.row = 3,   layout.pos.col = 2,   name = "delta")
    logo   <- viewport(  layout.pos.row = 4,   layout.pos.col = 2,   name = "logo")

    top.vp = viewport(layout=l)
    splot <- vpTree(top.vp, vpList(marginT, marginL, marginR, marginB, plot, legend, logo))
    pushViewport(splot)

    seekViewport("ECMAP")
    ComplexHeatmap::draw(ECM, newpage = FALSE) # heatmap_legend_side = "bottom" heatmap_legend_side='bottom',

    seekViewport("delta")
    ComplexHeatmap::draw(lgH, newpage = FALSE)

    seekViewport("logo")
    grid.text("GSECA 2019", just="center", gp=gpar(fontsize=6, col='black'))
    upViewport(0)

    dev.off()

  }
  return(ECM)
}

# EXECUTORs =========
# [*] Main Executor
# Wraps and execute all functions
# needed for the analysis
# M = Gene Expression matrix
# L = Sample label list
# geneset = gene set list
# symbol =  gene symbol (i.e. ensembl_gene_id)
# p_adj_th = adjusted p-value threshold
# analysis = analysis name
# outdir = outdir folder
# N.CORES = number of cores
# EMPIRICAL = true if empirical p-value is requested
# BOOTSTRP = true if success rate is requested
# nsim = number of bootstrapping
# AS = AS threshold
# PEMP = p.emp threshold
# SR = success rate threshold
# toprank = number of topranked gene sets
# iphen =  phenotype lables ("CASE","CNTR")
GSECA_executor <- function( M
                            , L
                            , geneset
                            , symbol="ensembl_gene_id"
                            , correction='fdr'
                            , p_adj_th = 0.1
                            , analysis = NULL
                            , outdir = NULL
                            , N.CORES = 3
                            , EMPIRICAL = F
                            , BOOTSTRP = F
                            , nsim = 10
                            , AS = 0.01
                            , PEMP    = 1
                            , SR   = 0.7
                            , toprank = 0
                            , iphen = c("CASE","CNTR")
){

  print.logo()
  if( N.CORES > nsim ) N.CORES <- nsim

  # 01. Prepare Data
  message("[*] loading datasets ...")
  expr <- get_expression_dataset(M, L, symbol)
  geneset <- prepare_pl(geneset, expr)

  message("[*] GSECA running ...")
  expr  = get_mixture_expr_class(expr, ne_value = 0.01 , nc  = N.CORES )
  gcr   = gene_class_representation(expr, pl = geneset, correction = correction)
  gseca = GSECA(  pl=geneset, gcr, correction = correction )

  gseca$p.emp = NA
  gseca$sr = NA

  if(EMPIRICAL){
    message(" [*] Empirical P-Value estimation ... ")
    gseca.emp <- GSECA.Bootstrap.empirical( expr, gseca, pl = geneset, NSIM = nsim, nc = N.CORES)
    gseca$p.emp <- gseca.emp$p.emp[match(gseca$gene_set, gseca.emp$gene_set)]

  }
  if(BOOTSTRP){
    message("[*] Success Rate for sample size ... ")
    gseca.boot <- GSECA.Bootstrap.sample_size( expr, gseca, pl = geneset, phen = iphen, NSIM = nsim, nc = N.CORES)
    if(!is.null(gseca.boot)) gseca$sr <- gseca.boot$success_rate[match(gseca$gene_set,gseca.boot$gene_set)]
  }
  gseca$gene_set <- fix_names(gseca$gene_set)

  if( !is.null(analysis) & !is.null(outdir)  ) {
    analysis <- create.output.folder(analysis, outdir)
    outfig <- paste0(analysis,"/ECmap.as_",AS, "_padj_",p_adj_th, "_pemp_",PEMP, "_srate_",SR,".pdf")
    message(paste0("[*] Saving files in folder: ", analysis))

    saveRDS(expr, file=paste0(analysis,"/expr.rds"))
    saveRDS(gseca, file=paste0(analysis,"/gseca.rds"))
    write.csv(gseca, file=paste0(analysis,"/gseca.csv"), row.names = F)


  } else {
    outfig <- NULL
  }

  message("[*] EC Map generation ...")
  ecmap <- tryCatch(GSECA.ECmap( gseca, filename=outfig , p_adj = p_adj_th, AS=AS, pemp = PEMP , SR=SR, toprank = toprank ), error=function(e) return(NULL))
  res = list( 'gseca' = gseca,'ECmap' = ecmap,'analysis' = analysis)
  return(res)
}
